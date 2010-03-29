! Module for parsing command line options and arguments similar to python
! package OptionParser.
! A command line is a blank-delimited list of arguments and options.
! An option consists of
!    - a long name (ex: '--long-option')
!    - a short name (ex: '-l')
!    - a description, printed if the program is called with -h or --help
!    - whether or not it is followed by a value
!    - an action (so far, only 'store true' and 'store false')
!    - a default value as character
!
! Author: P. Chanial

module module_optionparser

    use, intrinsic :: ISO_FORTRAN_ENV
    use            :: string, only : strlowcase, strinteger
    implicit none

    integer, private, parameter :: ARGUMENT_MAX = 128
    integer, private, parameter :: ARGUMENT_VALUE_LEN = 2048
    integer, private, parameter :: OPTION_MAX = 64
    integer, private, parameter :: OPTION_NAME_LEN = 32
    integer, private, parameter :: OPTION_DESCRIPTION_LEN = 80
    integer, private, parameter :: OPTION_VALUE_LEN = 256
    integer, private, parameter :: OPTION_ACTION_LEN = 16

    type, private :: option
        character(len=OPTION_NAME_LEN)        :: long
        character                             :: short
        character(len=OPTION_DESCRIPTION_LEN) :: description
        logical                               :: has_value
        character(len=OPTION_VALUE_LEN)       :: value
        character(len=OPTION_VALUE_LEN)       :: default
        character(len=OPTION_ACTION_LEN)      :: action
    end type option

    type, public :: optionparser

        character(len=80), private   :: usage
        character(len=2048), private :: command
        integer, private             :: narguments
        integer, private             :: noptions
        integer, private             :: narguments_min
        integer, private             :: narguments_max
        integer, private             :: argument(ARGUMENT_MAX)
        type(option), private        :: option(OPTION_MAX)

    contains

        procedure, public  :: init
        procedure, public  :: add_option
        procedure, public  :: parse
        procedure, public  :: reset
        procedure, public  :: get_option
        procedure, public  :: get_option_as_logical
        procedure, public  :: get_option_as_integer
        procedure, public  :: get_option_as_real
        procedure, public  :: get_option_count
        procedure, public  :: get_argument
        procedure, public  :: get_arguments
        procedure, public  :: get_argument_count
        procedure, public  :: print_options

        procedure, private :: parse_option
        procedure, private :: parse_argument
        procedure, private :: set_option
        procedure, private :: get_option_number
        procedure, private :: scan_command
        procedure, private :: print_usage

    end type optionparser


contains


    subroutine init(this, usage, narguments_min, narguments_max)
        class(optionparser), intent(inout) :: this
        character(len=*), intent(in)       :: usage
        integer, intent(in)                :: narguments_min
        integer, intent(in), optional      :: narguments_max

        this%narguments = 0
        this%noptions = 0
        
        ! command usage: 'mycommand [options] source target'
        this%usage = usage

        ! minimum number of arguments
        this%narguments_min = narguments_min

        ! maximum number of arguments -1: unlimited
        if (present(narguments_max)) then
            this%narguments_max = narguments_max
        else
            this%narguments_max = narguments_min
        end if

    end subroutine init


    !---------------------------------------------------------------------------


    subroutine add_option(this, long, short, description, has_value, default, action, status)
        class(optionparser), intent(inout)     :: this
        character(len=*), intent(in), optional :: long
        character(len=*), intent(in), optional :: short
        character(len=*), intent(in), optional :: description
        logical, intent(in), optional          :: has_value
        character(len=*), intent(in), optional :: default
        character(len=*), intent(in), optional :: action
        integer, intent(out), optional         :: status
        type(option)                           :: myoption
        character(len=256)                     :: buffer

        if (present(status)) status = 1

        ! check for overflow
        if (this%noptions == OPTION_MAX) then
            write (ERROR_UNIT,'(a,i0,a)') 'Maximum number of options is &
                  &reached (', OPTION_MAX, ')'
            if (.not. present(status)) stop 'Aborting.'
            return
        end if

        ! an option must have a short or a long name
        if (.not. present(long) .and. .not. present(short)) then
            write (ERROR_UNIT,'(a)') 'Option must have a short name or a long name.'
            if (.not. present(status)) stop 'Aborting.'
            return
        end if

        ! description is mandatory
        myoption%description = description
       
        ! option has a long name?
        if (present(long)) then
            myoption%long = long
        else
            myoption%long = ''
        end if

        ! option has a short name?
        myoption%short = ' '
        if (present(short)) then
            if (len_trim(short) > 1) then
                write (ERROR_UNIT,'(a)') "Short name '" // trim(short) // &
                      "' must not contain more than one character."
                if (.not. present(status)) stop 'Aborting.'
                return
            end if
            myoption%short = short
        end if

        ! option accept a value ?
        if (present(has_value)) then
            myoption%has_value = has_value
        else
            myoption%has_value = .false.
        end if

        ! option has action value ?
        if (present(action)) then
            buffer = strlowcase(action)
            if (buffer == 'store true' .or. buffer == 'store false') then
                myoption%action = buffer(1:len(myoption%action))
            else
                write (ERROR_UNIT, '(a)') "Invalid action '" // trim(buffer) // "'."
                if (.not. present(status)) stop 'Aborting.'
                return
            end if
        else
            if (myoption%has_value) then
                myoption%action = ''
            else
                myoption%action = 'store true'
            end if
        end if

        ! option has default value ?
        if (present(default)) then
            myoption%default = default
        else 
            myoption%default = ''
        end if

        ! store the option
        this%noptions = this%noptions + 1
        this%option(this%noptions) = myoption

        if (present(status)) status = 0

    end subroutine add_option


    !---------------------------------------------------------------------------


    subroutine parse(this, status, command)
        class(optionparser), intent(inout)     :: this
        integer, intent(out)                   :: status
        character(len=*), intent(in), optional :: command
        character(len=2048)                    :: buffer
        logical                                :: no_more_options
        integer                                :: count, current, length

        if (present(command)) then
            this%command = command
            write(OUTPUT_UNIT,'(a)') 'Command: ' // trim(this%command)
        else
            call get_command(buffer)
            count = index(buffer, ' ')
            if (count == 0) then
                this%command = ''
            else
                this%command = buffer(count+1:)
            end if
        end if

        ! set default values
        do current = 1, this%noptions
            if (this%option(current)%default /= '') then
                this%option(current)%value = this%option(current)%default
            else if (this%option(current)%action == 'store true') then
                this%option(current)%value = 'False'
            else if (this%option(current)%action == 'store false') then
                this%option(current)%value = 'True'
            else
                this%option(current)%value = ''
            end if
        end do

        status = 1
        no_more_options = .false.
        current = 1

        ! loop over the command arguments
        do
            
           call this%scan_command(current, buffer, length, status)
           if (status /= 0 .and. status /= -1) return
           if (status == -1) exit ! not found
           
           ! print usage
           if (buffer == '-h' .or. buffer == '--help') then
               call this%print_usage()
               call this%reset()
               status = -1
               return
           end if

           ! if '--' is encountered, all that follows is considered as arguments
           if (buffer == '--') then
               no_more_options = .true.
               current = current + 1
               cycle
           end if
           
           ! option
           if (buffer(1:1) == '-' .and. .not. no_more_options) then
               call this%parse_option(current, trim(buffer(2:)), status)
               if (status /= 0) then
                   call this%reset()
                   return
               end if

           ! argument
           else
               call this%parse_argument(current)
           end if

           current = current + 1

        end do

        ! parsing is successful
        if (this%narguments >= this%narguments_min .and. (&
            this%narguments <= this%narguments_max .or.   &
            this%narguments_max == -1)) then
            status = 0
            return
        end if

        ! not enough or two many arguments
        call this%print_usage()
        call this%reset()

    end subroutine parse
 

    !---------------------------------------------------------------------------


    subroutine print_usage(this)
        class(optionparser), intent(in) :: this
        character(len=OPTION_NAME_LEN+6):: buffer
        integer                         :: current, option_len

        ! not enough arguments or -h, --help option is specified
        write (OUTPUT_UNIT,'(a,2/,a)') 'Usage: ' // trim(this%usage), 'Options:'

        ! loop over the options to get the longest long name
        option_len = 0
        do current = 1, this%noptions
           option_len  = max(option_len, len_trim(this%option(current)%long))
        end do
        if (option_len /= 0) option_len = option_len + 2 ! --

        if (any(this%option(1:this%noptions)%short /= ' ')) then
            option_len = option_len + 4
        end if
        
        do current = 1, this%noptions
            if (this%option(current)%short /= ' ' .and. this%option(current)%long /= '') then
                buffer = '-' // this%option(current)%short // ', --' // trim(this%option(current)%long)
            else if (this%option(current)%short /= ' ') then
                buffer = '-' // this%option(current)%short
            else
                buffer = '--' // trim(this%option(current)%long)
            end if
            write (OUTPUT_UNIT,'(2x,a'//strinteger(option_len+1)//',x,a,$)') &
                  buffer, trim(this%option(current)%description)
            if (len_trim(this%option(current)%default) > 0) then
                write (OUTPUT_UNIT,'(a)')  &
                      ' [Default: ' // trim(this%option(current)%default) // ']'
            else
                write (OUTPUT_UNIT,*)
            end if
        end do

    end subroutine print_usage


    !---------------------------------------------------------------------------


    subroutine print_options(this)
        class(optionparser), intent(in) :: this
        character(len=2048)             :: buffer
        integer                         :: i, status

        ! print arguments
        if (this%narguments == 0) then
            write (OUTPUT_UNIT,'(a)') 'No arguments'
        end if
        do i = 1, this%narguments
            buffer = this%get_argument(i, status)
            write (OUTPUT_UNIT,'(a,i0,a)') 'Argument ', i, ': ' // trim(buffer)
        end do

        ! print options
        if (this%noptions == 0) then
            write (OUTPUT_UNIT,'(a)') 'No options'
        else
            write (OUTPUT_UNIT,'(a)') 'Options:'
            do i = 1, this%noptions
               if (this%option(i)%long /= '') then
                   write (OUTPUT_UNIT, '(2x,a,$)') trim(this%option(i)%long)
               else
                   write (OUTPUT_UNIT, '(2x,a,$)') trim(this%option(i)%short)
               end if
               write (OUTPUT_UNIT,'(a)') "='" // trim(this%option(i)%value)//"'"
            end do
        end if
    
    end subroutine print_options


    !---------------------------------------------------------------------------


    subroutine parse_option(this, current, buffer, status)
        class(optionparser), intent(inout) :: this
        integer, intent(inout)             :: current
        character(len=*), intent(in)       :: buffer
        integer, intent(out)               :: status
        integer                            :: i, number

        status = 1
        if (len(buffer) == 0) then
            write (ERROR_UNIT,'(a)') "Invalid option '-'."
            return
        endif

        ! test if the option has a long name '--option'
        if (buffer(1:1) == '-') then

            call this%set_option(current, buffer(2:), status)
            if (status /= 0) return

        ! the options are of the form, ex:  '-avz'
        else
        
            ! combined short options cannot accept a value
            if (len(buffer) > 1) then
                do i = 1, len(buffer)
                    number = this%get_option_number(buffer(i:i))
                    if (number == 0) cycle
                    if (this%option(number)%has_value) then
                        status = 1
                        write (ERROR_UNIT, '(a)') "Options '-"//buffer//"' cannot accept a value."
                        return
                    end if
                end do
            end if

            ! loop over the short options
            do i = 1, len(buffer)
                call this%set_option(current, buffer(i:i), status)
                if (status /= 0) return
            end do

        end if

    end subroutine parse_option


    !---------------------------------------------------------------------------


    subroutine set_option(this, current, buffer, status)
        class(optionparser), intent(inout) :: this
        integer, intent(inout)             :: current
        character(len=*), intent(in)       :: buffer
        integer, intent(out)               :: status
        integer                            :: pos, number, length

        pos = index(buffer, '=')
        if (pos == 0) then
            pos = len(buffer)+1
        end if

        number = this%get_option_number(buffer(:pos-1))
        if (number == 0) then
            status = 1
            write (ERROR_UNIT,'(a)') "Invalid option '" // buffer(:pos-1) // "'."
            return
        end if

        if (this%option(number)%has_value) then

            ! the option value can be '--option=value'...
            if (pos <= len(buffer)) then

                this%option(number)%value = adjustl(buffer(pos+1:))

            ! ...or '--option value'
            else

                current = current + 1
                call this%scan_command(current, this%option(number)%value, length, status)
                if (status /= 0 .and. status /= -1) return
                if (status == -1) then
                    write (*,*) "Option '" // buffer(:pos-1) // "' expects a value, but none is provided."
                    return
                endif

            end if
       
        else
            
            ! the option has no value, i.e. it is a boolean
            if (this%option(number)%action == 'store true') then

                this%option(number)%value = 'True'

            else

                this%option(number)%value = 'False'

            end if

        end if

        status = 0

    end subroutine set_option


    !---------------------------------------------------------------------------


    subroutine parse_argument(this, number)
        class(optionparser), intent(inout) :: this
        integer, intent(in)                :: number

        this%narguments = this%narguments + 1
        this%argument(this%narguments) = number

    end subroutine parse_argument


    !---------------------------------------------------------------------------


    subroutine reset(this)
        class(optionparser), intent(inout) :: this
        integer                            :: i
 
        do i = 1, this%noptions
            this%option(i)%value = ' '
        end do
        this%narguments = 0

    end subroutine reset


    !---------------------------------------------------------------------------


    function get_argument(this, number, status) result(value)
        class(optionparser), intent(in)   :: this
        integer, intent(in)               :: number
        integer, intent(out), optional    :: status
        character(len=ARGUMENT_VALUE_LEN) :: value
        integer                           :: status2, length

        if (present(status)) status = 1
        if (number <= 0) then
            write (ERROR_UNIT,'(a)') 'GET_ARGUMENT: invalid non-positive number.'
            if (.not. present(status)) stop 'Aborting.'
            return
        end if
        if (number > this%narguments) then
            write (ERROR_UNIT,'(a)') 'GET_ARGUMENT: argument number is greater than actual number of arguments.'
            if (.not. present(status)) stop 'Aborting.'
            return
        end if

        call this%scan_command(this%argument(number), value, length, status2)
        if (status2 == -1) then
            write (ERROR_UNIT,'(a)') 'GET_ARGUMENT: unknown error.'
        end if
        if (status2 /= 0) then
            if (.not. present(status)) stop 'Aborting.'
            return
        end if

        if (present(status)) status = 0

    end function get_argument


    !---------------------------------------------------------------------------


    subroutine get_arguments(this, arguments, status)
        class(optionparser), intent(in)                :: this
        character(len=ARGUMENT_VALUE_LEN), allocatable :: arguments(:)
        integer, intent(out), optional                 :: status
        integer :: narguments, iargument

        narguments = this%get_argument_count()
        allocate(arguments(narguments))
        do iargument = 1, narguments
            arguments(iargument) = this%get_argument(iargument, status)
            if (status /= 0) return
        end do

    end subroutine get_arguments


    !---------------------------------------------------------------------------


    function get_argument_count(this)
        class(optionparser), intent(in) :: this
        integer                         :: get_argument_count
        
        get_argument_count = this%narguments

    end function get_argument_count


    !---------------------------------------------------------------------------


    ! if number exceeds the actual number of arguments, status = -1, length=0
    subroutine scan_command(this, number, value, length, status)
        class(optionparser), intent(in) :: this
        integer, intent(in)             :: number
        integer, intent(out)            :: length
        character(len=*), intent(out)   :: value
        integer, intent(out)            :: status
        integer                         :: ipos, narg, ivalue
        character                       :: quote

        status = -1 ! not found
        length = 0
        ipos = 0
        narg = 0
        ivalue = 1
        value = ' '

        ! loop over the arguments
        do

            ! find first non-blank character
            do

                ipos = ipos + 1
                
                ! not found
                if (ipos > len(this%command)) return

                if (this%command(ipos:ipos) == ' ') cycle

                if (this%command(ipos:ipos) == "'") then
                    quote = "'"
                    ipos = ipos + 1
                else if (this%command(ipos:ipos) == '"') then
                    quote = '"'
                    ipos = ipos + 1
                else
                    quote = ' '
                end if
                exit

            end do
          
            ! find last character of argument
            narg = narg + 1
            do while (ipos <= len(this%command))

                if (this%command(ipos:ipos) == quote) go to 200

                ! look for character '\'=achar(92) (emacs fontifyer bug...)
                if (quote == ' ' .and. this%command(ipos:ipos) == achar(92)) then
                    if (ipos == len(this%command)) then
                        write (ERROR_UNIT,'(a)') "Command should not end with '\\'."
                        status = 1
                        return
                    end if
                    ipos = ipos + 1
                end if

                if (narg == number) then
                    if (ivalue > len(value)) then
                        write (ERROR_UNIT,'(a)') "Option value length is too &
                              &large: '" // value // "'."
                        status = 1
                        return
                    end if
                    value(ivalue:ivalue) = this%command(ipos:ipos)
                    ivalue = ivalue + 1
                end if

                ipos = ipos + 1
                     
            end do

            ! we haven't found the ending quote
            if (quote /= ' ') then
                write (ERROR_UNIT,'(a)') 'Command line does not terminate with proper quote ' // quote // '.'
                status = 1
                return
            end if

            ! check if we have found the requested argument
        200 if (narg == number) then
                length = ivalue - 1
                status = 0
                return
            end if
 
           ! check if we have reached the end of the command
            if (ipos > len(this%command)) return
                 
        end do

    end subroutine scan_command


    !---------------------------------------------------------------------------


    function get_option(this, name, status)
        class(optionparser), intent(in) :: this
        character(len=*), intent(in)    :: name
        integer, intent(out), optional  :: status
        character(len=OPTION_VALUE_LEN) :: get_option
        integer                         :: number
        
        number = this%get_option_number(name)
        if (number == 0) then
            write (ERROR_UNIT, '(a)') "Invalid option '" // trim(name) // "'."
            if (.not. present(status)) stop 'Aborting.'
            status = 1
            return
        end if

        if (present(status)) status = 0
        get_option = this%option(number)%value

    end function get_option


    !---------------------------------------------------------------------------


    function get_option_as_logical(this, name, status) result (value)
        class(optionparser), intent(in) :: this
        character(len=*), intent(in)    :: name
        integer, intent(out), optional  :: status
        logical                         :: value
        character(len=OPTION_VALUE_LEN) :: charvalue
        
        value = .false.
        charvalue = this%get_option(name, status)
        if (present(status)) then
            if (status /= 0) return
        end if

        if (charvalue == 'True') then
            value = .true.
        else if (charvalue /= 'False') then
            write (ERROR_UNIT,'(a)') "Invalid logical option value '" // &
                  trim(charvalue) // "'."
            if (.not. present(status)) stop 'Aborting.'
            status = 1
            return
        end if

    end function get_option_as_logical


    !---------------------------------------------------------------------------


    function get_option_as_integer(this, name, status) result (value)
        class(optionparser), intent(in) :: this
        character(len=*), intent(in)    :: name
        integer, intent(out), optional  :: status
        integer                         :: value
        character(len=OPTION_VALUE_LEN) :: charvalue
        integer                         :: iostatus
        
        value = 0
        charvalue = this%get_option(name, status)
        if (present(status)) then
            if (status /= 0) return
        end if

        read (charvalue,'(i20)',iostat=iostatus) value
        if (iostatus /= 0) then
            write (ERROR_UNIT,'(a)') "Invalid integer option value '" // &
                  trim(charvalue) // "'."
            if (.not. present(status)) stop 'Aborting.'
            status = 1
            return
        end if

    end function get_option_as_integer


    !---------------------------------------------------------------------------


    function get_option_as_real(this, name, status) result (value)
        class(optionparser), intent(in) :: this
        character(len=*), intent(in)    :: name
        integer, intent(out), optional  :: status
        real*8                          :: value
        character(len=OPTION_VALUE_LEN) :: charvalue
        integer                         :: iostatus
        
        !XXX should be NaN
        value = 0.d0
        charvalue = this%get_option(name, status)
        if (present(status)) then
            if (status /= 0) return
        end if

        read (charvalue,'(bn,f20.0)',iostat=iostatus) value
        if (iostatus /= 0) then
            write (ERROR_UNIT,'(a)') "Invalid real option value '" // &
                  trim(charvalue) // "'."
            status = 1
            if (.not. present(status)) stop 'Aborting.'
            return
        end if

    end function get_option_as_real


    !---------------------------------------------------------------------------


    function get_option_count(this)
        class(optionparser), intent(in) :: this
        integer                         :: get_option_count

        get_option_count = this%noptions

    end function get_option_count


    !---------------------------------------------------------------------------


    function get_option_number(this, name) result(number)
        class(optionparser), intent(in) :: this
        character(len=*), intent(in)    :: name
        integer                         :: number

        do number = this%noptions, 1, -1
            if (this%option(number)%long  == name .or. &
                this%option(number)%short == name) return
        end do
        
    end function get_option_number


end module module_optionparser
