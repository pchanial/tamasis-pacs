!-------------------------------------------------------------------------------
!
! This module defines the pacsobservation derived type.
! It contains all the information of a PACS observation relevant to its data
! processing. In particular:
!    - file name of the observation
!    - channel ('b', 'g' or 'r')
!    - transparent_mode (false or true)
!    - section to be processed
!
! Author: Pierre Chanial
!
!-------------------------------------------------------------------------------

module module_pacsobservation

    use iso_fortran_env,  only : ERROR_UNIT, OUTPUT_UNIT
    use module_fitstools, only : ft_open, ft_open_bintable, ft_read_column, ft_read_extension, ft_close
    use module_precision, only : p
    use module_string,    only : strinteger, strlowcase, strsection, strternary
    implicit none
    private

    public :: pacsobservation
    public :: pacsobsinfo

    type pacsobsinfo
        character(len=256) :: filename
        character          :: channel
        logical            :: transparent_mode
        integer            :: compression_factor
        integer*8          :: nsamples
        integer*8          :: first
        integer*8          :: last
    end type pacsobsinfo

    type pacsobservation

        private
        character, public :: channel
        logical, public   :: transparent_mode
        type(pacsobsinfo), allocatable, public :: info(:)

    contains

        private
        procedure, public :: init

    end type pacsobservation

contains


    subroutine init(this, filename, status, verbose)

        class(pacsobservation), intent(inout) :: this
        character(len=*), intent(in)          :: filename(:)
        integer, intent(out)                  :: status
        logical, intent(in), optional         :: verbose
        integer*8                             :: first, last
        integer                               :: iobs

        ! parameter checking
        if (size(filename) == 0) then 
            status = 1
            write (ERROR_UNIT,'(a)')'INIT_PACSOBSERVATION: Array has zero size.'
            return
        end if
        
        if (allocated(this%info)) deallocate(this%info)
        allocate(this%info(size(filename)))

        do iobs = 1, size(filename)
            call set_filename(this%info(iobs), filename(iobs), first, last,    &
                              status)
            if (status /= 0) return

            call set_valid_slice(this%info(iobs), first, last, status)
            if (status /= 0) return

            this%info(iobs)%nsamples = this%info(iobs)%last -                  &
                                       this%info(iobs)%first + 1

            call set_channel(this%info(iobs), status)
            if (status /= 0) return

            call set_compression_mode(this%info(iobs), status)
            if (status /= 0) return
        
            call set_transparent_mode(this%info(iobs), status)
            if (status /= 0) return

        end do

        ! make sure the channel is the same for all observations
        this%channel = this%info(1)%channel
        if (any(this%info%channel /= this%channel)) then
            status = 1
            write (ERROR_UNIT,'(a)') 'ERROR: Observations do not have the same channel.'
            return
        end if

        ! make sure the transparent mode is the same for all observations
        ! that could be relaxed, but we would need to properly handle the bad detector mask
        this%transparent_mode = this%info(1)%transparent_mode
        if (any(this%info%transparent_mode .neqv. this%transparent_mode)) then
            status = 1
            write (ERROR_UNIT,'(a)') 'ERROR: Observations do not have the same transparent mode.'
            return
        end if

        if (.not. present(verbose)) return
        if (.not. verbose) return
        
        ! print some info
        do iobs = 1, size(filename)

            ! observation number & file name
            write (OUTPUT_UNIT,'(a)') 'Info: Observation' // strternary(       &
                size(filename)>1, ' ' // strinteger(iobs), '') // ': ' //      &
                trim(this%info(iobs)%filename)
            write (OUTPUT_UNIT,'(a)') '      Section: [' //     &
                strsection(this%info(iobs)%first,this%info(iobs)%last) // ']'
            
            ! channel
            write (OUTPUT_UNIT,'(a,$)') "      Channel: "
            select case (this%info(iobs)%channel)
                case ('b')
                    write (OUTPUT_UNIT,'(a)') 'Blue'
                case ('g')
                    write (OUTPUT_UNIT,'(a)') 'Green'
                case ('r')
                    write (OUTPUT_UNIT,'(a)') 'Red'
                case default
                    write (OUTPUT_UNIT,'(a)') 'Unknown'
            end select

            ! compression mode
!            write (OUTPUT_UNIT,'(a)') "Info: Compression mode: '" // trim(compression) // "'"

            ! compression factor
            write (OUTPUT_UNIT,'(a,i0)') "      Compression factor: ",         &
                 this%info(iobs)%compression_factor

        end do

    end subroutine init


    !-------------------------------------------------------------------------------------------------------------------------------


    ! sets 'b', 'g' or 'r' for a blue, green or red channel observation
    subroutine set_channel(this, status)
        type(pacsobsinfo), intent(inout) :: this
        integer, intent(out)             :: status
        integer                          :: length, unit, nsamples
        integer                          :: status_close
        character(len=2), allocatable    :: channels(:)

        this%channel = ' '

        ! old style file format
        length = len_trim(this%filename)
        if (this%filename(length-4:length) /= '.fits') then
            status = 0
            if (strlowcase(this%filename(length-3:length)) == 'blue') then
               this%channel = 'b'
            else if (strlowcase(this%filename(length-4:length)) == 'green') then
               this%channel = 'g'
            else if (strlowcase(this%filename(length-2:length)) == 'red') then
               this%channel = 'r'
            else
               status = 1
               write (ERROR_UNIT,'(a)') 'File name does not contain the array channel identifier (blue, green, red).'
            end if
            return
        endif

        call ft_open_bintable(trim(this%filename) // '[Status]', unit,         &
                              nsamples, status)
        if (status /= 0) return

        allocate(channels(nsamples))
        call ft_read_column(unit, 'BAND', 1, nsamples, channels, status)
        if (status /= 0) return

        if (any(channels /= channels(1))) then
            status = 1
            write (ERROR_UNIT,'(a)') 'Observation is BANDSWITCH.'
        else
           select case (channels(1))
               case ('BS')
                   this%channel = 'b'
               case ('BL')
                   this%channel = 'g'
               case ('R ')
                   this%channel = 'r'
               case default
                   status = 1
                   write (ERROR_UNIT,'(a)') 'Invalid array BAND value: ' // channels(1)
           end select
        end if

        call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine set_channel


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine set_compression_mode(this, status)
        type(pacsobsinfo), intent(inout) :: this
        integer, intent(out)             :: status
        integer           :: unit, ikey, length, status_close
        character(len=72) :: algorithm, keyword, comment

        this%compression_factor = 0

        ! old style file format
        length = len_trim(this%filename)
        if (this%filename(length-4:length) /= '.fits') then
            status = 0
            this%compression_factor = 1
            write (*,'(a)') 'Warning: Assuming compression factor of one for ob&
                &solete file format.'
            return
        end if

        call ft_open(trim(this%filename), unit, status)
        if (status /= 0) return

        ikey = 1
        do
            call ftgkys(unit, 'key.META_'//strinteger(ikey), keyword, comment, status)
            if (status /= 0) then
                write (ERROR_UNIT,'(a)') "ERROR: FITS keyword 'algorithm' is &
                    &not found."
                go to 999
            end if
            if (keyword == 'algorithm') exit
            ikey = ikey + 1
        end do

        call ftgkys(unit, 'META_'//strinteger(ikey), algorithm, comment, status)

        if (algorithm(1:19) == 'Floating Average  :') then
            read (algorithm(20:),'(i3)', iostat=status) this%compression_factor
            if (status /= 0) then
                write (ERROR_UNIT, '(a)') "ERROR: The compression algorithm '" &
                    // trim(algorithm) // "' is not understood."
                go to 999
            end if
        else if (algorithm == 'None') then
            this%compression_factor = 1
        end if
        
    999 call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine set_compression_mode


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine set_transparent_mode(this, status)
        type(pacsobsinfo), intent(inout) :: this
        integer, intent(out)             :: status
        integer                          :: unit, ikey, length, status_close
        character(len=72)                :: compression, keyword, comment

        this%transparent_mode = .false.

        ! old style file format
        length = len_trim(this%filename)
        if (this%filename(length-4:length) /= '.fits') then
            status = 0
            this%transparent_mode = .true.
            write (*,'(a)') 'Warning: Transparent mode assumed for obsolete file f&
                &ormat.'
            return
        end if

        call ft_open(trim(this%filename), unit, status)
        if (status /= 0) return

        ikey = 1
        do
            call ftgkys(unit, 'key.META_'//strinteger(ikey), keyword, comment, status)
            if (status /= 0) then
                status = 0
                write (ERROR_UNIT,'(a)') "Warning: FITS keyword 'compMode' is n&
                    &ot found. Assuming normal mode."
                go to 999
            end if
            if (keyword == 'compMode') exit
            ikey = ikey + 1
        end do

        call ftgkys(unit, 'META_'//strinteger(ikey), compression, comment, status)
        if (compression == 'Photometry Default Mode') then
            this%transparent_mode = .false.
        else if (compression == 'Photometry Lossless Compression Mode') then
            this%transparent_mode = .true.
        else
            status = 1
            write (ERROR_UNIT,'(a)') "ERROR: Unknown compression mode: '" // trim(compression) // "'."
            go to 999
        end if

    999 call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine set_transparent_mode


    !-------------------------------------------------------------------------------------------------------------------------------


    ! if the first sample to consider is not specified, we search starting from the end for the first sample that satisfies
    ! abs(chopfpuangle) > 0.01.
    ! we then discard the first 10, that might be affected by relaxation.
    subroutine set_valid_slice(this, first, last, status)

        type(pacsobsinfo), intent(inout) :: this
        integer*8, intent(in)            :: first, last
        integer, intent(out)             :: status

        integer*8              :: nsamples
        integer                :: length
        integer*8              :: isample
        integer*8, allocatable :: chop(:)

        length = len_trim(this%filename)
        if (this%filename(length-4:length) /= '.fits') then
            write(*,'(a)') 'Info: Obsolete file format.'

            call ft_read_extension(trim(this%filename) // '_ChopFpuAngle.fits', chop, status)
        else 
           ! XXX we shouldn't need to read that if first is not 0
            call ft_read_column(trim(this%filename)//'[Status]', 'CHOPFPUANGLE', chop, status)
        end if
        if (status /= 0) return

        nsamples = size(chop)

        if (nsamples == 0) then
            status = 1
            write (ERROR_UNIT) 'ERROR: Status extension is empty.'
            return
        end if

        if (last > nsamples) then
           status = 1
           write (ERROR_UNIT,'(a,2(i0,a))') "ERROR: The last sample '", last, "' exceeds the size of the observation '",       &
                nsamples, "'."
           return
        end if

        ! set last valid sample
        if (last == 0) then
           this%last = nsamples
        else
           this%last = last
        end if

        ! set first valid sample
        if (first == 0) then
            if (abs(chop(this%last)) > 0.01d0) then
                write (OUTPUT_UNIT,'(a)') 'Warning: last sample is invalid. Automatic search for valid slice is disabled.'
                this%first = 1
            else
                do isample = this%last-1, 1, -1
                    if (abs(chop(isample)) > 1.d0) then
                        if (this%first + 10 > this%last) then
                            write (OUTPUT_UNIT,'(a)') 'Warning: It is not possible to discard the first 10 scan samples. The observ&
                                  &ation is too short or the specified last sample is not large enough.'                            
                        else
                            this%first = this%first + 9
                        end if
                        exit
                    end if
                end do
                this%first = isample + 1
            end if
        else
            this%first = first
        end if

    end subroutine set_valid_slice


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine set_filename(this, filename, first, last, status)
        type(pacsobsinfo), intent(inout) :: this
        character(len=*), intent(in)     :: filename
        integer*8, intent(out)           :: first, last
        integer, intent(out)             :: status
        integer                          :: pos, length, delim
 
        status = 0
        first = 0
        last  = 0

        length = len_trim(filename)
        if (length > len(this%filename)) then
            status = 1
            write (ERROR_UNIT,'(a,i0,a,i0,a)') 'ERROR: Input filename length is&
                 & too long ', length, '. Maximum is ', len(this%filename), '.'
            return
        end if

        this%filename = filename

        if (filename(length:length) /= ']') return

        do pos = length-1, 2, -1
            if (filename(pos:pos) == '[') exit
        end do

        if (pos == 1) then
            status = 1
            write (ERROR_UNIT,'(a)') "ERROR: Missing opening bracket in '" //  &
                 trim(filename) // "'."
            return
        end if
      
        pos = pos + 1
        length = length - 1

        ! find delimiter ':'
        delim = index(filename(pos:length), ':')
        if (delim == 0) then
            delim = length + 1
        else 
            delim = pos + delim - 1
        end if

        ! read last sample
        if (delim < length) then
            read (filename(delim+1:length), '(i20)', iostat=status) last
            if (status /= 0) then
                write (ERROR_UNIT,'(a)') "ERROR: Invalid last sample: '" //    &
                      filename(pos-1:length+1) // "'."
                return
            end if
        end if

        ! read first sample
        if (delim > pos) then
            read (filename(pos:delim-1), '(i20)', iostat=status) first
            if (status /= 0) then
                write (ERROR_UNIT,'(a)') "ERROR: Invalid first sample: '" //    &
                      filename(pos-1:length+1) // "'."
                return
            end if
        end if

        if (last /= 0 .and. last < first) then
            status = 1
            write (ERROR_UNIT,'(a,2(i0,a))')"ERROR: Last sample '", last,  &
                  "' is less than the first sample '", first, "'."
            return
        end if

        if (first < 0) then
            status = 1
            write (ERROR_UNIT,'(a,i0,a)')"ERROR: The first sample '", first,&
                  "' is less than 1."
            return
        end if

        ! range specifications could be read, remove it
        this%filename(pos-1:) = ' '

    end subroutine set_filename


 end module module_pacsobservation
