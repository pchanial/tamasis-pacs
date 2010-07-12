module module_tamasis

    use iso_fortran_env, only : ERROR_UNIT
    implicit none
    private

    public :: init_tamasis
    public :: get_tamasis_path
    public :: tamasis_path_len
    public :: POLICY_KEEP, POLICY_MASK, POLICY_REMOVE


    integer, parameter :: POLICY_KEEP = 0
    integer, parameter :: POLICY_MASK = 1
    integer, parameter :: POLICY_REMOVE = 2

    character(len=255) :: tamasis_path     ! Tamasis installation directory
    integer            :: tamasis_path_len ! length of the path to the tamasis installation directory

    interface init_tamasis
        module procedure  set_tamasis_path, set_tamasis_path_default
    end interface init_tamasis


contains


    subroutine set_tamasis_path(path)

        character(len=*), intent(in) :: path
        integer                      :: length
        
        length = len(path)

        if (length == 0) then
            tamasis_path = './'
            tamasis_path_len = 2
            return
        end if

        if (path(length:length) == '/') length = length - 1

        if (length >= len(tamasis_path)) then
            write (ERROR_UNIT,'(a)') 'The installation directory for TAMASIS is too long (>255 characters):'
            write (ERROR_UNIT,'(a)') path
            stop
        end if
        
        tamasis_path = path(1:length) // '/'
        tamasis_path_len = length+1
        
    end subroutine set_tamasis_path


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine set_tamasis_path_default()

        character(len=255) :: command
        integer            :: length, pos, status

        call get_command(command, length, status)
        if (status == 1) stop 'SET_TAMASIS_PATH_DEFAULT: GET_COMMAND failure.'
        
        if (status == -1) then
            write (ERROR_UNIT,'(a)') 'The command path is too long. It cannot be used to infer the TAMASIS installation directory.'
            stop
        end if

        pos = index(command, '/', .true.)
        if (pos == 0) then
            tamasis_path = ''
            tamasis_path_len = 0
        else
            tamasis_path = command(1:pos)
            tamasis_path_len = pos
        end if

    end subroutine


    !-------------------------------------------------------------------------------------------------------------------------------


    pure function get_tamasis_path() result (path)

        character(len=tamasis_path_len) :: path

        path = tamasis_path

    end function get_tamasis_path


end module module_tamasis
