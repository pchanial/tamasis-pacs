module module_tamasis

    use iso_fortran_env,  only : OUTPUT_UNIT
    use module_precision, only : sp, dp, qp
    implicit none
    private

    public :: p, POLICY_KEEP, POLICY_MASK, POLICY_REMOVE
    public :: tamasis_dir  ! Tamasis data directory
    public :: tamasis_version
    public :: info_time

    character(len=*), parameter :: tamasis_dir =                                                                                   &
TAMASIS_DIR
    character(len=*), parameter :: tamasis_version = TAMASIS_VERSION

#if PRECISION_REAL == 4
    integer, parameter :: p = sp
#elif PRECISION_REAL == 8
    integer, parameter :: p = dp
#elif PRECISION_REAL == 16
    integer, parameter :: p = qp
#else
#error Invalid or undefined macro PRECISION_REAL
#endif

    integer, parameter :: POLICY_KEEP   = 0
    integer, parameter :: POLICY_MASK   = 1
    integer, parameter :: POLICY_REMOVE = 2

contains


    subroutine info_time(message, count1, extra)

        character(len=*), intent(in)           :: message
        integer, intent(in)                    :: count1
        character(len=*), intent(in), optional :: extra

        character(len=100) :: hostname
        integer            :: count2, count_rate, status

#ifdef GFORTRAN
        status = hostnm(hostname)
#else
        call get_environment_variable(name="HOSTNAME", value=hostname)
#endif
        call system_clock(count2, count_rate)
        write (OUTPUT_UNIT, '(5a,f9.2,a,$)') 'Info ', trim(hostname), ': ', message, '... ', real(count2-count1)/count_rate, 's'
        if (present(extra)) then
            write (OUTPUT_UNIT, '(x,a)') extra
        else
            write (OUTPUT_UNIT,*)
        end if

    end subroutine info_time


end module module_tamasis
