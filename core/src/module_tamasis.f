module module_tamasis

    use iso_fortran_env,  only : ERROR_UNIT
    use module_precision, only : sp, dp, qp
    implicit none
    private

    public :: p, POLICY_KEEP, POLICY_MASK, POLICY_REMOVE
    public :: tamasis_dir  ! Tamasis data directory
    public :: tamasis_version

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


end module module_tamasis
