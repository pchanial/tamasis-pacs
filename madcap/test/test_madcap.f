program test_madcap

    use iso_fortran_env,  only : ERROR_UNIT
    use module_filtering, only : FilterUncorrelated
    use module_madcap
    use module_math,      only : neq_real
    use module_tamasis,   only : p
    implicit none

    character(len=*), parameter :: invnttfile1 = 'madcap/test/data/madmap1/invntt_'
    type(FilterUncorrelated), allocatable :: filter_le(:), filter_be(:)
    integer                     :: status
    integer*8, allocatable      :: nsamples(:)

    call read_filter(invnttfile1 // 'le', 'little_endian', 2, filter_le, nsamples, status)
    if (status /= 0) call failure('read_filter little_endian')

    if (size(filter_le) /= 3 .or. any(filter_le%ndetectors /= 2) .or.any(filter_le%ncorrelations/=100)) call failure('read_filter1')
    if (neq_real(filter_le(1)%data(1,1), 5597147.4155586753_p)) call failure('read_filter2')

    call read_filter(invnttfile1 // 'be', 'big_endian', 2, filter_be, nsamples, status)
    if (status /= 0) call failure('read_filter big_endian')

contains

    subroutine failure(errmsg)
        character(len=*), intent(in) :: errmsg
        write (ERROR_UNIT,'(a)'), 'FAILED: ' // errmsg
        stop 1
    end subroutine failure

end program test_madcap
