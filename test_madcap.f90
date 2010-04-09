program test_madcap

    use module_filtering,      only : filterset
    use module_madcap
    use module_math,           only : neq_real
    use module_pointingmatrix, only : pointingelement
    use module_precision,      only : p
    implicit none

    character(len=*), parameter :: invnttfile = 'tests/invntt_'
    character(len=*), parameter :: todfile = '/home/pchanial/work/test_madmap/rot=00/tod'
    type(filterset)             :: filter_le, filter_be
    integer                     :: status, nsamples, ndetectors, npixels_per_sample
    integer*8                   :: nsamples_tod
    real(p), allocatable        :: tod(:,:)
    type(pointingelement), allocatable :: pmatrix(:,:,:)

    call read_filter(invnttfile // 'le', 'little_endian', 2, filter_le, status)
    if (status /= 0) stop 'FAILED: read_filter little_endian'

    if (filter_le%nslices /= 3 .or. filter_le%ndetectors /= 2 .or. filter_le%ncorrelations /= 100) stop 'FAILED: read_filter'
    if (neq_real(filter_le%data(1,1,1), 5597147.4155586753_p, 15)) stop 'FAILED: read_filter'

    call read_filter(invnttfile // 'be', 'big_endian', 2, filter_be, status)
    if (status /= 0) stop 'FAILED: read_filter big_endian'

    if (filter_be%nslices /= 3 .or. filter_be%ndetectors /= 2 .or. filter_be%ncorrelations /= 100) stop 'FAILED: read_filter'
    if (any(neq_real(filter_le%data, filter_be%data, 15))) stop 'FAILED: read_filter data'

    nsamples = sum(filter_le%last - filter_le%first + 1)
    ndetectors = filter_le%ndetectors

    call read_tod_header(todfile, 'little_endian', nsamples_tod, npixels_per_sample, status)
    if (status /= 0 .or. npixels_per_sample /= 1) stop 'FAILED: read_tod_header'

    allocate (tod(nsamples,ndetectors))
    allocate (pmatrix(npixels_per_sample,nsamples,ndetectors))

    call read_tod(todfile, 'little_endian', filter_le%first, filter_le%last, tod, pmatrix, status)
    if (status /= 0) stop 'FAILED: read_tod'

    deallocate (tod)
    deallocate (pmatrix)
    stop 'OK.'

end program test_madcap
