program test_madcap

    use module_filtering,      only : filterset
    use module_fitstools,      only : ft_read_extension
    use module_madcap
    use module_math,           only : NaN, neq_real
    use module_pointingmatrix, only : backprojection_weighted, pointingelement
    use module_precision,      only : p
    implicit none

    character(len=*), parameter :: invnttfile1 = 'tests/madmap1/invntt_'
    character(len=*), parameter :: todfile = 'tests/madmap1/todSpirePsw_be'
    character(len=*), parameter :: invnttfile2 = 'tests/madmap1/invnttSpirePsw_be'
    type(filterset)             :: filter_le, filter_be, filter
    integer                     :: status, nsamples, ndetectors, npixels_per_sample, nx, ny
    real(p), allocatable        :: tod(:,:)
    type(pointingelement), allocatable :: pmatrix(:,:,:)
    real(p), allocatable        :: coverage(:,:), map1d(:), map(:,:), map_ref(:,:)

    call read_filter(invnttfile1 // 'le', 'little_endian', 2, filter_le, status)
    if (status /= 0) stop 'FAILED: read_filter little_endian'

    if (filter_le%nslices /= 3 .or. filter_le%ndetectors /= 2 .or. filter_le%ncorrelations /= 100) stop 'FAILED: read_filter'
    if (neq_real(filter_le%data(1,1,1), 5597147.4155586753_p, 15)) stop 'FAILED: read_filter'

    call read_filter(invnttfile1 // 'be', 'big_endian', 2, filter_be, status)
    if (status /= 0) stop 'FAILED: read_filter big_endian'

    if (filter_be%nslices /= 3 .or. filter_be%ndetectors /= 2 .or. filter_be%ncorrelations /= 100) stop 'FAILED: read_filter'
    if (any(neq_real(filter_le%data, filter_be%data, 15))) stop 'FAILED: read_filter data'

    ndetectors = 135
    call read_filter(invnttfile2, 'big_endian', ndetectors, filter, status)
    if (status /= 0) stop 'FAILED: read_filter spire'

    nsamples = sum(filter%last - filter%first + 1)
    npixels_per_sample = 1

    allocate (tod(nsamples,ndetectors))
    allocate (pmatrix(npixels_per_sample,nsamples,ndetectors))
    call read_tod(todfile, 'big_endian', filter%first, filter%last, tod, pmatrix, status)
    if (status /= 0) stop 'FAILED: read_tod spire'

    call ft_read_extension('tests/madmap1/naivemapSpirePsw.fits[image]', map_ref, status)
    if (status /= 0) stop 'FAILED: ft_read_extension image'
    nx = size(map_ref,1)
    ny = size(map_ref,2)

    call ft_read_extension('tests/madmap1/madmapSpirePsw.fits[coverage]', coverage, status)
    if (status /= 0) stop 'FAILED: ft_read_extension coverage'

#ifdef GFORTRAN
    where (coverage < 1e-15_p)
         coverage = NaN
         map_ref  = NaN
    end where
#endif

    allocate (map(nx, ny))
    allocate (map1d(0:maxval(pmatrix%pixel)))

    call backprojection_weighted(pmatrix, tod, map=map1d)
    
    map = unpack(map1d, coverage == coverage, NaN)

    if (any(neq_real(map, map_ref, 15))) stop 'FAILED: map_ref'

    deallocate (tod, pmatrix, map_ref, coverage, map, map1d)

    stop 'OK.'

end program test_madcap
