program test_madcap

    use module_filtering,      only : FilterUncorrelated
    use module_fitstools,      only : ft_read_image
    use module_madcap
    use module_math,           only : NaN, neq_real
    use module_pointingmatrix, only : backprojection_weighted, pointingelement
    use module_tamasis,        only : p
    implicit none

    character(len=*), parameter :: invnttfile1 = 'madcap/test/data/madmap1/invntt_'
    character(len=*), parameter :: todfile = 'madcap/test/data/madmap1/todSpirePsw_be'
    character(len=*), parameter :: invnttfile2 = 'madcap/test/data/madmap1/invnttSpirePsw_be'
    type(FilterUncorrelated), allocatable :: filter_le(:), filter_be(:), filter(:)
    type(pointingelement), allocatable    :: pmatrix(:,:,:)
    integer                     :: status, ndetectors, npixels_per_sample, nx, ny
    real(p), allocatable        :: tod(:,:)
    real(p), allocatable        :: coverage(:,:), map1d(:), weight1d(:), map(:,:), map_ref(:,:)
    integer, allocatable        :: nsamples(:)

    call read_filter(invnttfile1 // 'le', 'little_endian', 2, filter_le, nsamples, status)
    if (status /= 0) stop 'FAILED: read_filter little_endian'

    if (size(filter_le) /= 3 .or. any(filter_le%ndetectors /= 2) .or. any(filter_le%ncorrelations /= 100)) stop 'FAILED:read_filter'
    if (neq_real(filter_le(1)%data(1,1), 5597147.4155586753_p)) stop 'FAILED: read_filter'

    call read_filter(invnttfile1 // 'be', 'big_endian', 2, filter_be, nsamples, status)
    if (status /= 0) stop 'FAILED: read_filter big_endian'

    if (size(filter_be) /= 3 .or. any(filter_be%ndetectors /= 2) .or. any(filter_be%ncorrelations /= 100)) stop 'FAILED:read_filter'
    if (any(neq_real(filter_le(1)%data, filter_be(1)%data))) stop 'FAILED: read_filter data'
    if (any(neq_real(filter_le(2)%data, filter_be(2)%data))) stop 'FAILED: read_filter data'
    if (any(neq_real(filter_le(3)%data, filter_be(3)%data))) stop 'FAILED: read_filter data'

    ndetectors = 135
    call read_filter(invnttfile2, 'big_endian', ndetectors, filter, nsamples, status)
    if (status /= 0) stop 'FAILED: read_filter spire'

    npixels_per_sample = 1

    allocate (tod(sum(nsamples),ndetectors))
    allocate (pmatrix(npixels_per_sample,sum(nsamples),ndetectors))
    call read_tod(todfile, 'big_endian', nsamples, tod, pmatrix, status)
    if (status /= 0) stop 'FAILED: read_tod spire'

    call ft_read_image('madcap/test/data/madmap1/naivemapSpirePsw.fits[image]', map_ref, status)
    if (status /= 0) stop 'FAILED: ft_read_image image'
    nx = size(map_ref,1)
    ny = size(map_ref,2)

    call ft_read_image('madcap/test/data/madmap1/madmapSpirePsw.fits[coverage]', coverage, status)
    if (status /= 0) stop 'FAILED: ft_read_image coverage'

#ifdef GFORTRAN
    where (coverage < 1e-15_p .or. coverage > 5e+276_p)
         coverage = NaN
         map_ref  = NaN
    end where
#endif

    allocate (map(nx, ny))
    allocate (map1d(0:maxval(pmatrix%pixel)), weight1d(0:maxval(pmatrix%pixel)))

    call backprojection_weighted(pmatrix, tod, map=map1d, weight=weight1d)
    
    map = unpack(map1d, coverage == coverage, NaN)

    if (any(neq_real(map, map_ref))) stop 'FAILED: map_ref'

    deallocate (tod, pmatrix, map_ref, coverage, map, map1d, weight1d)

    stop 'OK.'

end program test_madcap
