program test_ngc6946_bpj

    use iso_fortran_env,        only : ERROR_UNIT, OUTPUT_UNIT
    use module_fitstools,       only : ft_read_keyword
    use module_math,            only : pInf, neq_real, sum_kahan
    use module_pacsinstrument,  only : NDIMS, NVERTICES, NEAREST_NEIGHBOUR, SHARP_EDGES, PacsInstrument
    use module_pacsobservation, only : PacsObservation, MaskPolicy
    use module_pointingmatrix,  only : PointingElement, pmatrix_direct, pmatrix_transpose
    use module_preprocessor,    only : subtract_meandim1, divide_vectordim2
    use module_projection,      only : surface_convex_polygon
    use module_string,          only : strinteger
    use module_tamasis,         only : p, POLICY_KEEP, POLICY_MASK, POLICY_REMOVE
    use module_wcs,             only : init_astrometry, ad2xy_gnomonic
    use omp_lib
    implicit none

    class(PacsInstrument), allocatable  :: pacs
    class(PacsObservation), allocatable :: obs
    type(MaskPolicy)                    :: policy
    character(len=*), parameter         :: dir = '/home/pchanial/data/pacs/transpScan/'
    character(len=*), parameter         :: filename(1) = dir // '1342184520_blue_PreparedFrames.fits[12001:16000]'
    integer, parameter                  :: npixels_per_sample = 6
    real(p), allocatable                :: signal(:,:)
    real(p), allocatable                :: coords(:,:), coords_yz(:,:)
    real(p)                             :: ra, dec, pa, chop, chop_old
    real(p), allocatable                :: surface1(:,:), surface2(:,:)
    logical*1, allocatable              :: mask(:,:)
    character(len=2880)                 :: header
    integer                             :: nx, ny
    integer                             :: status, count0, count1
    integer                             :: count2, count_rate, count_max
    integer                             :: idetector, isample, ipointing
    integer                             :: new_npixels_per_sample
    integer*8                           :: nsamples
    real(p), allocatable                :: map1d(:)
    type(pointingelement), allocatable  :: pmatrix(:,:,:)
    logical*1, allocatable              :: detector_mask(:,:)

    call system_clock(count0, count_rate, count_max)

    ! initialise observation
    allocate (obs)
    call obs%init(filename, policy, status, verbose=.true.)
    if (status == 104) then
       print *, 'Aborting test: file not found.'
       stop
    end if
    if (status /= 0) call failure('pacsobservation%init')
    nsamples = obs%slice(1)%nvalids

    call obs%print()

    ! initialise pacs instrument
    allocate (pacs)
    call pacs%read_detector_mask(obs%band, detector_mask, status,                                                      &
         transparent_mode=obs%slice(1)%observing_mode=='transparent')
    if (status /= 0) call failure('pacs%read_detector_mask')
    call pacs%init_with_calfiles(obs%band, detector_mask, 1, status)
    if (status /= 0) call failure('pacs%init_with_calfiles')

    ! get header map
    call pacs%compute_map_header(obs, .false., 3._p, header, status)
    if (status /= 0) call failure('compute_map_header.')

    call ft_read_keyword(header, 'naxis1', nx, status=status)
    if (status /= 0) call failure('compute_map_header 2.')
    call ft_read_keyword(header, 'naxis2', ny, status=status)
    if (status /= 0) call failure('compute_map_header 3.')

    ! allocate memory for the map
    allocate (map1d(0:nx*ny-1))

    ! read the signal file
    write (*,'(a)', advance='no') 'Reading tod file... '
    call system_clock(count1, count_rate, count_max)
    allocate (signal(nsamples, pacs%ndetectors))
    allocate (mask  (nsamples, pacs%ndetectors))
    call pacs%read(obs, [''], signal, mask, status)
    if (status /= 0) call failure('read_tod')
    call system_clock(count2, count_rate, count_max)
    write (*,'(f6.2,a)') real(count2-count1)/count_rate, 's'

    ! remove flat field
    write (*,'(a)') 'Flat-fielding... '
    call divide_vectordim2(signal, pacs%flatfield_detector)

    ! remove mean value in timeline
    write (*,'(a)') 'Removing mean value... '
    call subtract_meandim1(signal)

    ! compute the projector
    allocate (pmatrix(npixels_per_sample,nsamples,pacs%ndetectors))
    call pacs%compute_projection(SHARP_EDGES, obs, .false., header, nx, ny, pmatrix, new_npixels_per_sample, status)
    if (status /= 0) call failure('compute_projection_sharp_edges.')

    ! check flux conservation during backprojection
    write (*,'(a)', advance='no') 'Testing flux conservation...'
    call system_clock(count1, count_rate, count_max)
    allocate (surface1(nsamples, pacs%ndetectors))
    allocate (surface2(nsamples, pacs%ndetectors))
    allocate (coords(NDIMS, pacs%ndetectors*NVERTICES))
    allocate (coords_yz(NDIMS, pacs%ndetectors*NVERTICES))

    call init_astrometry(header, status=status)
    if (status /= 0) call failure('init_astrometry')

    chop_old = pInf

    write (*,*) 'surfaces:'
    !XXX IFORT bug
    !!$omp parallel do default(shared) firstprivate(chop_old) private(ipointing, isample, ra, dec, pa, chop, coords, coords_yz)
    isample = 1
    do ipointing = 1, obs%slice(1)%nsamples
        if (obs%slice(1)%p(ipointing)%removed) cycle
        call obs%get_position_index(1, ipointing, 1, ra, dec, pa, chop)
        if (abs(chop-chop_old) > 1.e-2_p) then
            coords_yz = pacs%uv2yz(pacs%detector_corner, pacs%distortion_yz, chop)
            chop_old = chop
        end if
        coords = pacs%yz2ad(coords_yz, ra, dec, pa)
        coords = ad2xy_gnomonic(coords)
        do idetector = 1, pacs%ndetectors
            surface1(isample,idetector) = abs(surface_convex_polygon(                                                              &
                                          real(coords(:,(idetector-1)*NVERTICES+1:idetector*NVERTICES),p)))
            surface2(isample,idetector) = sum_kahan(real(pmatrix(:,isample,idetector)%weight, p))
        end do
        isample = isample + 1
    end do
    !!$omp end parallel do

    call system_clock(count2, count_rate, count_max)
    write (*,'(f6.2,a)') real(count2-count1)/count_rate, 's'
    write (*,*) 'Difference: ', maxval(abs((surface1-surface2)/surface1))
    write (*,*) 'Sum surfaces: ', sum_kahan(surface1), sum_kahan(surface2), nx, ny
    write (*,*) 'Sum pmatrix weight: ', sum(pmatrix%weight), 'max:', maxval(pmatrix%weight)
    write (*,*) 'Sum pmatrix pixel: ', sum(int(pmatrix%pixel,kind=4)), 'max:', maxval(pmatrix%pixel)
    write (*,*) 'Sum signal: ', sum_kahan(signal), 'max:', maxval(signal)

    ! back project the timeline
    write (*,'(a)', advance='no') 'Computing the back projection... '
    signal = 1._p / surface1
    call system_clock(count1, count_rate, count_max)
    call pmatrix_transpose(pmatrix, signal, map1d)
    call system_clock(count2, count_rate, count_max)
    write (*,'(f6.2,a)') real(count2-count1)/count_rate, 's'

    ! test the back projected map
    write (OUTPUT_UNIT,*) 'Sum in map is ', sum_kahan(map1d), ' ...instead of ', strinteger(int(pacs%ndetectors*nsamples, kind=4))
    if (neq_real(sum_kahan(map1d), real(pacs%ndetectors*nsamples,p), 10._p * epsilon(1.0))) then
        stop 'FAILED:back projection'
    end if

    ! project the map
    write (*,'(a)', advance='no') 'Computing the forward projection... '
    map1d = 1._p
    call system_clock(count1, count_rate, count_max)
    call pmatrix_direct(pmatrix, map1d, signal)
    call system_clock(count2, count_rate, count_max)
    write (*,'(f6.2,a)') real(count2-count1)/count_rate, 's'

    ! test the back projected map
    if (any(neq_real(signal, surface1, 10._p * epsilon(1.0)))) then
        write (OUTPUT_UNIT,*), 'max N*epsilon:', maxval((signal-surface1)/signal/epsilon(1.0), signal==signal                      &
              .and. surface1 == surface1), ' instead of 6.07.'
        write (OUTPUT_UNIT,*), 'NaN:', count(signal/=signal), count(surface1/=surface1)
        call failure('wrong back projection.')
    end if

    call system_clock(count2, count_rate, count_max)
    write (*,'(a,f6.2,a)') 'Total elapsed time: ', real(count2-count0)/count_rate, 's'

contains

    subroutine failure(errmsg)
        character(len=*), intent(in) :: errmsg
        write (ERROR_UNIT,'(a)'), 'FAILED: ' // errmsg
        stop 1
    end subroutine failure

end program test_ngc6946_bpj
