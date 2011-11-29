program test_ngc6946_bpj

    use iso_fortran_env,        only : ERROR_UNIT, OUTPUT_UNIT
    use module_fitstools,       only : ft_read_keyword
    use module_math,            only : pInf, neq_real, sum_kahan
    use module_observation,     only : MaskPolicy, Observation
    use module_pacsinstrument,  only : NDIMS, NVERTICES, NEAREST_NEIGHBOUR, SHARP_EDGES, PacsInstrument
    use module_pacsobservation, only : PacsObservation
    use module_pointingmatrix,  only : PointingElement, pmatrix_direct, pmatrix_transpose
    use module_preprocessor,    only : subtract_meandim1, divide_vectordim2
    use module_projection,      only : surface_convex_polygon
    use module_string,          only : strinteger
    use module_tamasis,         only : p, POLICY_KEEP, POLICY_MASK, POLICY_REMOVE
    use module_wcs,             only : init_astrometry, ad2xy_gnomonic
    use omp_lib
    implicit none

    class(Observation), allocatable     :: obs
    class(PacsInstrument), allocatable  :: pacs
    class(PacsObservation), allocatable :: pacsobs
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
    integer                             :: status
    integer                             :: idetector, isample, ipointing
    integer                             :: new_npixels_per_sample
    real(p), allocatable                :: map1d(:)
    type(pointingelement), allocatable  :: pmatrix(:,:,:)
    logical*1                           :: detector_mask(32,64)

    ! initialise observation
    allocate (pacsobs)
    call pacsobs%init(filename, policy, status, verbose=.true.)
    if (status == 104) then
       print *, 'Aborting test: file not found.'
       stop
    end if
    if (status /= 0) call failure('pacsobservation%init')

    call pacsobs2obs(pacsobs, 1, obs, status)
    if (status /= 0) call failure('pacsobs2obs')
    obs%pointing%chop = 0

    ! initialise pacs instrument
    allocate (pacs)
    detector_mask = .true.
    detector_mask(1:16,17:32) = .false.
    call pacs%init_with_calfiles(pacsobs%band, detector_mask, 1, status)
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
    allocate (signal(obs%nvalids, pacs%ndetectors))
    allocate (mask  (obs%nvalids, pacs%ndetectors))
    call pacs%read(pacsobs%slice(1)%filename, obs%pointing%masked, obs%pointing%removed, [''], signal, mask, status)
    if (status /= 0) call failure('read_tod')

    ! remove flat field
    call divide_vectordim2(signal, pacs%flatfield_detector)
    where (mask) signal = 0

    ! remove mean value in timeline
    call subtract_meandim1(signal)

    ! compute the projector
    allocate (pmatrix(npixels_per_sample,obs%nvalids,pacs%ndetectors))
    call pacs%compute_projection(SHARP_EDGES, obs, .false., header, nx, ny, pmatrix, new_npixels_per_sample, status)
    if (status /= 0) call failure('compute_projection_sharp_edges.')

    ! check flux conservation during backprojection
    allocate (surface1(obs%nvalids, pacs%ndetectors))
    allocate (surface2(obs%nvalids, pacs%ndetectors))
    allocate (coords(NDIMS, pacs%ndetectors*NVERTICES))
    allocate (coords_yz(NDIMS, pacs%ndetectors*NVERTICES))

    call init_astrometry(header, status=status)
    if (status /= 0) call failure('init_astrometry')

    chop_old = pInf

    !XXX IFORT bug
    !!$omp parallel do default(shared) firstprivate(chop_old) private(ipointing, isample, ra, dec, pa, chop, coords, coords_yz)
    isample = 1
    do ipointing = 1, obs%nsamples
        if (obs%pointing(ipointing)%removed) cycle
        call obs%get_position_index(ipointing, 1, obs%offset, ra, dec, pa, chop)
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

    if (neq_real(sum(real(pmatrix%weight, kind=p)), 840698.47972_p, 10._p * epsilon(1.0))) call failure('sum weight')
    if (sum(int(pmatrix%pixel,kind=4)) /= 1974986011) call failure('sum pixel indices')

    ! back project the timeline
    signal = 1._p / surface1
    map1d = 0
    call pmatrix_transpose(pmatrix, signal, map1d)

    ! test the back projected map
    if (neq_real(sum_kahan(map1d), real(pacs%ndetectors*obs%nvalids,p), 10._p * epsilon(1.0))) then
        stop 'FAILED:back projection'
    end if

    ! project the map
    map1d = 1._p
    call pmatrix_direct(pmatrix, map1d, signal)

    ! test the back projected map
    if (any(neq_real(signal, surface1, 10._p * epsilon(1.0)))) then
        write (OUTPUT_UNIT,*), 'max N*epsilon:', maxval((signal-surface1)/signal/epsilon(1.0), signal==signal                      &
              .and. surface1 == surface1), ' instead of 6.07.'
        write (OUTPUT_UNIT,*), 'NaN:', count(signal/=signal), count(surface1/=surface1)
        call failure('wrong back projection.')
    end if

contains

    subroutine failure(errmsg)
        character(len=*), intent(in) :: errmsg
        write (ERROR_UNIT,'(a)'), 'FAILED: ' // errmsg
        stop 1
    end subroutine failure

    subroutine pacsobs2obs(pacsobs, islice, obs, status)
        class(PacsObservation), intent(in)           :: pacsobs
        integer, intent(in)                          :: islice
        class(Observation), allocatable, intent(out) :: obs
        integer, intent(out)                         :: status
        
        allocate (obs)
        call obs%init(pacsobs%slice(islice)%p%time, pacsobs%slice(islice)%p%ra, pacsobs%slice(islice)%p%dec,                       &
                      pacsobs%slice(islice)%p%pa, pacsobs%slice(islice)%p%chop, pacsobs%slice(islice)%p%masked,                    &
                      pacsobs%slice(islice)%p%removed, pacsobs%slice(islice)%compression_factor,                                   &
                      (pacsobs%slice(islice)%compression_factor - 1) / (2._p * pacsobs%slice(islice)%compression_factor), status)

    end subroutine pacsobs2obs

end program test_ngc6946_bpj
