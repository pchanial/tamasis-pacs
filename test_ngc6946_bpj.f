program test_ngc6946_bpj

    use iso_fortran_env,        only : ERROR_UNIT, OUTPUT_UNIT
    use module_fitstools,       only : ft_read_parameter, ft_write
    use module_math,            only : pInf, neq_real, sum_kahan
    use module_pacsinstrument,  only : ndims, nvertices, pacsinstrument
    use module_pacsobservation, only : pacsobservation
    use module_pacspointing,    only : pacspointing
    use module_pointingmatrix,  only : pointingelement, pmatrix_direct, pmatrix_transpose
    use module_preprocessor,    only : subtract_meandim1, divide_vectordim2
    use module_projection,      only : surface_convex_polygon
    use module_string,          only : strinteger
    use module_wcs,             only : init_astrometry, ad2xy_gnomonic
    use omp_lib
    implicit none

    class(pacsinstrument), allocatable  :: pacs
    class(pacsobservation), allocatable :: obs
    class(pacspointing), allocatable    :: pointing
    character(len=*), parameter         :: inputdir    = '/home/pchanial/work/pacs/data/transparent/NGC6946/'
    character(len=*), parameter         :: filename(1) = inputdir // '1342184520_blue[12001:86000]'
    integer, parameter                  :: npixels_per_sample = 6
    real*8, allocatable                 :: signal(:,:), coords(:,:), coords_yz(:,:)
    real*8                              :: ra, dec, pa, chop, chop_old
    real*8, allocatable                 :: surface1(:,:), surface2(:,:)
    logical*1, allocatable              :: mask(:,:)
    character(len=80)                   :: outfile
    character(len=2880)                 :: header
    integer                             :: nx, ny
    integer                             :: status, count, count0, count1
    integer                             :: count2, count_rate, count_max
    integer                             :: idetector, isample
    integer*8                           :: nsamples
    real*8, allocatable                 :: map1d(:)
    type(pointingelement), allocatable  :: pmatrix(:,:,:)

    call system_clock(count0, count_rate, count_max)

    ! initialise observation
    allocate(obs)
    call obs%init(filename, status)
    if (status /= 0) stop 'FAILED: pacsobservation%init'
    nsamples = obs%info(1)%nsamples

    ! initialise pacs instrument
    allocate(pacs)
    call pacs%init(obs%channel, obs%transparent_mode, 1, .false., status)
    if (status /= 0) stop 'FAILED: pacsinstrument%init'

    ! initialise pointing information
    allocate(pointing)
    call pointing%init(obs, status)
    if (status /= 0) stop 'FAILED: pacspointing%init'

    ! get header map
    call pacs%compute_map_header(pointing, .false., 3.d0, header, status)
    if (status /= 0) stop 'FAILED: compute_map_header.'

    call ft_read_parameter(header, 'naxis1', nx, count, status=status)
    if (status /= 0 .or. count == 0) stop 'FAILED: compute_map_header 2.'
    call ft_read_parameter(header, 'naxis2', ny, count, status=status)
    if (status /= 0 .or. count == 0) stop 'FAILED: compute_map_header 3.'

    ! allocate memory for the map
    allocate(map1d(0:nx*ny-1))

    ! read the signal file
    write(*,'(a)', advance='no') 'Reading tod file... '
    call system_clock(count1, count_rate, count_max)
    allocate(signal(nsamples, pacs%ndetectors))
    allocate(mask  (nsamples, pacs%ndetectors))
    call pacs%read(obs, signal, mask, status)
    if (status /= 0) stop 'FAILED: read_tod'
    call system_clock(count2, count_rate, count_max)
    write(*,'(f6.2,a)') real(count2-count1)/count_rate, 's'

    ! remove flat field
    write(*,'(a)') 'Flat-fielding... '
    call divide_vectordim2(signal, pacs%flatfield)

    ! remove mean value in timeline
    write(*,'(a)') 'Removing mean value... '
    call subtract_meandim1(signal)

    ! compute the projector
    write(*,'(a)', advance='no') 'Computing the projector... '
    allocate(pmatrix(npixels_per_sample,nsamples,pacs%ndetectors))
    call system_clock(count1, count_rate, count_max)
    call pacs%compute_projection_sharp_edges(pointing, .false., header, nx, ny, pmatrix, status)
    if (status /= 0) stop 'FAILED: compute_projection_sharp_edges.'
    call system_clock(count2, count_rate, count_max)
    write(*,'(f6.2,a)') real(count2-count1)/count_rate, 's'

    ! check flux conservation during backprojection
    write(*,'(a)', advance='no') 'Testing flux conservation...'
    call system_clock(count1, count_rate, count_max)
    allocate(surface1(nsamples, pacs%ndetectors))
    allocate(surface2(nsamples, pacs%ndetectors))
    allocate(coords(ndims, pacs%ndetectors*nvertices))
    allocate(coords_yz(ndims, pacs%ndetectors*nvertices))

    call init_astrometry(header, status=status)
    if (status /= 0) stop 'FAILED: init_astrometry'

    chop_old = pInf

    write(*,*) 'surfaces:'
    !XXX IFORT bug
    !!$omp parallel do default(shared) firstprivate(chop_old) private(isample, ra, dec, pa, chop, coords, coords_yz)
    do isample = 1, nsamples
        call pointing%get_position_index(1, isample, 1, ra, dec, pa, chop)
        if (abs(chop-chop_old) > 1.d-2) then
            coords_yz = pacs%uv2yz(pacs%corners_uv, pacs%distortion_yz, chop)
            chop_old = chop
        end if
        coords = pacs%yz2ad(coords_yz, ra, dec, pa)
        coords = ad2xy_gnomonic(coords)
        do idetector = 1, pacs%ndetectors
            surface1(isample,idetector) = abs(surface_convex_polygon(coords(:,(idetector-1)*nvertices+1:idetector*nvertices)))
            surface2(isample,idetector) = sum_kahan(real(pmatrix(:,isample,idetector)%weight, kind=8))
        end do
    end do
    !!$omp end parallel do

    call system_clock(count2, count_rate, count_max)
    write(*,'(f6.2,a)') real(count2-count1)/count_rate, 's'
    write(*,*) 'Difference: ', maxval(abs((surface1-surface2)/surface1))
    write(*,*) 'Sum surfaces: ', sum_kahan(surface1), sum_kahan(surface2), nx, ny
    write(*,*) 'Sum pmatrix weight: ', sum(pmatrix%weight), 'max:', maxval(pmatrix%weight)
    write(*,*) 'Sum pmatrix pixel: ', sum(int(pmatrix%pixel,kind=8)), 'max:', maxval(pmatrix%pixel)
    write(*,*) 'Sum signal: ', sum_kahan(signal), 'max:', maxval(signal)
    !XXX

    ! back project the timeline
    write(*,'(a)', advance='no') 'Computing the back projection... '
    signal = 1.d0 / surface1
    call system_clock(count1, count_rate, count_max)
    call pmatrix_transpose(pmatrix, signal, map1d)
    call system_clock(count2, count_rate, count_max)
    write(*,'(f6.2,a)') real(count2-count1)/count_rate, 's'

    ! write the map as fits file
    outfile = '/tmp/ngc6946_bpj_' // strinteger(omp_get_max_threads()) //'.fits'
    write(*,'(a)') 'Writing FITS file... ' // trim(outfile)
    call ft_write(trim(outfile), reshape(map1d, [nx,ny]), header, status)
    if (status /= 0) stop 'FAILED: ft_write.'

    ! test the back projected map
    write (OUTPUT_UNIT,*) 'Sum in map is ', sum_kahan(map1d), ' ...instead of ',&
                         strinteger(int(pacs%ndetectors*nsamples, kind=4))
    if (neq_real(sum_kahan(map1d), real(pacs%ndetectors*nsamples,kind=8), 6)) then
        stop 'FAILED.'
    end if

    ! project the map
    write(*,'(a)', advance='no') 'Computing the forward projection... '
    map1d = 1.d0
    call system_clock(count1, count_rate, count_max)
    call pmatrix_direct(pmatrix, map1d, signal)
    call system_clock(count2, count_rate, count_max)
    write(*,'(f6.2,a)') real(count2-count1)/count_rate, 's'

    ! test the back projected map
    if (any(neq_real(signal, surface1, 5))) then
        write (ERROR_UNIT,*) 'Invalid signal.'
        stop 'FAILED.'
    end if

    call system_clock(count2, count_rate, count_max)
    write(*,'(a,f6.2,a)') 'Total elapsed time: ', real(count2-count0)/count_rate, 's'

    flush(OUTPUT_UNIT)
    stop "OK."

end program test_ngc6946_bpj
