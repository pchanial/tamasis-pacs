program test_ngc6946_bpj

    use, intrinsic :: ISO_FORTRAN_ENV
    use module_fitstools
    use module_pacsinstrument
    use module_pacspointing
    use module_preprocessor
    use module_projection
    use module_wcs, only : init_wcslib, ad2xy_wcslib, free_wcslib
    use module_wcslib, only : WCSLEN, wcsfree
    use precision
    implicit none

    class(pacsinstrument), allocatable :: pacs
    class(pacspointing), allocatable   :: pointing
    character(len=*), parameter :: inputdir        = '/home/pchanial/work/pacs/data/transparent/NGC6946/'
    character(len=*), parameter :: outputdir       = '/home/pchanial/'
    character(len=*), parameter :: filename        = inputdir // '1342184520_blue'
    character(len=*), parameter :: filename_signal = inputdir // '1342184520_blue_Signal.fits'
    character(len=*), parameter :: filename_mask   = inputdir // '1342184520_blue_Mask.fits'
    character(len=*), parameter :: filename_time   = inputdir // '1342184520_blue_Time.fits'

    real*8, allocatable                :: signal(:,:), coords(:,:)
    real*8                             :: ra, dec, pa, chop
    real*4, allocatable                :: surface1(:,:), surface2(:,:)
    logical*1, allocatable             :: mask(:,:)
    real*8, allocatable                :: time(:)
    integer*8, allocatable             :: timeus(:)
    character(len=2880)                :: header
    integer                            :: wcs(WCSLEN), nx, ny, index
    integer                            :: status, count0, count1, count2, count_rate, count_max, idetector, isample
    integer*8                          :: first, last, nsamples
    real*8, allocatable                :: map1d(:)
    type(pointingelement), allocatable :: pmatrix(:,:,:)

    call system_clock(count0, count_rate, count_max)

    ! read pointing information
    allocate(pointing)
    call pointing%load(filename, status)
    if (status /= 0) stop 'pointing%load_filename: FAILED.'

    ! read the time file
    status = 0
    first = 12001
    first = 12001
    !last  = 86936
    last  = 86000
    nsamples = last - first + 1
    allocate(time(last-first+1))
    allocate(timeus(last-first+1))
    call ft_readslice(filename_time // '+1', first, last, timeus, status)
    if (status /= 0) stop 'FAILED: ft_readslice.'
    time = timeus * 1.0d-6

    ! get the pacs instance, read the calibration files
    allocate(pacs)
    call pacs%read_calibration_files(status)
    if (status /= 0) stop 'FAILED: read_calibration_files.'

    call pacs%filter_detectors('blue', transparent_mode=.true., status=status)
    if (status /= 0) stop 'FAILED: filter_detectors.'

    call pacs%compute_mapheader(pointing, time, 3.d0, header, status)
    if (status /= 0) stop 'FAILED: compute_mapheader.'

    call ft_header2wcs(header, wcs, nx, ny, status)
    if (status /= 0) stop 'FAILED: ft_header2wcs.'

    ! allocate memory for the map
    allocate(map1d(0:nx*ny-1))

    ! read the signal file
    write(*,'(a)', advance='no') 'Reading signal file... '
    call system_clock(count1, count_rate, count_max)
    allocate(signal(last-first+1, pacs%ndetectors))
    call pacs%read_signal_file(filename_signal, first, last, signal, status)
    if (status /= 0) stop 'FAILED: read_signal_file.'
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
    allocate(pmatrix(9,last-first+1,pacs%ndetectors))
    call system_clock(count1, count_rate, count_max)
    call pacs%compute_projection_sharp_edges(pointing, time, header, nx, ny, pmatrix, status)
    if (status /= 0) stop 'FAILED: compute_projection_sharp_edges.'
    call system_clock(count2, count_rate, count_max)
    write(*,'(f6.2,a)') real(count2-count1)/count_rate, 's'

    ! check flux conservation during backprojection
    write(*,'(a)', advance='no') 'Testing flux conservation...'
    call system_clock(count1, count_rate, count_max)
    allocate(surface1(nsamples, pacs%ndetectors))
    allocate(surface2(nsamples, pacs%ndetectors))
    allocate(coords(ndims, pacs%ndetectors*nvertices))

    call init_wcslib(header, status)
    if (status /= 0) stop 'FAILED: init_wcslib.'

    index = 2
    write(*,*) 'surfaces:'
    !$omp parallel do default(none) firstprivate(index) private(isample, ra, dec, pa, chop, coords, idetector) &
    !$omp shared(time, pointing, nx, pacs, pmatrix, wcs, surface1, surface2)
    do isample = 1, nsamples
        call pointing%get_position(time(isample), ra, dec, pa, chop, index)
        !write(*,*) isample, ', Time:',time(isample),'RA:',ra,'Dec',dec
        coords = pacs%uv2yz(pacs%corners_uv, pacs%distortion_yz_blue, chop)
        coords = pacs%yz2ad(coords, ra, dec, pa)
        coords = ad2xy_wcslib(coords)
        do idetector = 1, pacs%ndetectors
            surface1(isample,idetector) = abs(surface_convex_polygon(coords(:,(idetector-1)*nvertices+1:idetector*nvertices)))
            surface2(isample,idetector) = sum(pmatrix(:,isample,idetector)%weight)
        end do
        !write (*,*) isample, 1, surface1(isample,1), surface2(isample,1)
    end do
    !$omp end parallel do

    call system_clock(count2, count_rate, count_max)
    write(*,'(f6.2,a)') real(count2-count1)/count_rate, 's'
    write(*,*) 'Difference: ', maxval(abs((surface1-surface2)/surface1))
    !XXX

    ! back project the timeline
    write(*,'(a)', advance='no') 'Computing the back projection... '
    call system_clock(count1, count_rate, count_max)
    call pmatrix_transpose(pmatrix, signal, map1d)
    call system_clock(count2, count_rate, count_max)
    write(*,'(f6.2,a)') real(count2-count1)/count_rate, 's'

    ! test the back projected map
    if (.not. test_real_eq(sum(map1d), -0.19967872d0, 5)) then
        write (ERROR_UNIT,*) 'Sum in map is ', sum(map1d), ' instead of ', -0.19967872
    end if

    ! project the map
    write(*,'(a)', advance='no') 'Computing the forward projection... '
    call system_clock(count1, count_rate, count_max)
    call pmatrix_direct(pmatrix, map1d, signal)
    call system_clock(count2, count_rate, count_max)
    write(*,'(f6.2,a)') real(count2-count1)/count_rate, 's'
    write(*,*) 'total: ', sum(signal) ! -4883185.7722900705

    ! test the back projected map
    if (.not. test_real_eq(sum(signal), -4883185.772d0, 5)) then
        write (ERROR_UNIT,*) 'Sum in timeline is ', sum(signal), ' instead of ', -4883185.772d0
    end if

    ! write the map as fits file
    write(*,'(a)') 'Writing FITS file... '
    call ft_write(outputdir // 'ngc6946_bpj.fits', reshape(map1d, [nx,ny]), wcs, status)
    if (status /= 0) stop 'FAILED: ft_write.'

    call free_wcslib()

    call system_clock(count2, count_rate, count_max)
    write(*,'(a,f6.2,a)') 'Total elapsed time: ', real(count2-count0)/count_rate, 's'

    flush(OUTPUT_UNIT)
    stop "OK."

end program test_ngc6946_bpj
