program test_pacs

    use, intrinsic :: ISO_C_BINDING
    use module_fitstools
    use module_wcs, only : init_astrometry, ad2xy_gnomonic
    use module_pacsinstrument
    use module_pacspointing
    implicit none

    class(pacsinstrument), allocatable  :: pacs
    class(pacspointing), allocatable :: pointing
    character(len=*), parameter :: filename = '/home/pchanial/work/pacs/data/transparent/NGC6946/1342184520_blue'
    character(len=*), parameter :: filename_header = 'tests/csh_header_ngc6946.fits'

    real*8, allocatable    :: yz(:,:), ad(:,:), xy(:,:), time(:)
    integer*8, allocatable :: timeus(:)
    integer                :: status, i, j, nx, ny, index
    !integer                :: count, count2, count_rate, count_max
    integer*8              :: first, last
    real*8                 :: ra, dec, pa, chop, xmin, xmax, ymin, ymax, ra0, dec0
    character(len=2880)    :: header

    ! read observation

    allocate(pointing)
    call pointing%load(filename, status)
    if (status /= 0) stop 'pointing%load: FAILED.'
    call pointing%print()

    first = 12001
    last  = 86936
    allocate(time(last-first+1))
    allocate(timeus(last-first+1))
    call ft_readslice(filename // '_Time.fits+1', first, last, timeus, status)
    if (status /= 0) stop 'ft_readslice: FAILED.'

    time = timeus * 1.0d-6
    allocate(pacs)

    call pacs%read_calibration_files(status)
    if (status /= 0) stop 'ft_readslice: FAILED.'
    call pacs%filter_detectors(pacs%get_array_color(filename), transparent_mode=.true., status=status)
    if (status /= 0) stop 'filter_detectors: FAILED.'

    call ft_header2str(filename_header, header, status)
    if (status /= 0) stop 'ft_header2str: FAILED.'

    call init_astrometry(header, status=status)
    if (status /= 0) stop 'init_astrometry: FAILED.'

    write (*,*) 'Number of valid detectors:', pacs%ndetectors
    write (*,*) 'Number of samples:', size(time)

    write (*,*) 'shape corners_uv: ', shape(pacs%corners_uv)
    write (*,*) 'Detector 1 p,q,i,j: ', pacs%pq(:,1), pacs%ij(:,1)
    write (*,*) 'uv(1,1:4)', pacs%corners_uv(1,1:4)
    write (*,*) 'uv(2,1:4)', pacs%corners_uv(2,1:4)
    allocate(yz(ndims, size(pacs%corners_uv,2)))
    allocate(ad(ndims, size(pacs%corners_uv,2)))
    allocate(xy(ndims, size(pacs%corners_uv,2)))

    do i=1, 1
       write (*,*) 'U: ', i, (pacs%corners_uv(1, (i-1)*4+j), j = 1,4)
       write (*,*) 'V: ', i, (pacs%corners_uv(2, (i-1)*4+j), j = 1,4)
    end do

    index = 2
    call pointing%get_position(time(1), ra, dec, pa, chop, index)

    yz = pacs%uv2yz(pacs%corners_uv, pacs%distortion_yz_blue, chop)

    do i=1, 1
       write (*,*) 'Y: ', i, (yz(1, (i-1)*4+j), j = 1,4)
       write (*,*) 'Z: ', i, (yz(2, (i-1)*4+j), j = 1,4)
    end do

    ad = pacs%yz2ad(yz, ra, dec, pa)

    do i=1, 1
       write (*,*) 'RA:  ', i, (ad(1, (i-1)*4+j), j = 1,4)
       write (*,*) 'Dec: ', i, (ad(2, (i-1)*4+j), j = 1,4)
    end do

    xy = ad2xy_gnomonic(ad)

    do i=1, 1
        write (*,*) 'x: ', i, (xy(1, (i-1)*4+j), j = 1,4)
        write (*,*) 'y: ', i, (xy(2, (i-1)*4+j), j = 1,4)
    end do

    write(*,*) 'minmaxU', minval(pacs%corners_uv(1,:)), maxval(pacs%corners_uv(1,:))
    write(*,*) 'minmaxV', minval(pacs%corners_uv(2,:)), maxval(pacs%corners_uv(2,:))
    write(*,*) 'minmaxY', minval(yz(1,:)), maxval(yz(1,:))
    write(*,*) 'minmaxZ', minval(yz(2,:)), maxval(yz(2,:))
    write(*,*) 'minmaxX', minval(xy(1,:)), maxval(xy(1,:))
    write(*,*) 'minmaxY', minval(xy(2,:)), maxval(xy(2,:))

    ! minmax(pa) 252.34927       297.83408
    ! minmax(ra) 258.50445       309.04262
    ! minmax(dec) 59.980161       67.110813

    call pointing%compute_center(time, ra0, dec0)
    write(*,*) 'RA center: ', ra0, ', Dec center: ', dec0
    !print, format='(f17.12)', mean(ra.cube.__data__[12000:86936-1]), mean(dec.cube.__data__[12000:86936-1])
    !R.A.: 308.722596918714
    !Dec:   60.156005722598

    call pacs%find_minmax(pointing, time, xmin, xmax, ymin, ymax)
    write (*,*) 'X = ', xmin, xmax
    write (*,*) 'Y = ', ymin, ymax

    call pacs%compute_mapheader(pointing, time, 3.d0, header, status)
    if (status /= 0) stop 'compute_mapheader: FAILED.'

    !! check that all detector corners are inside the map
    !call system_clock(count, count_rate, count_max)
    !do i = 1, size(time)
    !    call pointing%get_position(time(i), ra, dec, pa, chop)
    !    yz = pacs%uv2yz(pacs%corners_uv, pacs%distortion_yz_blue, chop)
    !    ad = pacs%yz2ad(yz, ra, dec, pa)
    !    xy = pacs%ad2xy(ad, wcs)
    !    if (any(xy(1,:) < -0.5d0 .or. xy(1,:) > nx-0.5d0 .or. xy(2,:) < -0.5d0 .or. xy(2,:) > ny-0.5d0)) then
    !        write (*,*) 'XXX problem: ', i
    !    end if
    !end do
    !call system_clock(count2, count_rate, count_max)
    !write (*,'(a,f6.2,a)') 'elapsed time: ', real(count2-count)/count_rate, 's'
    !! get_position: 0.00s
    !!+uv2yz: 30.85s
    !!+yz2ad: 40.10s,40.72s
    !!+ad2xy: 71.54s,76.74s,69.16s
    !!+if: 70.52s

    stop "OK."

end program test_pacs
