program test_pacsinstrument

    use module_fitstools,       only : ft_read_parameter, ft_header2str
    use module_math,            only : mean, neq_real
    use module_pacsinstrument
    use module_pacsobservation, only : pacsobservation, maskarray
    use module_string,          only : strreal
    use module_wcs,             only : init_astrometry, ad2xy_gnomonic
    implicit none

    class(pacsinstrument), allocatable  :: pacs
    class(pacsobservation), allocatable :: obs
    type(maskarray)                     :: maskarray_policy
    character(len=*), parameter :: filename(1) = 'tests/frames_blue.fits'
    character(len=*), parameter :: filename_header = 'tests/csh_header_ngc6946.fits'

    real*8, allocatable    :: yz(:,:), ad(:,:), xy(:,:), time(:)
    integer                :: status, i, j, nx, ny, index
    integer                :: count1, count2, count_rate, count_max
    real*8                 :: ra, dec, pa, chop, xmin, xmax, ymin, ymax, ra0, dec0
    character(len=2880)    :: header

    real*8, allocatable    :: a_vect(:), d_vect(:), ad_vect(:,:)
    integer                :: n

    ! initialise observation
    allocate(obs)
    call obs%init(filename, maskarray_policy, status)
    if (status /= 0) stop 'FAILED: init_pacsobservation'

    ! initialise pacs instrument
    allocate(pacs)
    call pacs%init(obs%channel, obs%observing_mode == 'Transparent', 1, .false., status)
    if (status /= 0) stop 'FAILED: pacsinstrument%init'

    allocate(time(obs%nsamples_tot))
    time = obs%slice(1)%p%time

    if (pacs%ndetectors /= 1998) stop 'FAILED: ndetectors'
    if (size(time) /= 360) stop 'FAILED: nsamples'
    if (any(shape(pacs%corners_uv) /= [2,4*pacs%ndetectors])) stop 'FAILED: shape corners_uv'

    i = count(.not. pacs%mask_blue(1,:))

    if (any(pacs%pq(:,2)   /= [0,1])) stop 'FAILED pq1'
    if (any(pacs%pq(:,1+i) /= [1,0])) stop 'FAILED pq2'
    if (any(pacs%ij(:,2)   /= [0,1])) stop 'FAILED ij1'
    if (any(pacs%ij(:,1+i) /= [1,0])) stop 'FAILED ij2'

    if (any(pacs%corners_uv(:,5:8) /= pacs%corners_uv_blue(:,:,1,2)))          &
        stop 'FAILED uv1'
    if (any(pacs%corners_uv(:,i*4+1:i*4+4) /= pacs%corners_uv_blue(:,:,2,1)))  &
        stop 'FAILED uv2'

    allocate(yz(ndims, size(pacs%corners_uv,2)))
    allocate(ad(ndims, size(pacs%corners_uv,2)))
    allocate(xy(ndims, size(pacs%corners_uv,2)))
    call ft_header2str(filename_header, header, status)
    if (status /= 0) stop 'FAILED: ft_header2str.'
    call init_astrometry(header, status=status)
    if (status /= 0) stop 'FAILED: init_astrometry.'

    index = 0
    call obs%get_position_time(1, time(1), ra, dec, pa, chop, index)

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

    call obs%compute_center(ra0, dec0)
    if (neq_real(ra0,  mean(obs%slice(1)%p%ra), 15) .or. neq_real(dec0, mean(obs%slice(1)%p%dec), 15)) stop 'FAILED: compute_center'

    call pacs%find_minmax(obs, .false., xmin, xmax, ymin, ymax)
    write (*,*) 'Xmin: ', xmin, ', Xmax: ', xmax
    write (*,*) 'Ymin: ', ymin, ', Ymax: ', ymax

    call pacs%compute_map_header(obs, .false., 3.d0, header, status)
    if (status /= 0) stop 'FAILED: compute_map_header'

    call ft_read_parameter(header, 'naxis1', nx, count1, status=status)
    if (status /= 0 .or. count1 == 0) stop 'FAILED: read_param NAXIS1'

    call ft_read_parameter(header, 'naxis2', ny, count1, status=status)
    if (status /= 0 .or. count1 == 0) stop 'FAILED: read_param NAXIS2'

    call init_astrometry(header, status=status)
    if (status /= 0) stop 'FAILED: init_astrometry'

    ! check that all detector corners are inside the map
    index = 0
    call system_clock(count1, count_rate, count_max)
    do i = 1, size(time)
        call obs%get_position_time(1, time(i), ra, dec, pa, chop, index)
        yz = pacs%uv2yz(pacs%corners_uv, pacs%distortion_yz_blue, chop)
        ad = pacs%yz2ad(yz, ra, dec, pa)
        xy = ad2xy_gnomonic(ad)
        if (any(xy(1,:) < 0.5d0 .or. xy(1,:) > nx+0.5d0 .or. xy(2,:) < 0.5d0 .or. xy(2,:) > ny+0.5d0)) then
            write (*,'(3(a,i0))') 'sample:', i, ', nx: ', nx, ', ny: ', ny
            write (*,*) 'X: ', minval(xy(1,:)), maxval(xy(1,:))
            write (*,*) 'Y: ', minval(xy(2,:)), maxval(xy(2,:))
            stop 'FAILED: detector outside of map'
        end if
    end do
    call system_clock(count2, count_rate, count_max)
    write (*,'(a,f6.2,a)') 'elapsed time: ', real(count2-count1)/count_rate, 's'
    !! get_position: 0.00s
    !+uv2yz: 30.85s
    !+yz2ad: 40.10s,40.72s
    !+ad2xy: 71.54s,76.74s,69.16s
    !+if: 70.52s

    n = 10000000
    allocate(a_vect(n))
    allocate(d_vect(n))
    allocate(ad_vect(2,n))

    a_vect = [(2.3+(10.d0*i)/n, i=1,n)]
    d_vect = [(3.8+(10.d0*i)/n, i=1,n)]
    ad_vect(1,:) = a_vect
    ad_vect(2,:) = d_vect
    ra  = 12.3d0
    dec = 60.3d0
    pa  = 2.7d0

    call system_clock(count1, count_rate, count_max)
    ad_vect = pacs%yz2ad(ad_vect, ra, dec, pa)
    call system_clock(count2, count_rate, count_max)
    write(*,'(a)') 'yz2ad: ' // strreal(real(count2-count1,kind=8)/count_rate,4) // 's'
    call pacs%destroy()

    call pacs%init('r', .true., 1, .false., status)
    if (status /= 0) stop 'FAILED: pacs%init 1'
    call pacs%destroy()

    call pacs%init('g', .true., 1, .false., status)
    if (status /= 0) stop 'FAILED: pacs%init 2'
    call pacs%destroy()

    call pacs%init('b', .true., 1, .false., status)
    if (status /= 0) stop 'FAILED: pacs%init 3'
    call pacs%destroy()

    call pacs%init('x', .true., 1, .false., status)
    if (status == 0) stop 'FAILED: pacs%init 4'
    call pacs%init('b', .true., 3, .false., status)
    if (status == 0) stop 'FAILED: pacs%init 5'

    stop "OK."
    
end program test_pacsinstrument

! test before patch to handle several observations
!!$61) sapherschel1.extra.cea.fr:~pchanial/work/tamasis/tamasis> make test_pacs && ./test_pacs
!!$Warning: the pointing time is not evenly spaced or is drifting.
!!$ Time:    1634799159.4714119      ...   1634799427.3234940     
!!$ RA  :    245.97578268187323      ...   246.02256811506965     
!!$ Dec :    61.500650548286465      ...   61.529470882067791     
!!$ PA  :    219.88090872523742      ...   219.92216426552577     
!!$ Chop:  -6.50595714733942510E-004 ...  6.50599516571009895E-004
!!$Info: Channel: 'Blue'
!!$Info: Compression mode: 'Photometry Default Mode'
!!$ Number of valid detectors:        1998
!!$ Number of samples:         360
!!$ shape corners_uv:            2        7992
!!$ Detector 1 p,q,i,j:            0           0           0           0
!!$ uv(1,1:4)  -11.536449432373047       -11.556140899658203       -10.916444778442383       -10.896753311157227     
!!$ uv(2,1:4)  -25.931186676025391       -25.291488647460938       -25.271797180175781       -25.911495208740234     
!!$ U:            1  -11.536449432373047       -11.556140899658203       -10.916444778442383       -10.896753311157227     
!!$ V:            1  -25.931186676025391       -25.291488647460938       -25.271797180175781       -25.911495208740234     
!!$ Y:            1  3.19859511193604923E-002  3.12056263946425688E-002  3.11695109418158520E-002  3.19499492715473385E-002
!!$ Z:            1  1.33762097219791900E-002  1.34239330250266724E-002  1.26367501341196099E-002  1.25896466991417492E-002
!!$ RA:             1   245.93279096919014        245.93398291039887        245.93509840935835        245.93390544965087     
!!$ Dec:            1   61.510908778127934        61.510371582612237        61.510952273587293        61.511489066471981     
!!$ x:            1   33288.938397159887        33289.061945811991        33287.864028896096        33287.741427447902     
!!$ y:            1   19109.892561738641        19108.803862105047        19108.574881406410        19109.663865884180     
!!$ minmaxU  -13.276956558227539        13.294145584106445     
!!$ minmaxV  -25.931186676025391        25.915901184082031     
!!$ minmaxY -3.18624246565269373E-002  3.19859511193604923E-002
!!$ minmaxZ -1.64247369758722117E-002  1.70428953370297294E-002
!!$ minmaxX   33243.636219590051        33298.534260959212     
!!$ minmaxY   19014.092920793282        19109.892561738641     
!!$ RA center:    245.99842678324958      , Dec center:    61.514764712687644     
!!$ X =    33213.858197672103        33298.534260959212     
!!$ Y =    19006.140539825199        19143.214515560223     
!!$elapsed time:   1.28s
!!$STOP OK.
