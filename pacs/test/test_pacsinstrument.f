program test_pacsinstrument

    use iso_fortran_env,        only : ERROR_UNIT, OUTPUT_UNIT
    use module_fitstools,       only : ft_read_keyword, ft_header2str
    use module_math,            only : mean, neq_real
    use module_observation,     only : MaskPolicy, Observation
    use module_pacsinstrument
    use module_pacsobservation, only : PacsObservation
    use module_string,          only : strreal
    use module_tamasis,         only : p, POLICY_KEEP, POLICY_MASK, POLICY_REMOVE
    use module_wcs,             only : init_astrometry, ad2xy_gnomonic
    implicit none

    class(Observation), allocatable     :: obs
    class(PacsInstrument), allocatable  :: pacs
    class(PacsObservation), allocatable :: pacsobs
    type(MaskPolicy)                    :: policy
    character(len=*), parameter :: filename(1) = 'pacs/test/data/frames_blue.fits'
    character(len=*), parameter :: filename_header = 'core/test/data/header.fits'

    real(p), allocatable   :: yz(:,:), ad(:,:), xy(:,:), time(:)
    integer                :: status, i, j, nx, ny, index
    integer                :: count1, count2, count_rate, count_max
    real(p)                :: ra, dec, pa, chop, xmin, xmax, ymin, ymax, ra0, dec0
    character(len=2880*2)  :: header
    logical*1              :: detector_mask_blue(32,64)
    logical*1              :: detector_mask_red (16,32)

    logical*1, allocatable :: detector_bad_all(:,:)
    real(p), allocatable   :: detector_center_all(:,:,:)
    real(p), allocatable   :: detector_corner_all(:,:,:,:)
    real(p), allocatable   :: detector_area_all(:,:)
    real(p), allocatable   :: flatfield_optical_all(:,:)
    real(p), allocatable   :: flatfield_detector_all(:,:)
    real(p)                :: distortion_yz(2,3,3,3)
    real(p)                :: active_fraction
    real(p)                :: signal_ref(360)
    logical*1              :: mask_ref(360)
    real(p), allocatable   :: signal(:,:)
    logical*1, allocatable :: mask(:,:)
    integer                :: idetector, first, last

    real(p), allocatable   :: a_vect(:), d_vect(:), ad_vect(:,:)
    integer                :: n

    ! initialise pacs instrument
    allocate (pacs)

    ! read calibration files
    active_fraction = 0
    call pacs%read_calibration_files('blue', detector_bad_all, detector_center_all, detector_corner_all, detector_area_all,        &
                                     flatfield_optical_all, flatfield_detector_all, distortion_yz, active_fraction, status)
    if (status /= 0) call failure('read_calibration_files')

    detector_mask_blue = .false.
    call pacs%init_with_calfiles('blue', detector_mask_blue, 1, status)
    if (status /= 0) call failure('pacs%init_with_calfiles all')

    if (pacs%ndetectors /= 2048) call failure('ndetectors all')
    if (any(shape(pacs%detector_corner) /= [2,4*pacs%ndetectors])) call failure('shape detector_corner all')

    if (any(pacs%pq(:,2)   /= [0,1])) call failure('pq1 all')
    if (any(pacs%ij(:,2)   /= [0,1])) call failure('ij1 all')
    if (any(pacs%detector_corner(:,5:8) /= detector_corner_all(:,:,1,2))) call failure('uv1 all')

    i = count(.not. pacs%mask(1:8,1:16))
    if (any(pacs%pq(:,1+i) /= [8,0])) call failure('pq2 all')
    if (any(pacs%ij(:,1+i) /= [8,0])) call failure('ij2 all')
    if (any(pacs%detector_corner(:,i*4+1:i*4+4) /= detector_corner_all(:,:,9,1))) call failure('uv2 all')

    i = count(.not. pacs%mask(1:16,:)) + count(.not. pacs%mask(17:32,1:16))
    if (any(pacs%pq(:,1+i) /= [16,16])) call failure('pq3 all')
    if (any(pacs%ij(:,1+i) /= [0,0])) call failure('ij3 all')
    if (any(pacs%detector_corner(:,i*4+1:i*4+4) /= detector_corner_all(:,:,17,17))) call failure('uv3 all')

    call pacs%destroy

    ! test instrument setup, reading of signal & mask
    detector_mask_blue = detector_bad_all
    call pacs%init_with_calfiles('blue', detector_mask_blue, 1, status)
    if (status /= 0) call failure('pacs%init_with_calfiles rem')

    if (pacs%ndetectors /= 1997) call failure('ndetectors rem')
    if (any(shape(pacs%detector_corner) /= [2,4*pacs%ndetectors])) call failure('shape detector_corner rem')

    if (any(pacs%pq(:,2)   /= [0,1])) call failure('pq1 rem')
    if (any(pacs%ij(:,2)   /= [0,1])) call failure('ij1 rem')
    if (any(pacs%detector_corner(:,5:8) /= detector_corner_all(:,:,1,2))) call failure('uv1 rem')

    i = count(.not. pacs%mask(1:8,1:16))
    if (any(pacs%pq(:,1+i) /= [8,0])) call failure('pq2 rem')
    if (any(pacs%ij(:,1+i) /= [8,0])) call failure('ij2 rem')
    if (any(pacs%detector_corner(:,i*4+1:i*4+4) /= detector_corner_all(:,:,9,1))) call failure('uv2 rem')

    i = count(.not. pacs%mask(1:16,:)) + count(.not. pacs%mask(17:32,1:16))
    if (any(pacs%pq(:,1+i) /= [16,16])) call failure('pq3 rem')
    if (any(pacs%ij(:,1+i) /= [0,0])) call failure('ij3 rem')
    if (any(pacs%detector_corner(:,i*4+1:i*4+4) /= detector_corner_all(:,:,17,17))) call failure('uv3 rem')
    call pacs%destroy
    
    ! initialise observation
    detector_mask_blue = .false.
    call pacs%init_with_calfiles('blue', detector_mask_blue, 1, status)
    if (status /= 0) call failure('pacs%init_with_calfiles rem')
    allocate (pacsobs)
    call pacsobs%init(filename, policy, status)
    if (status /= 0) call failure('init_pacsobservation')
    allocate (obs)
    call pacsobs2obs(pacsobs, 1, obs, status)
    if (status /= 0) call failure('pacsobs2obs')

    allocate(time(obs%nsamples))
    time = obs%pointing%time
    if (size(time) /= 360) call failure('nsamples')

    if (pacsobs%band /= 'blue' .or. pacsobs%slice(1)%observing_mode /= 'prime') call failure('obs%init')

    allocate (signal(pacsobs%slice(1)%nvalids,pacs%ndetectors))
    allocate (mask  (pacsobs%slice(1)%nvalids,pacs%ndetectors))
    call pacs%read(pacsobs%slice(1)%filename, obs%pointing%masked, obs%pointing%removed, ['master'], signal, mask, status)
    if (status /= 0) call failure('pacs%read')

    ! get a detector that has been hit by a glitch
    do idetector = 1, pacs%ndetectors
        if (pacs%pq(1,idetector)+1 == 27.and.pacs%pq(2,idetector)+1 == 64) exit
    end do

    signal_ref = signal(:,idetector)
    mask_ref = .false.
    mask_ref([5,6,7,151]) = .true. ! (27,64)
    if (any(mask_ref .neqv. mask(:,idetector))) call failure('mask ref')
    deallocate (mask, signal)
    
    do first = 1, 360, 7
        do last = first, 360, 11
           obs%pointing%removed = .true.
           obs%pointing(first:last)%removed = .false.
           obs%nvalids = last - first + 1
           allocate(signal(obs%nvalids, pacs%ndetectors))
           allocate(mask  (obs%nvalids, pacs%ndetectors))
           call pacs%read(pacsobs%slice(1)%filename, obs%pointing%masked, obs%pointing%removed, ['master'], signal, mask, status)
           if (status /= 0) call failure('read_pacsobservation loop')
           if (any(signal(:,idetector) /= signal_ref(first:last))) call failure('read signal')
           if (any(mask(:,idetector) .neqv. mask_ref(first:last))) then
               print *,first, last
               print *,shape(mask(:,idetector))
               print *,shape(mask_ref)
               print *,mask(:,idetector)
               print *,mask_ref(first:last)
               call failure('read mask')
           end if
           deallocate (signal)
           deallocate (mask)
        end do
    end do
    obs%pointing%removed = .false.

    allocate(yz(ndims, size(pacs%detector_corner,2)))
    allocate(ad(ndims, size(pacs%detector_corner,2)))
    allocate(xy(ndims, size(pacs%detector_corner,2)))
    call ft_header2str(filename_header, header, status)
    if (status /= 0) call failure('ft_header2str.')
    call init_astrometry(header, status=status)
    if (status /= 0) call failure('init_astrometry.')

    index = 0
    call obs%get_position_time(time(1), ra, dec, pa, chop, index)

    yz = pacs%uv2yz(pacs%detector_corner, pacs%distortion_yz, chop)
    ad = pacs%yz2ad(yz, ra, dec, pa)
    xy = ad2xy_gnomonic(ad)

    if (any(neq_real([(yz(1,j), j = 1,4)], [3.06059075250975504E-002_p,  2.98567567627381451E-002_p, &
        2.98274651048868189E-002_p,  3.05767270814490018E-002_p]))) call failure('Y')
    if (any(neq_real([(yz(2,j), j = 1,4)], [1.30515691627312649E-002_p,  1.30921308057738196E-002_p, &
        1.23364621524597352E-002_p,  1.22964967432396752E-002_p]))) call failure('Z')
    if (any(neq_real([(ad(1,j), j = 1,4)], [245.93544838257338_p, 245.93659969792549_p, 245.93766190786474_p, &
        245.93650960929841_p]))) call failure('RA')
    if (any(neq_real([(ad(2,j), j = 1,4)], [61.510272502550855_p, 61.509760799263510_p, 61.510321692744220_p, &
        61.510833009997000_p]))) call failure('Dec')
    if (any(neq_real([(xy(1,j), j = 1,4)], [33288.533194678035_p, 33288.643767997528_p, 33287.493446513588_p, &
        33287.383783655569_p]))) call failure('x')
    if (any(neq_real([(xy(2,j), j = 1,4)], [19107.874613006607_p, 19106.828248107082_p, 19106.616041863876_p, &
        19107.662684113882_p]))) call failure('y')

    if (any(neq_real([minval(pacs%detector_corner(1,:)), maxval(pacs%detector_corner(1,:))], &
        [-13.276956558227539_p, 13.294145584106445_p]))) call failure('minmaxU')
    if (any(neq_real([minval(pacs%detector_corner(2,:)), maxval(pacs%detector_corner(2,:))], &
        [-25.931186676025391_p, 25.915901184082031_p]))) call failure('minmaxV')
    if (any(neq_real([minval(yz(1,:)), maxval(yz(1,:))], &
        [-3.04914370632606777E-002_p, 3.06059075250975504E-002_p]))) call failure('minmaxY')
    if (any(neq_real([minval(yz(2,:)), maxval(yz(2,:))], &
        [-1.55613991289617025E-002_p, 1.61528327291844132E-002_p]))) call failure('minmaxZ')
    if (any(neq_real([minval(ad(1,:)), maxval(ad(1,:))], &
        [245.93544838257338_p, 246.06808069685269_p]))) call failure('minmaxRA')
    if (any(neq_real([minval(ad(2,:)), maxval(ad(2,:))], &
        [61.469595978715510_p, 61.531092935763752_p]))) call failure('minmaxDec')
    if (any(neq_real([minval(xy(1,:)), maxval(xy(1,:))], [33245.023275401254_p, 33297.106116672716_p]))) call failure('minmaxx')
    if (any(neq_real([minval(xy(2,:)), maxval(xy(2,:))], [19016.089802392329_p, 19107.874613006607_p]))) call failure('minmaxy')

    call obs%compute_center(ra0, dec0)
    if (neq_real(ra0,  mean(obs%pointing%ra), 1.e-8_p) .or. &
        neq_real(dec0, mean(obs%pointing%dec), 1.e-8_p)) call failure('compute_center')

    call pacs%find_minmax(obs, .false., xmin, xmax, ymin, ymax)
    if (neq_real(xmin, 33215.245057167944_p)) call failure('xmin')
    if (neq_real(xmax, 33297.106116672716_p)) call failure('xmax')
    if (neq_real(ymin, 19008.137140577746_p)) call failure('ymin')
    if (neq_real(ymax, 19141.196499174162_p)) call failure('ymax')

    call pacs%compute_map_header(obs, .false., 3._p, header, status)
    if (status /= 0) call failure('compute_map_header')

    call ft_read_keyword(header, 'naxis1', nx, status=status)
    if (status /= 0) call failure('read_param NAXIS1')

    call ft_read_keyword(header, 'naxis2', ny, status=status)
    if (status /= 0) call failure('read_param NAXIS2')

    call init_astrometry(header, status=status)
    if (status /= 0) call failure('init_astrometry')

    ! check that all detector corners are inside the map
    index = 0
    call system_clock(count1, count_rate, count_max)
    do i = 1, size(time)
        call obs%get_position_time(time(i), ra, dec, pa, chop, index)
        yz = pacs%uv2yz(pacs%detector_corner, pacs%distortion_yz, chop)
        ad = pacs%yz2ad(yz, ra, dec, pa)
        xy = ad2xy_gnomonic(ad)
        if (any(xy(1,:) < 0.5_p .or. xy(1,:) > nx+0.5_p .or. xy(2,:) < 0.5_p .or. xy(2,:) > ny+0.5_p)) then
            write (*,'(3(a,i0))') 'sample:', i, ', nx: ', nx, ', ny: ', ny
            write (*,*) 'X: ', minval(xy(1,:)), maxval(xy(1,:))
            write (*,*) 'Y: ', minval(xy(2,:)), maxval(xy(2,:))
            call failure('detector outside of map')
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

    a_vect = [(2.3+(10._p*i)/n, i=1,n)]
    d_vect = [(3.8+(10._p*i)/n, i=1,n)]
    ad_vect(1,:) = a_vect
    ad_vect(2,:) = d_vect
    ra  = 12.3_p
    dec = 60.3_p
    pa  = 2.7_p

    call system_clock(count1, count_rate, count_max)
    ad_vect = pacs%yz2ad(ad_vect, ra, dec, pa)
    call system_clock(count2, count_rate, count_max)
    write(*,'(a)') 'yz2ad: ' // strreal(real(count2-count1,p)/count_rate,4) // 's'
    call pacs%destroy()

    detector_mask_red = .false.
    call pacs%init_with_calfiles('red', detector_mask_red, 1, status)
    if (status /= 0) call failure('pacs%init 1')
    call pacs%destroy()

    call pacs%init_with_calfiles('green', detector_mask_blue, 1, status)
    if (status /= 0) call failure('pacs%init 2')
    call pacs%destroy()

    call pacs%init_with_calfiles('blue', detector_mask_blue, 1, status)
    if (status /= 0) call failure('pacs%init 3')
    call pacs%destroy()

    write (OUTPUT_UNIT,'(a)') 'Testing error...'
    call pacs%init_with_calfiles('x', detector_mask_blue, 1, status)
    if (status == 0) call failure('pacs%init 4')
    write (OUTPUT_UNIT,'(a,/)') 'OK.'

    write (OUTPUT_UNIT,'(a)') 'Testing error...'
    call pacs%init_with_calfiles('blue', detector_mask_blue, 3, status)
    if (status == 0) call failure('pacs%init 5')
    write (OUTPUT_UNIT,'(a,/)') 'OK.'

contains

    subroutine failure(errmsg)
        character(len=*), intent(in) :: errmsg
        write (ERROR_UNIT,'(a)'), 'FAILED: ' // errmsg
        stop 1
    end subroutine failure


    subroutine pacsobs2obs(pacsobs, islice, obs, status)
        class(PacsObservation), intent(in) :: pacsobs
        integer, intent(in)                :: islice
        class(Observation), intent(out)    :: obs
        integer, intent(out)               :: status
        
        call obs%init(pacsobs%slice(islice)%p%time, pacsobs%slice(islice)%p%ra, pacsobs%slice(islice)%p%dec,                       &
                      pacsobs%slice(islice)%p%pa, pacsobs%slice(islice)%p%chop, pacsobs%slice(islice)%p%masked,                    &
                      pacsobs%slice(islice)%p%removed, pacsobs%slice(islice)%compression_factor,                                   &
                      (pacsobs%slice(islice)%compression_factor - 1) / (2._p * pacsobs%slice(islice)%compression_factor), status)

    end subroutine pacsobs2obs

end program test_pacsinstrument

! test before patch to handle several observations
!!$61) sapherschel1.extra.cea.fr:~pchanial/work/tamasis/tamasis> make test_pacs && ./test_pacs
!!$Warning: the pointing time is not evenly spaced or is drifting.
!!$ Time:    1634799159.4714119      ...   1634799427.3234940     
!!$ RA  :    245.97578268187323      ...   246.02256811506965     
!!$ Dec :    61.500650548286465      ...   61.529470882067791     
!!$ PA  :    219.88090872523742      ...   219.92216426552577     
!!$ Chop:  -6.50595714733942510E-004 ...  6.50599516571009895E-004
!!$Info: Band: 'Blue'
!!$Info: Compression mode: 'Photometry Default Mode'
!!$ Number of valid detectors:        1998
!!$ Number of samples:         360
!!$ shape detector_corner:            2        7992
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
