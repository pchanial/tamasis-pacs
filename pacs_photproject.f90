program pacs_photproject

    use, intrinsic :: ISO_FORTRAN_ENV
    use            :: module_fitstools
    use            :: module_deglitching, only : deglitch_l2b
    use            :: module_optionparser, only : optionparser
    use            :: module_pacsinstrument
    use            :: module_pacspointing
    use            :: module_pointingmatrix, only : backprojection_weighted
    use            :: module_preprocessor
    implicit none

    character(len=*), parameter :: filename = '/home/pchanial/work/pacs/data/transparent/NGC6946/1342184520_blue'
    type(pacsinstrument), allocatable  :: pacs
    type(pacspointing), allocatable    :: pointing

    character(len=2048)                :: infile, outfile, headerfile
    character(len=2880)                :: header
    real*8                             :: resolution
    integer                            :: npixels_per_sample
    logical                            :: do_flatfield, do_meansubtraction
    logical                            :: do_deglitching_l2b, do_deglitching_l2c
    integer*8                          :: first, last

    real*8, allocatable                :: signal(:,:)
    logical*1, allocatable             :: mask(:,:)
    real*8, allocatable                :: time(:)
    integer*8, allocatable             :: timeus(:)
    integer                            :: nx, ny
    integer                            :: status, count, count1
    integer                            :: count2, count_rate, count_max
    integer                            :: nsamples
    real*8                             :: nsigma_deglitching
    real*8, allocatable                :: map1d(:), weights(:)
    type(pointingelement), allocatable :: pmatrix(:,:,:)
    type(optionparser), allocatable    :: parser

    ! command parsing
    allocate(parser)
    call parser%init('pacs_photproject [options] fitsfile', 0, 1)
    call parser%add_option('', 'o', 'Filename of the output map (FITS format)',&
                           has_value=.true., default='photproject.fits')
    call parser%add_option('header', 'h', 'Input FITS header of the map',   &
                           has_value=.true.)
    call parser%add_option('resolution', '', 'Input pixel size of the map', &
                           has_value=.true., default='3.')
    call parser%add_option('npixels-per-sample', 'n', 'Maximum number of sky &
                           &pixels intercepted by a PACS detector',          &
                           has_value=.true., default='6')
    call parser%add_option('no-flatfield','','Do not divide by calibration flat field')
    call parser%add_option('filtering', 'f', 'Timeline filtering (mean|none)', &
                           has_value=.true., default='mean')
    call parser%add_option('deglitching', 'd', 'Timeline deglitching &
                           &(l2bl2bc|none)', has_value=.true., default='none')
    call parser%add_option('nsigma-deglitching', '', 'N-sigma for deglitching',&
                           has_value=.true., default='5.')
    call parser%add_option('first', '', 'First sample in timeline', &
                           has_value=.true., default='12001')
    call parser%add_option('last', '', 'Last sample in timeline',   &
                           has_value=.true., default='86000')

    call parser%parse(status)
    if (status == -1) stop
    if (status /=  0) stop 'Aborting.'

    call parser%print_options()

    if (parser%get_argument_count() == 1) then
       infile = parser%get_argument(1,status)
    else
       infile = filename
    end if
    outfile = parser%get_option('o', status) !XXX status gfortran bug
    headerfile = parser%get_option('header', status)
    resolution = parser%get_option_as_real('resolution', status)
    npixels_per_sample = parser%get_option_as_integer('npixels-per-sample', status)
    first = parser%get_option_as_integer('first', status)
    last  = parser%get_option_as_integer('last', status)
    do_flatfield = .not. parser%get_option_as_logical('no-flatfield', status)
    do_meansubtraction = parser%get_option('filtering', status) == 'mean'
    do_deglitching_l2b = parser%get_option('deglitching', status) == 'l2b'
    do_deglitching_l2c = parser%get_option('deglitching', status) == 'l2c'
    nsigma_deglitching = parser%get_option_as_real('nsigma-deglitching', status)

    ! read pointing information
    allocate(pointing)
    call pointing%load_filename(filename, status)
    if (status /= 0) stop 'FAILED: pointing%load_filename'

    ! read the time file
    status = 0

    nsamples = last - first + 1
    allocate(time(last-first+1))
    allocate(timeus(last-first+1))
    call ft_readslice(filename // '_Time.fits+1', first, last, timeus, status)
    if (status /= 0) stop 'FAILED: ft_readslice'
    time = timeus * 1.0d-6

    ! get the pacs instance, read the calibration files
    allocate(pacs)
    call pacs%read_calibration_files(status)
    if (status /= 0) stop 'FAILED: read_calibration_files'

    call pacs%filter_detectors('blue', transparent_mode=.true., status=status)
    if (status /= 0) stop 'FAILED: filter_detectors'

    call pacs%compute_mapheader(pointing, time, resolution, header, status)
    if (status /= 0) stop 'FAILED: compute_mapheader'

    call ft_readparam(header, 'naxis1', count, nx, status=status)
    if (status /= 0 .or. count == 0) stop 'FAILED: compute_mapheader 2'

    call ft_readparam(header, 'naxis2', count, ny, status=status)
    if (status /= 0 .or. count == 0) stop 'FAILED: compute_mapheader 3'

    ! allocate memory for the maps
    allocate(map1d(0:nx*ny-1))
    allocate(weights(0:nx*ny-1))

    ! compute the projector
    write(*,'(a)', advance='no') 'Computing the projector... '
    allocate(pmatrix(npixels_per_sample,last-first+1,pacs%ndetectors))
    call system_clock(count1, count_rate, count_max)

    call pacs%compute_projection_sharp_edges(pointing, time, header, nx, ny, pmatrix, status)
    if (status /= 0) stop 'FAILED: compute_projection_sharp_edges.'

    call system_clock(count2, count_rate, count_max)
    write(*,'(f6.2,a)') real(count2-count1)/count_rate, 's'

    ! read the signal file
    allocate(signal(last-first+1, pacs%ndetectors))
    allocate(mask(last-first+1, pacs%ndetectors))

    write(*,'(a)', advance='no') 'Reading timeline... '
    call system_clock(count1, count_rate, count_max)
    call pacs%read_signal_file(filename // '_Signal.fits', first, last, signal, status)
    if (status /= 0) stop 'FAILED: read_signal_file.'
    call pacs%read_mask_file(filename // '_Mask.fits', first, last, mask, status)
    if (status /= 0) stop 'FAILED: read_mask_file.'
    call system_clock(count2, count_rate, count_max)
    write(*,'(f6.2,a)') real(count2-count1)/count_rate, 's'

    ! remove flat field
    if (do_flatfield) then
        write(*,'(a)') 'Flat-fielding... '
        call divide_vectordim2(signal, pacs%flatfield)
    end if

    ! remove mean value in timeline
    if (do_meansubtraction) then
        write(*,'(a)') 'Removing mean value... '
        call subtract_meandim1(signal)
    end if

    ! remove mean value in timeline
    if (do_deglitching_l2b) then
        call deglitch_l2b(pmatrix, nx, ny, signal, mask, nsigma_deglitching, .false.)
    end if
    if (do_deglitching_l2c) then
        call deglitch_l2b(pmatrix, nx, ny, signal, mask, nsigma_deglitching, .true.)
    end if

    ! back project the timeline
    write(*,'(a)', advance='no') 'Computing the map... '
    call system_clock(count1, count_rate, count_max)
    call backprojection_weighted(pmatrix, signal, mask, map1d)
    call system_clock(count2, count_rate, count_max)
    write(*,'(f6.2,a)') real(count2-count1)/count_rate, 's'


    ! write the map as fits file
    write(*,'(a)') 'Writing FITS file... '
    call ft_write(outfile, reshape(map1d, [nx,ny]), header, status)
    if (status /= 0) stop 'FAILED: ft_write.'


end program pacs_photproject
