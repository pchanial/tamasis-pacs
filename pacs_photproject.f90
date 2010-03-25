program pacs_photproject

    use, intrinsic :: ISO_FORTRAN_ENV
    use            :: module_fitstools
    use            :: module_deglitching, only : deglitch_l2b
    use            :: module_optionparser, only : optionparser
    use            :: module_pacsinstrument
    use            :: module_pacspointing
    use            :: module_pointingmatrix, only : pointingelement,           &
                          backprojection_weighted
    use            :: module_preprocessor
    implicit none

    character(len=*), parameter :: defaultfile = '/home/pchanial/work/pacs/data/transparent/NGC6946/1342184520_blue'
    type(pacsinstrument)               :: pacs
    type(pacspointing)                 :: pointing

    character(len=2048)                :: infile, outfile, headerfile
    character(len=2880)                :: header
    real*8                             :: resolution
    integer                            :: npixels_per_sample
    logical                            :: do_flatfield, do_meansubtraction
    character(len=256)                 :: deglitching_method
    integer*8                          :: first, last

    real*8, allocatable                :: signal(:,:)
    logical*1, allocatable             :: mask(:,:)
    integer                            :: nx, ny
    integer                            :: status, count, count1
    integer                            :: count2, count_rate, count_max
    real*8                             :: deglitching_nsigma
    real*8, allocatable                :: map1d(:)
    type(pointingelement), allocatable :: pmatrix(:,:,:)
    type(optionparser)                 :: parser

    ! command parsing
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
    call parser%add_option('deglitching', 'd', 'Timeline deglitching (l2std|&
                           &l2mad|none)', has_value=.true., default='none')
    call parser%add_option('nsigma', '', 'N-sigma for deglitching',&
                           has_value=.true., default='5.')
    call parser%add_option('first', '', 'First sample in timeline', &
                           has_value=.true., default='12001')
    call parser%add_option('last', '', 'Last sample in timeline',   &
                           has_value=.true., default='86000')

    call parser%parse(status)
    if (status == -1) stop
    if (status /=  0) go to 999

    call parser%print_options()

    if (parser%get_argument_count() == 1) then
       infile = parser%get_argument(1,status)
    else
       infile = defaultfile
    end if
    outfile = parser%get_option('o', status) !XXX status gfortran bug
    headerfile = parser%get_option('header', status)
    resolution = parser%get_option_as_real('resolution', status)
    npixels_per_sample = parser%get_option_as_integer('npixels-per-sample', status)
    first = parser%get_option_as_integer('first', status)
    last  = parser%get_option_as_integer('last', status)
    do_flatfield = .not. parser%get_option_as_logical('no-flatfield', status)
    do_meansubtraction = parser%get_option('filtering', status) == 'mean'
    deglitching_method = parser%get_option('deglitching', status)
    deglitching_nsigma = parser%get_option_as_real('nsigma', status)

    ! read pointing information
    call pointing%load_filename(trim(infile), status)
    if (status /= 0) go to 999

    ! get the pacs instance, read the calibration files
    call pacs%init_filename(trim(infile), .false., status)
    if (status /= 0) go to 999

    ! get FITS header
    if (headerfile /= '') then
       call ft_header2str(headerfile, header, status)
    else
       call pacs%compute_mapheader(pointing, pointing%time(first:last),        &
                                   resolution, header, status)
    end if
    if (status /= 0) go to 999

    ! allocate memory for the maps
    call ft_readparam(header, 'naxis1', count, nx, status=status)
    if (status /= 0 .or. count == 0) go to 999

    call ft_readparam(header, 'naxis2', count, ny, status=status)
    if (status /= 0 .or. count == 0) go to 999

    allocate(map1d(0:nx*ny-1))

    ! compute the projector
    write(*,'(a)', advance='no') 'Computing the projector... '
    allocate(pmatrix(npixels_per_sample,last-first+1,pacs%ndetectors))

    call system_clock(count1, count_rate, count_max)
    call pacs%compute_projection_sharp_edges(pointing,                         &
                                             pointing%time(first:last), header,&
                                             nx, ny, pmatrix, status)
    if (status /= 0) go to 999

    call system_clock(count2, count_rate, count_max)
    write(*,'(f6.2,a)') real(count2-count1)/count_rate, 's'

    ! read the signal file
    allocate(signal(last-first+1, pacs%ndetectors))
    allocate(mask  (last-first+1, pacs%ndetectors))

    write(*,'(a)', advance='no') 'Reading timeline... '
    call system_clock(count1, count_rate, count_max)
    call pacs%read_tod_file(trim(infile), first, last, signal, mask, status)
    if (status /= 0) go to 999
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
    if (len_trim(deglitching_method) /= 0) then
        write (*,'(a,$)') 'Deglitching (' // trim(deglitching_method) // ')...'
        call system_clock(count1, count_rate, count_max)
        select case (deglitching_method)
            case ('l2std')
               call deglitch_l2b(pmatrix, nx, ny, signal, mask,                &
                                 deglitching_nsigma, .false.)
            case ('l2mad')
               call deglitch_l2b(pmatrix, nx, ny, signal, mask,                &
                                 deglitching_nsigma, .true.)
            case default
               write (ERROR_UNIT,'(/,a)') "ERROR: Invalid deglitching method '"&
                                          // trim(deglitching_method) // "'."
               go to 999
         end select
         call system_clock(count2, count_rate, count_max)
         write(*,'(f6.2,a)') real(count2-count1)/count_rate, 's'
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
    if (status /= 0) go to 999

    call exit(0)

999 stop 'Aborting.'

end program pacs_photproject
