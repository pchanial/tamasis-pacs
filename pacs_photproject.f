program pacs_photproject

    use iso_fortran_env,        only : ERROR_UNIT, OUTPUT_UNIT
    use module_fitstools,       only : ft_header2str, ft_read_parameter, ft_write
    use module_deglitching,     only : deglitch_l2b
    use module_optionparser,    only : optionparser
    use module_pacsinstrument,  only : pacsinstrument
    use module_pacsobservation, only : pacsobservation, maskarray
    use module_pointingmatrix,  only : backprojection_weighted, pointingelement
    use module_preprocessor,    only : divide_vectordim2, median_filtering_nocopy, subtract_meandim1
    use module_string,          only : strlowcase
    use module_tamasis,         only : init_tamasis
    implicit none

    class(pacsinstrument), allocatable :: pacs
    class(pacsobservation), allocatable:: obs
    type(maskarray)                    :: maskarray_policy
    type(pointingelement), allocatable :: pmatrix(:,:,:)
    class(optionparser), allocatable   :: parser
    real*8, allocatable                :: signal(:,:)
    logical*1, allocatable             :: mask(:,:)
    real*8, allocatable                :: map1d(:)
    character(len=2048), allocatable   :: infile(:)
    character(len=2048)                :: outfile, headerfile
    character(len=2880)                :: header
    real*8                             :: resolution
    integer                            :: npixels_per_sample
    logical                            :: do_flatfield
    character(len=256)                 :: filtering_method, deglitching_method
    integer*8                          :: nsamples, idetector
    integer                            :: nobs, iobs, nx, ny, start
    integer                            :: status, count1, count2, count_rate, count_max
    integer                            :: filtering_length
    real*8                             :: deglitching_nsigma

    ! command parsing
    allocate(parser)
    call parser%init('pacs_photproject [options] fitsfile...', 1, -1)
    call parser%add_option('', 'o', 'Filename of the output map (FITS format)', has_value=.true., default='photproject.fits')
    call parser%add_option('header', '', 'Input FITS header of the map',        has_value=.true.)
    call parser%add_option('resolution', '', 'Input pixel size of the map',     has_value=.true., default='3.2')
    call parser%add_option('npixels-per-sample', 'n', 'Maximum number of sky pixels intercepted by a PACS detector',               &
                           has_value=.true.)
    call parser%add_option('no-flatfield','','Do not divide by calibration flat field')
    call parser%add_option('filtering', 'f', 'Timeline filtering (median|none)',  has_value=.true., default='none')
    call parser%add_option('length', 'l', 'Filtering length',  has_value=.true., default='200')
    call parser%add_option('deglitching', 'd', 'Timeline deglitching (l2std|l2mad|none)', has_value=.true., default='none')
    call parser%add_option('nsigma', '', 'N-sigma for deglitching',             has_value=.true., default='5.')

    call parser%parse(status)
    if (status == -1) stop
    if (status /=  0) go to 999

    outfile = parser%get_option('o') !XXX status gfortran bug
    headerfile = parser%get_option('header')
    resolution = parser%get_option_as_real('resolution')
    npixels_per_sample = parser%get_option_as_integer('npixels-per-sample')
    do_flatfield       = .not. parser%get_option_as_logical('no-flatfield')
    filtering_method   = strlowcase(parser%get_option('filtering'))
    filtering_length   = parser%get_option_as_integer('length')
    deglitching_method = strlowcase(parser%get_option('deglitching'))
    deglitching_nsigma = parser%get_option_as_real('nsigma')

    if (filtering_method /= 'none' .and. filtering_method /= 'median') then
        write (ERROR_UNIT) "Invalid filtering method '" // filtering_method // "'. The only valid method is 'median'."
        go to 999
    end if
    
    if (deglitching_method /= 'none' .and. deglitching_method /= 'l2std' .and. deglitching_method /= 'l2mad') then
        write (ERROR_UNIT) "Invalid deglitching method '" // deglitching_method // "'. Valid methods are 'l2std' or 'l2mad'."
        go to 999
    end if
    
    write (OUTPUT_UNIT,*)

    ! initialise tamasis
    call init_tamasis

    ! initialise observations
    nobs = parser%get_argument_count()
    allocate(obs)
    call parser%get_arguments(infile, status)
    if (status /= 0) go to 999
    call obs%init(infile, maskarray_policy, status, verbose=.true.)
    if (status /= 0) go to 999
    nsamples = obs%nsamples_tot

    ! print observation info
    call obs%print()

    ! get npixels_per_sample
    if (parser%get_option('npixels-per-sample') == '') then
        if (obs%channel == 'r') then
            npixels_per_sample = 11
        else
            npixels_per_sample = 5
        endif
    else
        npixels_per_sample = parser%get_option_as_integer('npixels-per-sample')
    end if

    ! initialise pacs instrument
    allocate(pacs)
    call pacs%init(obs%channel, obs%observing_mode == 'Transparent', 1, .false., status)
    if (status /= 0) go to 999

    ! get FITS header
    if (headerfile /= '') then
       call ft_header2str(headerfile, header, status)
    else
       call pacs%compute_map_header(obs, .false., resolution, header, status)
    end if
    if (status /= 0) go to 999

    ! allocate memory for the maps
    call ft_read_parameter(header, 'naxis1', nx, count2, status=status)
    if (status /= 0 .or. count2 == 0) go to 999

    call ft_read_parameter(header, 'naxis2', ny, count2, status=status)
    if (status /= 0 .or. count2 == 0) go to 999

    allocate(map1d(0:nx*ny-1))

    ! compute the projector
    allocate (pmatrix(npixels_per_sample,nsamples,pacs%ndetectors))
    call pacs%compute_projection_sharp_edges(obs, .false., header,  nx, ny, pmatrix, status)
    if (status /= 0) go to 999

    ! read the signal file
    allocate (signal(nsamples, pacs%ndetectors))
    allocate (mask  (nsamples, pacs%ndetectors))
    call pacs%read(obs, signal, mask, status, verbose=.true.)
    if (status /= 0) go to 999

    ! remove flat field
    if (do_flatfield) then
        write (OUTPUT_UNIT,'(a)') 'Flat-fielding... '
        call divide_vectordim2(signal, pacs%flatfield)
    end if

    ! filtering
    select case (filtering_method)

        ! remove mean value in timeline
        case ('none')            
            write (OUTPUT_UNIT,'(a)') 'Removing mean value... '
            start = 1
            do iobs = 1, nobs
                call subtract_meandim1(signal(start:start+obs%slice(iobs)%nsamples-1,:))
                start = start + obs%slice(iobs)%nsamples
            end do


        case ('median')
            write (OUTPUT_UNIT,'(a,i0,a)', advance='no') 'Median filtering (length=', filtering_length, ')... '
            call system_clock(count1, count_rate, count_max)
            start = 1
            do iobs = 1, nobs

                !$omp parallel do
                do idetector = 1, pacs%ndetectors
                    call median_filtering_nocopy(signal(start:start+obs%slice(iobs)%nsamples-1,idetector), filtering_length)
                end do
                !$omp end parallel do
                start = start + obs%slice(iobs)%nsamples
                
            end do
            call system_clock(count2, count_rate, count_max)
            write (OUTPUT_UNIT,'(f7.2,a)') real(count2-count1)/count_rate, 's'

    end select

    ! deglitching
    select case (deglitching_method)

        case ('l2std')
            call deglitch_l2b(pmatrix, nx, ny, signal, mask, deglitching_nsigma, .false., verbose=.true.)

        case ('l2mad')
            call deglitch_l2b(pmatrix, nx, ny, signal, mask, deglitching_nsigma, .true., verbose=.true.)

    end select

    ! back project the timeline
    write(OUTPUT_UNIT,'(a)', advance='no') 'Computing the map... '
    call system_clock(count1, count_rate, count_max)
    call backprojection_weighted(pmatrix, signal, mask, map1d, 0.d0)
    call system_clock(count2, count_rate, count_max)
    write(OUTPUT_UNIT,'(f6.2,a)') real(count2-count1)/count_rate, 's'

    ! write the map as fits file
    write(OUTPUT_UNIT,'(a)') 'Writing FITS file... '
    call ft_write(outfile, reshape(map1d, [nx,ny]), header, status)
    if (status /= 0) go to 999

    call exit(0)

999 stop 'Aborting.'

end program pacs_photproject
