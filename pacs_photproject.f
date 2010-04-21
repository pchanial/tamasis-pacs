program pacs_photproject

    use iso_fortran_env,        only : ERROR_UNIT, OUTPUT_UNIT
    use module_fitstools,       only : ft_header2str, ft_read_parameter, ft_write
    use module_deglitching,     only : deglitch_l2b
    use module_optionparser,    only : optionparser
    use module_pacsinstrument,  only : pacsinstrument
    use module_pacsobservation, only : pacsobservation, maskarray
    use module_pointingmatrix,  only : backprojection_weighted, pointingelement
    use module_preprocessor,    only : divide_vectordim2, subtract_meandim1
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
    logical                            :: do_flatfield, do_meansubtraction
    character(len=256)                 :: deglitching_method
    integer*8                          :: nsamples
    integer                            :: nobs, iobs, nx, ny, dest
    integer                            :: status, count1, count2, count_rate, count_max
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
    call parser%add_option('filtering', 'f', 'Timeline filtering (mean|none)',  has_value=.true., default='mean')
    call parser%add_option('deglitching', 'd', 'Timeline deglitching (l2std|l2mad|none)', has_value=.true., default='none')
    call parser%add_option('nsigma', '', 'N-sigma for deglitching',             has_value=.true., default='5.')

    call parser%parse(status)
    if (status == -1) stop
    if (status /=  0) go to 999

    outfile = parser%get_option('o') !XXX status gfortran bug
    headerfile = parser%get_option('header')
    resolution = parser%get_option_as_real('resolution')
    npixels_per_sample = parser%get_option_as_integer('npixels-per-sample')
    do_flatfield = .not. parser%get_option_as_logical('no-flatfield')
    do_meansubtraction = parser%get_option('filtering') == 'mean'
    deglitching_method = parser%get_option('deglitching')
    deglitching_nsigma = parser%get_option_as_real('nsigma')
    
    write (OUTPUT_UNIT,*)

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
    allocate(pmatrix(npixels_per_sample,nsamples,pacs%ndetectors))
    call pacs%compute_projection_sharp_edges(obs, .false., header,  nx, ny, pmatrix, status)
    if (status /= 0) go to 999

    ! read the signal file
    allocate(signal(nsamples, pacs%ndetectors))
    allocate(mask  (nsamples, pacs%ndetectors))
    call pacs%read(obs, signal, mask, status)
    if (status /= 0) go to 999

    ! remove flat field
    if (do_flatfield) then
        write(OUTPUT_UNIT,'(a)') 'Flat-fielding... '
        call divide_vectordim2(signal, pacs%flatfield)
    end if

    ! remove mean value in timeline
    if (do_meansubtraction) then
        write(OUTPUT_UNIT,'(a)') 'Removing mean value... '
        dest = 1
        do iobs = 1, nobs
            call subtract_meandim1(signal(dest:dest+obs%slice(iobs)%nsamples-1,:))
            dest = dest + obs%slice(iobs)%nsamples
        end do
    end if

    ! deglitching
    if (deglitching_method /= 'none') then
        select case (deglitching_method)

            case ('l2std')
               call deglitch_l2b(pmatrix, nx, ny, signal, mask, deglitching_nsigma, .false.)
            case ('l2mad')
               call deglitch_l2b(pmatrix, nx, ny, signal, mask, deglitching_nsigma, .true.)
            case default
               write (ERROR_UNIT,'(/,a)') "ERROR: Invalid deglitching method '" // trim(deglitching_method) // "'."
               go to 999

         end select
    end if

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
