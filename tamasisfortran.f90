! Tamasis interface for f2py
!
! Author: P. Chanial

subroutine pacs_info_channel(filename, nfilenames, channel, status)

    use module_pacsobservation, only : pacsobservation, maskarray
    implicit none

    !f2py threadsafe
    !f2py intent(in)   filename
    !f2py intent(in)   nfilenames
    !f2py intent(out)  channel
    !f2py intent(out)  status

    character(len=*), intent(in) :: filename
    integer, intent(in)          :: nfilenames
    character, intent(out)       :: channel
    integer, intent(out)         :: status

    character(len=len(filename)/nfilenames), allocatable :: filename_(:)
    class(pacsobservation), allocatable :: obs
    type(maskarray)                     :: maskarray_policy
    integer                             :: iobs

    ! split input filename
    if (mod(len(filename), nfilenames) /= 0) then
        stop 'PACS_INFO: Invalid filename length.'
    end if

    allocate(filename_(nfilenames))

    do iobs = 1, nfilenames
        filename_(iobs) = filename((iobs-1)*len(filename_)+1:iobs*len(filename_))
    end do

    allocate(obs)
    call obs%init(filename_, maskarray_policy, status)
    if (status /= 0) return

    channel = obs%channel

end subroutine pacs_info_channel


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_info(tamasis_dir, filename, nfilenames, fine_sampling_factor, keep_bad_detectors, use_bad_detector_mask,           &
                     bad_detector_mask, nrows, ncolumns, ndetectors, output_mask, transparent_mode, compression_factor, nsamples,  &
                     status)

    use module_pacsinstrument,  only : pacsinstrument
    use module_pacsobservation, only : pacsobservation, maskarray
    use module_tamasis,         only : init_tamasis
    implicit none

    !f2py threadsafe
    !f2py intent(in)   tamasis_dir
    !f2py intent(in)   filename
    !f2py intent(in)   nfilenames
    !f2py intent(in)   fine_sampling_factor
    !f2py intent(in)   keep_bad_detectors
    !f2py intent(in)   use_bad_detector_mask
    !f2py intent(in)   bad_detector_mask
    !f2py intent(hide) nrows = shape(bad_detector_mask,0)
    !f2py intent(hide) ncolumns = shape(bad_detector_mask,1)
    !f2py intent(out)  ndetectors
    !f2py intent(out)  output_mask
    !f2py intent(out)  transparent_mode
    !f2py intent(out)  compression_factor
    !f2py intent(out)  nsamples
    !f2py intent(out)  status

    character(len=*), intent(in) :: tamasis_dir
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: nfilenames
    integer, intent(in)          :: fine_sampling_factor
    logical, intent(in)          :: keep_bad_detectors
    logical, intent(in)          :: use_bad_detector_mask
    logical*1, intent(in)        :: bad_detector_mask(nrows,ncolumns)
    integer, intent(in)          :: nrows, ncolumns
    integer, intent(out)         :: ndetectors
    logical*1, intent(out)       :: output_mask(nrows,ncolumns)
    logical, intent(out)         :: transparent_mode
    integer, intent(out)         :: compression_factor(nfilenames)
    integer, intent(out)         :: nsamples(nfilenames)
    integer, intent(out)         :: status

    character(len=len(filename)/nfilenames), allocatable :: filename_(:)
    class(pacsobservation), allocatable :: obs
    class(pacsinstrument), allocatable  :: pacs
    type(maskarray)                     :: maskarray_policy
    integer                             :: iobs

    ! initialise tamasis
    call init_tamasis(tamasis_dir)

    ! split input filename
    if (mod(len(filename), nfilenames) /= 0) then
        stop 'PACS_INFO: Invalid filename length.'
    end if

    allocate(filename_(nfilenames))

    do iobs = 1, nfilenames
        filename_(iobs) = filename((iobs-1)*len(filename_)+1:iobs*len(filename_))
    end do

    allocate(obs)
    call obs%init(filename_, maskarray_policy, status, verbose=.true.)
    if (status /= 0) return

    ! print some information about the observation
    call obs%print()

    transparent_mode = obs%slice(1)%observing_mode == 'Transparent'
    allocate(pacs)
    if (use_bad_detector_mask) then
        call pacs%init(obs%channel, transparent_mode, fine_sampling_factor, keep_bad_detectors, status, bad_detector_mask)
    else
        call pacs%init(obs%channel, transparent_mode, fine_sampling_factor, keep_bad_detectors, status)
    end if
    if (status /= 0) return

    ndetectors  = pacs%ndetectors
    output_mask = pacs%mask
    compression_factor = obs%slice%compression_factor
    nsamples = obs%slice%nsamples

end subroutine pacs_info


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_map_header(tamasis_dir, filename, nfilenames, finer_sampling, fine_sampling_factor, keep_bad_detectors,            &
                           use_bad_detector_mask, bad_detector_mask, nrows, ncolumns, resolution, header, status)

    use module_pacsinstrument,  only : pacsinstrument
    use module_pacsobservation, only : pacsobservation, maskarray
    use module_tamasis,         only : init_tamasis
    implicit none

    !f2py threadsafe
    !f2py intent(in)   tamasis_dir
    !f2py intent(in)   filename
    !f2py intent(in)   nfilenames
    !f2py intent(in)   finer_sampling
    !f2py intent(in)   fine_sampling_factor
    !f2py intent(in)   keep_bad_detectors
    !f2py intent(in)   use_bad_detector_mask
    !f2py intent(in)   bad_detector_mask
    !f2py intent(hide) nrows = shape(bad_detector_mask,0)
    !f2py intent(hide) ncolumns = shape(bad_detector_mask,1)
    !f2py intent(in)   resolution
    !f2py intent(out)  header
    !f2py intent(out)  status

    character(len=*), intent(in) :: tamasis_dir
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: nfilenames
    logical, intent(in)          :: finer_sampling
    integer, intent(in)          :: fine_sampling_factor
    logical, intent(in)          :: keep_bad_detectors
    logical, intent(in)          :: use_bad_detector_mask
    logical*1, intent(in)        :: bad_detector_mask(nrows,ncolumns)
    integer, intent(in)          :: nrows, ncolumns
    real*8, intent(in)           :: resolution
    character(len=2880), intent(out) :: header
    integer, intent(out)             :: status

    character(len=len(filename)/nfilenames), allocatable :: filename_(:)
    class(pacsobservation), allocatable :: obs
    class(pacsinstrument), allocatable  :: pacs
    type(maskarray)                     :: maskarray_policy
    integer                             :: iobs
    
    ! initialise tamasis
    call init_tamasis(tamasis_dir)

    ! split input filename
    if (mod(len(filename), nfilenames) /= 0) then
        stop 'PACS_INFO: Invalid filename length.'
    end if

    allocate(filename_(nfilenames))

    do iobs = 1, nfilenames
        filename_(iobs) = filename((iobs-1)*len(filename_)+1:iobs*len(filename_))
    end do

    allocate(obs)
    call obs%init(filename_, maskarray_policy, status)
    if (status /= 0) go to 999

    allocate(pacs)
    if (use_bad_detector_mask) then
        call pacs%init(obs%channel, obs%slice(1)%observing_mode == 'Transparent', fine_sampling_factor, keep_bad_detectors, status,&
             bad_detector_mask)
    else
        call pacs%init(obs%channel, obs%slice(1)%observing_mode == 'Transparent', fine_sampling_factor, keep_bad_detectors,status)
    end if
    if (status /= 0) go to 999
    
    call pacs%compute_map_header(obs, finer_sampling, resolution, header, status)
    
999 continue

end subroutine pacs_map_header


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_timeline(tamasis_dir, filename, nfilenames, nsamples, ndetectors, keep_bad_detectors, bad_detector_mask, nrow,     &
                         ncol, do_flatfielding, do_subtraction_mean, signal, mask, status)

    use iso_fortran_env,        only : ERROR_UNIT
    use module_pacsinstrument,  only : pacsinstrument
    use module_pacsobservation, only : pacsobservation, maskarray
    use module_preprocessor,    only : divide_vectordim2, subtract_meandim1
    use module_tamasis,         only : init_tamasis
    implicit none

    !f2py threadsafe
    !f2py intent(in)   :: tamasis_dir
    !f2py intent(in)   :: filename
    !f2py intent(in)   :: nfilenames
    !f2py intent(in)   :: nsamples
    !f2py intent(in)   :: ndetectors
    !f2py intent(in)   :: keep_bad_detectors
    !f2py intent(in)   :: bad_detector_mask
    !f2py intent(hide) :: nrow = shape(bad_detector_mask,0)
    !f2py intent(hide) :: ncol = shape(bad_detector_mask,1)
    !f2py intent(in)   :: do_flatfielding
    !f2py intent(in)   :: do_subtraction_mean
    !f2py intent(out)  :: signal
    !f2py intent(out)  :: mask
    !f2py intent(out)  :: status

    character(len=*), intent(in) :: tamasis_dir
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: nfilenames
    integer, intent(in)          :: nsamples
    integer, intent(in)          :: ndetectors
    logical, intent(in)          :: keep_bad_detectors
    logical*1, intent(in)        :: bad_detector_mask(nrow,ncol)
    integer, intent(in)          :: nrow, ncol
    logical, intent(in)          :: do_flatfielding, do_subtraction_mean
    real*8, intent(out)          :: signal(nsamples, ndetectors)
    logical*1, intent(out)       :: mask(nsamples, ndetectors)
    integer, intent(out)         :: status

    character(len=len(filename)/nfilenames), allocatable :: filename_(:)
    class(pacsobservation), allocatable :: obs
    class(pacsinstrument), allocatable  :: pacs
    type(maskarray)                     :: maskarray_policy
    integer                             :: iobs, destination

    ! initialise tamasis
    call init_tamasis(tamasis_dir)

    ! split input filename
    if (mod(len(filename), nfilenames) /= 0) then
        stop 'PACS_INFO: Invalid filename length.'
    end if

    allocate(filename_(nfilenames))

    do iobs = 1, nfilenames
        filename_(iobs) = filename((iobs-1)*len(filename_)+1:iobs*len(filename_))
    end do

    ! initialise observations
    allocate(obs)
    call obs%init(filename_, maskarray_policy, status)
    if (status /= 0) go to 999

    ! initialise pacs instrument
    allocate(pacs)
    call pacs%init(obs%channel, obs%observing_mode == 'Transparent', 1, keep_bad_detectors, status, bad_detector_mask)
    if (status /= 0) go to 999

    ! read timeline
    call pacs%read(obs, signal, mask, status, verbose=.true.)
    if (status /= 0) go to 999

    ! flat fielding
    if (do_flatfielding) then
        call divide_vectordim2(signal, pacs%flatfield)
    end if
    
    ! subtract the mean of each detector timeline
    if (do_subtraction_mean) then
        destination = 1
        do iobs=1, nfilenames
            call subtract_meandim1(signal(destination:destination+obs%slice(iobs)%nsamples-1,:))
            destination = destination + obs%slice(iobs)%nsamples
         end do
    end if

999 continue

end subroutine pacs_timeline


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_pointing_matrix_filename(tamasis_dir, filename, nfilenames, finer_sampling, fine_sampling_factor,                  &
                                         npixels_per_sample, nsamples, ndetectors, keep_bad_detectors, bad_detector_mask, nrow,    &
                                         ncol, header, pmatrix, status)

    use iso_fortran_env,        only : ERROR_UNIT
    use module_fitstools,       only : ft_read_parameter
    use module_pacsinstrument,  only : pacsinstrument
    use module_pacsobservation, only : pacsobservation, maskarray
    use module_pointingmatrix,  only : pointingelement
    use module_tamasis,         only : init_tamasis
    implicit none

    !f2py threadsafe
    !f2py intent(in)   :: tamasis_dir
    !f2py intent(in)   :: filename
    !f2py intent(in)   :: nfilenames
    !f2py intent(in)   :: finer_sampling
    !f2py intent(in)   :: fine_sampling_factor
    !f2py intent(in)   :: npixels_per_sample
    !f2py intent(in)   :: nsamples
    !f2py intent(in)   :: ndetectors
    !f2py intent(in)   :: keep_bad_detectors
    !f2py intent(in)   :: bad_detector_mask
    !f2py intent(hide) :: nrow = shape(bad_detector_mask,0)
    !f2py intent(hide) :: ncol = shape(bad_detector_mask,1)
    !f2py intent(in)   :: header
    !f2py integer*8, intent(inout) :: pmatrix(npixels_per_sample*nsamples*ndetectors)
    !f2py intent(out)  :: status

    character(len=*), intent(in) :: tamasis_dir
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: nfilenames
    logical, intent(in)          :: finer_sampling
    integer, intent(in)          :: fine_sampling_factor
    integer, intent(in)          :: npixels_per_sample
    integer, intent(in)          :: nsamples
    integer, intent(in)          :: ndetectors
    logical, intent(in)          :: keep_bad_detectors
    logical*1, intent(in)        :: bad_detector_mask(nrow,ncol)
    integer, intent(in)          :: nrow, ncol
    character(len=*), intent(in) :: header
    type(pointingelement), intent(inout) :: pmatrix(npixels_per_sample,nsamples,ndetectors)
    integer, intent(out)         :: status

    character(len=len(filename)/nfilenames), allocatable :: filename_(:)
    class(pacsobservation), allocatable :: obs
    class(pacsinstrument), allocatable  :: pacs
    type(maskarray)                     :: maskarray_policy
    integer                             :: iobs, nx, ny, count, nsamples_expected

    ! initialise tamasis
    call init_tamasis(tamasis_dir)

    ! split input filename
    if (mod(len(filename), nfilenames) /= 0) then
        stop 'PACS_INFO: Invalid filename length.'
    end if

    allocate(filename_(nfilenames))

    do iobs = 1, nfilenames
        filename_(iobs) = filename((iobs-1)*len(filename_)+1:iobs*len(filename_))
    end do

    ! initialise observations
    allocate(obs)
    call obs%init(filename_, maskarray_policy, status)
    if (status /= 0) return

    ! initialise pacs instrument
    allocate(pacs)
    call pacs%init(obs%channel, obs%slice(1)%observing_mode == 'Transparent', fine_sampling_factor, keep_bad_detectors, status,    &
                   bad_detector_mask)
    if (status /= 0) return

    ! check number of detectors
    if (pacs%ndetectors /= ndetectors) then
        status = 1
        write (ERROR_UNIT,'(a,2(i0,a))') "The specified number of detectors '", ndetectors, "' is incompatible with that from the i&
              &nput bad detector mask '", pacs%ndetectors, "'."
        return
    end if

    ! check number of fine samples
    if (finer_sampling) then
        nsamples_expected = sum(obs%slice%nsamples * obs%slice%compression_factor)*fine_sampling_factor
    else
        nsamples_expected = sum(obs%slice%nsamples)
    end if
    if (nsamples /= nsamples_expected) then
        status = 1
        write (ERROR_UNIT,'(a,2(i0,a))') "The specified total number of samples '", nsamples, "' is incompatible with that from the&
              & observations '", nsamples_expected, "'."
        return
    end if

    ! get the size of the map
    call ft_read_parameter(header, 'naxis1', nx, count, must_exist=.true., status=status)
    if (status /= 0) return
    call ft_read_parameter(header, 'naxis2', ny, count, must_exist=.true., status=status)
    if (status /= 0) return

    ! compute the projector
    call pacs%compute_projection_sharp_edges(obs, finer_sampling, header, nx, ny, pmatrix, status)

end subroutine pacs_pointing_matrix_filename


!-----------------------------------------------------------------------------------------------------------------------------------!!$
!!$
!!$subroutine pacs_map_header(array, time, ra, dec, pa, chop, npointings,         &
!!$                           finetime, nfinesamples, transparent_mode,           &
!!$                           keep_bad_detectors, bad_detector_mask, nrow, ncol,  &
!!$                           resolution, header)
!!$
!!$    use, intrinsic :: ISO_FORTRAN_ENV, only : ERROR_UNIT
!!$    use module_pacsinstrument, only : pacsinstrument
!!$    use module_pacspointing, only : pacspointing
!!$    implicit none
!!$
!!$    !f2py threadsafe
!!$    !f2py intent(in)   :: array
!!$    !f2py intent(in)   :: time, ra, dec, pa, chop
!!$    !f2py intent(hide) :: npointings = size(time)
!!$    !f2py intent(in)   :: finetime
!!$    !f2py intent(hide) :: nfinesamples = size(finetime)
!!$    !f2py intent(in)   :: transparent_mode
!!$    !f2py intent(in)   :: keep_bad_detectors
!!$    !f2py intent(in)   :: bad_detector_mask
!!$    !f2py intent(hide), depend(bad_detector_mask) :: nrow = shape(bad_detector_mask,0), ncol = shape(bad_detector_mask,1)
!!$    !f2py intent(in)   :: resolution
!!$    !f2py intent(out)  :: header
!!$
!!$    character(len=*), intent(in)     :: array
!!$    real*8, intent(in)               :: time(npointings), ra(npointings), dec(npointings), pa(npointings), chop(npointings)
!!$    integer, intent(in)              :: npointings
!!$    real*8, intent(in)               :: finetime(nfinesamples)
!!$    integer, intent(in)              :: nfinesamples
!!$    logical, intent(in)              :: transparent_mode
!!$    logical, intent(in)              :: keep_bad_detectors
!!$    logical*1, intent(in)            :: bad_detector_mask(nrow,ncol)
!!$    integer, intent(in)              :: nrow, ncol
!!$    real*8, intent(in)               :: resolution
!!$    character(len=2880), intent(out) :: header
!!$    integer, intent(out)             :: status
!!$
!!$    class(pacsinstrument), allocatable :: pacs
!!$    class(pacspointing), allocatable   :: pointing
!!$
!!$    ! read pointing information
!!$    allocate(pointing)
!!$    call pointing%init_sim(time, ra, dec, pa, chop, status)
!!$    if (status /= 0) go to 999
!!$
!!$    ! get the pacs instance, read the calibration files
!!$    allocate(pacs)
!!$    call pacs%init(array(1:1), transparent_mode, keep_bad_detectors, status,   &
!!$                   bad_detector_mask)
!!$    if (status /= 0) go to 999
!!$
!!$    if (any(pacs%mask .neqv. bad_detector_mask)) then
!!$        write (*,'(a)') "Info: using user's bad detector mask."
!!$    end if
!!$
!!$    call pacs%compute_map_header(pointing, .true., resolution, header, status)
!!$    if (status == 0) return
!!$
!!$999 deallocate(obs)
!!$    if (allocated(pointing)) deallocate(pointing)
!!$    if (allocate(pacs))      deallocate(pacs)
!!$
!!$end subroutine pacs_map_header
!!$
!!$
!-----------------------------------------------------------------------------------------------------------------------------------!!$
!!$
!!$subroutine pacs_pointing_matrix(array, time, ra, dec, pa, chop, npointings,    &
!!$                                finetime, nfinesamples, npixels_per_sample,    &
!!$                                ndetectors, transparent_mode,                  &
!!$                                keep_bad_detectors, bad_detector_mask,         &
!!$                                nrow, ncol, header, pmatrix)
!!$
!!$    use, intrinsic :: ISO_FORTRAN_ENV
!!$    use module_fitstools, only : ft_read_parameter
!!$    use module_pacsinstrument, only : pacsinstrument
!!$    use module_pacspointing, only : pacspointing
!!$    use module_pointingmatrix, only : pointingelement
!!$    implicit none
!!$
!!$    !f2py threadsafe
!!$    !f2py intent(in)   :: array
!!$    !f2py intent(in)   :: time, ra, dec, pa, chop
!!$    !f2py intent(hide) :: npointings = size(time)
!!$    !f2py intent(in)   :: finetime
!!$    !f2py intent(hide) :: nfinesamples = size(finetime)
!!$    !f2py intent(in)   :: npixels_per_sample
!!$    !f2py intent(in)   :: ndetectors
!!$    !f2py intent(in)   :: transparent_mode
!!$    !f2py intent(in)   :: keep_bad_detectors
!!$    !f2py intent(in)   :: bad_detector_mask
!!$    !f2py intent(hide), depend(bad_detector_mask) :: nrow = shape(bad_detector_mask,0), ncol = shape(bad_detector_mask,1)
!!$    !f2py intent(in)   :: header
!!$    !f2py integer*8, dimension(npixels_per_sample*nfinesamples*ndetectors), intent(inout) :: pmatrix
!!$
!!$    character(len=*), intent(in)         :: array
!!$    real*8, intent(in)                   :: time(npointings), ra(npointings), dec(npointings), pa(npointings), chop(npointings)
!!$    integer, intent(in)                  :: npointings
!!$    real*8, intent(in)                   :: finetime(nfinesamples)
!!$    integer, intent(in)                  :: nfinesamples
!!$    integer, intent(in)                  :: npixels_per_sample
!!$    integer, intent(in)                  :: ndetectors
!!$    logical, intent(in)                  :: transparent_mode
!!$    logical, intent(in)                  :: keep_bad_detectors
!!$    logical*1, intent(in)                :: bad_detector_mask(nrow,ncol)
!!$    integer, intent(in)                  :: nrow, ncol
!!$    character(len=*), intent(in)         :: header
!!$    type(pointingelement), intent(inout) :: pmatrix(npixels_per_sample, nfinesamples, ndetectors)
!!$
!!$    class(pacsinstrument), allocatable   :: pacs
!!$    class(pacspointing), allocatable     :: pointing
!!$    integer                              :: status, count,  count1, count2, count_rate, count_max, nx, ny
!!$
!!$    ! read pointing information
!!$    allocate(pointing)
!!$    call pointing%load_array(time, ra, dec, pa, chop, status)
!!$    if (status /= 0) go to 999
!!$
!!$    ! get the pacs instance, read the calibration files
!!$    allocate(pacs)
!!$    call pacs%init(array(1:1), transparent_mode, keep_bad_detectors, status,   &
!!$                   bad_detector_mask)
!!$    if (status /= 0) go to 999
!!$
!!$    if (any(pacs%mask .neqv. bad_detector_mask)) then
!!$        write (*,'(a)') "Info: using user's bad pixel mask."
!!$    end if
!!$
!!$    call ft_read_parameter(header, 'naxis1', nx, count, must_exist=.true., status=status)
!!$    if (status /= 0) return
!!$    call ft_read_parameter(header, 'naxis2', ny, count, must_exist=.true., status=status)
!!$    if (status /= 0) return
!!$
!!$    ! compute the projector
!!$    write(*,'(a)', advance='no') 'Info: computing the projector... '
!!$    call system_clock(count1, count_rate, count_max)
!!$    call pacs%compute_projection_sharp_edges(pointing, finetime, header, nx, ny, pmatrix, status)
!!$    call system_clock(count2, count_rate, count_max)
!!$    write(*,'(f6.2,a)') real(count2-count1)/count_rate, 's'
!!$
!!$    return
!!$
!!$999 write(ERROR_UNIT,'(a)') 'Aborting.'
!!$
!!$end subroutine pacs_pointing_matrix


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pointing_matrix_direct(pmatrix, map1d, signal, npixels_per_sample, nsamples, ndetectors, npixels)

    use module_pointingmatrix
    implicit none

    !f2py threadsafe
    !f2py integer*8, dimension(npixels_per_sample*nsamples*ndetectors), intent(inout) :: pmatrix
    !f2py intent(in)      :: map1d
    !f2py intent(inout)   :: signal
    !f2py intent(in)      :: npixels_per_sample
    !f2py intent(hide)    :: nsamples = shape(signal,0)
    !f2py intent(hide)    :: ndetectors = shape(signal,1)
    !f2py intent(hide)    :: npixels = size(map1d)

    type(pointingelement), intent(inout) :: pmatrix(npixels_per_sample, nsamples, ndetectors)
    real*8, intent(in)    :: map1d(npixels)
    real*8, intent(inout) :: signal(nsamples, ndetectors)
    integer, intent(in)   :: npixels_per_sample, nsamples, ndetectors, npixels

    call pmatrix_direct(pmatrix, map1d, signal)

end subroutine pointing_matrix_direct


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pointing_matrix_transpose(pmatrix, signal, map1d, npixels_per_sample, nsamples, ndetectors, npixels)

    use module_pointingmatrix
    implicit none

    !f2py threadsafe
    !f2py integer*8, dimension(npixels_per_sample*nsamples*ndetectors), intent(inout) :: pmatrix
    !f2py intent(in)    :: signal
    !f2py intent(inout) :: map1d
    !f2py intent(in)    :: npixels_per_sample
    !f2py intent(hide)  :: nsamples = shape(signal,0)
    !f2py intent(hide)  :: ndetectors = shape(signal,1)
    !f2py intent(hide)  :: npixels = size(map1d)

    type(pointingelement), intent(inout) :: pmatrix(npixels_per_sample, nsamples, ndetectors)
    real*8, intent(in)    :: signal(nsamples, ndetectors)
    real*8, intent(inout) :: map1d(npixels)
    integer, intent(in)   :: npixels_per_sample, nsamples, ndetectors, npixels

    call pmatrix_transpose(pmatrix, signal, map1d)

end subroutine pointing_matrix_transpose


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pointing_matrix_ptp(pmatrix, ptp, npixels_per_sample, nsamples, ndetectors, npixels)

    use module_pointingmatrix
    implicit none

    !f2py threadsafe
    !f2py integer*8, dimension(npixels_per_sample*nsamples*ndetectors), intent(inout) :: pmatrix
    !f2py intent(out) :: ptp
    !f2py intent(in)  :: npixels_per_sample
    !f2py intent(in)  :: nsamples
    !f2py intent(in)  :: ndetectors
    !f2py intent(in)  :: npixels

    type(pointingelement), intent(inout) :: pmatrix(npixels_per_sample, nsamples, ndetectors)
    real*8, intent(out)                  :: ptp(npixels, npixels)
    integer, intent(in)                  :: npixels_per_sample, nsamples, ndetectors, npixels

    call pmatrix_ptp(pmatrix, ptp)

end subroutine pointing_matrix_ptp


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_multiplexing_direct(signal, multiplexed, fine_sampling_factor, ij, nsamples, ndetectors)

    use module_pacsinstrument, only : multiplexing_direct
    implicit none

    !f2py threadsafe
    !f2py intent(in)      :: signal
    !f2py intent(inout)   :: multiplexed
    !f2py intent(in)      :: fine_sampling_factor
    !f2py intent(in)      :: ij
    !f2py intent(hide)    :: nsamples = shape(signal,0)
    !f2py intent(hide)    :: ndetectors = shape(signal,1)

    integer, intent(in)   :: nsamples, ndetectors
    real*8, intent(in)    :: signal(nsamples, ndetectors)
    integer, intent(in)   :: fine_sampling_factor, ij(2, ndetectors)
    real*8, intent(inout) :: multiplexed(nsamples/fine_sampling_factor, ndetectors)

    call multiplexing_direct(signal, multiplexed, fine_sampling_factor, ij)

end subroutine pacs_multiplexing_direct


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_multiplexing_transpose(multiplexed, signal, fine_sampling_factor, ij, nsamples, ndetectors)

    use module_pacsinstrument, only : multiplexing_transpose
    implicit none

    !f2py threadsafe
    !f2py intent(in)      :: multiplexed
    !f2py intent(inout)   :: signal
    !f2py intent(in)      :: fine_sampling_factor
    !f2py intent(in)      :: ij
    !f2py intent(hide)    :: nsamples = shape(signal,0)
    !f2py intent(hide)    :: ndetectors = shape(signal,1)
    
    integer, intent(in)   :: nsamples, ndetectors, fine_sampling_factor
    real*8, intent(in)    :: multiplexed(nsamples/fine_sampling_factor, ndetectors)
    integer, intent(in)   :: ij(2, ndetectors)
    real*8, intent(inout) :: signal(nsamples, ndetectors)

    call multiplexing_transpose(multiplexed, signal, fine_sampling_factor, ij)

end subroutine pacs_multiplexing_transpose


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine compression_average_direct(data, compressed, factor, nsamples, ndetectors)

    use module_compression, only : direct => compression_average_direct
    implicit none

    !f2py threadsafe
    !f2py intent(in)      :: data
    !f2py intent(inout)   :: compressed
    !f2py intent(in)      :: factor
    !f2py intent(hide)    :: nsamples = shape(data,0)/factor
    !f2py intent(hide)    :: ndetectors = shape(data,1)

    integer, intent(in)   :: factor, nsamples, ndetectors
    real*8, intent(in)    :: data(nsamples*factor,ndetectors)
    real*8, intent(out)   :: compressed(nsamples,ndetectors)

    call direct(data, compressed, factor)

end subroutine compression_average_direct


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine compression_average_transpose(compressed, data, factor, nsamples, ndetectors)

    use module_compression, only : transpos => compression_average_transpose
    implicit none

    !f2py threadsafe
    !f2py intent(inout)   :: compressed
    !f2py intent(in)      :: data
    !f2py intent(in)      :: factor
    !f2py intent(hide)    :: nsamples = shape(data,0)/factor
    !f2py intent(hide)    :: ndetectors = shape(data,1)

    integer, intent(in)   :: factor, nsamples, ndetectors
    real*8, intent(in)    :: compressed(nsamples,ndetectors)
    real*8, intent(out)   :: data(nsamples*factor,ndetectors)

    call transpos(compressed, data, factor)

end subroutine compression_average_transpose


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine downsampling_direct(data, compressed, factor, nsamples, ndetectors)

    use module_compression, only : direct => downsampling_direct
    implicit none

    !f2py threadsafe
    !f2py intent(in)      :: data
    !f2py intent(inout)   :: compressed
    !f2py intent(in)      :: factor
    !f2py intent(hide)    :: nsamples = shape(data,0)/factor
    !f2py intent(hide)    :: ndetectors = shape(data,1)

    integer, intent(in)   :: factor, nsamples, ndetectors
    real*8, intent(in)    :: data(nsamples*factor,ndetectors)
    real*8, intent(out)   :: compressed(nsamples,ndetectors)

    call direct(data, compressed, factor)

end subroutine downsampling_direct


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine downsampling_transpose(compressed, data, factor, nsamples, ndetectors)

    use module_compression, only : transpos => downsampling_transpose
    implicit none

    !f2py threadsafe
    !f2py intent(inout)   :: compressed
    !f2py intent(in)      :: data
    !f2py intent(in)      :: factor
    !f2py intent(hide)    :: nsamples = shape(data,0)/factor
    !f2py intent(hide)    :: ndetectors = shape(data,1)

    integer, intent(in)   :: factor, nsamples, ndetectors
    real*8, intent(in)    :: compressed(nsamples,ndetectors)
    real*8, intent(out)   :: data(nsamples*factor,ndetectors)

    call transpos(compressed, data, factor)

end subroutine downsampling_transpose


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine backprojection_weighted(pmatrix, data, mask, map1d, npixels_per_sample, nsamples, ndetectors, npixels)

    use module_pointingmatrix, only : bpw => backprojection_weighted, pointingelement
    implicit none

    !f2py threadsafe
    !f2py integer*8,intent(in):: pmatrix(npixels_per_sample*nsamples*ndetectors)
    !f2py intent(in)          :: data
    !f2py intent(in)          :: mask
    !f2py intent(inout)       :: map1d
    !f2py intent(in)          :: npixels_per_sample
    !f2py intent(hide)        :: nsamples = shape(data,0)
    !f2py intent(hide)        :: ndetectors = shape(data,1)
    !f2py intent(hide)        :: npixels = size(map1d)

    type(pointingelement), intent(in) :: pmatrix(npixels_per_sample,nsamples,ndetectors)
    real*8, intent(in)                :: data(nsamples,ndetectors)
    logical*1, intent(in)             :: mask(nsamples,ndetectors)
    real*8, intent(inout)             :: map1d(npixels)
    integer, intent(in)               :: npixels_per_sample
    integer, intent(in)               :: nsamples
    integer, intent(in)               :: ndetectors
    integer, intent(in)               :: npixels

    call bpw(pmatrix, data, mask, map1d)

end subroutine backprojection_weighted


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine deglitch_l2b_std(pmatrix, nx, ny, data, mask, nsigma, npixels_per_sample, nsamples, ndetectors)

    use module_pointingmatrix, only : pointingelement
    use module_deglitching, only : deglitch_l2b
    implicit none

    !f2py threadsafe
    !f2py integer*8, intent(in) :: pmatrix(npixels_per_sample*nsamples*ndetectors)
    !f2py intent(in)    :: nx, ny
    !f2py intent(in)    :: data
    !f2py intent(inout) :: mask
    !f2py intent(in)    :: nsigma
    !f2py intent(in)    :: npixels_per_sample
    !f2py intent(hide)  :: nsamples = shape(data,0)
    !f2py intent(hide)  :: ndetectors = shape(data,1)

    type(pointingelement), intent(in) :: pmatrix(npixels_per_sample,nsamples,ndetectors)
    integer, intent(in)               :: nx, ny
    real*8, intent(in)                :: data(nsamples,ndetectors)
    logical*1, intent(inout)          :: mask(nsamples,ndetectors)
    real*8, intent(in)                :: nsigma
    integer, intent(in)               :: npixels_per_sample, nsamples, ndetectors

    call deglitch_l2b(pmatrix, nx, ny, data, mask, nsigma, .false., verbose=.true.)

end subroutine deglitch_l2b_std


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine deglitch_l2b_mad(pmatrix, nx, ny, data, mask, nsigma, npixels_per_sample, nsamples, ndetectors)

    use module_pointingmatrix, only : pointingelement
    use module_deglitching,    only : deglitch_l2b
    implicit none

    !f2py threadsafe
    !f2py integer*8, intent(in) :: pmatrix(npixels_per_sample*nsamples*ndetectors)
    !f2py intent(in)    :: nx, ny
    !f2py intent(in)    :: data
    !f2py intent(inout) :: mask
    !f2py intent(in)    :: nsigma
    !f2py intent(in)    :: npixels_per_sample
    !f2py intent(hide)  :: nsamples = shape(data,0)
    !f2py intent(hide)  :: ndetectors = shape(data,1)

    type(pointingelement), intent(in) :: pmatrix(npixels_per_sample,nsamples,ndetectors)
    integer, intent(in)               :: nx, ny
    real*8, intent(in)                :: data(nsamples,ndetectors)
    logical*1, intent(inout)          :: mask(nsamples,ndetectors)
    real*8, intent(in)                :: nsigma
    integer, intent(in)               :: npixels_per_sample, nsamples, ndetectors

    call deglitch_l2b(pmatrix, nx, ny, data, mask, nsigma, .true., verbose=.true.)

end subroutine deglitch_l2b_mad


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine filter_median(data, length, nsamples, nsamples_tot, nslices, ndetectors, status)

    use iso_fortran_env,     only : ERROR_UNIT, OUTPUT_UNIT
    use module_preprocessor, only : median_filtering_nocopy
    implicit none

    !f2py threadsafe
    !f2py intent(in)   :: data
    !f2py intent(in)   :: length
    !f2py intent(in)   :: nsamples
    !f2py intent(hide) :: nsamples_tot=shape(data,0)
    !f2py intent(hide) :: nslices=size(nsamples)
    !f2py intent(hide) :: ndetectors=shape(data,1)
    !f2py intent(out)  :: status

    real*8, intent(inout) :: data(nsamples_tot,ndetectors)
    integer, intent(in)   :: length
    integer, intent(in)   :: nsamples(nslices)
    integer, intent(in)   :: nsamples_tot
    integer, intent(in)   :: nslices
    integer, intent(in)   :: ndetectors
    integer, intent(out)  :: status

    integer :: islice, start
    integer :: count1, count2, count_rate, count_max

    if (sum(nsamples) /= nsamples_tot) then
        status = 1
        write (ERROR_UNIT,'(a)') 'ERROR: The total number of samples is not the sum of the size of the slices.'
        return
    end if

    if (length <= 0) then
        status = 1
        write (ERROR_UNIT,'(a)') 'ERROR: The filter length must not be negative.'
        return
    end if

    status = 0

    write (OUTPUT_UNIT,'(a,i0,a)', advance='no') 'Median filtering (length=', length, ')... '
    call system_clock(count1, count_rate, count_max)
    start = 1
    do islice = 1, nslices
        call median_filtering_nocopy(data(start:start+nsamples(islice)-1,:), length)
        start = start + nsamples(islice)
    end do
    call system_clock(count2, count_rate, count_max)
    write (OUTPUT_UNIT,'(f7.2,a)') real(count2-count1)/count_rate, 's'

end subroutine filter_median


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine read_madmap1_nslices(invnttfile, ndetectors, nslices, status)

    use module_madcap, only : get_nfilters
    implicit none

    !f2py threadsafe
    !f2py intent(in)  :: invnttfile
    !f2py intent(in)  :: ndetectors
    !f2py intent(out) :: nslices
    !f2py intent(out) :: status

    character(len=*), intent(in) :: invnttfile
    integer, intent(in)          :: ndetectors
    integer, intent(out)         :: nslices
    integer, intent(out)         :: status

    integer                      :: nfilterfiles

    nfilterfiles = get_nfilters(invnttfile, ndetectors, status)
    if (status /= 0) return
    nslices = nfilterfiles / ndetectors

end subroutine read_madmap1_nslices


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine read_madmap1_info(todfile, invnttfile, convert, ndetectors, nslices, npixels_per_sample, nsamples, status)

    use iso_fortran_env, only : ERROR_UNIT
    use module_madcap, only   : read_filter_headers, read_tod_header
    use module_string, only   : strinteger
    implicit none

    !f2py threadsafe
    !f2py intent(in)  :: todfile
    !f2py intent(in)  :: invnttfile
    !f2py intent(in)  :: convert
    !f2py intent(in)  :: ndetectors
    !f2py intent(in)  :: nslices
    !f2py intent(out) :: npixels_per_sample
    !f2py intent(out) :: nsamples
    !f2py intent(out) :: status

    character(len=*), intent(in) :: todfile
    character(len=*), intent(in) :: invnttfile
    character(len=*), intent(in) :: convert
    integer, intent(in)          :: ndetectors
    integer, intent(in)          :: nslices
    integer, intent(out)         :: npixels_per_sample
    integer, intent(out)         :: nsamples(nslices)
    integer, intent(out)         :: status

    integer*8, allocatable       :: first(:), last(:)
    integer, allocatable         :: ncorrelations(:)
    integer                      :: idetector, ifilter, islice
    integer*8                    :: nsamples_tot

    call read_tod_header(todfile, convert, ndetectors, nsamples_tot, npixels_per_sample, status)
    if (status /= 0) return

    call read_filter_headers(invnttfile, convert, ndetectors, first, last, ncorrelations, status)
    if (status /= 0) return

    ! check number of correlations
    if (any(ncorrelations /= ncorrelations(1))) then
        status = 1
        write (ERROR_UNIT,'(a)') 'Error: Filter files do not have the same correlation length.'
        return
    end if

    ! check number of samples per slice
    do islice=1, nslices
        
        nsamples(islice) = last((islice-1)*ndetectors+1) - first((islice-1)*ndetectors+1) + 1
        
        do idetector = 1, ndetectors
            
            ifilter = (islice-1)*ndetectors + idetector
            if (last(ifilter) - first(ifilter) + 1 /= nsamples(islice)) then
                write (ERROR_UNIT,'(a)') "Error: Filter '" // invnttfile // '.' // strinteger(ifilter-1)// "' does not apply to th&
                     &e same number of samples than filter '" // invnttfile // '.' // strinteger((islice-1)*ndetectors) // "'."
                status = 1
                return
            end if
            
        end do

    end do

    ! check number of samples
    if (nsamples_tot /= sum(nsamples)) then
        status = 1
        write (ERROR_UNIT,'(a,2(i0,a))') "Error: Invalid number of samples in tod ('", nsamples_tot, "' instead of '",    &
             sum(nsamples), "')."
        return
    end if

end subroutine read_madmap1_info


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine read_madmap1(todfile, invnttfile, convert, npixels_per_sample, nsamples, ndetectors, tod, pmatrix, status)

    use module_filtering,      only : filterset
    use module_madcap,         only : read_filter, read_tod
    use module_pointingmatrix, only : pointingelement
    implicit none

    !f2py threadsafe
    !f2py intent(in)               :: todfile
    !f2py intent(in)               :: invnttfile
    !f2py intent(in)               :: convert
    !f2py intent(in)               :: npixels_per_sample
    !f2py intent(hide)             :: nsamples = shape(tod,0)
    !f2py intent(hide)             :: ndetectors = shape(tod,1)
    !f2py intent(inout)            :: tod(nsamples,ndetectors)
    !f2py integer*8, intent(inout) :: pmatrix(npixels_per_sample*nsamples*ndetectors)
    !f2py intent(out)              :: status

    character(len=*), intent(in)         :: todfile
    character(len=*), intent(in)         :: invnttfile
    character(len=*), intent(in)         :: convert
    integer, intent(in)                  :: npixels_per_sample
    integer*8, intent(in)                :: nsamples
    integer, intent(in)                  :: ndetectors
    real*8, intent(inout)                :: tod(nsamples,ndetectors)
    type(pointingelement), intent(inout) :: pmatrix(npixels_per_sample,nsamples,ndetectors)
    integer, intent(out)                 :: status

    type(filterset)                      :: filter
    

    call read_filter(invnttfile, convert, ndetectors, filter, status)
    if (status /= 0) return

    call read_tod(todfile, convert, filter%first, filter%last, tod, pmatrix, status)

end subroutine read_madmap1


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine invntt_madmap1(filename, convert, nslices, nsamples, nsamples_tot, ndetectors, tod_filter, ncorrelations,      &
                                   status)

    !f2py threadsafe
    !f2py intent(in)   :: filename
    !f2py intent(in)   :: convert
    !f2py intent(hide) :: nslices = size(nsamples)
    !f2py intent(in)   :: nsamples
    !f2py intent(in)   :: nsamples_tot
    !f2py intent(in)   :: ndetectors
    !f2py intent(out)  :: tod(nsamples_tot,ndetectors)
    !f2py intent(out)  :: ncorrelations
    !f2py intent(out)  :: status

    use module_filtering, only : create_filter_uncorrelated, filterset
    use module_madcap,    only : read_filter
    implicit none

    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: convert
    integer, intent(in)          :: nslices
    integer, intent(in)          :: nsamples(nslices)
    integer, intent(in)          :: nsamples_tot
    integer, intent(in)          :: ndetectors
    real*8, intent(out)          :: tod_filter(nsamples_tot,ndetectors)
    integer, intent(out)         :: ncorrelations
    integer, intent(out)         :: status

    type(filterset)              :: filter

    call read_filter(filename, convert, ndetectors, filter, status)
    if (status /= 0) return

    ncorrelations = filter%ncorrelations

    call create_filter_uncorrelated(filter, nsamples, nsamples_tot, tod_filter, status)
    if (status /= 0) return

end subroutine invntt_madmap1


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine sqrt_invntt_madmap1(filename, convert, nslices, nsamples, nsamples_tot, ndetectors, tod_filter, ncorrelations,      &
                                   status)

    !f2py threadsafe
    !f2py intent(in)   :: filename
    !f2py intent(in)   :: convert
    !f2py intent(hide) :: nslices = size(nsamples)
    !f2py intent(in)   :: nsamples
    !f2py intent(in)   :: nsamples_tot
    !f2py intent(in)   :: ndetectors
    !f2py intent(out)  :: tod(nsamples_tot,ndetectors)
    !f2py intent(out)  :: ncorrelations
    !f2py intent(out)  :: status

    use iso_fortran_env,  only : OUTPUT_UNIT
    use module_filtering, only : create_filter_uncorrelated, filterset
    use module_madcap,    only : read_filter
    implicit none

    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: convert
    integer, intent(in)          :: nslices
    integer, intent(in)          :: nsamples(nslices)
    integer, intent(in)          :: nsamples_tot
    integer, intent(in)          :: ndetectors
    real*8, intent(out)          :: tod_filter(nsamples_tot,ndetectors)
    integer, intent(out)         :: ncorrelations
    integer, intent(out)         :: status

    type(filterset)              :: filter

    call read_filter(filename, convert, ndetectors, filter, status)
    if (status /= 0) return

    ncorrelations = filter%ncorrelations

    call create_filter_uncorrelated(filter, nsamples, nsamples_tot, tod_filter, status)
    if (status /= 0) return

    if (any(tod_filter < 0)) then
        write (OUTPUT_UNIT,'(a)') 'SQRT_INVNTT_MADMAP1: some values are negative.'
    end if

    tod_filter = sqrt(tod_filter)

end subroutine sqrt_invntt_madmap1


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine unpack_direct(input, nvalids, mask, nx, ny, output, field)

    use iso_fortran_env, only : ERROR_UNIT
    implicit none

    !f2py threadsafe
    !f2py intent(in)    :: input
    !f2py intent(hide)  :: nvalids = size(input)
    !f2py intent(in)    :: mask
    !f2py intent(hide)  :: nx=shape(mask,0)
    !f2py intent(hide)  :: ny=shape(mask,1)
    !f2py intent(inout) :: output(nx,ny)
    !f2py intent(in)    :: field

    real*8, intent(in)    :: input(nvalids)
    integer, intent(in)   :: nvalids
    logical*1, intent(in) :: mask(nx,ny)
    integer, intent(in)   :: nx, ny
    real*8, intent(out)   :: output(nx,ny)
    real*8, intent(in)    :: field

    if (count(.not. mask) /= nvalids) then
        write (ERROR_UNIT,'(a)') 'UNPACK_DIRECT: The mask is not compatible with the input size.'
        return
    endif
    output = unpack(input, .not. mask, field)

end subroutine unpack_direct


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine unpack_transpose(input, mask, nx, ny, nvalids, output)

    use iso_fortran_env, only : ERROR_UNIT
    implicit none

    !f2py threadsafe
    !f2py intent(in)    :: input
    !f2py intent(in)    :: mask(nx,ny)
    !f2py intent(hide)  :: nx=shape(input,0)
    !f2py intent(hide)  :: ny=shape(input,1)
    !f2py intent(hide)  :: nvalids = size(output)
    !f2py intent(inout) :: output(nvalids)

    real*8, intent(in)    :: input(nx,ny)
    logical*1, intent(in) :: mask(nx,ny)
    integer, intent(in)   :: nx, ny
    integer, intent(in)   :: nvalids
    real*8, intent(out)   :: output(nvalids)

    if (count(.not. mask) /= nvalids) then
        write (ERROR_UNIT,'(a)') 'UNPACK_TRANSPOSE: The mask is not compatible with the output size.'
        return
    endif
    output = pack(input, .not. mask)

end subroutine unpack_transpose


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine masking(input, ninputs, mask, nmasks, status)

    use iso_fortran_env, only : ERROR_UNIT
    implicit none

    !f2py threadsafe
    !f2py intent(inout) :: input
    !f2py intent(hide)  :: ninputs=size(input)
    !f2py intent(in)    :: mask(nx)
    !f2py intent(hide)  :: nmasks=size(mask)
    !f2py intent(out)   :: status

    real*8, intent(inout) :: input(ninputs)
    logical*1, intent(in) :: mask(nmasks)
    integer, intent(in)   :: ninputs, nmasks
    integer, intent(out)  :: status

    if (ninputs /= nmasks) then
        write (ERROR_UNIT,'(a,2(i0,a))') "The data array has a size incompatible with the mask ('", ninputs, "' instead of '",     &
              nmasks, "')."
        status = 1
    end if

    !$omp parallel workshare
    where (mask)
        input = 0.d0
    end where
    !$omp end parallel workshare

    status = 0

end subroutine masking
