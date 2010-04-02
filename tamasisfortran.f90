! Tamasis interface for f2py
!
! Author: P. Chanial

subroutine pacs_info_channel(filename, nfilenames, channel, status)
    use module_pacsobservation, only : pacsobservation
    implicit none

    !f2py intent(in)   filename
    !f2py intent(in)   nfilenames
    !f2py intent(out)  channel
    !f2py intent(out)  status
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: nfilenames
    character, intent(out)       :: channel
    integer, intent(out)         :: status

    class(pacsobservation), allocatable :: obs
    character(len=len(filename)/nfilenames), allocatable :: filename_(:)
    integer :: iobs

    ! split input filename
    if (mod(len(filename), nfilenames) /= 0) then
        stop 'PACS_INFO: Invalid filename length.'
    end if

    allocate(filename_(nfilenames))

    do iobs = 1, nfilenames
        filename_(iobs) = filename((iobs-1)*len(filename_)+1:iobs*len(filename_))
    end do

    allocate(obs)
    call obs%init(filename_, status)
    if (status /= 0) go to 999

    channel = obs%info(1)%channel

!XXX GFORTRAN: segfaults with allocatable scalars... claims that obs is unallocated
!999 deallocate(obs)
!    if (allocated(pacs)) deallocate(pacs)
! deallocate filenames too
999 continue

end subroutine pacs_info_channel


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_info(filename, nfilenames, fine_sampling_factor, keep_bad_detectors, use_bad_detector_mask, bad_detector_mask,     &
                     nrows, ncolumns, ndetectors, output_mask, transparent_mode, compression_factor, nsamples, status)
    use module_pacsinstrument, only : pacsinstrument
    use module_pacsobservation, only : pacsobservation
    implicit none

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

    class(pacsobservation), allocatable :: obs
    class(pacsinstrument), allocatable  :: pacs
    character(len=len(filename)/nfilenames), allocatable :: filename_(:)
    integer :: iobs

    
    ! split input filename
    if (mod(len(filename), nfilenames) /= 0) then
        stop 'PACS_INFO: Invalid filename length.'
    end if

    allocate(filename_(nfilenames))

    do iobs = 1, nfilenames
        filename_(iobs) = filename((iobs-1)*len(filename_)+1:iobs*len(filename_))
    end do

    allocate(obs)
    call obs%init(filename_, status, verbose=.true.)
    if (status /= 0) go to 999

    allocate(pacs)
    if (use_bad_detector_mask) then
        call pacs%init(obs%channel, obs%transparent_mode, fine_sampling_factor, keep_bad_detectors, status, bad_detector_mask)
    else
        call pacs%init(obs%channel, obs%transparent_mode, fine_sampling_factor, keep_bad_detectors, status)
    end if
    if (status /= 0) go to 999

    ndetectors  = pacs%ndetectors
    output_mask = pacs%mask
    transparent_mode = obs%info(1)%transparent_mode
    compression_factor = obs%info%compression_factor
    nsamples = obs%info%nsamples

!XXX GFORTRAN: segfaults... claims that obs is unallocated
!999 deallocate(obs)
!    if (allocated(pacs)) deallocate(pacs)
999 continue

end subroutine pacs_info


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_map_header(filename, nfilenames, finer_sampling, fine_sampling_factor, keep_bad_detectors, use_bad_detector_mask,  &
                           bad_detector_mask, nrows, ncolumns, resolution, header, status)
    use module_pacsinstrument,  only : pacsinstrument
    use module_pacsobservation, only : pacsobservation
    use module_pacspointing,    only : pacspointing
    implicit none

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

    class(pacsobservation), allocatable :: obs
    class(pacsinstrument), allocatable  :: pacs
    class(pacspointing), allocatable    :: pointing
    character(len=len(filename)/nfilenames), allocatable :: filename_(:)
    integer :: iobs
    
    ! split input filename
    if (mod(len(filename), nfilenames) /= 0) then
        stop 'PACS_INFO: Invalid filename length.'
    end if

    allocate(filename_(nfilenames))

    do iobs = 1, nfilenames
        filename_(iobs) = filename((iobs-1)*len(filename_)+1:iobs*len(filename_))
    end do

    allocate(obs)
    call obs%init(filename_, status)
    if (status /= 0) go to 999

    allocate(pacs)
    if (use_bad_detector_mask) then
        call pacs%init(obs%channel, obs%transparent_mode, fine_sampling_factor, keep_bad_detectors, status, bad_detector_mask)
    else
        call pacs%init(obs%channel, obs%transparent_mode, fine_sampling_factor, keep_bad_detectors, status)
    end if
    if (status /= 0) go to 999
    
    allocate(pointing)
    call pointing%init(obs, status)
    if (status /= 0) go to 999

    call pacs%compute_mapheader(pointing, finer_sampling, resolution, header, status)
    
999 continue

end subroutine pacs_map_header


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_timeline(filename, nfilenames, nsamples, ndetectors, keep_bad_detectors, bad_detector_mask, nrow, ncol,            &
                         do_flatfielding, do_subtraction_mean, signal, mask, status)
    use iso_fortran_env,        only : ERROR_UNIT
    use module_pacsinstrument,  only : pacsinstrument
    use module_pacsobservation, only : pacsobservation
    use module_preprocessor,    only : divide_vectordim2, subtract_meandim1
    implicit none

    !f2py threadsafe
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

    class(pacsobservation), allocatable :: obs
    class(pacsinstrument), allocatable  :: pacs
    character(len=len(filename)/nfilenames), allocatable :: filename_(:)
    integer                                              :: iobs, destination

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
    call obs%init(filename_, status)
    if (status /= 0) go to 999

    ! initialise pacs instrument
    allocate(pacs)
    call pacs%init(obs%channel, obs%transparent_mode, 1, keep_bad_detectors, status, bad_detector_mask)
    if (status /= 0) go to 999

    ! read timeline
    call pacs%read(obs, signal, mask, status)
    if (status /= 0) go to 999

    ! flat fielding
    if (do_flatfielding) then
        call divide_vectordim2(signal, pacs%flatfield)
    end if
    
    ! subtract the mean of each detector timeline
    if (do_subtraction_mean) then
        destination = 1
        do iobs=1, nfilenames
            call subtract_meandim1(signal(destination:destination+obs%info(iobs)%nsamples-1,:))
            destination = destination + obs%info(iobs)%nsamples
         end do
    end if

999 continue

end subroutine pacs_timeline


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_pointing_matrix_filename(filename, nfilenames, finer_sampling, fine_sampling_factor, npixels_per_sample, nsamples, &
                                         ndetectors, keep_bad_detectors, bad_detector_mask, nrow, ncol, header, pmatrix, status)
    use iso_fortran_env,        only : ERROR_UNIT
    use module_fitstools,       only : ft_readparam
    use module_pacsinstrument,  only : pacsinstrument
    use module_pacsobservation, only : pacsobservation
    use module_pacspointing,    only : pacspointing
    use module_pointingmatrix,  only : pointingelement
    implicit none

    !f2py threadsafe
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

    class(pacsobservation), allocatable :: obs
    class(pacsinstrument), allocatable  :: pacs
    class(pacspointing), allocatable    :: pointing
    character(len=len(filename)/nfilenames), allocatable :: filename_(:)
    integer                                              :: iobs, nx, ny, count, nsamples_expected

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
    call obs%init(filename_, status)
    if (status /= 0) go to 999

    ! initialise pacs instrument
    allocate(pacs)
    call pacs%init(obs%channel, obs%transparent_mode, fine_sampling_factor, keep_bad_detectors, status, bad_detector_mask)
    if (status /= 0) go to 999

    ! check number of detectors
    if (pacs%ndetectors /= ndetectors) then
        status = 1
        write (ERROR_UNIT,'(a,2(i0,a))') "The specified number of detectors '", ndetectors, "' is incompatible with that from the i&
              &nput bad detector mask '", pacs%ndetectors, "'."
        return
    end if

    ! check number of fine samples
    if (finer_sampling) then
        nsamples_expected = sum(obs%info%nsamples * obs%info%compression_factor)*fine_sampling_factor
    else
        nsamples_expected = sum(obs%info%nsamples)
    end if
    if (nsamples /= nsamples_expected) then
        status = 1
        write (ERROR_UNIT,'(a,2(i0,a))') "The specified total number of samples '", nsamples, "' is incompatible with that from the&
              & observations '", nsamples_expected, "'."
        return
    end if

    ! initialise pointing
    allocate(pointing)
    call pointing%init(obs, status)
    if (status /= 0) go to 999

    ! get the size of the map
    call ft_readparam(header, 'naxis1', count, nx, status=status)
    if (status /= 0 .or. count == 0) go to 999
    call ft_readparam(header, 'naxis2', count, ny, status=status)
    if (status /= 0 .or. count == 0) go to 999

    ! compute the projector
    call pacs%compute_projection_sharp_edges(pointing, finer_sampling, header,  nx, ny, pmatrix, status)

999 continue

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
!!$    call pacs%compute_mapheader(pointing, .true., resolution, header, status)
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
!!$    use module_fitstools, only : ft_readparam
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
!!$    call ft_readparam(header, 'naxis1', count, nx, status=status)
!!$    if (status /= 0 .or. count == 0) go to 999
!!$    call ft_readparam(header, 'naxis2', count, ny, status=status)
!!$    if (status /= 0 .or. count == 0) go to 999
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


subroutine pacs_projection_sharp_edges_direct(pmatrix, map1d, signal, npixels_per_sample, nsamples, ndetectors, npixels)

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

end subroutine pacs_projection_sharp_edges_direct


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_projection_sharp_edges_transpose(pmatrix, signal, map1d, npixels_per_sample, nsamples, ndetectors, npixels)

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

end subroutine pacs_projection_sharp_edges_transpose


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


subroutine backprojection_weighted(pmatrix, data, mask, map1d, npixels_per_sample, nsamples, ndetectors, npixels)
    use module_pointingmatrix, only : bpw => backprojection_weighted, pointingelement
    implicit none
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


subroutine deglitch_l2b_std(pmatrix, nx, ny, data, mask, nsigma, outmask, npixels_per_sample, nsamples, ndetectors)
    use module_pointingmatrix, only : pointingelement
    use module_deglitching, only : deglitch_l2b
    implicit none
    !f2py integer*8, intent(in) :: pmatrix(npixels_per_sample*nsamples*ndetectors)
    !f2py intent(in)    :: nx, ny
    !f2py intent(in)    :: data
    !f2py intent(in)    :: mask
    !f2py intent(in)    :: nsigma
    !f2py intent(in)    :: npixels_per_sample
    !f2py intent(hide)  :: nsamples = shape(data,0)
    !f2py intent(hide)  :: ndetectors = shape(data,1)
    !f2py intent(out)   :: outmask
    type(pointingelement), intent(in) :: pmatrix(npixels_per_sample,nsamples,ndetectors)
    integer, intent(in)               :: nx, ny
    real*8, intent(in)                :: data(nsamples,ndetectors)
    logical*1, intent(in)             :: mask(nsamples,ndetectors)
    logical*1, intent(out)            :: outmask(nsamples,ndetectors)
    real*8, intent(in)                :: nsigma
    integer, intent(in)               :: npixels_per_sample, nsamples, ndetectors
    integer                           :: count1, count2, count_rate, count_max

    write(*,'(a)', advance='no') 'Info: deglitching (std)... '
    call system_clock(count1, count_rate, count_max)
    outmask = mask
    call deglitch_l2b(pmatrix, nx, ny, data, outmask, nsigma, .false.)
    call system_clock(count2, count_rate, count_max)
    write(*,'(f6.2,a)') real(count2-count1)/count_rate, 's'

end subroutine deglitch_l2b_std


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine deglitch_l2b_mad(pmatrix, nx, ny, data, mask, nsigma, outmask, npixels_per_sample, nsamples, ndetectors)
    use module_pointingmatrix, only : pointingelement
    use module_deglitching, only : deglitch_l2b
    implicit none
    !f2py integer*8, intent(in) :: pmatrix(npixels_per_sample*nsamples*ndetectors)
    !f2py intent(in)    :: nx, ny
    !f2py intent(in)    :: data
    !f2py intent(in)    :: mask
    !f2py intent(in)    :: nsigma
    !f2py intent(in)    :: npixels_per_sample
    !f2py intent(hide)  :: nsamples = shape(data,0)
    !f2py intent(hide)  :: ndetectors = shape(data,1)
    !f2py intent(out)   :: outmask
    type(pointingelement), intent(in) :: pmatrix(npixels_per_sample,nsamples,ndetectors)
    integer, intent(in)               :: nx, ny
    real*8, intent(in)                :: data(nsamples,ndetectors)
    logical*1, intent(in)             :: mask(nsamples,ndetectors)
    logical*1, intent(out)            :: outmask(nsamples,ndetectors)
    real*8, intent(in)                :: nsigma
    integer, intent(in)               :: npixels_per_sample, nsamples, ndetectors
    integer                           :: count1, count2, count_rate, count_max

    write(*,'(a)', advance='no') 'Info: deglitching (mad)... '
    call system_clock(count1, count_rate, count_max)
    outmask = mask
    call deglitch_l2b(pmatrix, nx, ny, data, outmask, nsigma, .true.)
    call system_clock(count2, count_rate, count_max)
    write(*,'(f6.2,a)') real(count2-count1)/count_rate, 's'

end subroutine deglitch_l2b_mad

