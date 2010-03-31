
    ! check that we're not using the wrong calib file for blue observation


subroutine pacs_info(filename, keep_bad_detectors, bad_detector_mask, channel, transparent_mode, compression_factor, first, last, ndetectors, status)
    use, intrinsic :: ISO_FORTRAN_ENV
    use module_pacsinstrument
    implicit none

    !f2py intent(in)  filename
    !f2py intent(in)  keep_bad_detectors
    !f2py intent(out) channel
    !f2py intent(out) transparent_mode
    !f2py intent(out) compression_factor
    !f2py intent(out) ndetectors
    !f2py intent(out) first, last
    !f2py intent(out) status
    character(len=*), intent(in) :: filename
    logical, intent(in)          :: keep_bad_detectors
    character, intent(out)       :: channel
    logical, intent(out)         :: transparent_mode
    integer, intent(out)         :: compression_factor
    integer, intent(out)         :: ndetectors
    integer, intent(out)         :: first, last
    integer, intent(out)         :: status

    class(pacsobservation), allocatable :: obs
    class(pacsinstrument), allocatable  :: pacs
    integer                             :: status

    allocate(obs)
    call obs%init([filename], status)
    if (status /= 0) return

    allocate(pacs)
    call pacs%init(obs, .false., status, bad_detector_mask)
    if (status /= 0) return

    channel = obs%info(1)%channel
    transparent_mode = obs%info(1)%transparent_mode
    compression_factor = obs%info(1)%compression_factor
    ndetectors = pacs%ndetectors
    first = obs%info(1)%first
    last = obs%info(1)%last

999 deallocate(obs)
    if allocated(pacs) deallocate(pacs)

end subroutine pacs_info_ndetectors


!-------------------------------------------------------------------------------


subroutine pacs_timeline(filename, nsamples, ndetectors, bad_detector_mask, &
                         nrow, ncol, keep_bad_detectors, signal, mask, status)
    use, intrinsic :: ISO_FORTRAN_ENV, only : ERROR_UNIT
    use module_pacsinstrument, only : pacsinstrument
    use module_preprocessor, only : divide_vectordim2, subtract_meandim1
    implicit none

    !f2py threadsafe
    !f2py intent(in)   :: filename
    !f2py intent(in)   :: first
    !f2py intent(in)   :: last
    !f2py intent(in)   :: ndetectors
    !f2py intent(in)   :: bad_detector_mask
    !f2py intent(hide) :: nrow = shape(bad_detector_mask,0)
    !f2py intent(hide) :: ncol = shape(bad_detector_mask,1)
    !f2py intent(in)   :: keep_bad_detectors
    !f2py intent(out)  :: signal
    !f2py intent(out)  :: mask
    !f2py intent(out)  :: status
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: nsamples
    integer, intent(in)          :: ndetectors
    logical*1, intent(in)        :: bad_detector_mask(nrow,ncol)
    integer, intent(in)          :: nrow, ncol
    logical, intent(in)          :: keep_bad_detectors
    real*8, intent(out)          :: signal(nsamples, ndetectors)
    logical*1, intent(out)       :: mask(nsamples, ndetectors)
    integer, intent(out)         :: status

    class(pacsinstrument), allocatable :: obs
    class(pacsinstrument), allocatable :: pacs

    allocate(obs)
    call obs%init([filename], status)
    if (status /= 0) go to 999

    allocate(pacs)
    call pacs%init(obs, fine_sampling_factor, keep_bad_detectors, status,      &
                   bad_detector_mask)
    if (status /= 0) goto 999

    call pacs%read(obs, pacs%pq, signal, mask, status)
    if (status /= 0) goto 999

    call divide_vectordim2(signal, pacs%flatfield)
    
    do iobs=1, size(obs%info)
        call subtract_meandim1(signal(obs%info(iobs)%first:                    &
            obs%info(iobs)%last,:))
    end do

999 deallocate(obs)
    if (allocated(pacs)) deallocate(pacs)

end subroutine pacs_timeline


!-------------------------------------------------------------------------------


subroutine pacs_map_header(array, time, ra, dec, pa, chop, npointings,         &
                           finetime, nfinesamples, transparent_mode,           &
                           keep_bad_detectors, bad_detector_mask, nrow, ncol,  &
                           resolution, header)

    use, intrinsic :: ISO_FORTRAN_ENV, only : ERROR_UNIT
    use module_pacsinstrument, only : pacsinstrument
    use module_pacspointing, only : pacspointing
    implicit none

    !f2py threadsafe
    !f2py intent(in)   :: array
    !f2py intent(in)   :: time, ra, dec, pa, chop
    !f2py intent(hide) :: npointings = size(time)
    !f2py intent(in)   :: finetime
    !f2py intent(hide) :: nfinesamples = size(finetime)
    !f2py intent(in)   :: transparent_mode
    !f2py intent(in)   :: keep_bad_detectors
    !f2py intent(in)   :: bad_detector_mask
    !f2py intent(hide), depend(bad_detector_mask) :: nrow = shape(bad_detector_mask,0), ncol = shape(bad_detector_mask,1)
    !f2py intent(in)   :: resolution
    !f2py intent(out)  :: header

    character(len=*), intent(in)     :: array
    real*8, intent(in)               :: time(npointings), ra(npointings), dec(npointings), pa(npointings), chop(npointings)
    integer, intent(in)              :: npointings
    real*8, intent(in)               :: finetime(nfinesamples)
    integer, intent(in)              :: nfinesamples
    logical, intent(in)              :: transparent_mode
    logical, intent(in)              :: keep_bad_detectors
    logical*1, intent(in)            :: bad_detector_mask(nrow,ncol)
    integer, intent(in)              :: nrow, ncol
    real*8, intent(in)               :: resolution
    character(len=2880), intent(out) :: header
    integer, intent(out)             :: status

    class(pacsinstrument), allocatable :: pacs
    class(pacspointing), allocatable   :: pointing

    ! read pointing information
    allocate(pointing)
    call pointing%init_sim(time, ra, dec, pa, chop, status)
    if (status /= 0) goto 999

    ! get the pacs instance, read the calibration files
    allocate(pacs)
    call pacs%init(array(1:1), transparent_mode, keep_bad_detectors, status,   &
                   bad_detector_mask)
    if (status /= 0) goto 999

    if (any(pacs%mask .neqv. bad_detector_mask)) then
        write (*,'(a)') "Info: using user's bad detector mask."
    end if

    call pacs%compute_mapheader(pointing, .true., resolution, header, status)
    if (status == 0) return

999 deallocate(obs)
    if (allocated(pointing)) deallocate(pointing)
    if (allocate(pacs))      deallocate(pacs)

end subroutine pacs_map_header


!-------------------------------------------------------------------------------


subroutine pacs_pointing_matrix(array, time, ra, dec, pa, chop, npointings,    &
                                finetime, nfinesamples, npixels_per_sample,    &
                                ndetectors, transparent_mode,                  &
                                keep_bad_detectors, bad_detector_mask,         &
                                nrow, ncol, header, pmatrix)

    use, intrinsic :: ISO_FORTRAN_ENV
    use module_fitstools, only : ft_readparam
    use module_pacsinstrument, only : pacsinstrument
    use module_pacspointing, only : pacspointing
    use module_pointingmatrix, only : pointingelement
    implicit none

    !f2py threadsafe
    !f2py intent(in)   :: array
    !f2py intent(in)   :: time, ra, dec, pa, chop
    !f2py intent(hide) :: npointings = size(time)
    !f2py intent(in)   :: finetime
    !f2py intent(hide) :: nfinesamples = size(finetime)
    !f2py intent(in)   :: npixels_per_sample
    !f2py intent(in)   :: ndetectors
    !f2py intent(in)   :: transparent_mode
    !f2py intent(in)   :: keep_bad_detectors
    !f2py intent(in)   :: bad_detector_mask
    !f2py intent(hide), depend(bad_detector_mask) :: nrow = shape(bad_detector_mask,0), ncol = shape(bad_detector_mask,1)
    !f2py intent(in)   :: header
    !f2py integer*8, dimension(npixels_per_sample*nfinesamples*ndetectors), intent(inout) :: pmatrix

    character(len=*), intent(in)         :: array
    real*8, intent(in)                   :: time(npointings), ra(npointings), dec(npointings), pa(npointings), chop(npointings)
    integer, intent(in)                  :: npointings
    real*8, intent(in)                   :: finetime(nfinesamples)
    integer, intent(in)                  :: nfinesamples
    integer, intent(in)                  :: npixels_per_sample
    integer, intent(in)                  :: ndetectors
    logical, intent(in)                  :: transparent_mode
    logical, intent(in)                  :: keep_bad_detectors
    logical*1, intent(in)                :: bad_detector_mask(nrow,ncol)
    integer, intent(in)                  :: nrow, ncol
    character(len=*), intent(in)         :: header
    type(pointingelement), intent(inout) :: pmatrix(npixels_per_sample, nfinesamples, ndetectors)

    class(pacsinstrument), allocatable   :: pacs
    class(pacspointing), allocatable     :: pointing
    integer                              :: status, count,  count1, count2, count_rate, count_max, nx, ny

    ! read pointing information
    allocate(pointing)
    call pointing%load_array(time, ra, dec, pa, chop, status)
    if (status /= 0) goto 999

    ! get the pacs instance, read the calibration files
    allocate(pacs)
    call pacs%init(array(1:1), transparent_mode, keep_bad_detectors, status,   &
                   bad_detector_mask)
    if (status /= 0) goto 999

    if (any(pacs%mask .neqv. bad_detector_mask)) then
        write (*,'(a)') "Info: using user's bad pixel mask."
    end if

    call ft_readparam(header, 'naxis1', count, nx, status=status)
    if (status /= 0 .or. count == 0) goto 999
    call ft_readparam(header, 'naxis2', count, ny, status=status)
    if (status /= 0 .or. count == 0) goto 999

    ! compute the projector
    write(*,'(a)', advance='no') 'Info: computing the projector... '
    call system_clock(count1, count_rate, count_max)
    call pacs%compute_projection_sharp_edges(pointing, finetime, header, nx, ny, pmatrix, status)
    call system_clock(count2, count_rate, count_max)
    write(*,'(f6.2,a)') real(count2-count1)/count_rate, 's'

    return

999 write(ERROR_UNIT,'(a)') 'Aborting.'

end subroutine pacs_pointing_matrix


!-------------------------------------------------------------------------------


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


!-------------------------------------------------------------------------------


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


!-------------------------------------------------------------------------------


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


!-------------------------------------------------------------------------------


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


!-------------------------------------------------------------------------------


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


!-------------------------------------------------------------------------------


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


!-------------------------------------------------------------------------------


subroutine backprojection_weighted(pmatrix, data, mask, map1d, npixels_per_sample, nsamples, ndetectors, npixels)
    use module_pointingmatrix, only : bpw => backprojection_weighted,           &
                                      pointingelement
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


!-------------------------------------------------------------------------------


subroutine deglitch_l2b_std(pmatrix, nx, ny, data, mask, nsigma, outmask, npixels_per_sample, nsamples, ndetectors)
    use module_pointingmatrix, only : pointingelement
    use module_deglitching, only : deglitch_l2b
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


!-------------------------------------------------------------------------------


subroutine deglitch_l2b_mad(pmatrix, nx, ny, data, mask, nsigma, outmask, npixels_per_sample, nsamples, ndetectors)
    use module_pointingmatrix, only : pointingelement
    use module_deglitching, only : deglitch_l2b
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

