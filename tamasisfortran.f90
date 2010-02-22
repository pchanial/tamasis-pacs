
    ! check that we're not using the wrong calib file for blue observation


subroutine pacs_info_ndetectors(filename, transparent_mode, ndetectors)
    use module_pacsinstrument
    implicit none

    !f2py intent(in) filename
    !f2py intent(in) transparent_mode
    !f2py intent(out) ndetectors
    character(len=*), intent(in) :: filename
    logical, intent(in)          :: transparent_mode
    integer, intent(out)         :: ndetectors

    class(pacsinstrument), allocatable :: pacs

    allocate(pacs)
    call pacs%read_calibration_files()
    call pacs%filter_detectors(pacs%get_array_color(filename), transparent_mode=transparent_mode)
    ndetectors = pacs%ndetectors

end subroutine pacs_info_ndetectors


!-------------------------------------------------------------------------------


subroutine pacs_info_ij(filename, transparent_mode, ndetectors, ij)
    use module_pacsinstrument
    implicit none

    !f2py intent(in) filename
    !f2py intent(in) transparent_mode
    !f2py intent(in) ndetectors
    !f2py intent(out) ij(2,ndetectors)
    character(len=*), intent(in) :: filename
    logical, intent(in)          :: transparent_mode
    integer, intent(in)          :: ndetectors
    integer, intent(out)         :: ij(2,ndetectors)

    class(pacsinstrument), allocatable :: pacs

    allocate(pacs)
    call pacs%read_calibration_files()
    call pacs%filter_detectors(pacs%get_array_color(filename), transparent_mode=transparent_mode)
    ij = pacs%ij

end subroutine pacs_info_ij


!-------------------------------------------------------------------------------


subroutine pacs_info_nsamples(filename, nsamples)
    use, intrinsic :: ISO_FORTRAN_ENV
    use module_fitstools
    use module_pacsinstrument
    implicit none

    !f2py intent(in) filename
    !f2py intent(in) transparent_mode
    !f2py intent(out) nsamples
    character(len=*), intent(in) :: filename
    integer, intent(out)         :: nsamples

    integer, allocatable               :: imageshape(:)
    integer                            :: unit, status
    status = 0
    call ft_openimage(filename // "_Signal.fits", unit, 3, imageshape, status)
    call ft_close(unit, status)
    call ft_printerror(status, filename)
    if (status /= 0) return

    if (size(imageshape) /= 3) then
        write(ERROR_UNIT, *) 'Number of dimensions is not 3 in ' // filename // "_Signal.fits"
        return
    endif
    nsamples = imageshape(1)

end subroutine pacs_info_nsamples


!-------------------------------------------------------------------------------


subroutine pacs_timeline(filename, transparent_mode, first, last, ndetectors, signal, mask)
    use module_pacsinstrument
    use module_preprocessor
    implicit none

    !f2py threadsafe
    !f2py intent(in)  :: filename
    !f2py intent(in)  :: transparent_mode
    !f2py intent(in)  :: first
    !f2py intent(in)  :: last
    !f2py intent(in)  :: ndetectors
    !f2py intent(out) :: signal
    !f2py intent(out) :: mask
    character(len=*), intent(in) :: filename
    logical, intent(in)          :: transparent_mode
    integer*8, intent(in)        :: first, last
    integer, intent(in)          :: ndetectors
    real*8, intent(out)          :: signal(last-first+1, ndetectors)
    logical*1, intent(out)       :: mask(last-first+1, ndetectors)

    class(pacsinstrument), allocatable :: pacs

    allocate(pacs)
    call pacs%read_calibration_files()
    call pacs%filter_detectors(get_array_color(filename), transparent_mode=transparent_mode)
    call pacs%read_signal_file(filename // '_Signal.fits', first, last, signal)
    call pacs%read_mask_file(filename // '_Mask.fits', first, last, mask)
    call divide_vectordim2(signal, pacs%flatfield)
    call subtract_meandim1(signal)

end subroutine pacs_timeline


!-------------------------------------------------------------------------------


subroutine pacs_map_header_file(filename, transparent_mode, first, last, fine_sampling_factor, compression_factor, &
                                resolution, header)
    use module_fitstools
    use module_pacsinstrument
    use module_pacspointing
    implicit none

    !f2py intent(in) filename
    !f2py intent(in) transparent_mode
    !f2py intent(in) first
    !f2py intent(in) last
    !f2py intent(in) fine_sampling_factor
    !f2py intent(in) compression_factor
    !f2py intent(in) resolution
    !f2py intent(out) header
    character(len=*), intent(in)       :: filename
    logical, intent(in)                :: transparent_mode
    integer*8, intent(in)              :: first, last
    integer, intent(in)                :: fine_sampling_factor, compression_factor
    real*8, intent(in)                 :: resolution
    character(len=2880), intent(out)   :: header

    class(pacsinstrument), allocatable :: pacs
    class(pacspointing), allocatable   :: pointing
    integer                            :: nsamples, status, i, j, sampling_factor
    real*8                             :: delta_time
    real*8, allocatable                :: time(:)
    integer*8, allocatable             :: timeus(:)

    nsamples = last - first + 1
    sampling_factor = fine_sampling_factor * compression_factor
    allocate(time(nsamples * sampling_factor))
    allocate(timeus(nsamples))
    status = 0
    call ft_readslice(filename // '_Time.fits+1', first, last, timeus, status)
    call ft_printerror(status, filename // '_Time.fits+1')
    if (status /= 0) return

    delta_time = (timeus(2)-timeus(1)) * 1.0d-6
    do i = 1, nsamples
       do j = 1, sampling_factor
           time((i-1)*sampling_factor+j) = timeus(i) * 1.0d-6 - (sampling_factor-j) * delta_time / sampling_factor
       end do
    end do

    allocate(pointing)
    call pointing%load(filename)

    ! get the pacs instance, read the calibration files
    allocate(pacs)
    call pacs%read_calibration_files()
    call pacs%filter_detectors(get_array_color(filename), transparent_mode=transparent_mode)
    call pacs%compute_mapheader(pointing, time, resolution, header)

end subroutine pacs_map_header_file


!-------------------------------------------------------------------------------


subroutine pacs_pointing_matrix_file(filename, transparent_mode, npixels_per_sample, first, last, ndetectors, &
                                     fine_sampling_factor, compression_factor, header, pmatrix)
    use module_fitstools
    use module_pacsinstrument
    use module_pacspointing
    use module_wcslib, only : WCSLEN, wcsfree
    implicit none

    !f2py threadsafe
    !f2py intent(in) :: filename
    !f2py intent(in) :: transparent_mode
    !f2py intent(in) :: npixels_per_sample
    !f2py intent(in) :: first
    !f2py intent(in) :: last
    !f2py intent(in) :: ndetectors
    !f2py intent(in) :: fine_sampling_factor
    !f2py intent(in) :: compression_factor
    !f2py intent(in) :: header
    !f2py integer*8, dimension(npixels_per_sample*(last-first+1)*ndetectors*fine_sampling_factor*compression_factor), intent(inout) :: pmatrix

    character(len=*), intent(in)         :: filename
    logical, intent(in)                  :: transparent_mode
    integer, intent(in)                  :: npixels_per_sample, ndetectors, fine_sampling_factor, compression_factor
    integer*8, intent(in)                :: first, last
    character(len=*), intent(in)         :: header
    type(pointingelement), intent(inout) :: pmatrix(npixels_per_sample, (last-first+1)*fine_sampling_factor*compression_factor, &
                                                    ndetectors)
    class(pacsinstrument), allocatable   :: pacs
    class(pacspointing), allocatable     :: pointing
    real*8, allocatable                  :: finetime(:)
    integer*8, allocatable               :: timeus(:)
    real*8                               :: delta_time
    integer                              :: wcs(WCSLEN), nx, ny
    integer                              :: i, j, status, count1, count2, count_rate, count_max, sampling_factor
    integer*8                            :: nsamples

    ! read pointing information
    allocate(pointing)
    call pointing%load(filename)

    ! read the time file
    status = 0
    nsamples = last - first + 1
    sampling_factor = fine_sampling_factor * compression_factor
    allocate(finetime(nsamples * sampling_factor))
    allocate(timeus(nsamples))
    call ft_readslice(filename // '_Time.fits+1', first, last, timeus, status)
    call ft_printerror(status, filename // '_Time.fits+1')

    delta_time = (timeus(2)-timeus(1)) * 1.0d-6 / sampling_factor
    do i = 1, nsamples
       do j = 1, sampling_factor
           finetime((i-1)*sampling_factor+j) = timeus(i) * 1.0d-6 + (j - 1) * delta_time
       end do
    end do

    ! get the pacs instance, read the calibration files
    allocate(pacs)
    call pacs%read_calibration_files()
    call pacs%filter_detectors(get_array_color(filename), transparent_mode=transparent_mode)
    call ft_header2wcs(header, wcs, nx, ny)

    ! compute the projector
    write(*,'(a)', advance='no') 'Info: computing the projector... '
    call system_clock(count1, count_rate, count_max)
    call pacs%compute_projection_sharp_edges(pointing, finetime, wcs, nx, pmatrix)
    call system_clock(count2, count_rate, count_max)
    write(*,'(f6.2,a)') real(count2-count1)/count_rate, 's'

    ! free the wcs
    status = wcsfree(wcs)

    ! fortran2008: sizeof_pmatrix = storage_size(pmatrix)
    !sizeof_pmatrix = 8 * npixels_per_sample * nsamples * ndetectors

end subroutine pacs_pointing_matrix_file


!-------------------------------------------------------------------------------


subroutine pacs_pointing_matrix_array(time, ra, dec, pa, chop, npointings, finetime, nfinesamples, npixels_per_sample, ndetectors, &
                                      array, transparent_mode, header, pmatrix)

    use module_fitstools
    use module_pacsinstrument
    use module_pacspointing
    use module_wcslib, only : WCSLEN, wcsfree
    implicit none

    !f2py threadsafe
    !f2py intent(in)   :: time, ra, dec, pa, chop
    !f2py intent(hide) :: npointings = size(time)
    !f2py intent(in)   :: finetime
    !f2py intent(hide) :: nfinesamples = size(finetime)
    !f2py intent(in)   :: npixels_per_sample
    !f2py intent(in)   :: ndetectors
    !f2py intent(in)   :: array
    !f2py intent(in)   :: transparent_mode
    !f2PY intent(in)   :: nsamples
    !f2py intent(in)   :: header
    !f2py integer*8, dimension(npixels_per_sample*nfinesamples*ndetectors), intent(inout) :: pmatrix

    real*8, intent(in)                   :: time(npointings), ra(npointings), dec(npointings), pa(npointings), chop(npointings)
    integer, intent(in)                  :: npointings
    real*8, intent(in)                   :: finetime(nfinesamples)
    integer, intent(in)                  :: nfinesamples
    integer, intent(in)                  :: npixels_per_sample
    integer, intent(in)                  :: ndetectors
    character(len=*), intent(in)         :: array
    logical, intent(in)                  :: transparent_mode
    character(len=*), intent(in)         :: header
    type(pointingelement), intent(inout) :: pmatrix(npixels_per_sample, nfinesamples, ndetectors)

    class(pacsinstrument), allocatable   :: pacs
    class(pacspointing), allocatable     :: pointing
    integer                              :: wcs(WCSLEN), nx, ny, status, count1, count2, count_rate, count_max

    ! read pointing information
    allocate(pointing)
    call pointing%load_array(time, ra, dec, pa, chop)

    ! get the pacs instance, read the calibration files
    allocate(pacs)
    call pacs%read_calibration_files()
    call pacs%filter_detectors(array, transparent_mode=transparent_mode)
    call ft_header2wcs(header, wcs, nx, ny)

    ! compute the projector
    write(*,'(a)', advance='no') 'Info: computing the projector... '
    call system_clock(count1, count_rate, count_max)
    call pacs%compute_projection_sharp_edges(pointing, finetime, wcs, nx, pmatrix)
    call system_clock(count2, count_rate, count_max)
    write(*,'(f6.2,a)') real(count2-count1)/count_rate, 's'

    ! free the wcs
    status = wcsfree(wcs)

end subroutine pacs_pointing_matrix_array


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
    real*8, intent(in)    :: signal(nsamples, ndetectors)
    real*8, intent(inout) :: multiplexed(nsamples/fine_sampling_factor, ndetectors)
    integer, intent(in)   :: fine_sampling_factor, ij(2, ndetectors), nsamples, ndetectors

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
    real*8, intent(in)    :: multiplexed(nsamples/fine_sampling_factor, ndetectors)
    real*8, intent(inout) :: signal(nsamples, ndetectors)
    integer, intent(in)   :: fine_sampling_factor, ij(2, ndetectors), nsamples, ndetectors

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
    real*8, intent(in)    :: data(nsamples*factor,ndetectors)
    real*8, intent(out)   :: compressed(nsamples,ndetectors)
    integer, intent(in)   :: factor, nsamples, ndetectors

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
    real*8, intent(in)    :: compressed(nsamples,ndetectors)
    real*8, intent(out)   :: data(nsamples*factor,ndetectors)
    integer, intent(in)   :: factor, nsamples, ndetectors

    call transpos(compressed, data, factor)

end subroutine compression_average_transpose


!-------------------------------------------------------------------------------


subroutine apply_mask(signal, mask, nsamples, ndetectors)

    use module_preprocessor, only : apply => apply_mask
    implicit none

    !f2py threadsafe
    !f2py intent(inout) :: signal
    !f2py intent(in)    :: mask
    !f2py intent(hide)  :: nsamples = shape(signal,0)
    !f2py intent(hide)  :: ndetectors = shape(signal,1)
    real*8, intent(inout) :: signal(nsamples, ndetectors)
    logical*1, intent(in) :: mask(nsamples, ndetectors)
    integer, intent(in)   :: nsamples, ndetectors

    call apply(signal, mask)

end subroutine apply_mask
