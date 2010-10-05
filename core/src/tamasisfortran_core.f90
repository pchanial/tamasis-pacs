! Tamasis interface for f2py
!
! Author: P. Chanial

subroutine info_nbytes_real(nbytes)

    implicit none

    !f2py intent(out)    :: nbytes

    integer, intent(out) :: nbytes
    
    nbytes = PRECISION_REAL

end subroutine info_nbytes_real


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine info_nthreads(nthreads)

    use omp_lib, only : omp_get_max_threads
    implicit none

    !f2py intent(out)    :: nthreads

    integer, intent(out) :: nthreads
    
    nthreads = omp_get_max_threads()

end subroutine info_nthreads


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pointing_matrix_direct(pmatrix, map1d, signal, npixels_per_sample, nsamples, ndetectors, npixels)

    use module_pointingmatrix, only : PointingElement, pmatrix_direct
    use module_tamasis,        only : p
    implicit none

    !f2py threadsafe
    !f2py integer*8, dimension(npixels_per_sample*nsamples*ndetectors), intent(inout) :: pmatrix
    !f2py intent(in)      :: map1d
    !f2py intent(inout)   :: signal
    !f2py intent(in)      :: npixels_per_sample
    !f2py intent(hide)    :: nsamples = shape(signal,0)
    !f2py intent(hide)    :: ndetectors = shape(signal,1)
    !f2py intent(hide)    :: npixels = size(map1d)

    type(PointingElement), intent(inout) :: pmatrix(npixels_per_sample, nsamples, ndetectors)
    real(p), intent(in)    :: map1d(npixels)
    real(p), intent(inout) :: signal(nsamples, ndetectors)
    integer, intent(in)    :: npixels_per_sample
    integer*8, intent(in)  :: nsamples
    integer, intent(in)    :: ndetectors
    integer, intent(in)    :: npixels

    call pmatrix_direct(pmatrix, map1d, signal)

end subroutine pointing_matrix_direct


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pointing_matrix_transpose(pmatrix, signal, map1d, npixels_per_sample, nsamples, ndetectors, npixels)

    use module_pointingmatrix, only : PointingElement, pmatrix_transpose
    use module_tamasis,        only : p
    implicit none

    !f2py threadsafe
    !f2py integer*8, dimension(npixels_per_sample*nsamples*ndetectors), intent(inout) :: pmatrix
    !f2py intent(in)    :: signal
    !f2py intent(inout) :: map1d
    !f2py intent(in)    :: npixels_per_sample
    !f2py intent(hide)  :: nsamples = shape(signal,0)
    !f2py intent(hide)  :: ndetectors = shape(signal,1)
    !f2py intent(hide)  :: npixels = size(map1d)

    type(PointingElement), intent(inout) :: pmatrix(npixels_per_sample, nsamples, ndetectors)
    real(p), intent(in)    :: signal(nsamples, ndetectors)
    real(p), intent(inout) :: map1d(npixels)
    integer, intent(in)    :: npixels_per_sample
    integer*8, intent(in)  :: nsamples
    integer, intent(in)    :: ndetectors
    integer, intent(in)    :: npixels

    call pmatrix_transpose(pmatrix, signal, map1d)

end subroutine pointing_matrix_transpose


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pointing_matrix_ptp(pmatrix, ptp, npixels_per_sample, nsamples, ndetectors, npixels)

    use module_pointingmatrix, only : PointingElement, pmatrix_ptp
    use module_tamasis,        only : p
    implicit none

    !f2py threadsafe
    !f2py integer*8, dimension(npixels_per_sample*nsamples*ndetectors), intent(inout) :: pmatrix
    !f2py intent(out) :: ptp
    !f2py intent(in)  :: npixels_per_sample
    !f2py intent(in)  :: nsamples
    !f2py intent(in)  :: ndetectors
    !f2py intent(in)  :: npixels

    type(PointingElement), intent(inout) :: pmatrix(npixels_per_sample, nsamples, ndetectors)
    real(p), intent(out)                 :: ptp(npixels, npixels)
    integer, intent(in)                  :: npixels_per_sample
    integer*8, intent(in)                :: nsamples
    integer, intent(in)                  :: ndetectors
    integer, intent(in)                  :: npixels

    call pmatrix_ptp(pmatrix, ptp)

end subroutine pointing_matrix_ptp


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine compression_average_direct(data, compressed, factor, nsamples, ndetectors)

    use module_compression, only : direct => compression_average_direct
    use module_tamasis,     only : p
    implicit none

    !f2py threadsafe
    !f2py intent(in)      :: data
    !f2py intent(inout)   :: compressed
    !f2py intent(in)      :: factor
    !f2py intent(hide)    :: nsamples = shape(data,0)/factor
    !f2py intent(hide)    :: ndetectors = shape(data,1)

    integer, intent(in)   :: factor, nsamples, ndetectors
    real(p), intent(in)   :: data(nsamples*factor,ndetectors)
    real(p), intent(out)  :: compressed(nsamples,ndetectors)

    call direct(data, compressed, factor)

end subroutine compression_average_direct


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine compression_average_transpose(compressed, data, factor, nsamples, ndetectors)

    use module_compression, only : transpos => compression_average_transpose
    use module_tamasis,     only : p
    implicit none

    !f2py threadsafe
    !f2py intent(inout)   :: compressed
    !f2py intent(in)      :: data
    !f2py intent(in)      :: factor
    !f2py intent(hide)    :: nsamples = shape(data,0)/factor
    !f2py intent(hide)    :: ndetectors = shape(data,1)

    integer, intent(in)   :: factor, nsamples, ndetectors
    real(p), intent(in)   :: compressed(nsamples,ndetectors)
    real(p), intent(out)  :: data(nsamples*factor,ndetectors)

    call transpos(compressed, data, factor)

end subroutine compression_average_transpose


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine downsampling_direct(data, compressed, factor, nsamples, ndetectors)

    use module_compression, only : direct => downsampling_direct
    use module_tamasis,     only : p
    implicit none

    !f2py threadsafe
    !f2py intent(in)      :: data
    !f2py intent(inout)   :: compressed
    !f2py intent(in)      :: factor
    !f2py intent(hide)    :: nsamples = shape(data,0)/factor
    !f2py intent(hide)    :: ndetectors = shape(data,1)

    integer, intent(in)   :: factor, nsamples, ndetectors
    real(p), intent(in)   :: data(nsamples*factor,ndetectors)
    real(p), intent(out)  :: compressed(nsamples,ndetectors)

    call direct(data, compressed, factor)

end subroutine downsampling_direct


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine downsampling_transpose(compressed, data, factor, nsamples, ndetectors)

    use module_compression, only : transpos => downsampling_transpose
    use module_tamasis,     only : p
    implicit none

    !f2py threadsafe
    !f2py intent(inout)   :: compressed
    !f2py intent(in)      :: data
    !f2py intent(in)      :: factor
    !f2py intent(hide)    :: nsamples = shape(data,0)/factor
    !f2py intent(hide)    :: ndetectors = shape(data,1)

    integer, intent(in)   :: factor, nsamples, ndetectors
    real(p), intent(in)   :: compressed(nsamples,ndetectors)
    real(p), intent(out)  :: data(nsamples*factor,ndetectors)

    call transpos(compressed, data, factor)

end subroutine downsampling_transpose


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine backprojection_weighted(pmatrix, data, mask, map1d, weight1d, npixels_per_sample, nsamples, ndetectors, npixels)

    use module_pointingmatrix, only : bpw => backprojection_weighted, PointingElement
    use module_tamasis,        only : p
    implicit none

    !f2py threadsafe
    !f2py integer*8,intent(in) :: pmatrix(npixels_per_sample*nsamples*ndetectors)
    !f2py intent(in)           :: data
    !f2py intent(in)           :: mask
    !f2py intent(inout)        :: map1d
    !f2py intent(inout)        :: weight1d
    !f2py intent(in)           :: npixels_per_sample
    !f2py intent(hide)         :: nsamples = shape(data,0)
    !f2py intent(hide)         :: ndetectors = shape(data,1)
    !f2py intent(hide)         :: npixels = size(map1d)

    type(PointingElement), intent(in) :: pmatrix(npixels_per_sample,nsamples,ndetectors)
    real(p), intent(in)               :: data(nsamples,ndetectors)
    logical*1, intent(in)             :: mask(nsamples,ndetectors)
    real(p), intent(inout)            :: map1d(npixels)
    real(p), intent(inout)            :: weight1d(npixels)
    integer, intent(in)               :: npixels_per_sample
    integer*8, intent(in)             :: nsamples
    integer, intent(in)               :: ndetectors
    integer, intent(in)               :: npixels

    call bpw(pmatrix, data, mask, map1d, weight1d)

end subroutine backprojection_weighted


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine deglitch_l2b_std(pmatrix, nx, ny, data, mask, nsigma, npixels_per_sample, nsamples, ndetectors)

    use module_pointingmatrix, only : PointingElement
    use module_deglitching,    only : deglitch_l2b
    use module_tamasis,        only : p
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

    type(PointingElement), intent(in) :: pmatrix(npixels_per_sample,nsamples,ndetectors)
    integer, intent(in)               :: nx, ny
    real(p), intent(in)               :: data(nsamples,ndetectors)
    logical*1, intent(inout)          :: mask(nsamples,ndetectors)
    real(p), intent(in)               :: nsigma
    integer, intent(in)               :: npixels_per_sample
    integer*8, intent(in)             :: nsamples
    integer, intent(in)               :: ndetectors

    call deglitch_l2b(pmatrix, nx, ny, data, mask, nsigma, .false., verbose=.true.)

end subroutine deglitch_l2b_std


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine deglitch_l2b_mad(pmatrix, nx, ny, data, mask, nsigma, npixels_per_sample, nsamples, ndetectors)

    use module_pointingmatrix, only : PointingElement
    use module_deglitching,    only : deglitch_l2b
    use module_tamasis,        only : p
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

    type(PointingElement), intent(in) :: pmatrix(npixels_per_sample,nsamples,ndetectors)
    integer, intent(in)               :: nx, ny
    real(p), intent(in)               :: data(nsamples,ndetectors)
    logical*1, intent(inout)          :: mask(nsamples,ndetectors)
    real(p), intent(in)               :: nsigma
    integer, intent(in)               :: npixels_per_sample
    integer*8, intent(in)             :: nsamples
    integer, intent(in)               :: ndetectors

    call deglitch_l2b(pmatrix, nx, ny, data, mask, nsigma, .true., verbose=.true.)

end subroutine deglitch_l2b_mad


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine filter_median(data, mask, length, nsamples, nsamples_tot, nslices, ndetectors, status)

    use iso_fortran_env,     only : ERROR_UNIT, OUTPUT_UNIT
    use module_preprocessor, only : median_filtering
    use module_tamasis,      only : p
    implicit none

    !f2py threadsafe
    !f2py intent(in)   :: data
    !f2py intent(in)   :: mask
    !f2py intent(in)   :: length
    !f2py intent(in)   :: nsamples
    !f2py intent(hide) :: nsamples_tot=shape(data,0)
    !f2py intent(hide) :: nslices=size(nsamples)
    !f2py intent(hide) :: ndetectors=shape(data,1)
    !f2py intent(out)  :: status

    real(p), intent(inout) :: data(nsamples_tot,ndetectors)
    logical*1, intent(in)  :: mask(nsamples_tot,ndetectors)
    integer, intent(in)    :: length
    integer, intent(in)    :: nsamples(nslices)
    integer, intent(in)    :: nsamples_tot
    integer, intent(in)    :: nslices
    integer, intent(in)    :: ndetectors
    integer, intent(out)   :: status

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
        call median_filtering(data(start:start+nsamples(islice)-1,:), mask(start:start+nsamples(islice)-1,:), length)
        start = start + nsamples(islice)
    end do
    call system_clock(count2, count_rate, count_max)
    write (OUTPUT_UNIT,'(f7.2,a)') real(count2-count1)/count_rate, 's'

end subroutine filter_median


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine fft_filter_uncorrelated(data, nsamples, nsamples_tot, ncorrelations, ndetectors, nslices, tod_filter, status)

    use iso_fortran_env,  only : ERROR_UNIT
    use module_filtering, only : FilterUncorrelated, fft_filter => create_filter_uncorrelated
    use module_tamasis,   only : p
    implicit none

    !f2py threadsafe
    !f2py intent(in)   :: data
    !f2py intent(in)   :: nsamples
    !f2py intent(in)   :: nsamples_tot
    !f2py intent(hide) :: ncorrelations = size(data,0) - 1
    !f2py intent(hide) :: ndetectors = size(data,1)
    !f2py intent(hide) :: nslices = size(data,2)
    !f2py intent(out)  :: tod_filter
    !f2py intent(out)  :: status

    integer, intent(in)  :: nsamples_tot
    integer, intent(in)  :: ncorrelations
    integer, intent(in)  :: ndetectors
    integer, intent(in)  :: nslices
    real(p), intent(in)  :: data(ncorrelations+1,ndetectors,nslices)
    integer, intent(in)  :: nsamples(nslices)
    real(p), intent(out) :: tod_filter(nsamples_tot,ndetectors)
    integer, intent(out) :: status

    type(FilterUncorrelated), allocatable :: filter(:)
    integer                               :: islice

    allocate (filter(nslices))
    filter%ncorrelations = ncorrelations
    filter%bandwidth = 2 * ncorrelations + 1
    filter%ndetectors = ndetectors
    do islice = 1, nslices
        allocate (filter(islice)%data(ncorrelations+1,ndetectors))
        filter(islice)%data = data(:,:,islice)
    end do

    call fft_filter(filter, nsamples, ndetectors, tod_filter, status)
    if (status /= 0) return

end subroutine fft_filter_uncorrelated


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine fft_plan(data, nsamples, nslices, plan, nsamples_tot, ndetectors)

    use module_filtering, only : fft_tod
    use module_tamasis,   only : p
    implicit none

    !f2py threadsafe
    !f2py intent(inout) :: data
    !f2py intent(in)    :: nsamples
    !f2py intent(hide)  :: nslices = size(nsamples)
    !f2py intent(in)    :: plan
    !f2py intent(hide)  :: nsamples_tot = size(data,2)
    !f2py intent(hide)  :: ndetectors = size(data,1)

    real(p), intent(inout) :: data(nsamples_tot,ndetectors)
    integer*8, intent(in)  :: nsamples(nslices)
    integer, intent(in)    :: nslices
    integer*8, intent(in)  :: plan(nslices)
    integer, intent(in)    :: nsamples_tot
    integer, intent(in)    :: ndetectors

    integer                :: islice, dest

    dest = 1
    do islice = 1, nslices
       
        call fft_tod(plan(islice), data(dest:dest+nsamples(islice)-1,:))
        dest = dest + nsamples(islice)

    end do

end subroutine fft_plan


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine unpack_direct(input, nvalids, mask, nx, ny, output, field)

    use iso_fortran_env, only : ERROR_UNIT
    use module_tamasis,  only : p
    implicit none

    !f2py threadsafe
    !f2py intent(in)    :: input
    !f2py intent(hide)  :: nvalids = size(input)
    !f2py intent(in)    :: mask
    !f2py intent(hide)  :: nx=shape(mask,0)
    !f2py intent(hide)  :: ny=shape(mask,1)
    !f2py intent(inout) :: output(nx,ny)
    !f2py intent(in)    :: field

    real(p), intent(in)    :: input(nvalids)
    integer, intent(in)    :: nvalids
    logical*1, intent(in)  :: mask(nx,ny)
    integer, intent(in)    :: nx, ny
    real(p), intent(out)   :: output(nx,ny)
    real(p), intent(in)    :: field

    if (count(.not. mask) /= nvalids) then
        write (ERROR_UNIT,'(a)') 'UNPACK_DIRECT: The mask is not compatible with the input size.'
        return
    endif
    output = unpack(input, .not. mask, field)

end subroutine unpack_direct


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine unpack_transpose(input, mask, nx, ny, nvalids, output)

    use iso_fortran_env, only : ERROR_UNIT
    use module_tamasis,  only : p
    implicit none

    !f2py threadsafe
    !f2py intent(in)    :: input
    !f2py intent(in)    :: mask(nx,ny)
    !f2py intent(hide)  :: nx=shape(input,0)
    !f2py intent(hide)  :: ny=shape(input,1)
    !f2py intent(hide)  :: nvalids = size(output)
    !f2py intent(inout) :: output(nvalids)

    real(p), intent(in)   :: input(nx,ny)
    logical*1, intent(in) :: mask(nx,ny)
    integer, intent(in)   :: nx, ny
    integer, intent(in)   :: nvalids
    real(p), intent(out)  :: output(nvalids)

    if (count(.not. mask) /= nvalids) then
        write (ERROR_UNIT,'(a)') 'UNPACK_TRANSPOSE: The mask is not compatible with the output size.'
        return
    endif
    output = pack(input, .not. mask)

end subroutine unpack_transpose


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine masking(input, ninputs, mask, nmasks, status)

    use iso_fortran_env, only : ERROR_UNIT
    use module_tamasis,  only : p
    implicit none

    !f2py threadsafe
    !f2py intent(inout) :: input
    !f2py intent(hide)  :: ninputs=size(input)
    !f2py intent(in)    :: mask(nx)
    !f2py intent(hide)  :: nmasks=size(mask)
    !f2py intent(out)   :: status

    real(p), intent(inout) :: input(ninputs)
    logical*1, intent(in)  :: mask(nmasks)
    integer, intent(in)    :: ninputs, nmasks
    integer, intent(out)   :: status

    if (ninputs /= nmasks) then
        write (ERROR_UNIT,'(a,2(i0,a))') "The data array has a size incompatible with the mask ('", ninputs, "' instead of '",     &
              nmasks, "')."
        status = 1
        return
    end if

    !$omp parallel workshare
    where (mask)
        input = 0
    end where
    !$omp end parallel workshare

    status = 0

end subroutine masking


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine mean_degrees(array, n, mean)

    use module_math,    only : mean_degrees_ => mean_degrees
    use module_tamasis, only : p
    implicit none

    !f2py threadsafe
    !f2py intent(in)   :: array
    !f2py intent(hide) :: n = size(array)
    !f2py intent(out)  :: mean

    real(p), intent(in)  :: array(n)
    integer, intent(in)  :: n
    real(p), intent(out) :: mean

    mean = mean_degrees_(array)

end subroutine mean_degrees

