! Tamasis interface for f2py
!
! Author: P. Chanial

subroutine test_broken_locale(ok)

    implicit none

    !f2py intent(out)      :: ok

    logical*1, intent(out) :: ok
    character(len=70)      :: svalue
    double precision       :: value

    svalue="0.5"
    read (svalue,*) value
    ok = value == 0.5

end subroutine test_broken_locale


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine info_version(version)

    use module_tamasis, only : tamasis_version
    implicit none

    !f2py intent(out)              :: version

    character(len=80), intent(out) :: version

    version = tamasis_version

end subroutine info_version


!-----------------------------------------------------------------------------------------------------------------------------------


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
    !f2py intent(in)       :: map1d
    !f2py intent(inout)    :: signal
    !f2py intent(in)       :: npixels_per_sample
    !f2py intent(hide)     :: nsamples = shape(signal,0)
    !f2py intent(hide)     :: ndetectors = shape(signal,1)
    !f2py intent(hide)     :: npixels = size(map1d)

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
    !f2py intent(in)       :: signal
    !f2py intent(inout)    :: map1d
    !f2py intent(in)       :: npixels_per_sample
    !f2py intent(hide)     :: nsamples = shape(signal,0)
    !f2py intent(hide)     :: ndetectors = shape(signal,1)
    !f2py intent(hide)     :: npixels = size(map1d)

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
    !f2py intent(out)                    :: ptp
    !f2py intent(in)                     :: npixels_per_sample
    !f2py intent(in)                     :: nsamples
    !f2py intent(in)                     :: ndetectors
    !f2py intent(in)                     :: npixels

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
    !f2py intent(in)     :: data
    !f2py intent(inout)  :: compressed
    !f2py intent(in)     :: factor
    !f2py intent(hide)   :: nsamples = shape(data,0)
    !f2py intent(hide)   :: ndetectors = shape(data,1)

    integer, intent(in)  :: factor, nsamples, ndetectors
    real(p), intent(in)  :: data(nsamples,ndetectors)
    real(p), intent(out) :: compressed(nsamples/factor,ndetectors)

    call direct(data, compressed, factor)

end subroutine compression_average_direct


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine compression_average_transpose(compressed, data, factor, nsamples, ndetectors)

    use module_compression, only : transpos => compression_average_transpose
    use module_tamasis,     only : p
    implicit none

    !f2py threadsafe
    !f2py intent(in)     :: compressed
    !f2py intent(inout)  :: data
    !f2py intent(in)     :: factor
    !f2py intent(hide)   :: nsamples = shape(compressed,0)
    !f2py intent(hide)   :: ndetectors = shape(compressed,1)

    integer, intent(in)  :: factor, nsamples, ndetectors
    real(p), intent(in)  :: compressed(nsamples,ndetectors)
    real(p), intent(out) :: data(nsamples*factor,ndetectors)

    call transpos(compressed, data, factor)

end subroutine compression_average_transpose


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine downsampling_direct(data, compressed, factor, nsamples, ndetectors)

    use module_compression, only : direct => downsampling_direct
    use module_tamasis,     only : p
    implicit none

    !f2py threadsafe
    !f2py intent(in)     :: data
    !f2py intent(inout)  :: compressed
    !f2py intent(in)     :: factor
    !f2py intent(hide)   :: nsamples = shape(data,0)/factor
    !f2py intent(hide)   :: ndetectors = shape(data,1)

    integer, intent(in)  :: factor, nsamples, ndetectors
    real(p), intent(in)  :: data(nsamples*factor,ndetectors)
    real(p), intent(out) :: compressed(nsamples,ndetectors)

    call direct(data, compressed, factor)

end subroutine downsampling_direct


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine downsampling_transpose(compressed, data, factor, nsamples, ndetectors)

    use module_compression, only : transpos => downsampling_transpose
    use module_tamasis,     only : p
    implicit none

    !f2py threadsafe
    !f2py intent(inout)  :: compressed
    !f2py intent(in)     :: data
    !f2py intent(in)     :: factor
    !f2py intent(hide)   :: nsamples = shape(data,0)/factor
    !f2py intent(hide)   :: ndetectors = shape(data,1)

    integer, intent(in)  :: factor, nsamples, ndetectors
    real(p), intent(in)  :: compressed(nsamples,ndetectors)
    real(p), intent(out) :: data(nsamples*factor,ndetectors)

    call transpos(compressed, data, factor)

end subroutine downsampling_transpose


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine backprojection_weighted(pmatrix, data, mask, map1d, weight1d, npixels_per_sample, nsamples, ndetectors, npixels)

    use module_pointingmatrix, only : bpw => backprojection_weighted, PointingElement
    use module_tamasis,        only : p
    implicit none

    !f2py threadsafe
    !f2py integer*8,intent(in)        :: pmatrix(npixels_per_sample*nsamples*ndetectors)
    !f2py intent(in)                  :: data
    !f2py intent(in)                  :: mask
    !f2py intent(inout)               :: map1d
    !f2py intent(inout)               :: weight1d
    !f2py intent(in)                  :: npixels_per_sample
    !f2py intent(hide)                :: nsamples = shape(data,0)
    !f2py intent(hide)                :: ndetectors = shape(data,1)
    !f2py intent(hide)                :: npixels = size(map1d)

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
    !f2py integer*8, intent(in)       :: pmatrix(npixels_per_sample*nsamples*ndetectors)
    !f2py intent(in)                  :: nx, ny
    !f2py intent(in)                  :: data
    !f2py intent(inout)               :: mask
    !f2py intent(in)                  :: nsigma
    !f2py intent(in)                  :: npixels_per_sample
    !f2py intent(hide)                :: nsamples = shape(data,0)
    !f2py intent(hide)                :: ndetectors = shape(data,1)

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
    !f2py integer*8, intent(in)       :: pmatrix(npixels_per_sample*nsamples*ndetectors)
    !f2py intent(in)                  :: nx, ny
    !f2py intent(in)                  :: data
    !f2py intent(inout)               :: mask
    !f2py intent(in)                  :: nsigma
    !f2py intent(in)                  :: npixels_per_sample
    !f2py intent(hide)                :: nsamples = shape(data,0)
    !f2py intent(hide)                :: ndetectors = shape(data,1)

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


subroutine filter_median(data, mask, length, nsamples, nsamples_tot, ndetectors, nslices, status)

    use iso_fortran_env,     only : ERROR_UNIT, OUTPUT_UNIT
    use module_preprocessor, only : median_filtering
    use module_tamasis,      only : p, info_time
    use module_string,       only : strinteger
    implicit none

    !f2py threadsafe
    !f2py intent(in)       :: data
    !f2py intent(in)       :: mask
    !f2py intent(in)       :: length
    !f2py intent(in)       :: nsamples
    !f2py intent(hide)     :: nsamples_tot=shape(data,0)
    !f2py intent(hide)     :: ndetectors=shape(data,1)
    !f2py intent(hide)     :: nslices=size(nsamples)
    !f2py intent(out)      :: status

    real(p), intent(inout) :: data(nsamples_tot,ndetectors)
    logical*1, intent(in)  :: mask(nsamples_tot,ndetectors)
    integer, intent(in)    :: length
    integer, intent(in)    :: nsamples(nslices)
    integer, intent(in)    :: nsamples_tot
    integer, intent(in)    :: nslices
    integer, intent(in)    :: ndetectors
    integer, intent(out)   :: status

    integer :: islice, idetector, start
    integer :: count_start

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

    call system_clock(count_start)
    start = 1
    do islice = 1, nslices
        call median_filtering(data(start:start+nsamples(islice)-1,:), mask(start:start+nsamples(islice)-1,:), length)
        start = start + nsamples(islice)
    end do

    ! the filtered timeline has NaN only if it is completely masked
    do idetector = 1, ndetectors
        if (data(1,idetector) /= data(1,idetector)) then
            data(:,idetector) = 0
        end if
    end do

    call info_time('Median filtering (length=' // strinteger(length) // ')', count_start)

end subroutine filter_median


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine fft_filter_uncorrelated(data, nsamples, nsamples_tot, ncorrelations, ndetectors, nslices, tod_filter, status)

    use iso_fortran_env,  only : ERROR_UNIT
    use module_filtering, only : FilterUncorrelated, fft_filter => create_filter_uncorrelated
    use module_tamasis,   only : p
    implicit none

    !f2py threadsafe
    !f2py intent(in)     :: data
    !f2py intent(in)     :: nsamples
    !f2py intent(in)     :: nsamples_tot
    !f2py intent(hide)   :: ncorrelations = shape(data,0) - 1
    !f2py intent(hide)   :: ndetectors = shape(data,1)
    !f2py intent(hide)   :: nslices = shape(data,2)
    !f2py intent(out)    :: tod_filter
    !f2py intent(out)    :: status

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
    !f2py intent(inout)    :: data
    !f2py intent(in)       :: nsamples
    !f2py intent(hide)     :: nslices = size(nsamples)
    !f2py intent(in)       :: plan
    !f2py intent(hide)     :: nsamples_tot = shape(data,0)
    !f2py intent(hide)     :: ndetectors = shape(data,1)

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
    !f2py intent(in)       :: input
    !f2py intent(hide)     :: nvalids = size(input)
    !f2py intent(in)       :: mask
    !f2py intent(hide)     :: nx=shape(mask,0)
    !f2py intent(hide)     :: ny=shape(mask,1)
    !f2py intent(inout)    :: output(nx,ny)
    !f2py intent(in)       :: field

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
    !f2py intent(in)      :: input
    !f2py intent(in)      :: mask(nx,ny)
    !f2py intent(hide)    :: nx=shape(input,0)
    !f2py intent(hide)    :: ny=shape(input,1)
    !f2py intent(hide)    :: nvalids = size(output)
    !f2py intent(inout)   :: output(nvalids)

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


subroutine add_inplace(a, b, n)

    use module_tamasis,  only : p
    use module_math, only : add
    implicit none

    !f2py threadsafe
    !f2py intent(inout)    :: a
    !f2py intent(hide)     :: n = size(a)
    !f2py intent(in)       :: b

    real(p), intent(inout) :: a(n)
    real(p), intent(in)    :: b(n)
    integer, intent(in)    :: n

    call add(a, b, n)

end subroutine add_inplace


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine subtract_inplace(a, b, n)

    use module_tamasis,  only : p
    implicit none

    !f2py threadsafe
    !f2py intent(inout)    :: a
    !f2py intent(hide)     :: n = size(a)
    !f2py intent(in)       :: b

    real(p), intent(inout) :: a(n)
    real(p), intent(in)    :: b(n)
    integer, intent(in)    :: n
    integer                :: i

    !$omp parallel do
    do i = 1, n
        a(i) = a(i) - b(i)
    end do
    !$omp end parallel do

end subroutine subtract_inplace


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine multiply_inplace(a, b, n)

    use module_tamasis,  only : p
    implicit none

    !f2py threadsafe
    !f2py intent(inout)    :: a
    !f2py intent(hide)     :: n = size(a)
    !f2py intent(in)       :: b

    real(p), intent(inout) :: a(n)
    real(p), intent(in)    :: b(n)
    integer, intent(in)    :: n
    integer                :: i

    !$omp parallel do
    do i = 1, n
        a(i) = a(i) * b(i)
    end do
    !$omp end parallel do

end subroutine multiply_inplace


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine divide_inplace(a, b, n)

    use module_tamasis,  only : p
    implicit none

    !f2py threadsafe
    !f2py intent(inout)    :: a
    !f2py intent(hide)     :: n = size(a)
    !f2py intent(in)       :: b

    real(p), intent(inout) :: a(n)
    real(p), intent(in)    :: b(n)
    integer, intent(in)    :: n
    integer                :: i

    !$omp parallel do
    do i = 1, n
        a(i) = a(i) / b(i)
    end do
    !$omp end parallel do

end subroutine divide_inplace


!!$!-----------------------------------------------------------------------------------------------------------------------------------
!!$
!!$
!!$subroutine add_inplace_blas(a, b, n)
!!$
!!$    use module_tamasis,  only : p
!!$    implicit none
!!$
!!$    !f2py threadsafe
!!$    !f2py intent(inout)    :: a
!!$    !f2py intent(hide)     :: n = size(a)
!!$    !f2py intent(in)       :: b
!!$
!!$    real(p), intent(inout) :: a(n)
!!$    real(p), intent(in)    :: b(n)
!!$    integer, intent(in)    :: n
!!$
!!$    call daxpy(n, 1._p, b, 1, a, 1)
!!$
!!$end subroutine add_inplace_blas
!!$
!!$
!!$!-----------------------------------------------------------------------------------------------------------------------------------
!!$
!!$
!!$subroutine subtract_inplace_blas(a, b, n)
!!$
!!$    use module_tamasis,  only : p
!!$    implicit none
!!$
!!$    !f2py threadsafe
!!$    !f2py intent(inout)    :: a
!!$    !f2py intent(hide)     :: n = size(a)
!!$    !f2py intent(in)       :: b
!!$
!!$    real(p), intent(inout) :: a(n)
!!$    real(p), intent(in)    :: b(n)
!!$    integer, intent(in)    :: n
!!$
!!$    call daxpy(n, -1._p, b, 1, a, 1)
!!$
!!$end subroutine subtract_inplace_blas
!!$
!!$
!!$!-----------------------------------------------------------------------------------------------------------------------------------


subroutine diff(array, asize, dim, ashape, arank)

    use module_math,    only : diff_fast, diff_slow
    use module_tamasis, only : p
    implicit none

    !f2py intent(inout)    :: array(asize)
    !f2py intent(hide)     :: asize
    !f2py intent(in)       :: dim
    !f2py intent(in)       :: ashape
    !f2py intent(hide)     :: arank=size(ashape)

    integer*8, intent(in)  :: asize
    real(p), intent(inout) :: array(asize)
    integer, intent(in)    :: dim
    integer, intent(in)    :: arank
    integer*8, intent(in)  :: ashape(arank)

    integer :: idim, nfast, ndiff, nslow

    nfast = 1
    nslow = 1
    do idim = 1, dim-1
        nfast = nfast * ashape(idim)
    end do
    ndiff = ashape(dim)
    do idim = dim+1, arank
        nslow = nslow * ashape(idim)
    end do

    if (dim == 1) then
        call diff_fast(array(1), ndiff, nslow)
    else
        call diff_slow(array(1), nfast, ndiff, nslow)
    end if

end subroutine diff


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine diffT(array, asize, dim, ashape, arank)

    use module_math,    only : diffT_fast, diffT_slow
    use module_tamasis, only : p
    implicit none

    !f2py intent(inout)    :: array(size)
    !f2py intent(hide)     :: asize
    !f2py intent(in)       :: dim
    !f2py intent(in)       :: ashape
    !f2py intent(hide)     :: arank=size(ashape)

    integer*8, intent(in)  :: asize
    real(p), intent(inout) :: array(asize)
    integer, intent(in)    :: dim
    integer, intent(in)    :: arank
    integer*8, intent(in)  :: ashape(arank)

    integer :: idim, nfast, ndiff, nslow

    nfast = 1
    nslow = 1
    do idim = 1, dim-1
        nfast = nfast * ashape(idim)
    end do
    ndiff = ashape(dim)
    do idim = dim+1, arank
        nslow = nslow * ashape(idim)
    end do

    if (dim == 1) then
        call diffT_fast(array(1), ndiff, nslow)
    else
        call diffT_slow(array(1), nfast, ndiff, nslow)
    end if

end subroutine diffT


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine diffTdiff(array, asize, dim, ashape, arank, scalar)

    use module_math,    only : diffTdiff_fast, diffTdiff_slow
    use module_tamasis, only : p
    implicit none

    !f2py intent(inout)    :: array(asize)
    !f2py intent(hide)     :: asize
    !f2py intent(in)       :: dim
    !f2py intent(in)       :: ashape
    !f2py intent(hide)     :: arank=size(ashape)
    !f2py intent(in)       :: scalar

    integer*8, intent(in)  :: asize
    real(p), intent(inout) :: array(asize)
    integer, intent(in)    :: dim
    integer, intent(in)    :: arank
    integer*8, intent(in)  :: ashape(arank)
    real(p), intent(in)    :: scalar

    integer :: idim, nfast, ndiff, nslow

    nfast = 1
    nslow = 1
    do idim = 1, dim-1
        nfast = nfast * ashape(idim)
    end do
    ndiff = ashape(dim)
    do idim = dim+1, arank
        nslow = nslow * ashape(idim)
    end do

    if (dim == 1) then
        call diffTdiff_fast(array(1), ndiff, nslow, scalar)
    else
        call diffTdiff_slow(array(1), nfast, ndiff, nslow, scalar)
    end if

end subroutine diffTdiff


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine shift(array, asize, dim, ashape, arank, offset, offsetsize)

    use module_math,    only : shift_fast, shift_slow
    use module_tamasis, only : p
    implicit none

    !f2py intent(inout)    :: array(asize)
    !f2py intent(hide)     :: asize
    !f2py intent(in)       :: dim
    !f2py intent(in)       :: ashape
    !f2py intent(hide)     :: arank=size(ashape)
    !f2py intent(in)       :: offset(offsetsize)
    !f2py intent(hide)     :: offsetsise=size(offset)

    integer*8, intent(in)  :: asize
    real(p), intent(inout) :: array(asize)
    integer, intent(in)    :: dim
    integer, intent(in)    :: arank
    integer*8, intent(in)  :: ashape(arank)
    integer, intent(in)    :: offsetsize
    integer, intent(in)    :: offset(offsetsize)

    integer :: idim, nfast, nshift, nslow

    nfast = 1
    nslow = 1
    do idim = 1, dim-1
        nfast = nfast * ashape(idim)
    end do
    nshift = ashape(dim)
    do idim = dim+1, arank
        nslow = nslow * ashape(idim)
    end do

    if (dim == 1) then
        call shift_fast(array(1), nshift, nslow, offset)
    else
        call shift_slow(array(1), nfast, nshift, nslow, offset)
    end if

end subroutine shift


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine masking(input, ninputs, mask, nmasks, status)

    use iso_fortran_env, only : ERROR_UNIT
    use module_tamasis,  only : p
    implicit none

    !f2py threadsafe
    !f2py intent(inout)    :: input
    !f2py intent(hide)     :: ninputs=size(input)
    !f2py intent(in)       :: mask(nx)
    !f2py intent(hide)     :: nmasks=size(mask)
    !f2py intent(out)      :: status

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


subroutine remove_nan(signal, mask, nsamples, ndetectors)

    use module_preprocessor, only : remove_nan_ => remove_nan
    use module_tamasis,      only : p
    implicit none

    !f2py threadsafe
    !f2py intent(inout)      :: signal
    !f2py intent(in)         :: mask
    !f2py intent(hide)       :: nsamples = shape(signal,0)
    !f2py intent(hide)       :: ndetectors = shape(signal,1)

    integer, intent(in)      :: nsamples
    integer, intent(in)      :: ndetectors
    real(p), intent(inout)   :: signal(nsamples, ndetectors)
    logical*1, intent(inout) :: mask(nsamples, ndetectors)

    call remove_nan_(signal, mask)

end subroutine remove_nan


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine mean_degrees(array, n, mean)

    use module_math,    only : mean_degrees_ => mean_degrees
    use module_tamasis, only : p
    implicit none

    !f2py threadsafe
    !f2py intent(in)     :: array
    !f2py intent(hide)   :: n = size(array)
    !f2py intent(out)    :: mean

    real(p), intent(in)  :: array(n)
    integer, intent(in)  :: n
    real(p), intent(out) :: mean

    mean = mean_degrees_(array)

end subroutine mean_degrees


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine minmax_degrees(array, n, min, max)

    use module_math,    only : minmax_degrees_ => minmax_degrees
    use module_tamasis, only : p
    implicit none

    !f2py threadsafe
    !f2py intent(in)     :: array
    !f2py intent(hide)   :: n = size(array)
    !f2py intent(out)    :: min, max

    real(p), intent(in)  :: array(n)
    integer, intent(in)  :: n
    real(p), intent(out) :: min, max

    call minmax_degrees_(array, min, max)

end subroutine minmax_degrees


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine barycenter_lonlat(lon, lat, n, lon0, lat0)

    use module_math,    only : barycenter_lonlat_ => barycenter_lonlat
    use module_tamasis, only : p
    implicit none

    !f2py threadsafe
    !f2py intent(in)     :: lon, lat
    !f2py intent(hide)   :: n = size(lon)
    !f2py intent(out)    :: lon0, lat0

    real(p), intent(in)  :: lon(n), lat(n)
    integer, intent(in)  :: n
    real(p), intent(out) :: lon0, lat0

    call barycenter_lonlat_(lon, lat, lon0, lat0)

end subroutine barycenter_lonlat


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine angle_lonlat(lon1, lat1, m, lon2, lat2, n, angle)

    use module_math,    only : angle_lonlat_ => angle_lonlat
    use module_tamasis, only : p
    implicit none

    !f2py threadsafe
    !f2py intent(in)     :: lon1, lat1
    !f2py intent(hide)   :: m = size(lon1)
    !f2py intent(in)     :: lon2, lat2
    !f2py intent(hide)   :: n = size(lon2)
    !f2py intent(out)    :: angle(max(m,n))

    real(p), intent(in)  :: lon1(m), lat1(m)
    integer, intent(in)  :: m
    real(p), intent(in)  :: lon2(n), lat2(n)
    integer, intent(in)  :: n
    real(p), intent(out) :: angle(max(m,n))

    if (m == 1) then
        call angle_lonlat_(lon1(1), lat1(1), lon2, lat2, angle)
    else if (n == 1) then
        call angle_lonlat_(lon1, lat1, lon2(1), lat2(1), angle)
    else
        call angle_lonlat_(lon1, lat1, lon2, lat2, angle)
    end if

end subroutine angle_lonlat


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine distance_1d(nx, origin, resolution, array)

    use module_math,    only : distance
    use module_tamasis, only : p
    implicit none

    !f2py threadsafe
    !f2py intent(in)     :: nx
    !f2py intent(in)     :: origin
    !f2py intent(in)     :: resolution
    !f2py intent(out)    :: array(nx)

    integer, intent(in)  :: nx
    real(p), intent(in)  :: origin
    real(p), intent(in)  :: resolution
    real(p), intent(out) :: array(nx)

    call distance(array, origin, resolution)

end subroutine distance_1d


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine distance_2d(nx, ny, origin, resolution, array)

    use module_math,    only : distance
    use module_tamasis, only : p
    implicit none

    !f2py threadsafe
    !f2py intent(in)     :: nx
    !f2py intent(in)     :: ny
    !f2py intent(in)     :: origin(2)
    !f2py intent(in)     :: resolution(2)
    !f2py intent(out)    :: array(nx,ny)

    integer, intent(in)  :: nx, ny
    real(p), intent(in)  :: origin(2)
    real(p), intent(in)  :: resolution(2)
    real(p), intent(out) :: array(nx,ny)

    call distance(array, origin, resolution)

end subroutine distance_2d


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine distance_3d(nx, ny, nz, origin, resolution, array)

    use module_math,    only : distance
    use module_tamasis, only : p
    implicit none

    !f2py threadsafe
    !f2py intent(in)     :: nx
    !f2py intent(in)     :: ny
    !f2py intent(in)     :: nz
    !f2py intent(in)     :: origin(3)
    !f2py intent(in)     :: resolution(3)
    !f2py intent(out)    :: array(nx,ny,nz)

    integer, intent(in)  :: nx, ny, nz
    real(p), intent(in)  :: origin(3)
    real(p), intent(in)  :: resolution(3)
    real(p), intent(out) :: array(nx,ny,nz)

    call distance(array, origin, resolution)

end subroutine distance_3d


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine profile_axisymmetric_2d(array, nx, ny, origin, bin, nbins, x, y, n)

    use module_math,    only : profile
    use module_tamasis, only : p
    implicit none

    !f2py threadsafe
    !f2py intent(in)     :: array
    !f2py intent(hide)   :: nx = shape(array,0)
    !f2py intent(hide)   :: ny = shape(array,1)
    !f2py intent(in)     :: origin(2)
    !f2py intent(in)     :: bin
    !f2py intent(in)     :: nbins
    !f2py intent(out)    :: x, y, n

    real(p), intent(in)  :: array(nx,ny)
    integer, intent(in)  :: nx, ny
    real(p), intent(in)  :: origin(2)
    real(p), intent(in)  :: bin
    integer, intent(in)  :: nbins
    real(p), intent(out) :: x(nbins)
    real(p), intent(out) :: y(nbins)
    integer, intent(out) :: n(nbins)

    call profile(array, origin, bin, nbins, x, y, n)

end subroutine profile_axisymmetric_2d

!-----------------------------------------------------------------------------------------------------------------------------------


subroutine projection_scale(header, nx, ny, array, status)

    use module_tamasis, only : p
    use module_wcs,     only : Astrometry, init_astrometry, projection_scale_ => projection_scale
    implicit none

    !f2py intent(in)             :: header
    !f2py intent(in)             :: nx
    !f2py intent(in)             :: ny
    !f2py intent(out)            :: array
    !f2py intent(out)            :: status

    character(len=*), intent(in) :: header
    integer, intent(in)          :: nx, ny
    real(p), intent(out)         :: array(nx,ny)
    integer, intent(out)         :: status

    type(Astrometry) :: astr

    call init_astrometry(header, astr, status)
    if (status /= 0) return

    call projection_scale_(astr, array, status)
    if (status /= 0) return

end subroutine projection_scale


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine convolution_trexp_direct(data, nsamples, nslices, tau, nsamples_tot, ndetectors)

    use module_filtering, only : trexp_direct => convolution_trexp_direct
    use module_tamasis,   only : p
    implicit none

    !f2py threadsafe
    !f2py intent(inout)    :: data
    !f2py intent(in)       :: nsamples
    !f2py intent(hide)     :: nslices = size(nsamples)
    !f2py intent(in)       :: tau
    !f2py intent(hide)     :: nsamples_tot = shape(data,0)
    !f2py intent(hide)     :: ndetectors = shape(data,1)

    real(p), intent(inout) :: data(nsamples_tot,ndetectors)
    integer*8, intent(in)  :: nsamples(nslices)
    integer, intent(in)    :: nslices
    real(p), intent(in)    :: tau(ndetectors)
    integer, intent(in)    :: nsamples_tot
    integer, intent(in)    :: ndetectors

    integer                :: islice, dest

    dest = 1
    do islice = 1, nslices

        call trexp_direct(data(dest:dest+nsamples(islice)-1,:), tau)
        dest = dest + nsamples(islice)

    end do

end subroutine convolution_trexp_direct


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine convolution_trexp_transpose(data, nsamples, nslices, tau, nsamples_tot, ndetectors)

    use module_filtering, only : trexp_transpose => convolution_trexp_transpose
    use module_tamasis,   only : p
    implicit none

    !f2py threadsafe
    !f2py intent(inout)    :: data
    !f2py intent(in)       :: nsamples
    !f2py intent(hide)     :: nslices = size(nsamples)
    !f2py intent(in)       :: tau
    !f2py intent(hide)     :: nsamples_tot = shape(data,0)
    !f2py intent(hide)     :: ndetectors = shape(data,1)

    real(p), intent(inout) :: data(nsamples_tot,ndetectors)
    integer*8, intent(in)  :: nsamples(nslices)
    integer, intent(in)    :: nslices
    real(p), intent(in)    :: tau(ndetectors)
    integer, intent(in)    :: nsamples_tot
    integer, intent(in)    :: ndetectors

    integer                :: islice, dest

    dest = 1
    do islice = 1, nslices

        call trexp_transpose(data(dest:dest+nsamples(islice)-1,:), tau)
        dest = dest + nsamples(islice)

    end do

end subroutine convolution_trexp_transpose
