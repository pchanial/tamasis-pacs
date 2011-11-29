! Copyright 2010-2011 Pierre Chanial
! All rights reserved
!
! Tamasis interface for f2py
!
! Author: P. Chanial

subroutine test_broken_locale(ok)

    implicit none

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

    character(len=80), intent(out) :: version

    version = tamasis_version

end subroutine info_version


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine info_nbytes_real(nbytes)

    implicit none

    integer, intent(out) :: nbytes

    nbytes = PRECISION_REAL

end subroutine info_nbytes_real


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pointing_matrix_direct(pmatrix, map1d, signal, npixels_per_sample, nsamples, ndetectors, npixels)

    use module_pointingmatrix, only : PointingElement, pmatrix_direct, pmatrix_direct_one_pixel_per_sample
    use module_tamasis,        only : p
    implicit none

    !f2py integer*8, dimension(npixels_per_sample*nsamples*ndetectors), intent(inout) :: pmatrix

    type(PointingElement), intent(inout) :: pmatrix(npixels_per_sample, nsamples, ndetectors)
    real(p), intent(in)    :: map1d(npixels)
    real(p), intent(inout) :: signal(nsamples, ndetectors)
    integer, intent(in)    :: npixels_per_sample
    integer*8, intent(in)  :: nsamples
    integer, intent(in)    :: ndetectors
    integer, intent(in)    :: npixels

    if (npixels_per_sample == 1) then
        call pmatrix_direct_one_pixel_per_sample(pmatrix(1,:,:), map1d, signal)
    else
        call pmatrix_direct(pmatrix, map1d, signal)
    end if

end subroutine pointing_matrix_direct


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pointing_matrix_transpose(pmatrix, signal, map1d, npixels_per_sample, nsamples, ndetectors, npixels)

    use module_pointingmatrix, only : PointingElement, pmatrix_transpose, pmatrix_transpose_one_pixel_per_sample
    use module_tamasis,        only : p
    implicit none

    !f2py integer*8, dimension(npixels_per_sample*nsamples*ndetectors), intent(inout) :: pmatrix

    type(PointingElement), intent(inout) :: pmatrix(npixels_per_sample, nsamples, ndetectors)
    real(p), intent(in)    :: signal(nsamples, ndetectors)
    real(p), intent(inout) :: map1d(npixels)
    integer, intent(in)    :: npixels_per_sample
    integer*8, intent(in)  :: nsamples
    integer, intent(in)    :: ndetectors
    integer, intent(in)    :: npixels

    if (npixels_per_sample == 1) then
        call pmatrix_transpose_one_pixel_per_sample(pmatrix(1,:,:), signal, map1d)
    else
        call pmatrix_transpose(pmatrix, signal, map1d)
    end if

end subroutine pointing_matrix_transpose


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pointing_matrix_ptp(pmatrix, ptp, npixels_per_sample, nsamples, ndetectors, npixels)

    use module_pointingmatrix, only : PointingElement, pmatrix_ptp
    use module_tamasis,        only : p
    implicit none

    !f2py integer*8, dimension(npixels_per_sample*nsamples*ndetectors), intent(inout) :: pmatrix

    type(PointingElement), intent(inout) :: pmatrix(npixels_per_sample, nsamples, ndetectors)
    real(p), intent(out)                 :: ptp(npixels, npixels)
    integer, intent(in)                  :: npixels_per_sample
    integer*8, intent(in)                :: nsamples
    integer, intent(in)                  :: ndetectors
    integer, intent(in)                  :: npixels

    call pmatrix_ptp(pmatrix, ptp)

end subroutine pointing_matrix_ptp


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pointing_matrix_mask(pmatrix, mask1d, npixels_per_sample, nsamples, ndetectors, npixels)

    use module_pointingmatrix, only : PointingElement, pmatrix_mask
    use module_tamasis,        only : p
    implicit none

    !f2py integer*8, dimension(npixels_per_sample*nsamples*ndetectors), intent(in) :: pmatrix

    type(PointingElement), intent(in) :: pmatrix(npixels_per_sample, nsamples, ndetectors)
    logical*1, intent(inout)          :: mask1d(npixels)
    integer, intent(in)               :: npixels_per_sample
    integer*8, intent(in)             :: nsamples
    integer, intent(in)               :: ndetectors
    integer, intent(in)               :: npixels

    call pmatrix_mask(pmatrix, mask1d)

end subroutine pointing_matrix_mask


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pointing_matrix_pack(pmatrix, mask1d, npixels_per_sample, nsamples, ndetectors, npixels)

    use module_pointingmatrix, only : PointingElement, pmatrix_pack
    use module_tamasis,        only : p
    implicit none

    !f2py integer*8, dimension(npixels_per_sample*nsamples*ndetectors), intent(inout) :: pmatrix

    type(PointingElement), intent(inout) :: pmatrix(npixels_per_sample, nsamples, ndetectors)
    logical*1, intent(in)                :: mask1d(npixels)
    integer, intent(in)                  :: npixels_per_sample
    integer*8, intent(in)                :: nsamples
    integer, intent(in)                  :: ndetectors
    integer, intent(in)                  :: npixels

    call pmatrix_pack(pmatrix, mask1d)

end subroutine pointing_matrix_pack


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine compression_average_direct(input, ninputs, ndetectors, isize, istride, output, noutputs, osize, ostride, factor)

    use module_compression, only : func => compression_average_direct
    use module_tamasis,     only : p
    implicit none

    real(p), intent(in)    :: input(ninputs)
    real(p), intent(inout) :: output(noutputs)
    integer, intent(in)    :: ninputs, ndetectors, isize, istride
    integer, intent(in)    :: noutputs, osize, ostride
    integer, intent(in)    :: factor

    integer :: i

    !$omp parallel do
    do i = 1, ndetectors
        call func(input((i-1)*istride+1:(i-1)*istride+isize), output((i-1)*ostride+1:(i-1)*ostride+osize), factor)
    end do
    !$omp end parallel do

end subroutine compression_average_direct


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine compression_average_transpose(input, ninputs, ndetectors, isize, istride, output, noutputs, osize, ostride, factor)

    use module_compression, only : func => compression_average_transpose
    use module_tamasis,     only : p
    implicit none

    real(p), intent(in)    :: input(ninputs)
    real(p), intent(inout) :: output(noutputs)
    integer, intent(in)    :: ninputs, ndetectors, isize, istride
    integer, intent(in)    :: noutputs, osize, ostride
    integer, intent(in)    :: factor

    integer :: i

    !$omp parallel do
    do i = 1, ndetectors
        call func(input((i-1)*istride+1:(i-1)*istride+isize), output((i-1)*ostride+1:(i-1)*ostride+osize), factor)
    end do
    !$omp end parallel do

end subroutine compression_average_transpose


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine downsampling_direct(input, ninputs, ndetectors, isize, istride, output, noutputs, osize, ostride, factor)

    use module_compression, only : func => downsampling_direct
    use module_tamasis,     only : p
    implicit none

    real(p), intent(in)    :: input(ninputs)
    real(p), intent(inout) :: output(noutputs)
    integer, intent(in)    :: ninputs, ndetectors, isize, istride
    integer, intent(in)    :: noutputs, osize, ostride
    integer, intent(in)    :: factor

    integer :: i

    !$omp parallel do
    do i = 1, ndetectors
        call func(input((i-1)*istride+1:(i-1)*istride+isize), output((i-1)*ostride+1:(i-1)*ostride+osize), factor)
    end do
    !$omp end parallel do

end subroutine downsampling_direct


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine downsampling_transpose(input, ninputs, ndetectors, isize, istride, output, noutputs, osize, ostride, factor)

    use module_compression, only : func => downsampling_transpose
    use module_tamasis,     only : p
    implicit none

    real(p), intent(in)    :: input(ninputs)
    real(p), intent(inout) :: output(noutputs)
    integer, intent(in)    :: ninputs, ndetectors, isize, istride
    integer, intent(in)    :: noutputs, osize, ostride
    integer, intent(in)    :: factor

    integer :: i

    !$omp parallel do
    do i = 1, ndetectors
        call func(input((i-1)*istride+1:(i-1)*istride+isize), output((i-1)*ostride+1:(i-1)*ostride+osize), factor)
    end do
    !$omp end parallel do

end subroutine downsampling_transpose


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine backprojection_weighted(pmatrix, data, mask, map1d, weight1d, npixels_per_sample, nsamples, ndetectors, npixels)

    use module_pointingmatrix, only : bpw => backprojection_weighted, PointingElement
    use module_tamasis,        only : p
    implicit none

    !f2py integer*8,intent(in)        :: pmatrix(npixels_per_sample*nsamples*ndetectors)

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

    !f2py integer*8, intent(in)       :: pmatrix(npixels_per_sample*nsamples*ndetectors)

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

    !f2py integer*8, intent(in)       :: pmatrix(npixels_per_sample*nsamples*ndetectors)

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


subroutine fft_filter_uncorrelated(data, nsamples, ncorrelations, ndetectors, fft_filter, status)

    use iso_fortran_env,  only : ERROR_UNIT
    use module_filtering, only : FilterUncorrelated, create_filter_uncorrelated
    use module_tamasis,   only : p
    implicit none

    real(p), intent(in)  :: data(ncorrelations+1,ndetectors)
    integer, intent(in)  :: nsamples
    integer, intent(in)  :: ncorrelations
    integer, intent(in)  :: ndetectors
    real(p), intent(out) :: fft_filter(nsamples,ndetectors)
    integer, intent(out) :: status

    type(FilterUncorrelated) :: filter(1)

    filter%ncorrelations = ncorrelations
    filter%bandwidth = 2 * ncorrelations + 1
    filter%ndetectors = ndetectors
    allocate (filter(1)%data(ncorrelations+1,ndetectors))
    filter(1)%data = data(:,:)

    call create_filter_uncorrelated(filter, [nsamples], ndetectors, fft_filter, status)
    if (status /= 0) return

end subroutine fft_filter_uncorrelated


!-----------------------------------------------------------------------------------------------------------------------------------

subroutine fft_plan_inplace(input, ninputs, ndetectors, isize, istride, plan)

    use module_tamasis, only : p
    implicit none

    real(p), intent(in)    :: input(ninputs)
    integer, intent(in)    :: ninputs, ndetectors, isize, istride
    integer*8, intent(in)  :: plan

    integer :: i

    !$omp parallel do
    do i = 1, ndetectors
        call dfftw_execute_r2r(plan, input((i-1)*istride+1:(i-1)*istride+isize), input((i-1)*istride+1:(i-1)*istride+isize))
    end do
    !$omp end parallel do

end subroutine fft_plan_inplace


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine invntt_uncorrelated(input, ninputs, nsamples, istride, output, noutputs, ostride, fft_filter, filter_length, ndetectors,&
                               fplan, bplan, left, right)

    use module_tamasis, only : p
    implicit none

    real(p), intent(in)    :: input(ninputs)
    integer, intent(in)    :: ninputs
    integer, intent(in)    :: nsamples, istride
    real(p), intent(inout) :: output(noutputs)
    integer, intent(in)    :: noutputs
    integer, intent(in)    :: ostride
    integer, intent(in)    :: filter_length, ndetectors
    real(p), intent(in)    :: fft_filter(filter_length, ndetectors)
    integer*8, intent(in)  :: fplan, bplan
    integer, intent(in)    :: left, right

    real(p) :: buffer(filter_length)

    integer :: i

    !$omp parallel do private(buffer)
    do i = 1, ndetectors
        buffer(:left) = 0
        buffer(left+1:filter_length-right) = input((i-1)*istride+1:(i-1)*istride+nsamples)
        buffer(filter_length-right+1:) = 0
        call dfftw_execute_r2r(fplan, buffer, buffer)
        buffer = buffer * fft_filter(:,i)
        call dfftw_execute_r2r(bplan, buffer, buffer)
        output((i-1)*ostride+1:(i-1)*ostride+nsamples) = buffer(left+1:filter_length-right)
    end do
    !$omp end parallel do

end subroutine invntt_uncorrelated


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine fft_plan_outplace(input, ninputs, ndetectors, isize, istride, output, noutputs, osize, ostride, plan)

    use module_tamasis, only : p
    implicit none

    real(p), intent(in)    :: input(ninputs)
    real(p), intent(inout) :: output(noutputs)
    integer, intent(in)    :: ninputs, ndetectors, isize, istride
    integer, intent(in)    :: noutputs, osize, ostride
    integer*8, intent(in)  :: plan

    integer :: i

    !$omp parallel do
    do i = 1, ndetectors
        call dfftw_execute_r2r(plan, input((i-1)*istride+1:(i-1)*istride+isize), output((i-1)*ostride+1:(i-1)*ostride+osize))
    end do
    !$omp end parallel do

end subroutine fft_plan_outplace


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine packing(input, mask, ninputs, output, noutputs)

    use module_tamasis,  only : p
    implicit none

    real(p), intent(in)    :: input(ninputs)
    logical*1, intent(in)  :: mask(ninputs)
    integer, intent(in)    :: ninputs
    real(p), intent(inout) :: output(noutputs)
    integer, intent(in)    :: noutputs

    output = pack(input, .not. mask)

end subroutine packing


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine unpacking(input, ninputs, mask, nmasks, output, noutputs)

    use module_tamasis,  only : p
    implicit none

    real(p), intent(in)    :: input(ninputs)
    logical*1, intent(in)  :: mask(nmasks)
    real(p), intent(inout) :: output(noutputs)
    integer, intent(in)    :: ninputs, nmasks, noutputs

    integer :: ii, io

    ! unlike unpack, the following works in-place
    ii = ninputs
    do io = noutputs, 1, -1
        if (mask(io)) then
            output(io) = 0
        else
            output(io) = input(ii)
            ii = ii - 1
        end if
    end do

end subroutine unpacking


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine multiply_conjugate_complex(a, m, b, n, c)

    use module_tamasis, only : p
    implicit none

    complex(p), intent(in)    :: a(m)
    complex(p), intent(in)    :: b(n)
    complex(p), intent(inout) :: c(m)
    integer, intent(in)       :: m, n
    integer                   :: i, o

    if (m == n) then
        !$omp parallel workshare
        c = a * conjg(b)
        !$omp end parallel workshare
        return
    end if

    o = m / n
    !$omp parallel do
    do i = 1, n
        c((i-1)*o+1:i*o) = a((i-1)*o+1:i*o) * conjg(b(i))
    end do
    !$omp end parallel do

end subroutine multiply_conjugate_complex


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine round_rtz(input, output, n)

    use module_tamasis, only : p
    implicit none

    real(p), intent(in)    :: input(n)
    real(p), intent(inout) :: output(n)
    integer, intent(in)    :: n

    !$omp parallel workshare
    output = aint(input)
    !$omp end parallel workshare

end subroutine round_rtz


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine round_rti(input, output, n)

    use module_tamasis, only : p
    implicit none

    real(p), intent(in)    :: input(n)
    real(p), intent(inout) :: output(n)
    integer, intent(in)    :: n

    !$omp parallel workshare
    output = sign(real(ceiling(abs(input)), kind=p), input)
    !$omp end parallel workshare

end subroutine round_rti


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine round_rtmi(input, output, n)

    use module_tamasis, only : p
    implicit none

    real(p), intent(in)    :: input(n)
    real(p), intent(inout) :: output(n)
    integer, intent(in)    :: n

    !$omp parallel workshare
    output = floor(input)
    !$omp end parallel workshare

end subroutine round_rtmi


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine round_rtpi(input, output, n)

    use module_tamasis, only : p
    implicit none

    real(p), intent(in)    :: input(n)
    real(p), intent(inout) :: output(n)
    integer, intent(in)    :: n

    !$omp parallel workshare
    output = ceiling(input)
    !$omp end parallel workshare

end subroutine round_rtpi


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine round_rhtz(input, output, n)

    use module_tamasis, only : p
    implicit none

    real(p), intent(in)    :: input(n)
    real(p), intent(inout) :: output(n)
    integer, intent(in)    :: n

    integer :: i
    real(p) :: x

    !$omp parallel do private(x)
    do i = 1, n
        x = anint(input(i))
        if (abs(x-input(i)) ==  0.5_p) then
            x = x - sign(1._p, input(i))
        end if
        output(i) = x
    end do
    !$omp end parallel do

end subroutine round_rhtz


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine round_rhti(input, output, n)

    use module_tamasis, only : p
    implicit none

    real(p), intent(in)    :: input(n)
    real(p), intent(inout) :: output(n)
    integer, intent(in)    :: n

    !$omp parallel workshare
    output = anint(input)
    !$omp end parallel workshare

end subroutine round_rhti


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine round_rhtmi(input, output, n)

    use module_math,    only : nint_down
    use module_tamasis, only : p
    implicit none

    real(p), intent(in)    :: input(n)
    real(p), intent(inout) :: output(n)
    integer, intent(in)    :: n

    !$omp parallel workshare
    output = nint_down(input)
    !$omp end parallel workshare

end subroutine round_rhtmi


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine round_rhtpi(input, output, n)

    use module_math,    only : nint_up
    use module_tamasis, only : p
    implicit none

    real(p), intent(in)    :: input(n)
    real(p), intent(inout) :: output(n)
    integer, intent(in)    :: n

    !$omp parallel workshare
    output = nint_up(input)
    !$omp end parallel workshare

end subroutine round_rhtpi


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine round_rhs(input, output, n)

    use module_tamasis, only : p
    implicit none

    real(p), intent(in)    :: input(n)
    real(p), intent(inout) :: output(n)
    integer, intent(in)    :: n

    integer                :: i, s, clock
    integer, allocatable   :: seed(:)
    real(p)                :: x, rnd

    ! initialise the random seed
    call random_seed(size=s)
    allocate(seed(s))
    call system_clock(count=clock)
    seed = clock + 37 * [ (i - 1, i = 1, s) ]
    call random_seed(put=seed)
    deallocate(seed)

    !$omp parallel do private(x)
    do i = 1, n
        x = anint(input(i))
        if (abs(x-input(i)) == 0.5_p) then
            call random_number(rnd)
            if (rnd >= 0.5_p) then
                x = x - sign(1._p, input(i))
            end if
        end if
        output(i) = x
    end do
    !$omp end parallel do

end subroutine round_rhs


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine sum(array, n, output)

    use omp_lib,        only : omp_get_max_threads
    use module_math,    only : sum_kahan
    use module_tamasis, only : p
    implicit none

    real(p), intent(in)  :: array(n)
    integer, intent(in)  :: n
    real(p), intent(out) :: output

    integer :: nthreads, work, i, a, z

    if (n <= 1024) then
        output = sum_kahan(array)
        return
    end if

    output = 0
    nthreads = omp_get_max_threads()
    work = max(n / nthreads, 1)
    !$omp parallel do reduction(+:output) private(a,z)
    do i = 1, min(n, nthreads)
        a = min((i-1)*work+1,n+1)
        z = min(i*work,n)
        output = output + sum_kahan(array(a:z))
    end do
    !$omp end parallel do

end subroutine sum


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine dot(array1, array2, n, m, output)

    use omp_lib,        only : omp_get_max_threads
    use module_math,    only : d => dot
    use module_tamasis, only : p
    implicit none

    real(p), intent(in)  :: array1(n)
    real(p), intent(in)  :: array2(m)
    integer, intent(in)  :: n, m
    real(p), intent(out) :: output

    integer :: nthreads, work, i, a, z

    if (n <= 1024) then
        output = d(array1, array2)
        return
    end if

    output = 0
    nthreads = omp_get_max_threads()
    work = max(n / nthreads, 1)
    !$omp parallel do reduction(+:output) private(a,z)
    do i = 1, min(n, nthreads)
        a = min((i-1)*work+1,n+1)
        z = min(i*work,n)
        output = output + d(array1(a:z), array2(a:z))
    end do
    !$omp end parallel do

end subroutine dot


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine norm_l1(array, n, comm, output)

    use omp_lib,        only : omp_get_max_threads
    use module_math,    only : norm => norm_l1
    use module_tamasis, only : p
#ifdef HAVE_MPI_MODULE
    use mpi
#endif
    implicit none

    real(p), intent(in)  :: array(n)
    integer, intent(in)  :: n
    integer, intent(in)  :: comm
    real(p), intent(out) :: output

    integer :: nthreads, work, i, a, z, status

#ifndef HAVE_MPI_MODULE
    include 'mpif.h'
#endif

    if (n <= 1024) then
        output = norm(array)
    else
        output = 0
        nthreads = omp_get_max_threads()
        work = max(n / nthreads, 1)
        !$omp parallel do reduction(+:output) private(a,z)
        do i = 1, min(n, nthreads)
            a = min((i-1)*work+1,n+1)
            z = min(i*work,n)
            output = output + norm(array(a:z))
        end do
        !$omp end parallel do
    end if

    call MPI_Allreduce(MPI_IN_PLACE, output, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, status)

end subroutine norm_l1


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine dnorm_l1(array, n, output)

    use module_tamasis, only : p
    implicit none

    real(p), intent(in)    :: array(n)
    integer, intent(in)    :: n
    real(p), intent(inout) :: output(n)

    !$omp parallel workshare
    output = sign(1._p, array)
    !$omp end parallel workshare

end subroutine dnorm_l1

    
!-----------------------------------------------------------------------------------------------------------------------------------


subroutine dnorm_l2(array, n, comm, output)

    use module_tamasis, only : p
    implicit none

    real(p), intent(in)    :: array(n)
    integer, intent(in)    :: n
    integer, intent(in)    :: comm
    real(p), intent(inout) :: output(n)

    external :: norm2
    real(p)  :: a

    call norm2(array, n, comm, a)
    a = a ** (-0.5_p)
    !$omp parallel workshare
    output = sign(abs(array) * a, array)
    !$omp end parallel workshare

end subroutine dnorm_l2

    
!-----------------------------------------------------------------------------------------------------------------------------------


subroutine dnorm_lp(array, n, lp, comm, output)

    use module_tamasis, only : p
    implicit none

    real(p), intent(in)    :: array(n)
    integer, intent(in)    :: n
    real(p), intent(in)    :: lp
    integer, intent(in)    :: comm
    real(p), intent(inout) :: output(n)

    external :: normp
    real(p)  :: a

    call normp(array, n, lp, comm, a)
    a = a ** (1 / lp - 1)
    !$omp parallel workshare
    output = sign(abs(array) ** (lp - 1) * a, array)
    !$omp end parallel workshare

end subroutine dnorm_lp


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine norm_linf(array, n, comm, output)

    use module_tamasis, only : p
#ifdef HAVE_MPI_MODULE
    use mpi
#endif
    implicit none

    real(p), intent(in)  :: array(n)
    integer, intent(in)  :: n
    integer, intent(in)  :: comm
    real(p), intent(out) :: output

    integer :: status

#ifndef HAVE_MPI_MODULE
    include 'mpif.h'
#endif

    output = maxval(abs(array))
    call MPI_Allreduce(MPI_IN_PLACE, output, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, status)

end subroutine norm_linf


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine norm_huber(array, n, delta, comm, output)

    use omp_lib,        only : omp_get_max_threads
    use module_math,    only : norm => norm_huber
    use module_tamasis, only : p
#ifdef HAVE_MPI_MODULE
    use mpi
#endif
    implicit none

    real(p), intent(in)  :: array(n)
    integer, intent(in)  :: n
    real(p), intent(in)  :: delta
    integer, intent(in)  :: comm
    real(p), intent(out) :: output

    integer :: nthreads, work, i, a, z, status

#ifndef HAVE_MPI_MODULE
    include 'mpif.h'
#endif

    if (n <= 1024) then
        output = norm(array, delta)
    else
        output = 0
        nthreads = omp_get_max_threads()
        work = max(n / nthreads, 1)
        !$omp parallel do reduction(+:output) private(a,z)
        do i = 1, min(n, nthreads)
            a = min((i-1)*work+1,n+1)
            z = min(i*work,n)
            output = output + norm(array(a:z), delta)
        end do
        !$omp end parallel do
    end if

    call MPI_Allreduce(MPI_IN_PLACE, output, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, status)

end subroutine norm_huber


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine dnorm_huber(array, n, delta, output)

    use module_tamasis, only : p
    implicit none

    real(p), intent(in)    :: array(n)
    integer, intent(in)    :: n
    real(p), intent(in)    :: delta
    real(p), intent(inout) :: output(n)

    !$omp parallel workshare
    where (array >= delta)
        output = 2 * delta
    elsewhere (array <= -delta)
        output = - 2 * delta
    elsewhere
        output = 2 * array
    end where
    !$omp end parallel workshare

end subroutine dnorm_huber


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine norm2(array, n, comm, output)

    use omp_lib,        only : omp_get_max_threads
    use module_math,    only : norm => norm2
    use module_tamasis, only : p
#ifdef HAVE_MPI_MODULE
    use mpi
#endif
    implicit none

    real(p), intent(in)  :: array(n)
    integer, intent(in)  :: n
    integer, intent(in)  :: comm
    real(p), intent(out) :: output

    integer :: nthreads, work, i, a, z, status

#ifndef HAVE_MPI_MODULE
    include 'mpif.h'
#endif

    if (n <= 1024) then
        output = norm(array)
    else
        output = 0
        nthreads = omp_get_max_threads()
        work = max(n / nthreads, 1)
        !$omp parallel do reduction(+:output) private(a,z)
        do i = 1, min(n, nthreads)
            a = min((i-1)*work+1,n+1)
            z = min(i*work,n)
            output = output + norm(array(a:z))
        end do
        !$omp end parallel do
    end if

    call MPI_Allreduce(MPI_IN_PLACE, output, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, status)

end subroutine norm2


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine dnorm2(array, n, output)

    use module_tamasis, only : p
    implicit none

    real(p), intent(in)    :: array(n)
    integer, intent(in)    :: n
    real(p), intent(inout) :: output(n)

    !$omp parallel workshare
    output = 2 * array
    !$omp end parallel workshare

end subroutine dnorm2


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine normp(array, n, lp, comm, output)

    use omp_lib,        only : omp_get_max_threads
    use module_math,    only : norm => normp
    use module_tamasis, only : p
#ifdef HAVE_MPI_MODULE
    use mpi
#endif
    implicit none

    real(p), intent(in)  :: array(n)
    integer, intent(in)  :: n
    real(p), intent(in)  :: lp
    integer, intent(in)  :: comm
    real(p), intent(out) :: output

    integer :: nthreads, work, i, a, z, status

#ifndef HAVE_MPI_MODULE
    include 'mpif.h'
#endif

    if (n <= 1024) then
        output = norm(array, lp)
    else
        output = 0
        nthreads = omp_get_max_threads()
        work = max(n / nthreads, 1)
        !$omp parallel do reduction(+:output) private(a,z)
        do i = 1, min(n, nthreads)
            a = min((i-1)*work+1,n+1)
            z = min(i*work,n)
            output = output + norm(array(a:z), lp)
        end do
        !$omp end parallel do
    end if

    call MPI_Allreduce(MPI_IN_PLACE, output, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, status)

end subroutine normp


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine dnormp(array, n, lp, output)

    use module_tamasis, only : p
    implicit none

    real(p), intent(in)    :: array(n)
    integer, intent(in)    :: n
    real(p), intent(in)    :: lp
    real(p), intent(inout) :: output(n)

    !$omp parallel workshare
    output = lp * sign(abs(array) ** (lp - 1._p), array)
    !$omp end parallel workshare

end subroutine dnormp


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine diff(input, output, asize, dim, ashape, arank)

    use module_math,    only : diff_fast, diff_medium
    use module_tamasis, only : p
    implicit none

    integer*8, intent(in)  :: asize
    real(p), intent(in)    :: input(asize)
    real(p), intent(inout) :: output(asize)
    integer, intent(in)    :: dim
    integer, intent(in)    :: arank
    integer*8, intent(in)  :: ashape(arank)

    integer :: nfast, ndiff, nslow

    nfast = product(ashape(1:dim-1))
    ndiff = ashape(dim)
    nslow = product(ashape(dim+1:arank))

    if (dim == 1) then
        call diff_fast(input(1), output(1), ndiff, nslow)
    else
        call diff_medium(input(1), output(1), nfast, ndiff, nslow)
    end if

end subroutine diff


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine diff_mpi(input, output, asize, dim, ashape, arank, comm, status)

    use module_math,    only : diff_fast, diff_medium, diff_slow
    use module_tamasis, only : p
#ifdef HAVE_MPI_MODULE
    use mpi
#endif
    implicit none

    integer*8, intent(in)  :: asize
    real(p), intent(in)    :: input(asize)
    real(p), intent(inout) :: output(asize)
    integer, intent(in)    :: dim
    integer, intent(in)    :: arank
    integer*8, intent(in)  :: ashape(arank)
    integer, intent(in)    :: comm
    integer, intent(out)   :: status

    integer                :: rank, size, nfast, ndiff, nslow, mpistatus(6)
    real(p), allocatable   :: boundary(:)

#ifndef HAVE_MPI_MODULE
    include 'mpif.h'
#endif

    nfast = product(ashape(1:dim-1))
    ndiff = ashape(dim)
    nslow = product(ashape(dim+1:arank))

    if (dim == arank) then
        call MPI_Comm_rank(comm, rank, status)
        if (status /= 0) return
        call MPI_Comm_size(comm, size, status)
        if (status /= 0) return
        allocate (boundary(nfast))
        call MPI_Sendrecv(input(1:nfast), nfast, MPI_DOUBLE_PRECISION, modulo(rank-1, size), 0,                                    &
                          boundary, nfast, MPI_DOUBLE_PRECISION, modulo(rank+1, size), 0, comm, mpistatus, status)
        if (status /= 0) return
        call diff_slow(input(1), output(1), nfast, ndiff, boundary, rank == size - 1)
    else if (dim > 1) then
        call diff_medium(input(1), output(1), nfast, ndiff, nslow)
    else
        call diff_fast(input(1), output(1), ndiff, nslow)
    end if

end subroutine diff_mpi


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine diffT(input, output, asize, dim, ashape, arank)

    use module_math,    only : diffT_fast, diffT_medium
    use module_tamasis, only : p
    implicit none

    integer*8, intent(in)  :: asize
    real(p), intent(in)    :: input(asize)
    real(p), intent(inout) :: output(asize)
    integer, intent(in)    :: dim
    integer, intent(in)    :: arank
    integer*8, intent(in)  :: ashape(arank)

    integer :: nfast, ndiff, nslow

    nfast = product(ashape(1:dim-1))
    ndiff = ashape(dim)
    nslow = product(ashape(dim+1:arank))

    if (dim == 1) then
        call diffT_fast(input(1), output(1), ndiff, nslow)
    else
        call diffT_medium(input(1), output(1), nfast, ndiff, nslow)
    end if

end subroutine diffT


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine diffT_mpi(input, output, asize, dim, ashape, arank, comm, status)

    use module_math,    only : diffT_fast, diffT_medium, diffT_slow
    use module_tamasis, only : p
#ifdef HAVE_MPI_MODULE
    use mpi
#endif
    implicit none

    integer*8, intent(in)  :: asize
    real(p), intent(in)    :: input(asize)
    real(p), intent(inout) :: output(asize)
    integer, intent(in)    :: dim
    integer, intent(in)    :: arank
    integer*8, intent(in)  :: ashape(arank)
    integer, intent(in)    :: comm
    integer, intent(out)   :: status

    integer                :: rank, size, nfast, ndiff, nslow, mpistatus(6)
    real(p), allocatable   :: boundary(:)

#ifndef HAVE_MPI_MODULE
    include 'mpif.h'
#endif

    nfast = product(ashape(1:dim-1))
    ndiff = ashape(dim)
    nslow = product(ashape(dim+1:arank))

    if (dim == arank) then
        call MPI_Comm_rank(comm, rank, status)
        if (status /= 0) return
        call MPI_Comm_size(comm, size, status)
        if (status /= 0) return
        allocate (boundary(nfast))
        call MPI_Sendrecv(input(asize-nfast+1:asize), nfast, MPI_DOUBLE_PRECISION, modulo(rank+1, size), 0,                        &
                          boundary, nfast, MPI_DOUBLE_PRECISION, modulo(rank-1, size), 0, comm, mpistatus, status)
        if (status /= 0) return
        call diffT_slow(input(1), output(1), nfast, ndiff, boundary, rank == 0, rank == size - 1)
    else if (dim > 1) then
        call diffT_medium(input(1), output(1), nfast, ndiff, nslow)
    else
        call diffT_fast(input(1), output(1), ndiff, nslow)
    end if

end subroutine difft_mpi


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine diffTdiff(input, output, asize, dim, ashape, arank, scalar)

    use module_math,    only : diffTdiff_fast, diffTdiff_medium
    use module_tamasis, only : p
    implicit none

    integer*8, intent(in)  :: asize
    real(p), intent(in)    :: input(asize)
    real(p), intent(inout) :: output(asize)
    integer, intent(in)    :: dim
    integer, intent(in)    :: arank
    integer*8, intent(in)  :: ashape(arank)
    real(p), intent(in)    :: scalar

    integer :: nfast, ndiff, nslow

    nfast = product(ashape(1:dim-1))
    ndiff = ashape(dim)
    nslow = product(ashape(dim+1:arank))

    if (dim == 1) then
        call diffTdiff_fast(input(1), output(1), ndiff, nslow, scalar)
    else
        call diffTdiff_medium(input(1), output(1), nfast, ndiff, nslow, scalar)
    end if

end subroutine diffTdiff


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine diffTdiff_mpi(input, output, asize, dim, ashape, arank, scalar, comm, status)

    use module_math,    only : diffTdiff_fast, diffTdiff_medium, diffTdiff_slow
    use module_tamasis, only : p
#ifdef HAVE_MPI_MODULE
    use mpi
#endif
    implicit none

    integer*8, intent(in)  :: asize
    real(p), intent(in)    :: input(asize)
    real(p), intent(inout) :: output(asize)
    integer, intent(in)    :: dim
    integer, intent(in)    :: arank
    integer*8, intent(in)  :: ashape(arank)
    real(p), intent(in)    :: scalar
    integer, intent(in)    :: comm
    integer, intent(out)   :: status

    integer                :: rank, size, nfast, ndiff, nslow, mpistatus(6)
    real(p), allocatable   :: boundary1(:), boundary2(:)

#ifndef HAVE_MPI_MODULE
    include 'mpif.h'
#endif

    nfast = product(ashape(1:dim-1))
    ndiff = ashape(dim)
    nslow = product(ashape(dim+1:arank))

    if (dim == arank) then
        call MPI_Comm_rank(comm, rank, status)
        if (status /= 0) return
        call MPI_Comm_size(comm, size, status)
        if (status /= 0) return
        allocate (boundary1(nfast), boundary2(nfast))
        call MPI_Sendrecv(input(asize-nfast+1:asize), nfast, MPI_DOUBLE_PRECISION, modulo(rank+1, size), 0,                        &
                          boundary1, nfast, MPI_DOUBLE_PRECISION, modulo(rank-1, size), 0, comm, mpistatus, status)
        if (status /= 0) return
        call MPI_Sendrecv(input(1:nfast), nfast, MPI_DOUBLE_PRECISION, modulo(rank-1, size), 0,                                    &
                          boundary2, nfast, MPI_DOUBLE_PRECISION, modulo(rank+1, size), 0, comm, mpistatus, status)
        if (status /= 0) return
        call diffTdiff_slow(input(1), output(1), nfast, ndiff, scalar, boundary1, boundary2, rank == 0, rank == size - 1)
    else if (dim > 1) then
        call diffTdiff_medium(input(1), output(1), nfast, ndiff, nslow, scalar)
    else
        call diffTdiff_fast(input(1), output(1), ndiff, nslow, scalar)
    end if

end subroutine difftdiff_mpi


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine shift(input, output, asize, dim, ashape, arank, offset, offsetsize)

    use module_math,    only : shift_fast, shift_medium
    use module_tamasis, only : p
    implicit none

    integer*8, intent(in)  :: asize
    real(p), intent(in)    :: input(asize)
    real(p), intent(inout) :: output(asize)
    integer, intent(in)    :: dim
    integer, intent(in)    :: arank
    integer*8, intent(in)  :: ashape(arank)
    integer, intent(in)    :: offsetsize
    integer, intent(in)    :: offset(offsetsize)

    integer :: nfast, nshift, nslow

    nfast = product(ashape(1:dim-1))
    nshift = ashape(dim)
    nslow = product(ashape(dim+1:arank))

    if (dim == 1) then
        call shift_fast(input(1), output(1), nshift, nslow, offset)
    else
        call shift_medium(input(1), output(1), nfast, nshift, nslow, offset)
    end if

end subroutine shift


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine masking_inplace(input, m, mask, n)

    use iso_fortran_env, only : ERROR_UNIT
    use module_tamasis,  only : p
    implicit none

    real(p), intent(inout) :: input(m)
    logical*1, intent(in)  :: mask(n)
    integer, intent(in)    :: m, n

    integer :: i, o

    if (m == n) then
        !$omp parallel workshare
        where (mask)
            input = 0
        end where
        !$omp end parallel workshare
    else if (n == 1) then
        if (mask(1)) then
            !$omp parallel workshare
            input = 0
            !$omp end parallel workshare
        end if
    else
        o = m / n
        !$omp parallel do
        do i = 1, n
            if (mask(i)) input((i-1)*o+1:i*o) = 0
        end do
        !$omp end parallel do
    end if

end subroutine masking_inplace


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine masking_outplace(input, m, mask, n, output)

    use iso_fortran_env, only : ERROR_UNIT
    use module_tamasis,  only : p
    implicit none

    real(p), intent(in)    :: input(m)
    logical*1, intent(in)  :: mask(n)
    integer, intent(in)    :: m, n
    real(p), intent(inout) :: output(m)

    integer :: i, o

    if (m == n) then
        !$omp parallel workshare
        where (mask)
            output = 0
        elsewhere
            output = input
        end where
        !$omp end parallel workshare
    else if (n == 1) then
        if (mask(1)) then
            !$omp parallel workshare
            output = 0
            !$omp end parallel workshare
        else
            !$omp parallel workshare
            output = input
            !$omp end parallel workshare
        end if
    else
        o = m / n
        !$omp parallel do
        do i = 1, n
            if (mask(i)) then
                output((i-1)*o+1:i*o) = 0
            else
                output((i-1)*o+1:i*o) = input((i-1)*o+1:i*o)
            end if
        end do
        !$omp end parallel do
    end if

end subroutine masking_outplace


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine remove_nonfinite(signal, nsamples)

    use module_preprocessor, only : remove_nonfinite_ => remove_nonfinite
    use module_tamasis,      only : p
    implicit none

    integer, intent(in)      :: nsamples
    real(p), intent(inout)   :: signal(nsamples)

    call remove_nonfinite_(signal)

end subroutine remove_nonfinite


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine remove_nonfinite_mask(signal, mask, nsamples)

    use module_preprocessor, only : remove_nonfinite
    use module_tamasis,      only : p
    implicit none

    integer, intent(in)      :: nsamples
    real(p), intent(inout)   :: signal(nsamples)
    logical*1, intent(inout) :: mask(nsamples)

    call remove_nonfinite(signal, mask)

end subroutine remove_nonfinite_mask


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine mean_degrees(array, n, mean)

    use module_math,    only : mean_degrees_ => mean_degrees
    use module_tamasis, only : p
    implicit none

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


subroutine convolution_trexp_direct(input, ninputs, ndetectors, isize, istride, output, noutputs, osize, ostride, tau, ntaus)

    use module_filtering, only : func => convolution_trexp_direct
    use module_tamasis,   only : p
    implicit none

    real(p), intent(in)    :: input(ninputs)
    real(p), intent(inout) :: output(noutputs)
    integer, intent(in)    :: ninputs, ndetectors, isize, istride
    integer, intent(in)    :: noutputs, osize, ostride
    real(p), intent(in)    :: tau(ntaus)
    integer, intent(in)    :: ntaus

    integer :: i

    !$omp parallel do
    do i = 1, ndetectors
        call func(input((i-1)*istride+1:(i-1)*istride+isize), output((i-1)*ostride+1:(i-1)*ostride+osize), tau(i))
    end do
    !$omp end parallel do

end subroutine convolution_trexp_direct


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine convolution_trexp_transpose(input, ninputs, ndetectors, isize, istride, output, noutputs, osize, ostride, tau, ntaus)

    use module_filtering, only : func => convolution_trexp_transpose
    use module_tamasis,   only : p
    implicit none

    real(p), intent(in)    :: input(ninputs)
    real(p), intent(inout) :: output(noutputs)
    integer, intent(in)    :: ninputs, ndetectors, isize, istride
    integer, intent(in)    :: noutputs, osize, ostride
    real(p), intent(in)    :: tau(ntaus)
    integer, intent(in)    :: ntaus

    integer :: i

    !$omp parallel do
    do i = 1, ndetectors
        call func(input((i-1)*istride+1:(i-1)*istride+isize), output((i-1)*ostride+1:(i-1)*ostride+osize), tau(i))
    end do
    !$omp end parallel do

end subroutine convolution_trexp_transpose


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine mpi_allreducelocal(input, ninputs, mask, nmasks, output, noutputs, op, comm, status)
    
    use module_math,    only : mInf, pInf
    use module_tamasis, only : p
#ifdef HAVE_MPI_MODULE 
    use mpi
#endif
    implicit none
 
    real(p), intent(in)    :: input(ninputs)
    logical*1, intent(in)  :: mask(nmasks)
    real(p), intent(inout) :: output(noutputs)
    integer, intent(in)    :: ninputs, nmasks, noutputs, op, comm
    integer, intent(out)   :: status

    integer :: size, root, a, z, iinput, imask
    real(p) :: input_(noutputs), field

#ifndef HAVE_MPI_MODULE
    include 'mpif.h'
#endif

    call MPI_Comm_size(comm, size, status)
    if (status /= 0) return

    select case (op)
        case (MPI_PROD)
            field = 1
        case (MPI_MIN)
           field = pInf
        case (MPI_MAX)
           field = mInf
        case default
           field = 0
    end select
 
    output = 0
    iinput = 1
    do root = 0, size - 1

        ! unpack the input under the control of mask
        a = root * noutputs + 1
        z = (root + 1) * noutputs
        do imask = a, z
            if (imask > nmasks) then
                input_(imask-a+1:) = field
                exit
            end if
            if (mask(imask)) then
                input_(imask-a+1) = field
            else
                input_(imask-a+1) = input(iinput)
                iinput = iinput + 1
            end if
        end do

        ! reduce on node root
        call MPI_Reduce(input_, output, noutputs, MPI_DOUBLE_PRECISION, op, root, comm, status)
        if (status /= 0) return

    end do

end subroutine mpi_allreducelocal


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine mpi_allscatterlocal(input, ninputs, mask, nmasks, output, noutputs, comm, status)
    
    use module_tamasis, only : p
#ifdef HAVE_MPI_MODULE
    use mpi
#endif
    implicit none

    real(p), intent(in)    :: input(ninputs) ! local image of the distributed map
    logical*1, intent(in)  :: mask(nmasks)   ! global mask for the locally observed pixels
    real(p), intent(inout) :: output(noutputs)
    integer, intent(in)    :: ninputs, nmasks, noutputs, comm
    integer, intent(out)   :: status

    integer :: a_unpacked, z_unpacked, tag, nrecvs, nsends, mpistatus(6), ioutput, imask
    integer :: source, dest, dp, rank, size
    integer, allocatable :: a_packed(:), z_packed(:)

    logical*1 :: maskbuffer(ninputs)
    real(p)   :: databuffer(min(ninputs,noutputs))

#ifndef HAVE_MPI_MODULE
    include 'mpif.h'
#endif

    call MPI_Comm_size(comm, size, status)
    if (status /= 0) return
    call MPI_Comm_rank(comm, rank, status)
    if (status /= 0) return

    ! before exchanging data, compute the local map bounds in the packed output
    allocate (a_packed(0:size-1), z_packed(0:size-1))
    a_packed(0) = 1
    z_packed(0) = 0
    source = 0
    do
        do imask = source * ninputs + 1, min((source+1) * ninputs, nmasks)
            if (mask(imask)) cycle
            z_packed(source) = z_packed(source) + 1
        end do
        if (source == size - 1) exit
        source = source + 1
        a_packed(source) = min(z_packed(source-1) + 1, noutputs)
        z_packed(source) = a_packed(source) - 1
    end do

    ! loop over the distances between mpi nodes, logically organised in a ring
    ioutput = 1
    tag = 99
    do dp = 0, size - 1

        dest = modulo(rank + dp, size)
        source = modulo(rank - dp, size)

        ! send mask to 'dest' and receive it from source
        a_unpacked = min(dest * ninputs + 1, nmasks+1)
        z_unpacked = min((dest+1) * ninputs, nmasks)
        call MPI_Sendrecv(mask(a_unpacked:z_unpacked), z_unpacked-a_unpacked+1, MPI_BYTE, dest, tag,                               &
                          maskbuffer, ninputs, MPI_BYTE, source, tag, comm, mpistatus, status)
        if (status /= 0) return
        call MPI_Get_count(mpistatus, MPI_BYTE, nrecvs, status)
        if (status /= 0) return

        ! pack pixels observed by 'source' from local map of the current processor
        nsends = 0
        do imask = 1, nrecvs
            if (maskbuffer(imask)) cycle
            nsends = nsends + 1
            databuffer(nsends) = input(imask)
        end do

        ! send data to 'source' and receive them from 'dest'
        call MPI_Sendrecv(databuffer, nsends, MPI_DOUBLE_PRECISION, source, tag + 1, output(a_packed(dest)),                       &
             z_packed(dest) - a_packed(dest) + 1, MPI_DOUBLE_PRECISION, dest, tag + 1, comm, mpistatus, status)
        if (status /= 0) return
        call MPI_Get_count(mpistatus, MPI_DOUBLE_PRECISION, nrecvs, status)
        if (status /= 0) return

    end do

end subroutine mpi_allscatterlocal
