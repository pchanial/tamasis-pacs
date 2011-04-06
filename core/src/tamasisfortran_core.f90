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


subroutine info_nthreads(nthreads)

    use omp_lib, only : omp_get_max_threads
    implicit none

    integer, intent(out) :: nthreads

    nthreads = omp_get_max_threads()

end subroutine info_nthreads


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pointing_matrix_direct(pmatrix, map1d, signal, npixels_per_sample, nsamples, ndetectors, npixels)

    use module_pointingmatrix, only : PointingElement, pmatrix_direct
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

    call pmatrix_direct(pmatrix, map1d, signal)

end subroutine pointing_matrix_direct


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pointing_matrix_transpose(pmatrix, signal, map1d, npixels_per_sample, nsamples, ndetectors, npixels)

    use module_pointingmatrix, only : PointingElement, pmatrix_transpose
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

    call pmatrix_transpose(pmatrix, signal, map1d)

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
    integer, intent(in)               :: nsamples
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
    integer, intent(in)                  :: nsamples
    integer, intent(in)                  :: ndetectors
    integer, intent(in)                  :: npixels

    call pmatrix_pack(pmatrix, mask1d)

end subroutine pointing_matrix_pack


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine compression_average_direct(data, compressed, nsamples, nslices, factor, nfactors, nsamples_data, nsamples_compressed,   &
                                      ndetectors)

    use module_compression, only : direct => compression_average_direct
    use module_tamasis,     only : p
    implicit none

    integer, intent(in)    :: nfactors, nslices
    integer, intent(in)    :: factor(nfactors), nsamples(nslices), nsamples_data, nsamples_compressed, ndetectors
    real(p), intent(in)    :: data(nsamples_data,ndetectors)
    real(p), intent(inout) :: compressed(nsamples_compressed,ndetectors)

    integer :: i, i2, j, j2, islice, f

    i = 1
    j = 1
    islice = 1
    f = factor(1)
    do
        i2 = i
        j2 = j
        do
            i2 = i2 + nsamples(islice)
            j2 = j2 + nsamples(islice) / f
            if (islice == nslices) exit
            if (factor(min(islice+1, nfactors)) /= f) exit
            islice = islice + 1
        end do
        call direct(data(i:i2-1,:), compressed(j:j2-1,:), f)
        if (islice == nslices) exit
        i = i2
        j = j2
        islice = islice + 1
        f = factor(islice)
    end do

end subroutine compression_average_direct


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine compression_average_transpose(compressed, data, nsamples, nslices, factor, nfactors, nsamples_compressed, nsamples_data,&
                                         ndetectors)

    use module_compression, only : transpose => compression_average_transpose
    use module_tamasis,     only : p
    implicit none

    integer, intent(in)    :: nfactors, nslices
    integer, intent(in)    :: factor(nfactors), nsamples(nslices), nsamples_data, nsamples_compressed, ndetectors
    real(p), intent(in)    :: compressed(nsamples_compressed,ndetectors)
    real(p), intent(inout) :: data(nsamples_data,ndetectors)

    integer :: i, i2, j, j2, islice, f

    i = 1
    j = 1
    islice = 1
    f = factor(1)
    do
        i2 = i
        j2 = j
        do
            i2 = i2 + nsamples(islice)
            j2 = j2 + nsamples(islice) * f
            if (islice == nslices) exit
            if (factor(min(islice+1, nfactors)) /= f) exit
            islice = islice + 1
        end do
        call transpose(compressed(i:i2-1,:), data(j:j2-1,:), f)
        if (islice == nslices) exit
        i = i2
        j = j2
        islice = islice + 1
        f = factor(islice)
    end do

end subroutine compression_average_transpose


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine downsampling_direct(data, compressed, nsamples, nslices, factor, nfactors, nsamples_data, nsamples_compressed,          &
                               ndetectors)

    use module_compression, only : direct => downsampling_direct
    use module_tamasis,     only : p
    implicit none

    integer, intent(in)    :: nfactors, nslices
    integer, intent(in)    :: factor(nfactors), nsamples(nslices), nsamples_data, nsamples_compressed, ndetectors
    real(p), intent(in)    :: data(nsamples_data,ndetectors)
    real(p), intent(inout) :: compressed(nsamples_compressed,ndetectors)

    integer :: i, i2, j, j2, islice, f

    i = 1
    j = 1
    islice = 1
    f = factor(1)
    do
        i2 = i
        j2 = j
        do
            i2 = i2 + nsamples(islice)
            j2 = j2 + nsamples(islice) / f
            if (islice == nslices) exit
            if (factor(min(islice+1, nfactors)) /= f) exit
            islice = islice + 1
        end do
        call direct(data(i:i2-1,:), compressed(j:j2-1,:), f)
        if (islice == nslices) exit
        i = i2
        j = j2
        islice = islice + 1
        f = factor(islice)
    end do

end subroutine downsampling_direct


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine downsampling_transpose(compressed, data, nsamples, nslices, factor, nfactors, nsamples_compressed, nsamples_data,       &
                                  ndetectors)

    use module_compression, only : transpose => downsampling_transpose
    use module_tamasis,     only : p
    implicit none

    integer, intent(in)    :: nfactors, nslices
    integer, intent(in)    :: factor(nfactors), nsamples(nslices), nsamples_data, nsamples_compressed, ndetectors
    real(p), intent(in)    :: compressed(nsamples_compressed,ndetectors)
    real(p), intent(inout) :: data(nsamples_data,ndetectors)

    integer :: i, i2, j, j2, islice, f

    i = 1
    j = 1
    islice = 1
    f = factor(1)
    do
        i2 = i
        j2 = j
        do
            i2 = i2 + nsamples(islice)
            j2 = j2 + nsamples(islice) * f
            if (islice == nslices) exit
            if (factor(min(islice+1, nfactors)) /= f) exit
            islice = islice + 1
        end do
        call transpose(compressed(i:i2-1,:), data(j:j2-1,:), f)
        if (islice == nslices) exit
        i = i2
        j = j2
        islice = islice + 1
        f = factor(islice)
    end do

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


subroutine fft_filter_uncorrelated(data, nsamples, nsamples_tot, ncorrelations, ndetectors, nslices, tod_filter, status)

    use iso_fortran_env,  only : ERROR_UNIT
    use module_filtering, only : FilterUncorrelated, fft_filter => create_filter_uncorrelated
    use module_tamasis,   only : p
    implicit none

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

    real(p), intent(in)    :: input(nvalids)
    integer, intent(in)    :: nvalids
    logical*1, intent(in)  :: mask(nx,ny)
    integer, intent(in)    :: nx, ny
    real(p), intent(inout) :: output(nx,ny)
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

    real(p), intent(in)    :: input(nx,ny)
    logical*1, intent(in)  :: mask(nx,ny)
    integer, intent(in)    :: nx, ny
    integer, intent(in)    :: nvalids
    real(p), intent(inout) :: output(nvalids)

    if (count(.not. mask) /= nvalids) then
        write (ERROR_UNIT,'(a)') 'UNPACK_TRANSPOSE: The mask is not compatible with the output size.'
        return
    endif
    output = pack(input, .not. mask)

end subroutine unpack_transpose


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine add_inplace(a, m, b, n)

    use module_tamasis,  only : p
    implicit none

    real(p), intent(inout) :: a(m)
    real(p), intent(in)    :: b(n)
    integer, intent(in)    :: m, n
    integer                :: i, o

    if (m == n) then
        !$omp parallel do
        do i = 1, m
            a(i) = a(i) + b(i)
        end do
        !$omp end parallel do
        return
    end if

    o = m / n
    !$omp parallel do
    do i = 1, n
        a((i-1)*o+1:i*o) = a((i-1)*o+1:i*o) + b(i)
    end do
    !$omp end parallel do

end subroutine add_inplace


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine subtract_inplace(a, m, b, n)

    use module_tamasis,  only : p
    implicit none

    real(p), intent(inout) :: a(m)
    real(p), intent(in)    :: b(n)
    integer, intent(in)    :: m, n
    integer                :: i, o

    if (m == n) then
        !$omp parallel do
        do i = 1, m
            a(i) = a(i) - b(i)
        end do
        !$omp end parallel do
        return
    end if

    o = m / n
    !$omp parallel do
    do i = 1, n
        a((i-1)*o+1:i*o) = a((i-1)*o+1:i*o) - b(i)
    end do
    !$omp end parallel do

end subroutine subtract_inplace


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine multiply_inplace(a, m, b, n)

    use module_tamasis,  only : p
    implicit none

    real(p), intent(inout) :: a(m)
    real(p), intent(in)    :: b(n)
    integer, intent(in)    :: m, n
    integer                :: i, o

    if (m == n) then
        !$omp parallel do
        do i = 1, m
            a(i) = a(i) * b(i)
        end do
        !$omp end parallel do
        return
    end if

    o = m / n
    !$omp parallel do
    do i = 1, n
        a((i-1)*o+1:i*o) = a((i-1)*o+1:i*o) * b(i)
    end do
    !$omp end parallel do

end subroutine multiply_inplace


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine divide_inplace(a, m, b, n)

    use module_tamasis,  only : p
    implicit none

    real(p), intent(inout) :: a(m)
    real(p), intent(in)    :: b(n)
    integer, intent(in)    :: m, n
    integer                :: i, o

    if (m == n) then
        !$omp parallel do
        do i = 1, m
            a(i) = a(i) / b(i)
        end do
        !$omp end parallel do
        return
    end if

    o = m / n
    !$omp parallel do
    do i = 1, n
        a((i-1)*o+1:i*o) = a((i-1)*o+1:i*o) / b(i)
    end do
    !$omp end parallel do

end subroutine divide_inplace


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine sum(array, n, output)

    use omp_lib,        only : omp_get_max_threads
    use module_math,    only : sum_kahan
    use module_tamasis, only : p
    implicit none

    real(p), intent(in)  :: array(n)
    integer, intent(in)  :: n
    real(p), intent(out) :: output

    integer :: nthreads, work, i

    if (n <= 1024) then
        output = sum_kahan(array)
        return
    end if

    output = 0

    nthreads = omp_get_max_threads()
    work = max(n / nthreads, 1)
    !$omp parallel do reduction(+:output)
    do i = 1, min(n, nthreads)
        output = output + sum_kahan(array(min((i-1)*work+1,n+1):min(i*work,n)))
    end do
    !$omp end parallel do

end subroutine sum


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine norm2(array, n, output)

    use omp_lib,        only : omp_get_max_threads
    use module_math,    only : norm => norm2
    use module_tamasis, only : p
    implicit none

    real(p), intent(in)  :: array(n)
    integer, intent(in)  :: n
    real(p), intent(out) :: output

    integer :: nthreads, work, i

    if (n <= 1024) then
        output = norm(array)
        return
    end if

    output = 0

    nthreads = omp_get_max_threads()
    work = max(n / nthreads, 1)
    !$omp parallel do reduction(+:output)
    do i = 1, min(n, nthreads)
        output = output + norm(array(min((i-1)*work+1,n+1):min(i*work,n)))
    end do
    !$omp end parallel do

end subroutine norm2


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine dot(array1, array2, n, output)

    use omp_lib,        only : omp_get_max_threads
    use module_math,    only : d => dot
    use module_tamasis, only : p
    implicit none

    real(p), intent(in)  :: array1(n)
    real(p), intent(in)  :: array2(n)
    integer, intent(in)  :: n
    real(p), intent(out) :: output

    integer :: nthreads, work, i, a, z

    if (n <= 1024) then
        output = d(array1, array2)
        return
    end if

    output = 0

    nthreads = omp_get_max_threads()
    work = max(n / nthreads, 1)
    !$omp parallel do reduction(+:output)
    do i = 1, min(n, nthreads)
        a = min((i-1)*work+1,n+1)
        z = min(i*work,n)
        output = output + d(array1(a:z), array2(a:z))
    end do
    !$omp end parallel do

end subroutine dot


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine diff(array, asize, dim, ashape, arank)

    use module_math,    only : diff_fast, diff_slow
    use module_tamasis, only : p
    implicit none

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


subroutine masking(input, m, mask, n, status)

    use iso_fortran_env, only : ERROR_UNIT
    use module_tamasis,  only : p
    implicit none

    real(p), intent(inout) :: input(m)
    logical*1, intent(in)  :: mask(n)
    integer, intent(in)    :: m, n
    integer, intent(out)   :: status

    integer :: i, o

    if (modulo(m, n) /= 0) then
        write (ERROR_UNIT,'(a,2(i0,a))') "The data array has a size incompatible with the mask ('", m, "' instead of '", n, "')."
        status = 1
        return
    end if

    status = 0

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

end subroutine masking


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


subroutine convolution_trexp_direct(data, nsamples, nslices, tau, nsamples_tot, ndetectors)

    use module_filtering, only : trexp_direct => convolution_trexp_direct
    use module_tamasis,   only : p
    implicit none

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

    ! before exchanging data, compute the local map bounds of the packed output
    allocate (a_packed(0:size-1), z_packed(0:size-1))
    a_packed(0) = 1
    z_packed(0) = 0
    source = 0
    do
        do imask = 1, ninputs
            if (imask > nmasks) exit
            if (mask(imask + source * ninputs)) cycle
            z_packed(source) = z_packed(source) + 1
        end do
        if (source == size - 1) exit
        source = source + 1
        a_packed(source) = z_packed(source-1) + 1
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
