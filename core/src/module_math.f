module module_math

    use iso_fortran_env, only : ERROR_UNIT
    use module_tamasis,  only : p

    implicit none
    private

    public :: PI
    public :: DEG2RAD
    public :: RAD2DEG
    public :: mInf
    public :: pInf
    public :: NaN

    public :: distance
    public :: profile
    public :: linspace
    public :: logspace
    public :: mad
    public :: mean
    public :: mean_degrees
    public :: minmax_degrees
    public :: angle_lonlat
    public :: barycenter_lonlat
    public :: median
    public :: median_nocopy
    public :: moment
    public :: nint_down
    public :: nint_up
    public :: sigma_clipping
    public :: stddev
    public :: sum_kahan
    public :: swap
    public :: eq_real
    public :: neq_real
    public :: add
    public :: add_blas
    public :: diff_fast
    public :: diff_slow
    public :: diffT_fast
    public :: diffT_slow
    public :: diffTdiff_fast
    public :: diffTdiff_slow
    public :: shift_fast
    public :: shift_slow

    interface distance
        module procedure distance_1d, distance_2d, distance_3d
    end interface

    interface profile
        module procedure profile_axisymmetric_2d
    end interface

    interface sum_kahan
        module procedure sum_kahan_1d, sum_kahan_2d, sum_kahan_3d
    end interface

    interface median
        module procedure median_nomask, median_mask
    end interface

    real(p), parameter :: PI = 4._p * atan(1._p)
    real(p), parameter :: DEG2RAD = PI / 180._p
    real(p), parameter :: RAD2DEG = 180._p / PI
    !XXX should use ieee_arithmetic instead when gfortran implements it. So far, gfortran doesn't allow NaN, mInf, pInf conversion
    ! between different real kinds...
#if PRECISION_REAL == 4
    real(p), parameter ::                                                                                                          &
        NaN  = transfer('11111111110000000000000000000000'b, 0._p),                                                                &
        mInf = transfer('11111111100000000000000000000000'b, 0._p),                                                                &
        pInf = transfer('01111111100000000000000000000000'b, 0._p)
#elif PRECISION_REAL == 8
    real(p), parameter ::                                                                                                          &
        NaN  = transfer('1111111111111000000000000000000000000000000000000000000000000000'b, 0._p),                                &
        mInf = transfer('1111111111110000000000000000000000000000000000000000000000000000'b, 0._p),                                &
        pInf = transfer('0111111111110000000000000000000000000000000000000000000000000000'b, 0._p)
#elif PRECISION_REAL == 16
    real(16), parameter :: NaN  = 'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF'z
    real(16), parameter :: mInf = 'FFFF0000000000000000000000000000'z
    real(16), parameter :: pInf = '7FFF0000000000000000000000000000'z
#endif


contains


    subroutine moment(input, mean, variance, skewness, kurtosis, stddev, meandev, mask)

        real(p), intent(in)            :: input(:)
        real(p), intent(out), optional :: mean, variance, skewness
        real(p), intent(out), optional :: kurtosis, stddev, meandev
        logical, intent(in), optional  :: mask(:)

        integer              :: nsamples
        real(p)              :: m, var, sdev
        real(p), allocatable :: residuals(:)

        nsamples = size(input)

        if (present(mask)) then
            if (size(mask) /= nsamples) then
                write (ERROR_UNIT,'(a,2(i0,a))'), 'Error: the mask has an incompatible dimension (', size(mask), ' instead of ',   &
                     size(input), ').'
                stop 1
            end if
            nsamples = nsamples - count(mask)
        end if

        ! compute the mean
        if (nsamples == 0) then
            m = NaN
        else
            m = sum_kahan(input, mask) / nsamples
        end if
        if (present(mean)) mean = m

        if (.not. present(meandev)  .and. .not. present(stddev)   .and.                                                            &
            .not. present(variance) .and. .not. present(skewness) .and.                                                            &
            .not. present(kurtosis)) return

        allocate(residuals(size(input)))
        residuals = input - m

        ! compute mean deviation
        if (present(meandev)) meandev = sum_kahan(abs(residuals), mask) / nsamples

        ! compute variance
        if (nsamples <= 1) then
            var = NaN
        else
            var = sum_kahan(residuals**2, mask) / (nsamples-1)
        end if
        if (present(variance)) variance = var

        ! compute standard deviation
        sdev = sqrt(var)
        if (present(stddev)) stddev = sdev

        ! compute skewness
        if (present(skewness)) then
            if (sdev == 0) then
                skewness = NaN
            else
                skewness = sum_kahan(residuals**3, mask) / (nsamples * sdev**3)
            end if
        end if

        ! compute kurtosis
        if (present(kurtosis)) then
            if (sdev == 0) then
                kurtosis = NaN
            else
                kurtosis = sum_kahan(residuals**4, mask) / (nsamples * sdev**4) - 3
            end if
        end if

        deallocate(residuals)

    end subroutine moment


    !-------------------------------------------------------------------------------------------------------------------------------


    function mean(input, mask)

        real(p)                      :: mean
        real(p), intent(in)          :: input(:)
        logical, intent(in),optional :: mask(:)

        call moment(input, mean, mask=mask)

    end function mean


    !-------------------------------------------------------------------------------------------------------------------------------


    function stddev(input, mask)

        real(p)                       :: stddev
        real(p), intent(in)           :: input(:)
        logical, intent(in), optional :: mask(:)

        call moment(input, stddev=stddev, mask=mask)

    end function stddev


    !-------------------------------------------------------------------------------------------------------------------------------


    function sum_kahan_1d(input, mask) result (sum)

        real(p), intent(in)           :: input(:)
        logical, intent(in), optional :: mask(:)
        real(p)                       :: sum, c, t, y

        integer                       :: i

        if (size(input) == 0) then
            sum = NaN
            return
        end if

        if (present(mask)) then
            if (size(mask) /= size(input)) then
                write (ERROR_UNIT,'(a,2(i0,a))'), 'Error: the mask has an incompatible dimension (', size(mask), ' instead of ',   &
                     size(input), ').'
                stop 1
            end if
            do i = 1, size(input)
                if (.not. mask(i)) exit
            end do
            if (i > size(input)) then
                sum = NaN
                return
            end if
        else
            i = 1
        end if

        sum = input(i)
        c = 0
        do i = i+1, size(input)
            if (present(mask)) then
                if (mask(i)) cycle
            end if
            y = input(i) - c
            t = sum + y
            c = (t - sum) - y
            sum = t
        end do

    end function sum_kahan_1d


    !-------------------------------------------------------------------------------------------------------------------------------


    function sum_kahan_2d(input) result (sum)

        real(p), intent(in) :: input(:,:)
        real(p)             :: sum, c, t, y

        integer             :: i

        if (size(input) == 0) then
            sum = NaN
            return
        end if

        sum = sum_kahan_1d(input(:,1))
        c = 0
        do i = 2, size(input,2)
            y = sum_kahan_1d(input(:,i)) - c
            t = sum + y
            c = (t - sum) - y
            sum = t
        end do
    end function sum_kahan_2d


    !-------------------------------------------------------------------------------------------------------------------------------


    function sum_kahan_3d(input) result (sum)

        real(p), intent(in) :: input(:,:,:)
        real(p)             :: sum, c, t, y

        integer             :: i

        if (size(input) == 0) then
            sum = NaN
            return
        end if

        sum = sum_kahan_2d(input(:,:,1))
        c = 0
        do i = 2, size(input,3)
            y = sum_kahan_2d(input(:,:,i)) - c
            t = sum + y
            c = (t - sum) - y
            sum = t
        end do
    end function sum_kahan_3d


    !-------------------------------------------------------------------------------------------------------------------------------


    function mean_degrees(array, mask)

        real(p), intent(in)           :: array(:)
        logical, optional, intent(in) :: mask(:)
        real(p)                       :: mean_degrees

        real(p) :: value
        integer :: isample, n180, nvalids
        logical :: zero_minus, zero_plus

        mean_degrees = 0
        zero_minus = .false.
        zero_plus  = .false.
        n180 = 0
        nvalids = 0

        !$omp parallel do default(shared) reduction(+:n180, nvalids,mean_degrees)           &
        !$omp reduction(.or.:zero_minus,zero_plus) private(isample, value)
        do isample=1, size(array)
            if (present(mask)) then
                if (mask(isample)) cycle
            end if
            value = modulo(array(isample), 360._p)
            if (value /= value) cycle
            zero_minus = zero_minus .or. value > 270._p
            zero_plus  = zero_plus  .or. value <= 90._p
            if (value >= 180._p) n180 = n180 + 1
            mean_degrees = mean_degrees + value
            nvalids = nvalids + 1
        end do
        !$omp end parallel do

        if (zero_minus .and. zero_plus) mean_degrees = mean_degrees - 360._p * n180

        if (nvalids > 0) then
            mean_degrees = modulo(mean_degrees / nvalids, 360._p)
        else
            mean_degrees = NaN
        end if

    end function mean_degrees


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine minmax_degrees(array, minv, maxv)

        real(p), intent(in)  :: array(:)
        real(p), intent(out) :: minv, maxv

        real(p) :: value, meanv
        integer :: isample

        meanv = mean_degrees(array)
        if (meanv /= meanv) then
            minv = NaN
            maxv = NaN
            return
        end if

        minv = pInf
        maxv = mInf
        !$omp parallel do default(shared) reduction(min:minv) reduction(max:maxv) private(value)
        do isample=1, size(array)
            value = array(isample)
            if (value /= value) cycle
            value = modulo(value, 360._p)
            if (value > meanv) then
                if (abs(value-360._p-meanv) < abs(value-meanv)) value = value - 360._p
            else
                if (abs(value+360._p-meanv) < abs(value-meanv)) value = value + 360._p
            end if
            minv = min(minv, value)
            maxv = max(maxv, value)
        end do
        !$omp end parallel do

        minv = modulo(minv, 360._p)
        maxv = modulo(maxv, 360._p)

    end subroutine minmax_degrees


    !-------------------------------------------------------------------------------------------------------------------------------


    pure elemental subroutine angle_lonlat(lon1, lat1, lon2, lat2, angle)

        real(p), intent(in)  :: lon1, lat1, lon2, lat2
        real(p), intent(out) :: angle
        
        angle = acos(cos(lat1*DEG2RAD) * cos(lat2*DEG2RAD) * cos((lon2 - lon1)*DEG2RAD) + sin(lat1*DEG2RAD) * sin(lat2*DEG2RAD)) * &
                RAD2DEG

    end subroutine angle_lonlat


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine barycenter_lonlat(lon, lat, lon0, lat0)

        real(p), intent(in)  :: lon(:), lat(:)
        real(p), intent(out) :: lon0, lat0

        real(p) :: x, y, z, phi, cotheta, cocotheta
        integer :: i

        x = 0
        y = 0
        z = 0
        !$omp parallel do reduction(+:x,y,z) private(phi,cotheta,cocotheta)
        do i = 1, size(lon)
            phi = lon(i) * DEG2RAD
            cotheta = lat(i) * DEG2RAD
            if (cotheta /= cotheta .or. phi /= phi) cycle
            cocotheta = cos(cotheta)
            x = x + cocotheta * cos(phi)
            y = y + cocotheta * sin(phi)
            z = z + sin(cotheta)
        end do
        !$omp end parallel do
        
        lon0 = modulo(atan2(y,x) * RAD2DEG + 360._p, 360._p)
        lat0 = acos(sqrt((x**2 + y**2)/(x**2 + y**2 + z**2))) * RAD2DEG
        lat0 = sign(lat0, z)
        
    end subroutine barycenter_lonlat


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine distance_1d(distance, origin, resolution)

        real(p), intent(out) :: distance(:)
        real(p), intent(in)  :: origin
        real(p), intent(in)  :: resolution

        integer :: i
        real(p) :: resolution_

        resolution_ = abs(resolution)

        !$omp parallel do
        do i = 1, size(distance,1)
            distance(i) = resolution_ * abs(i - origin)
        end do
        !$omp end parallel do

    end subroutine distance_1d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine distance_2d(distance, origin, resolution)

        real(p), intent(out) :: distance(:,:)
        real(p), intent(in)  :: origin(2)
        real(p), intent(in)  :: resolution(2)

        integer :: i, j
        real(p) :: tmpj

        !$omp parallel do private(tmpj)
        do j = 1, size(distance,2)
            tmpj = (resolution(2) * (j - origin(2)))**2
            do i = 1, size(distance,1)
                distance(i,j) = sqrt((resolution(1) * (i - origin(1)))**2 + tmpj)
            end do
        end do
        !$omp end parallel do

    end subroutine distance_2d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine distance_3d(distance, origin, resolution)

        real(p), intent(out) :: distance(:,:,:)
        real(p), intent(in)  :: origin(3)
        real(p), intent(in)  :: resolution(3)

        integer :: i, j, k
        real(p) :: tmpj, tmpk

        !$omp parallel do private(tmpj, tmpk)
        do k = 1, size(distance,3)
            tmpk = (resolution(3) * (k - origin(3)))**2
            do j = 1, size(distance,2)
                tmpj = (resolution(2) * (j - origin(2)))**2
                do i = 1, size(distance,1)
                    distance(i,j,k) = sqrt((resolution(1) * (i - origin(1)))**2 + tmpj + tmpk)
                end do
            end do
        end do
        !$omp end parallel do

    end subroutine distance_3d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine profile_axisymmetric_2d(array, origin, bin, nbins, x, y, n)

        real(p), intent(in)  :: array(:,:)
        real(p), intent(in)  :: origin(2)
        real(p), intent(in)  :: bin
        integer, intent(in)  :: nbins
        real(p), intent(out) :: x(nbins)
        real(p), intent(out) :: y(nbins)
        integer, intent(out) :: n(nbins)

        integer :: i, j, ibin
        real(p) :: distance, value, tmpj

        x = 0
        y = 0
        n = 0
        !$omp parallel do private(distance, ibin, tmpj, value) reduction(+:x,y,n)
        do j = 1, size(array,2)
            tmpj = (j - origin(2))**2
            do i = 1, size(array,1)
                value = array(i,j)
                if (value /= value) cycle
                distance = sqrt((i - origin(1))**2 + tmpj)
                ibin = int(distance / bin) + 1
                if (ibin > nbins) cycle
                x(ibin) = x(ibin) + distance
                y(ibin) = y(ibin) + value
                n(ibin) = n(ibin) + 1
            end do
        end do
        !$omp end parallel do

        do i = 1, nbins
            if (n(i) /= 0) then
                x(i) = x(i) / n(i)
                y(i) = y(i) / n(i)
            else
                x(i) = bin * (i - 0.5_p)
                y(i) = NaN
            end if
        end do

    end subroutine profile_axisymmetric_2d


    !-------------------------------------------------------------------------------------------------------------------------------


    function linspace(min, max, n)

        real(p)             :: linspace(n)
        real(p), intent(in) :: min, max
        integer, intent(in) :: n

        integer             :: i

        linspace = min + (max - min) / (n-1) * [(i, i=0, n-1)]

    end function linspace


    !-------------------------------------------------------------------------------------------------------------------------------


    function logspace(min, max, n)
        real(p)             :: logspace(n)
        real(p), intent(in) :: min, max
        integer, intent(in) :: n
        integer             :: i

        logspace = exp(log(min)+(log(max)-log(min)) / (n-1) * [(i, i=0,n-1)])

    end function logspace


    !-------------------------------------------------------------------------------------------------------------------------------


    elemental function nint_down(x)

        integer             :: nint_down
        real(p), intent(in) :: x

        nint_down = nint(x)
        if (x > 0 .and. abs(x-nint_down) == 0.5_p) then
            nint_down = nint_down - 1
        end if

    end function nint_down


    !-------------------------------------------------------------------------------------------------------------------------------


    elemental function nint_up(x)

        integer             :: nint_up
        real(p), intent(in) :: x

        nint_up = nint(x)
        if (x < 0 .and. abs(x-nint_up) == 0.5_p) then
            nint_up = nint_up + 1
        end if

    end function nint_up


    !-------------------------------------------------------------------------------------------------------------------------------


    elemental subroutine swap(a,b)

        real(p), intent(inout) :: a, b

        real(p)                :: tmp

        tmp = a
        a   = b
        b   = tmp

    end subroutine swap


    !-------------------------------------------------------------------------------------------------------------------------------


    function median_nomask(array) result (median) 

        real(p)             :: median
        real(p), intent(in) :: array(:)

        real(p), allocatable :: array_copy(:)

#ifdef IFORT
        allocate (array_copy(count(array == array)))
#endif
        array_copy = pack(array, array == array)
        median = median_nocopy(array_copy)

    end function median_nomask


    !-------------------------------------------------------------------------------------------------------------------------------


    function median_mask(array, mask) result (median) 

        real(p)               :: median
        real(p), intent(in)   :: array(:)
        logical*1, intent(in) :: mask(size(array))

        real(p), allocatable  :: array_copy(:)

#ifdef IFORT
        allocate (array_copy(count(array == array .and. .not. mask)))
#endif
        array_copy = pack(array, array == array .and. .not. mask)
        median = median_nocopy(array_copy)

    end function median_mask


    !-------------------------------------------------------------------------------------------------------------------------------


    ! This Quickselect routine is based on the algorithm described in
    ! "Numerical recipes in C", Second Edition,
    ! Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
    ! input array may be reordered
    ! This code by Nicolas Devillard - 1998. Public domain.
    function median_nocopy(arr) result (median)

        real(p)                :: median
        real(p), intent(inout) :: arr(0:)

        integer                :: low, high, imedian, middle, ll, hh

        if (size(arr) == 0) then
            median = NaN
            return
        end if

        low = 0
        high = size(arr)-1
        imedian = (low + high) / 2
        do
            if (high <= low) then
                median = arr(imedian)
                return
            end if

            if (high == low + 1) then  ! Two elements only
                if (arr(low) > arr(high)) call swap(arr(low), arr(high))
                median = arr(imedian)
                return
            end if

            ! Find imedian of low, middle and high items swap into position low
            middle = (low + high) / 2
            if (arr(middle) > arr(high)) call swap(arr(middle), arr(high))
            if (arr(low)    > arr(high)) call swap(arr(low),    arr(high))
            if (arr(middle) > arr(low))  call swap(arr(middle), arr(low))

            ! Swap low item (now in position middle) into position (low+1)
            call swap(arr(middle), arr(low+1)) 

            ! Nibble from each end towards middle, swapping items when stuck
            ll = low + 1
            hh = high
            do
                do 
                    ll = ll + 1
                    if (arr(low) <= arr(ll)) exit
                end do

                do 
                    hh = hh - 1
                    if (arr(hh)  <= arr(low)) exit
                end do

                if (hh < ll) exit

                call swap(arr(ll), arr(hh)) 

            end do

            ! Swap middle item (in position low) back into correct position
            call swap(arr(low), arr(hh)) 

            ! Re-set active partition
            if (hh <= imedian) low = ll
            if (hh >= imedian) high = hh - 1

        end do

    end function median_nocopy


    !-------------------------------------------------------------------------------------------------------------------------------


    ! returns the median absolute deviation
    function mad(x, m)

        real(p)                        :: mad
        real(p), intent(in)            :: x(:)
        real(p), intent(out), optional :: m

        real(p)                        :: x_(size(x)), med

        x_ = x
        med = median_nocopy(x_)
        x_ = abs(x_ - med)
        mad = median_nocopy(x_)

        if (present(m)) m = med

    end function mad


    !-------------------------------------------------------------------------------------------------------------------------------


    ! sigma clip an input vector. In the output mask, .true. means rejected
    subroutine sigma_clipping(input, mask, nsigma, nitermax)

        real(p), intent(in)           :: input(:)
        logical, intent(out)          :: mask(size(input))
        real(p), intent(in)           :: nsigma
        integer, intent(in), optional :: nitermax

        integer :: iter, nitermax_
        real(p) :: mean, stddev
        logical :: newmask(size(input))

        if (present(nitermax)) then
            nitermax_ = nitermax
        else
            nitermax_ = 0
        end if

        mask = .false.
        iter = 1
        do

            call moment(input, mean, stddev=stddev, mask=mask)
            newmask = mask .or. abs(input - mean) > nsigma * stddev
            if (count(newmask) == count(mask)) exit
            mask = newmask
            if (iter == nitermax_) exit

        end do

    end subroutine sigma_clipping


    !-------------------------------------------------------------------------------------------------------------------------------


    elemental function eq_real(a, b, rtol)

        logical                       :: eq_real
        real(p), intent(in)           :: a, b
        real(p), intent(in), optional :: rtol

        real(p) :: rtol_

        ! check for NaN values
        if (a /= a) then
            eq_real = b /= b
            return
        end if
        if (b /= b) then
            eq_real = .false.
            return
        end if

        if (present(rtol)) then
            rtol_ = rtol
        else
            rtol_ = 10._p * epsilon(1._p)
        end if

        eq_real = abs(a-b) <= rtol_ * max(abs(a),abs(b))

    end function eq_real


    !-------------------------------------------------------------------------------------------------------------------------------


    elemental function neq_real(a, b, rtol)

        logical                       :: neq_real
        real(p), intent(in)           :: a, b
        real(p), intent(in), optional :: rtol

        neq_real = .not. eq_real(a, b, rtol)

    end function neq_real


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine add(a, b, n)
        real(p), intent(inout) :: a(n)
        real(p), intent(in)    :: b(n)
        integer, intent(in)    :: n
        integer                :: i

        !$omp parallel do
        do i = 1, n
            a(i) = a(i) + b(i)
        end do
        !$omp end parallel do

    end subroutine add


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine add_blas(a, b)

        real(p), intent(inout) :: a(:)
        real(p), intent(in)    :: b(:)

        call daxpy(size(a), 1._p, b, 1, a, 1)

    end subroutine add_blas


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine diff_fast(array, m, n)

        real(p), intent(inout) :: array(m,n)
        integer, intent(in)    :: m, n

        integer :: i, j

        !$omp parallel do private(i,j)
        do j = 1, n
            do i = 1, m-1
                array(i,j) = array(i,j) - array(i+1,j)
            end do
            array(m,j) = 0
        end do
        !$omp end parallel do

    end subroutine diff_fast


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine diff_slow(array, m, n, o)

        integer, intent(in)    :: m, n, o
        real(p), intent(inout) :: array(m,n,o)

        integer, parameter     :: block = 4096
        integer                :: i, j, k, a, z

        !$omp parallel do private(i,j,k,a,z)
        do k = 1, o
            ! block computation to avoid cache exhaustion
            do i = 1, (m-1) / block + 1
                a = (i-1) * block + 1
                z = min(m, i * block)
                do j = 1, n - 1
                    array(a:z,j,k) = array(a:z,j,k) - array(a:z,j+1,k)
                end do
            end do
            array(:,n,k) = 0
        end do
        !$omp end parallel do

    end subroutine diff_slow


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine diffT_fast(array, m, n)

        real(p), intent(inout) :: array(m,n)
        integer, intent(in)    :: m, n

        integer :: i, j

        !$omp parallel do private(i,j)
        do j = 1, n
            array(m,j) = -array(m-1,j)
            do i = m-1, 2, -1
                array(i,j) = array(i,j) - array(i-1,j)
            end do
        end do
        !$omp end parallel do

    end subroutine diffT_fast


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine diffT_slow(array, m, n, o)

        integer, intent(in)    :: m, n, o
        real(p), intent(inout) :: array(m,n,o)

        integer, parameter     :: block = 4096
        integer                :: i, j, k, a, z

        !$omp parallel do private(i,j,k,a,z)
        do k = 1, o
            do i = 1, (m-1) / block + 1
                a = (i-1) * block + 1
                z = min(m, i * block)
                array(a:z,n,k) = -array(a:z,n-1,k)
                do j = n-1, 2, -1
                    array(a:z,j,k) = array(a:z,j,k) - array(a:z,j-1,k)
                end do
            end do
        end do
        !$omp end parallel do

    end subroutine diffT_slow


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine diffTdiff_fast(array, m, n, scalar)

        use module_tamasis,  only : p
        implicit none

        real(p), intent(inout)        :: array(m,n)
        integer, intent(in)           :: m, n
        real(p), intent(in), optional :: scalar

        integer :: i, j
        real(p) :: v, w, s

        if (m == 1) then
            array = 0
            return
        end if

        if (present(scalar)) then
            s = scalar
        else
            s = 1._p
        end if

        !$omp parallel do private(i,j,v,w)
        do j = 1, n
            v = array(1,j)
            array(1,j) = s * (array(1,j) - array(2,j))
            do i = 2, m-1
                w = array(i,j)
                array(i,j) = s * (2 * w - v - array(i+1,j))
                v = w
            end do
            array(m,j) = s * (array(m,j) - v)
        end do
        !$omp end parallel do

    end subroutine diffTdiff_fast


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine diffTdiff_slow(array, m, n, o, scalar)

        use module_tamasis,  only : p
        implicit none

        integer, intent(in)           :: m, n, o
        real(p), intent(inout)        :: array(m,n,o)
        real(p), intent(in), optional :: scalar

        integer, parameter :: block = 4096
        integer :: h, i, j, k, a, z
        real(p) :: w, v(block), s
        
        if (n == 1) then
            array = 0
            return
        end if

        if (present(scalar)) then
            s = scalar
        else
            s = 1._p
        end if

        !$omp parallel do private(i,j,k,a,z,v,w)
        do k = 1, o
            do i = 1, (m-1) / block + 1
                a = (i-1) * block + 1
                z = min(i * block, m)
                v(1:z-a+1) = array(a:z,1,k)
                array(a:z,1,k) = s * (v(1:z-a+1) - array(a:z,2,k))
                do j = 2, n-1
                    do h = a, z
                        w = array(h,j,k)
                        array(h,j,k) = s * (2 * w - v(h-a+1) - array(h,j+1,k))
                        v(h-a+1) = w
                    end  do
                end do
                array(a:z,n,k) = s * (array(a:z,n,k) - v(1:z-a+1))
            end do
        end do
        !$omp end parallel do

    end subroutine diffTdiff_slow


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine shift_fast(array, m, n, offset)

        real(p), intent(inout) :: array(m,n)
        integer, intent(in)    :: m, n
        integer, intent(in)    :: offset(:)

        integer :: i, j, d

        !$omp parallel do private(d)
        do j = 1, n
            d = offset((j-1) / (n / size(offset)) + 1)
            if (d == 0) cycle
            if (d > 0) then
                do i = m, d+1, -1
                    array(i,j) = array(i-d,j)
                end do
                array(1:min(d,m),j) = 0
            else
                do i = 1, m+d
                    array(i,j) = array(i-d,j)
                end do
                array(max(m+d+1,1):m,j) = 0             
            end if
        end do
        !$omp end parallel do

    end subroutine shift_fast


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine shift_slow(array, m, n, o, offset)

        real(p), intent(inout) :: array(m,n,o)
        integer, intent(in)    :: m, n, o
        integer, intent(in)    :: offset(:)

        integer :: j, k, d

        !$omp parallel do private(d)
        do k = 1, o
            d = offset((k-1) / (o / size(offset)) + 1)
            if (d == 0) cycle
            if (d > 0) then
                do j = n, d+1, -1
                    array(:,j,k) = array(:,j-d,k)
                end do
                array(:,1:min(d,n),k) = 0
            else
                do j = 1, n+d
                    array(:,j,k) = array(:,j-d,k)
                end do
                array(:,max(n+d+1,1):n,k) = 0             
            end if
        end do
        !$omp end parallel do

    end subroutine shift_slow


end module module_math
