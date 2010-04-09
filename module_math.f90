module module_math

    use module_precision, only : p, dp
    implicit none
    private

    public :: PI
    public :: DEG2RAD
    public :: RAD2DEG
    public :: mInf
    public :: pInf
    public :: NaN

    public :: linspace
    public :: logspace
    public :: mad
    public :: mean
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

    interface sum_kahan
        module procedure sum_kahan_1d, sum_kahan_2d, sum_kahan_3d
    end interface sum_kahan

    real(p), parameter :: PI = 4.0_p * atan(1.0_p)
    real(p), parameter :: DEG2RAD = PI / 180._p
    real(p), parameter :: RAD2DEG = 180._p / PI
    !XXX should use ieee_arithmetic instead when gfortran implements it
    real(p), parameter ::                                                                       &
        NaN  = transfer('1111111111111000000000000000000000000000000000000000000000000000'b, 0._dp), &
        mInf = transfer('1111111111110000000000000000000000000000000000000000000000000000'b,0._dp),  &
        pInf = transfer('0111111111110000000000000000000000000000000000000000000000000000'b,0._dp)


contains


    subroutine moment(input, mean, variance, skewness, kurtosis, stddev, meandev, mask)

        real(kind=p), intent(in)            :: input(:)
        real(kind=p), intent(out), optional :: mean, variance, skewness
        real(kind=p), intent(out), optional :: kurtosis, stddev, meandev
        logical, intent(in), optional       :: mask(:)

        integer                   :: nsamples
        real(kind=p)              :: m, var, sdev
        real(kind=p), allocatable :: residuals(:)

        nsamples = size(input)
        if (present(mask)) then
            nsamples = nsamples - count(mask)
        end if

        ! compute the mean
        if (nsamples == 0) then
            m = NaN
        else
            m = sum_kahan(input, mask) / nsamples
        end if
        if (present(mean)) mean = m

        if (.not. present(meandev)  .and. .not. present(stddev)   .and.&
            .not. present(variance) .and. .not. present(skewness) .and.&
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
                kurtosis = sum_kahan(residuals**4, mask) / (nsamples * sdev**4)-3.0_p
            end if
        end if

        deallocate(residuals)

    end subroutine moment


    !-------------------------------------------------------------------------------------------------------------------------------


    function mean(input, mask)

        real(kind=p)                 :: mean
        real(kind=p), intent(in)     :: input(:)
        logical, intent(in),optional :: mask(:)
        
        call moment(input, mean, mask=mask)

    end function mean


    !-------------------------------------------------------------------------------------------------------------------------------


    function stddev(input, mask)

        real(kind=p)                  :: stddev
        real(kind=p), intent(in)      :: input(:)
        logical, intent(in), optional :: mask(:)
        
        call moment(input, stddev=stddev, mask=mask)

    end function stddev


    !-------------------------------------------------------------------------------------------------------------------------------


    function sum_kahan_1d(input, mask) result(sum)

        real(kind=p), intent(in)      :: input(:)
        logical, intent(in), optional :: mask(:)
        real(kind=p)                  :: sum, c, t, y

        integer                       :: i

        if (size(input) == 0) then
            sum = NaN
            return
        end if

        if (present(mask)) then
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
        c = 0.0_p
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


    function sum_kahan_2d(input) result(sum)

        real(kind=p), intent(in) :: input(:,:)
        real(kind=p)             :: sum, c, t, y

        integer                  :: i

        if (size(input) == 0) then
            sum = NaN
            return
        end if

        sum = sum_kahan_1d(input(:,1))
        c = 0.0_p
        do i = 2, size(input,2)
            y = sum_kahan_1d(input(:,i)) - c
            t = sum + y
            c = (t - sum) - y
            sum = t
        end do
    end function sum_kahan_2d


    !-------------------------------------------------------------------------------------------------------------------------------


    function sum_kahan_3d(input) result(sum)

        real(kind=p), intent(in) :: input(:,:,:)
        real(kind=p)             :: sum, c, t, y

        integer                  :: i

        if (size(input) == 0) then
            sum = NaN
            return
        end if

        sum = sum_kahan_2d(input(:,:,1))
        c = 0.0_p
        do i = 2, size(input,3)
            y = sum_kahan_2d(input(:,:,i)) - c
            t = sum + y
            c = (t - sum) - y
            sum = t
        end do
    end function sum_kahan_3d


    !-------------------------------------------------------------------------------------------------------------------------------


    function linspace(min, max, n)

        real(kind=p)              :: linspace(n)
        real(kind=p), intent(in)  :: min, max
        integer, intent(in)       :: n

        integer                   :: i

        linspace = min + (max - min) / (n-1) * [(i, i=0, n-1)]

    end function linspace


    !-------------------------------------------------------------------------------------------------------------------------------


    function logspace(min, max, n)
        real(kind=p)              :: logspace(n)
        real(kind=p), intent(in)  :: min, max
        integer, intent(in)       :: n
        integer                   :: i

        logspace = exp(log(min)+(log(max)-log(min)) / (n-1) * [(i, i=0,n-1)])

    end function logspace


    !-------------------------------------------------------------------------------------------------------------------------------


    elemental function nint_down(x)

        integer                  :: nint_down
        real(kind=p), intent(in) :: x
        
        nint_down = nint(x)
        if (x > 0 .and. abs(x-nint_down) == 0.5_p) then
            nint_down = nint_down - 1
        end if

    end function nint_down


    !-------------------------------------------------------------------------------------------------------------------------------


    elemental function nint_up(x)

        integer                  :: nint_up
        real(kind=p), intent(in) :: x
        
        nint_up = nint(x)
        if (x < 0 .and. abs(x-nint_up) == 0.5_p) then
            nint_up = nint_up + 1
        end if

    end function nint_up


    !-------------------------------------------------------------------------------------------------------------------------------


    elemental subroutine swap(a,b)
        real(kind=p), intent(inout) :: a, b
        real(kind=p)                :: tmp
        tmp = a
        a   = b
        b   = tmp
    end subroutine swap


    !-------------------------------------------------------------------------------------------------------------------------------


    function median(arr) 

        real(p)             :: median
        real(p), intent(in) :: arr(:)

        real(p)             :: arr_copy(size(arr))

        arr_copy = arr
        median = median_nocopy(arr_copy)

    end function median


    !-------------------------------------------------------------------------------------------------------------------------------


    ! This Quickselect routine is based on the algorithm described in
    ! "Numerical recipes in C", Second Edition,
    ! Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
    ! input array may be reordered
    ! This code by Nicolas Devillard - 1998. Public domain.
    function median_nocopy(arr) result(median)

        real(kind=p)                :: median
        real(kind=p), intent(inout) :: arr(0:)

        integer                     :: low, high, imedian, middle, ll, hh

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

        real(kind=p)                        :: mad
        real(kind=p), intent(in)            :: x(:)
        real(kind=p), intent(out), optional :: m

        real(kind=p)                        :: x_(size(x)), med
 
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


    elemental function eq_real(a, b, n)

        logical                  :: eq_real
        real(kind=p), intent(in) :: a, b
        integer, intent(in)      :: n

        real(kind=p)             :: epsilon

        ! check for NaN values
        if (a /= a) then
            eq_real = b /= b
            return
        end if
        if (b /= b) then
            eq_real = .false.
            return
        end if

        epsilon = 10_p**(-real(n, kind=p))
        eq_real = abs(a-b) <= epsilon * abs(max(a,b))

    end function eq_real


    !-------------------------------------------------------------------------------------------------------------------------------


    elemental function neq_real(a, b, n)

        logical                  :: neq_real
        real(kind=p), intent(in) :: a, b

        integer, intent(in)      :: n
        
        neq_real = .not. eq_real(a, b, n)

    end function neq_real


end module module_math
