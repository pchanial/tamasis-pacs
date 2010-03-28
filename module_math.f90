module module_math

    use precision, only : p, dp
    implicit none
    private

    public :: linspace
    public :: logspace
    public :: mad
    public :: mean
    public :: median
    public :: moment
    public :: nint_down
    public :: nint_up
    public :: stddev
    public :: sum_kahan
    public :: swap
    public :: test_real_eq
    public :: test_real_neq
    public :: NaN

    interface sum_kahan
        module procedure sum_kahan_1d, sum_kahan_2d, sum_kahan_3d
    end interface sum_kahan

    !XXX should use ieee_arithmetic instead when gfortran implements it
    real(kind=p), parameter ::                                                                       &
        NaN = transfer('1111111111111000000000000000000000000000000000000000000000000000'b, 0._dp),  &
        Infinity = transfer('0111111111110000000000000000000000000000000000000000000000000000'b,0._dp)


contains


    subroutine moment(input, mean, variance, skewness, kurtosis, stddev, meandev)
        real(kind=p), intent(in)            :: input(:)
        real(kind=p), intent(out), optional :: mean, variance, skewness
        real(kind=p), intent(out), optional :: kurtosis, stddev, meandev
        integer                   :: nsamples
        real(kind=p)              :: m, var, sdev
        real(kind=p), allocatable :: residuals(:)

        nsamples = size(input)

        ! compute the mean
        if (nsamples == 0) then
            m = NaN
        else
            m = sum_kahan(input) / nsamples
        end if
        if (present(mean)) mean = m

        if (.not. present(meandev)  .and. .not. present(stddev)   .and.&
            .not. present(variance) .and. .not. present(skewness) .and.&
            .not. present(kurtosis)) return

        allocate(residuals(nsamples))
        residuals = input - m

        ! compute mean deviation
        if (present(meandev)) meandev = sum_kahan(abs(residuals)) / nsamples

        ! compute variance
        if (nsamples <= 1) then
            var = NaN
        else
            var = sum_kahan(residuals**2) / (nsamples-1)
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
                skewness = sum_kahan(residuals**3) / (nsamples * sdev**3)
            end if
        end if

        ! compute kurtosis
        if (present(kurtosis)) then
            if (sdev == 0) then
                kurtosis = NaN
            else
                kurtosis = sum_kahan(residuals**4) / (nsamples * sdev**4)-3.0_p
            end if
        end if

        deallocate(residuals)

    end subroutine moment


    !---------------------------------------------------------------------------


    function mean(input)
        real(kind=p)             :: mean
        real(kind=p), intent(in) :: input(:)
        
        call moment(input, mean)

    end function mean


    !---------------------------------------------------------------------------


    function stddev(input)
        real(kind=p)             :: stddev
        real(kind=p), intent(in) :: input(:)
        
        call moment(input, stddev=stddev)

    end function stddev


    !---------------------------------------------------------------------------


    function sum_kahan_1d(input) result(sum)
        real(kind=p), intent(in) :: input(:)
        real(kind=p)             :: sum, c, t, y
        integer                  :: i

        if (size(input) == 0) then
            sum = NaN
            return
        end if

        sum = input(1)
        c = 0.0_p
        do i = 2, size(input)
            y = input(i) - c
            t = sum + y
            c = (t - sum) - y
            sum = t
        end do
    end function sum_kahan_1d


    !---------------------------------------------------------------------------


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


    !---------------------------------------------------------------------------


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


    !---------------------------------------------------------------------------


    function linspace(min, max, n)
        real(kind=p), allocatable :: linspace(:)
        real(kind=p), intent(in)  :: min, max
        integer, intent(in)       :: n
        integer                   :: i

        allocate(linspace(n))
        linspace = min + (max - min) / (n-1) * [(i, i=0,n-1)]

    end function linspace


    !---------------------------------------------------------------------------


    function logspace(min, max, n)
        real(kind=p), allocatable :: logspace(:)
        real(kind=p), intent(in)  :: min, max
        integer, intent(in)       :: n
        integer                   :: i

        allocate(logspace(n))
        logspace = exp(log(min)+(log(max)-log(min)) / (n-1) * [(i, i=0,n-1)])

    end function logspace


    !---------------------------------------------------------------------------


    elemental function nint_down(x)
        integer                  :: nint_down
        real(kind=p), intent(in) :: x
        
        nint_down = nint(x)
        if (x > 0 .and. abs(x-nint_down) == 0.5_p) then
            nint_down = nint_down - 1
        end if

    end function nint_down


    !---------------------------------------------------------------------------


    elemental function nint_up(x)
        integer                  :: nint_up
        real(kind=p), intent(in) :: x
        
        nint_up = nint(x)
        if (x < 0 .and. abs(x-nint_up) == 0.5_p) then
            nint_up = nint_up + 1
        end if

    end function nint_up


    !---------------------------------------------------------------------------


    elemental subroutine swap(a,b)
        real(kind=p), intent(inout) :: a, b
        real(kind=p)                :: tmp
        tmp = a
        a   = b
        b   = tmp
    end subroutine swap


    !---------------------------------------------------------------------------


    ! This Quickselect routine is based on the algorithm described in
    ! "Numerical recipes in C", Second Edition,
    ! Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
    ! input array may be reordered
    ! This code by Nicolas Devillard - 1998. Public domain.
    function median(arr) 
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

    end function median 


    !---------------------------------------------------------------------------


    ! returns the median absolute deviation
    function mad(x, m)
        real(kind=p)                        :: mad
        real(kind=p), intent(in)            :: x(:)
        real(kind=p), intent(out), optional :: m
        real(kind=p)                        :: x_(size(x)), med
 
        x_ = x
        med = median(x_)
        x_ = abs(x_ - med)
        mad = median(x_)

        if (present(m)) m = med

    end function mad


    !---------------------------------------------------------------------------


    elemental function test_real_eq(a, b, n)
        logical                  :: test_real_eq
        real(kind=p), intent(in) :: a, b
        integer, intent(in)      :: n
        real(kind=p)             :: epsilon

        ! check for NaN values
        if (a /= a) then
            test_real_eq = b /= b
            return
        end if
        if (b /= b) then
            test_real_eq = .false.
            return
        end if

        epsilon = 10_p**(-real(n, kind=p))
        test_real_eq = abs(a-b) <= epsilon * abs(max(a,b))

    end function test_real_eq


    !---------------------------------------------------------------------------


    elemental function test_real_neq(a, b, n)
        logical                  :: test_real_neq
        real(kind=p), intent(in) :: a, b
        integer, intent(in)      :: n
        
        test_real_neq = .not. test_real_eq(a, b, n)

    end function test_real_neq


end module module_math
