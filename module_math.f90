module module_math

    use precision, only : p
    implicit none
    private

    public :: mean
    public :: stddev
    public :: moment
    public :: sum_kahan
    public :: test_real_eq

    interface sum_kahan
        module procedure sum_kahan_1d, sum_kahan_2d, sum_kahan_3d
    end interface sum_kahan


contains


    subroutine moment(input, mean, variance, skewness, kurtosis, stddev, meandev)
        real(kind=p), intent(in)            :: input(:)
        real(kind=p), intent(out), optional :: mean, variance, skewness
        real(kind=p), intent(out), optional :: kurtosis, stddev, meandev
        integer                   :: nsamples
        real(kind=p)              :: m, var, sdev, NaN, zero
        real(kind=p), allocatable :: residuals(:)

        !XXX NaN should come from intrinsic module
        zero = 0_p
        NaN = 0_p / zero

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


    recursive function sum_kahan_1d(input) result(sum)
        real(kind=p), intent(in) :: input(:)
        real(kind=p)             :: sum, c, t, y
        integer                  :: i

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


    recursive function sum_kahan_2d(input) result(sum)
        real(kind=p), intent(in) :: input(:,:)
        real(kind=p)             :: sum, c, t, y
        integer                  :: i

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


    recursive function sum_kahan_3d(input) result(sum)
        real(kind=p), intent(in) :: input(:,:,:)
        real(kind=p)             :: sum, c, t, y
        integer                  :: i

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


    elemental function test_real_eq(a, b, n)
        real(kind=p), intent(in) :: a, b
        integer, intent(in)      :: n
        logical                  :: test_real_eq
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
        test_real_eq = abs(a-b) <= epsilon * abs(a)

   end function test_real_eq

end module module_math
