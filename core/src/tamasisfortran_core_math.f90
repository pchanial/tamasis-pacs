module math

    use module_tamasis, only : p
    use module_math,    only : median_nocopy
    implicit none

contains

    subroutine median(input, output, n)

        integer*8, intent(in)  :: n
        real(p), intent(in)    :: input(n)
        real(p), intent(inout) :: output

        real(p) :: input_(n)

        input_ = input
        output = median_nocopy(input_, .true.)

    end subroutine median

    subroutine median_mask(input, mask, output, n)

        integer*8, intent(in)  :: n
        real(p), intent(in)    :: input(n)
        logical*1, intent(in)  :: mask(n)
        real(p), intent(inout) :: output

        real(p)   :: input_(n)
        logical*1 :: mask_(n)

        input_ = input
        mask_ = mask
        output = median_nocopy(input_, .true., mask_)

    end subroutine median_mask

    subroutine median_axis(input, output, m, n, o)

        integer*8, intent(in)  :: m, n, o
        real(p), intent(in)    :: input(m, n, o)
        real(p), intent(inout) :: output(m, o)

        real(p)   :: input_(n)
        integer*8 :: i, j

#ifdef IFORT
        integer*8 :: k
        !$omp parallel do private(i, j, input_)
        do k = 1, size(input, 1) * size(input, 3)
            i = mod(k - 1, size(input, 1)) + 1
            j = (k - 1) / size(input, 1) + 1
#else
        !$omp parallel do collapse(2) private(input_)
        do j = 1, size(input, 3)
            do i = 1, size(input, 1)
#endif
                input_ = input(i,:,j)
                output(i,j) = median_nocopy(input_, .true.)
            end do
#ifndef IFORT
        end do
#endif
        !$omp end parallel do

    end subroutine median_axis

    subroutine median_mask_axis(input, mask, output, m, n, o)

        integer*8, intent(in)  :: m, n, o
        real(p), intent(in)    :: input(m, n, o)
        logical*1, intent(in)  :: mask (m, n, o)
        real(p), intent(inout) :: output(m, o)

        real(p)   :: input_(n)
        logical*1 :: mask_ (n)
        integer*8 :: i, j

#ifdef IFORT
        integer*8 :: k
        !$omp parallel do private(i, j, input_, mask_)
        do k = 1, size(input, 1) * size(input, 3)
            i = mod(k - 1, size(input, 1)) + 1
            j = (k - 1) / size(input, 1) + 1
#else
        !$omp parallel do collapse(2) private(input_, mask_)
        do j = 1, size(input, 3)
            do i = 1, size(input, 1)
#endif
                input_ = input(i,:,j)
                mask_  = mask (i,:,j)
                output(i,j) = median_nocopy(input_, .true., mask_)
            end do
#ifndef IFORT
        end do
#endif
        !$omp end parallel do
    end subroutine median_mask_axis

end module math
