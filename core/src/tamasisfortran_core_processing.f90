module processing

    use module_tamasis, only : p
    use module_math,    only : mInf, pInf
    implicit none

contains

    subroutine filter_nonfinite_inplace(signal, nsamples)

        integer*8, intent(in)  :: nsamples
        real(p), intent(inout) :: signal(nsamples)

        integer*8              :: i
        real(p)                :: value

        !$omp parallel do private(value)
        do i = 1, nsamples
            value = signal(i)
            if (value /= value .or. value == mInf .or. value == pInf) then
                signal(i) = 0._p
            end if
        end do
        !$omp end parallel do

    end subroutine filter_nonfinite_inplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine filter_nonfinite_outplace(signal, output, nsamples)

        integer*8, intent(in)    :: nsamples
        real(p), intent(in)      :: signal(nsamples)
        real(p), intent(inout)   :: output(nsamples)

        integer*8                :: i
        real(p)                  :: value

        !$omp parallel do private(value)
        do i = 1, nsamples
            value = signal(i)
            if (value /= value .or. value == mInf .or. value == pInf) then
                output(i) = 0._p
            else
                output(i) = value
            end if
        end do
        !$omp end parallel do

    end subroutine filter_nonfinite_outplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine filter_nonfinite_mask_inplace(signal, mask, nsamples)

        integer*8, intent(in)    :: nsamples
        real(p), intent(inout)   :: signal(nsamples)
        logical*1, intent(inout) :: mask(nsamples)

        integer*8                :: i
        real(p)                  :: value

        !$omp parallel do private(value)
        do i = 1, nsamples
            value = signal(i)
            if (value /= value .or. value == mInf .or. value == pInf) then
                signal(i) = 0._p
                mask(i)   = .true.
            end if
        end do
        !$omp end parallel do

    end subroutine filter_nonfinite_mask_inplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine filter_nonfinite_mask_outplace(signal, output, mask, nsamples)

        integer*8, intent(in)    :: nsamples
        real(p), intent(in)      :: signal(nsamples)
        real(p), intent(inout)   :: output(nsamples)
        logical*1, intent(inout) :: mask(nsamples)

        integer*8                :: i
        real(p)                  :: value

        !$omp parallel do private(value)
        do i = 1, nsamples
            value = signal(i)
            if (value /= value .or. value == mInf .or. value == pInf) then
                output(i) = 0._p
                mask(i)   = .true.
            else
                output(i) = value
            end if
        end do
        !$omp end parallel do

    end subroutine filter_nonfinite_mask_outplace


end module processing
