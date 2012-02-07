module processing

    use module_tamasis, only : p
    use module_math,    only : mInf, pInf
    implicit none

contains

    subroutine filter_nonfinite_inplace(signal, nsamples)

        integer, intent(in)    :: nsamples
        real(p), intent(inout) :: signal(nsamples)

        !$omp parallel workshare
        where (signal /= signal .or. signal == mInf .or. signal == pInf)
            signal = 0._p
        end where
        !$omp end parallel workshare

    end subroutine filter_nonfinite_inplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine filter_nonfinite_outplace(signal, output, nsamples)

        integer, intent(in)      :: nsamples
        real(p), intent(in)      :: signal(nsamples)
        real(p), intent(inout)   :: output(nsamples)

        !$omp parallel workshare
        where (signal /= signal .or. signal == mInf .or. signal == pInf)
            output = 0._p
        else where
            output = signal
        end where
        !$omp end parallel workshare

    end subroutine filter_nonfinite_outplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine filter_nonfinite_mask_inplace(signal, mask, nsamples)

        integer, intent(in)      :: nsamples
        real(p), intent(inout)   :: signal(nsamples)
        logical*1, intent(inout) :: mask(nsamples)

        !$omp parallel workshare
        where (signal /= signal .or. signal == mInf .or. signal == pInf)
            signal = 0._p
            mask   = .true.
        end where
        !$omp end parallel workshare

    end subroutine filter_nonfinite_mask_inplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine filter_nonfinite_mask_outplace(signal, output, mask, nsamples)

        integer, intent(in)      :: nsamples
        real(p), intent(in)      :: signal(nsamples)
        real(p), intent(inout)   :: output(nsamples)
        logical*1, intent(inout) :: mask(nsamples)

        !$omp parallel workshare
        where (signal /= signal .or. signal == mInf .or. signal == pInf)
            output = 0._p
            mask   = .true.
        else where
            output = signal
        end where
        !$omp end parallel workshare

    end subroutine filter_nonfinite_mask_outplace


end module processing
