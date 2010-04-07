module module_preprocessor

    use module_math, only : sum_kahan
    implicit none
    private

    public :: subtract_meandim1
    public :: add_vectordim2
    public :: subtract_vectordim2
    public :: multiply_vectordim2
    public :: divide_vectordim2
    public :: apply_mask


contains


    ! Kahan sum
    subroutine subtract_meandim1(data)

        real*8, intent(inout) :: data(:,:)
        integer               :: nsamples, ndetectors, idetector

        nsamples   = size(data,1)
        ndetectors = size(data,2)
        !$omp parallel do
        do idetector=1, ndetectors
            data(:,idetector) = data(:,idetector) - sum_kahan(data(:,idetector)) / nsamples
        end do
        !$omp end parallel do

    end subroutine subtract_meandim1


    !-------------------------------------------------------------------------------------------------------------------------------


    ! add a vector v(n) to a 2d array data(m,n)
    subroutine add_vectordim2(data, vector)

        real*8, intent(inout) :: data(:,:)
        real*8, intent(in)    :: vector(size(data,2))
        integer               :: nsamples, ndetectors, idetector

        nsamples   = size(data,1)
        ndetectors = size(data,2)
        !$omp parallel do
        do idetector=1, ndetectors
            data(:,idetector) = data(:,idetector) + vector(idetector)
        end do
        !$omp end parallel do

    end subroutine add_vectordim2


    !-------------------------------------------------------------------------------------------------------------------------------


    ! subtract a vector v(n) to a 2d array data(m,n)
    subroutine subtract_vectordim2(data, vector)

        real*8, intent(inout) :: data(:,:)
        real*8, intent(in)    :: vector(size(data,2))
        integer               :: nsamples, ndetectors, idetector

        nsamples   = size(data,1)
        ndetectors = size(data,2)
        !$omp parallel do
        do idetector=1, ndetectors
            data(:,idetector) = data(:,idetector) - vector(idetector)
        end do
        !$omp end parallel do

    end subroutine subtract_vectordim2


    !-------------------------------------------------------------------------------------------------------------------------------


    ! multiply a vector v(n) to a 2d array data(m,n)
    subroutine multiply_vectordim2(data, vector)

        real*8, intent(inout) :: data(:,:)
        real*8, intent(in)    :: vector(size(data,2))
        integer               :: nsamples, ndetectors, idetector

        nsamples   = size(data,1)
        ndetectors = size(data,2)
        !$omp parallel do
        do idetector=1, ndetectors
            data(:,idetector) = data(:,idetector) * vector(idetector)
        end do
        !$omp end parallel do

    end subroutine multiply_vectordim2


    !-------------------------------------------------------------------------------------------------------------------------------


    ! divide a vector v(n) to a 2d array data(m,n)
    subroutine divide_vectordim2(data, vector)

        real*8, intent(inout) :: data(:,:)
        real*8, intent(in)    :: vector(size(data,2))
        integer               :: nsamples, ndetectors, idetector

        nsamples   = size(data,1)
        ndetectors = size(data,2)
        !$omp parallel do
        do idetector=1, ndetectors
            data(:,idetector) = data(:,idetector) / vector(idetector)
        end do
        !$omp end parallel do

    end subroutine divide_vectordim2


    !-------------------------------------------------------------------------------------------------------------------------------


    ! set to zero masked values
    subroutine apply_mask(data, mask)
        real*8, intent(inout) :: data(:,:)
        logical*1, intent(in) :: mask(:,:)

        !$omp parallel workshare
        where(mask)
            data = 0.d0
        end where
        !$omp end parallel workshare

    end subroutine apply_mask


end module module_preprocessor
