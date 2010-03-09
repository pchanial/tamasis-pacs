module module_preprocessor

    implicit none

    public :: subtract_meandim1
    public :: add_vectordim2
    public :: subtract_vectordim2
    public :: multiply_vectordim2
    public :: divide_vectordim2
    public :: apply_mask
    public :: sum_kahan

    interface sum_kahan
        module procedure sum_kahan_1d, sum_kahan_2d, sum_kahan_3d
    end interface sum_kahan


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


    !---------------------------------------------------------------------------


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


    !---------------------------------------------------------------------------


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


    !---------------------------------------------------------------------------


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


    !---------------------------------------------------------------------------


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


    !---------------------------------------------------------------------------


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


    !---------------------------------------------------------------------------


    function sum_kahan_1d(input) result(sum)
        real*8, intent(in) :: input(:)
        real*8             :: sum, c, t, y
        integer            :: i

        sum = input(1)
        c = 0.0d0
        do i = 2, size(input)
            y = input(i) - c
            t = sum + y
            c = (t - sum) - y
            sum = t
        end do
    end function sum_kahan_1d


    !---------------------------------------------------------------------------


    function sum_kahan_2d(input) result(sum)
        real*8, intent(in) :: input(:,:)
        real*8             :: sum, c, t, y
        integer            :: i

        sum = sum_kahan_1d(input(:,1))
        c = 0.0d0
        do i = 2, size(input,2)
            y = sum_kahan_1d(input(:,i)) - c
            t = sum + y
            c = (t - sum) - y
            sum = t
        end do
    end function sum_kahan_2d


    !---------------------------------------------------------------------------


    function sum_kahan_3d(input) result(sum)
        real*8, intent(in) :: input(:,:,:)
        real*8             :: sum, c, t, y
        integer            :: i

        sum = sum_kahan_2d(input(:,:,1))
        c = 0.0d0
        do i = 2, size(input,3)
            y = sum_kahan_2d(input(:,:,i)) - c
            t = sum + y
            c = (t - sum) - y
            sum = t
        end do
    end function sum_kahan_3d

end module module_preprocessor
