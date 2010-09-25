module module_preprocessor

    use module_math,    only : median, sum_kahan
    use module_sort,    only : histogram, reorder
    use module_tamasis, only : p
    implicit none
    private

    public :: subtract_meandim1
    public :: add_vectordim2
    public :: subtract_vectordim2
    public :: multiply_vectordim2
    public :: divide_vectordim2
    public :: apply_mask
    public :: median_filtering

    interface median_filtering
        module procedure median_filtering_1d_1d, median_filtering_1d_2d
    end interface median_filtering
    interface subtract_meandim1
        module procedure subtract_meandim1, subtract_meandim1_mask
    end interface subtract_meandim1

contains


    ! Kahan sum
    subroutine subtract_meandim1(data)

        real(p), intent(inout) :: data(:,:)
        integer                :: nsamples, ndetectors, idetector

        nsamples   = size(data,1)
        ndetectors = size(data,2)
        !$omp parallel do
        do idetector=1, ndetectors
            data(:,idetector) = data(:,idetector) - sum_kahan(data(:,idetector)) / nsamples
        end do
        !$omp end parallel do

    end subroutine subtract_meandim1


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine subtract_meandim1_mask(data, mask)

        real(p), intent(inout) :: data(:,:)
        logical*1, intent(in)  :: mask(:,:)
        integer                :: nsamples, ndetectors, idetector

        ndetectors = size(data,2)
        !$omp parallel do
        do idetector=1, ndetectors
            nsamples = count(.not. mask(:,idetector))
            if (nsamples == 0) cycle
            data(:,idetector) = data(:,idetector) - sum_kahan(data(:,idetector), logical(mask(:,idetector), kind(.true.)))/nsamples
        end do
        !$omp end parallel do

    end subroutine subtract_meandim1_mask


    !-------------------------------------------------------------------------------------------------------------------------------


    ! add a vector v(n) to a 2d array data(m,n)
    subroutine add_vectordim2(data, vector)

        real(p), intent(inout) :: data(:,:)
        real(p), intent(in)    :: vector(size(data,2))
        integer                :: nsamples, ndetectors, idetector

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

        real(p), intent(inout) :: data(:,:)
        real(p), intent(in)    :: vector(size(data,2))
        integer                :: nsamples, ndetectors, idetector

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

        real(p), intent(inout) :: data(:,:)
        real(p), intent(in)    :: vector(size(data,2))
        integer                :: nsamples, ndetectors, idetector

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

        real(p), intent(inout) :: data(:,:)
        real(p), intent(in)    :: vector(size(data,2))
        integer                :: nsamples, ndetectors, idetector

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

        real(p), intent(inout) :: data(:,:)
        logical*1, intent(in)  :: mask(:,:)

        !$omp parallel workshare
        where(mask)
            data = 0
        end where
        !$omp end parallel workshare

    end subroutine apply_mask


    !-------------------------------------------------------------------------------------------------------------------------------


    ! median filtering in O(1) for the window length
    ! the samples inside the window form an histogram which is updated as the window slides.
    subroutine median_filtering_1d_1d(data, length)

        real(p), intent(inout) :: data(:)
        integer, intent(in)    :: length

        integer                :: nbins, ndata, order(size(data)), i, old, new, half_minus, half_plus, ibin, irank
        integer, allocatable   :: hist(:)
        real(p), allocatable   :: table(:)
        logical                :: even

        ndata = size(data)

        if (ndata <= length) then
           data = data - median(data)
           return
        end if

        call reorder(data, order, nbins, table, 1e-12_p)

        allocate (hist(nbins))
        half_minus = min((length-1) / 2, ndata-1)
        half_plus  = min(length / 2, ndata-1)
        hist = histogram(order(1:1+half_plus), nbins)
        even = mod(1 + half_plus, 2) == 0

        call median_filtering_find_rank(hist, 1+half_plus, ibin, irank)

        data(1) = data(1) - table(ibin)
        do i = 2, ndata

            ! there are missing values on the left hand side, we'll add them one by one until we reach a sample of size 'length' 
            if (i <= half_minus + 1) then

                new = order(i+half_plus)
                hist(new) = hist(new) + 1

                if (new >= ibin) then
                    if (even) then
                        irank = irank + 1
                        if (irank > hist(ibin)) then
                            call median_filtering_next(hist, ibin, irank)
                        end if
                    end if
                else if (.not. even) then
                    irank = irank - 1
                    if (irank == 0) then
                        call median_filtering_previous(hist, ibin, irank)
                    end if
                end if
                even = .not. even

            ! there are missing values on the right hand side, we'll subtract them one by one until we reach the end
            else if (i > ndata - half_plus) then
                
                old = order(i-half_minus - 1)
                hist(old) = hist(old) - 1

                if (old >= ibin) then
                    if (.not. even) then
                        irank = irank - 1
                        if (irank == 0) then
                            call median_filtering_previous(hist, ibin, irank)
                        end if
                    end if
                else if (even) then
                    irank = irank + 1
                    if (irank > hist(ibin)) then
                        call median_filtering_next(hist, ibin, irank)
                    end if
                end if
                even = .not. even
               
            ! the full sample of size 'length' is available
            else

                old = order(i - half_minus - 1)
                new = order(i + half_plus)
                hist(old) = hist(old) - 1
                hist(new) = hist(new) + 1

                if (new > ibin .and. old < ibin) then ! (1)
                    irank = irank + 1
                    if (irank > hist(ibin)) then
                        call median_filtering_next(hist, ibin, irank)
                    end if
                else if (new < ibin .and. old > ibin) then ! (2)
                    irank = irank - 1
                    if (irank == 0) then
                        call median_filtering_previous(hist, ibin, irank)
                    end if
                else if (new == ibin .and. old < ibin) then ! (3)
                    irank = irank + 1
                else if (new < ibin .and. old == ibin) then ! (4)
                    irank = irank - 1
                    if (irank == 0) then
                        call median_filtering_previous(hist, ibin, irank)
                    end if
                else if (new > ibin .and. old == ibin .and. irank > hist(ibin)) then ! (5)
                    call median_filtering_next(hist, ibin, irank)
                end if

            end if

            data(i) = data(i) - table(ibin)

        end do

    end subroutine median_filtering_1d_1d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine median_filtering_1d_2d(data, length)

        real(p), intent(inout) :: data(:,:)
        integer, intent(in)    :: length

        integer                :: i

        !$omp parallel do
        do i = 1, size(data, 2)
            call median_filtering(data(:,i), length)
        end do
        !$omp end parallel do

    end subroutine median_filtering_1d_2d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine median_filtering_find_rank(hist, length, ibin, irank)
        
        integer, intent(in)  :: hist(:), length
        integer, intent(out) :: ibin, irank

        integer :: nbins

        nbins = size(hist)
        irank = (length + 1) / 2
        do ibin=1, nbins
            if (irank - hist(ibin) <= 0) exit
            irank = irank - hist(ibin)
        end do
        
    end subroutine median_filtering_find_rank
        

    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine median_filtering_previous(hist, ibin, irank)
        
        integer, intent(in)    :: hist(:)
        integer, intent(inout) :: ibin, irank
        
        do
            ibin = ibin - 1
            irank = hist(ibin)
            if (irank /= 0) exit
        end do

    end subroutine median_filtering_previous


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine median_filtering_next(hist, ibin, irank)
        
        integer, intent(in)    :: hist(:)
        integer, intent(inout) :: ibin, irank
        
        do
            ibin = ibin + 1
            if (hist(ibin) /= 0) exit
        end do
        irank = 1

    end subroutine median_filtering_next


end module module_preprocessor
