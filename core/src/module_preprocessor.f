module module_preprocessor

    use module_math,    only : NaN, median, sum_kahan
    use module_sort,    only : histogram, reorder
    use module_tamasis, only : p
    implicit none
    private

    public :: add_vectordim2
    public :: apply_mask
    public :: divide_vectordim2
    public :: interpolate_linear
    public :: median_filtering
    public :: multiply_vectordim2
    public :: subtract_meandim1
    public :: subtract_vectordim2

    interface interpolate_linear
        module procedure interpolate_linear_equally_spaced_1d_1d, interpolate_linear_equally_spaced_1d_2d
    end interface interpolate_linear

    interface median_filtering
        module procedure median_filtering_mask_1d_1d, median_filtering_mask_1d_2d
        module procedure median_filtering_nomask_1d_1d, median_filtering_nomask_1d_2d
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


    subroutine interpolate_linear_equally_spaced_1d_1d(data)

        real(p), intent(inout) :: data(:)

        real(p) :: val, delta
        logical :: hasNaN
        integer :: i, j, iNaN, ival

        if (size(data) <= 1) then
            return
        end if

        val = NaN
        iNaN = 1
        hasNaN = .false.

        do i = 1, size(data)

            ! we hit a NaN, update index of first NaN if necessary
            if (data(i) /= data(i)) then
                if (hasNaN) cycle
                hasNaN = .true.
                iNaN = i
                cycle
            endif

            if (hasNaN .and. val == val) then
                ! some values need to be interpolated loop over the values which might be NaN
                delta = (data(i) - val) / (i - ival)
                do j = iNaN, i-1
                    data(j) = val + delta * (j - ival)
                end do
                hasNaN = .false.
            end if

            ! update last defined value
            val = data(i)
            ival = i

        end do

        ! deal with trailing NaN
        if (hasNaN .and. val == val .and. data(ival-1) == data(ival-1)) then
            delta = val - data(ival-1)
            do j = ival+1, size(data)
                data(j) = val + delta * (j - ival)
            end do
        end if

    end subroutine interpolate_linear_equally_spaced_1d_1d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine interpolate_linear_equally_spaced_1d_2d(data)

        real(p), intent(inout) :: data(:,:)

        integer :: i

        !$omp parallel do
        do i = 1, size(data,2)
            call interpolate_linear(data(:,i))
        end do
        !$omp end parallel do

    end subroutine interpolate_linear_equally_spaced_1d_2d


    !-------------------------------------------------------------------------------------------------------------------------------


    ! median filtering in O(1) for the window length
    ! the samples inside the window form an histogram which is updated as the window slides.
    ! the number of elements in the histogram is length minus the number of NaN value in the window.
    subroutine median_filtering_mask_1d_1d(data, mask, length)

        real(p), intent(inout) :: data(:)
        logical*1, intent(in)  :: mask(size(data))
        integer, intent(in)    :: length

        integer                :: nbins, ndata, nvalids, order(size(data)), i, old, new, half_minus, half_plus, ibin, irank
        real(p)                :: filter(size(data))
        integer, allocatable   :: hist(:)
        real(p), allocatable   :: table(:)
        logical                :: even
        integer, parameter     :: iNaN = huge(order)

        integer, allocatable :: iw(:)

        ndata = size(data)

        if (ndata <= length) then
           data = data - median(data, mask)
           return
        end if

        call reorder(data, mask, order, nbins, table, 1000._p*epsilon(data))

        allocate (hist(nbins))
        half_minus = (length-1) / 2
        half_plus  = length / 2
        hist = histogram(order(1:1+half_plus), nbins)
        nvalids = count(order(1:1+half_plus) /= iNaN)
        even = modulo(nvalids, 2) == 0
        call median_filtering_find_rank(hist, nvalids, ibin, irank)

        filter(1) = table(ibin)
        do i = 2, ndata

            if (i <= half_minus + 1) then
                old = iNaN
            else
                old = order(i-half_minus-1)
            end if
            
            if (i > ndata - half_plus) then
                new = iNaN
            else
                new = order(i+half_plus)
            end if

            ! no old value to be removed
            if (old == iNaN) then

                ! nothing to do
                if (new == iNaN) then
                    if (nvalids == 0) then
                        filter(i) = NaN
                    else
                        filter(i) = table(ibin)
                    end if
                    cycle
                end if

                hist(new) = hist(new) + 1
                nvalids = nvalids + 1

                if (nvalids == 1) then
                    ibin = new
                    irank = 1
                else if (new >= ibin) then
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

            ! no new value to be added
            else if (new == iNaN) then
                
                hist(old) = hist(old) - 1
                nvalids = nvalids - 1

                if (old >= ibin) then
                    if (hist(ibin) == 0) then
                        if (nvalids == 0) then
                            filter(i) = NaN
                            even = .true.
                            cycle
                        else if (even) then
                            call median_filtering_next(hist, ibin, irank)
                        else
                            call median_filtering_previous(hist, ibin, irank)
                        end if
                    else if (.not. even) then
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
               
            ! add new value and remove old one
            else

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

            filter(i) = table(ibin)

        end do

        call interpolate_linear(filter)
        data = data - filter

    end subroutine median_filtering_mask_1d_1d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine median_filtering_nomask_1d_1d(data, length)

        real(p), intent(inout) :: data(:)
        integer, intent(in)    :: length

        integer                :: i

        call median_filtering(data, [(logical(.false.,1), i=1, size(data))], length)

    end subroutine median_filtering_nomask_1d_1d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine median_filtering_mask_1d_2d(data, mask, length)

        real(p), intent(inout) :: data(:,:)
        logical*1, intent(in)  :: mask(size(data,1),size(data,2))
        integer, intent(in)    :: length

        integer                :: i

        !$omp parallel do
        do i = 1, size(data, 2)
            call median_filtering(data(:,i), mask(:,i), length)
        end do
        !$omp end parallel do

    end subroutine median_filtering_mask_1d_2d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine median_filtering_nomask_1d_2d(data, length)

        real(p), intent(inout) :: data(:,:)
        integer, intent(in)    :: length

        integer                :: i

        !$omp parallel do
        do i = 1, size(data, 2)
            call median_filtering(data(:,i), length)
        end do
        !$omp end parallel do

    end subroutine median_filtering_nomask_1d_2d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine median_filtering_find_rank(hist, length, ibin, irank)
        
        integer, intent(in)  :: hist(:), length
        integer, intent(out) :: ibin, irank

        integer :: nbins

        irank = (length + 1) / 2
        if (irank == 0) then
            ibin = 0
            return
        end if

        nbins = size(hist)
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
            if (irank /= 0) return
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
