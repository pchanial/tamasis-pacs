module module_pointingmatrix

    use precision, only : p, sp
    use module_pointingelement, only : pointingelement
    implicit none
    private

    public :: pointingelement
    public :: pmatrix_direct
    public :: pmatrix_transpose
    public :: backprojection_weighted
    public :: backprojection_weighted_roi


contains


    subroutine pmatrix_direct(pmatrix, map, timeline)

        type(pointingelement), intent(in) :: pmatrix(:,:,:)
        real(kind=p), intent(in)          :: map(0:)
        real(kind=p), intent(inout)       :: timeline(:,:)
        integer                           :: idetector, isample, ipixel, npixels, nsamples, ndetectors

        npixels    = size(pmatrix, 1)
        nsamples   = size(pmatrix, 2)
        ndetectors = size(pmatrix, 3)

        !$omp parallel do private(idetector, isample, ipixel)
        do idetector = 1, ndetectors
            do isample = 1, nsamples
                timeline(isample,idetector) = 0
                do ipixel = 1, npixels
                    timeline(isample,idetector) = timeline(isample,idetector) + map(pmatrix(ipixel,isample,idetector)%pixel) * &
                        pmatrix(ipixel,isample,idetector)%weight
                end do
                !faster with ifort:
                !timeline(isample,idetector) = sum(map(pmatrix(:,isample,idetector)%pixel)* &
                !                              pmatrix(:,isample,idetector)%weight)
            end do
        end do
        !$omp end parallel do

    end subroutine pmatrix_direct


    !---------------------------------------------------------------------------


    subroutine pmatrix_transpose(pmatrix, timeline, map)
        type(pointingelement), intent(in) :: pmatrix(:,:,:)
        real(kind=p), intent(in)          :: timeline(:,:)
        real(kind=p), intent(out)         :: map(0:)
        integer                           :: idetector, isample, ipixel, npixels, nsamples, ndetectors

        npixels    = size(pmatrix, 1)
        nsamples   = size(pmatrix, 2)
        ndetectors = size(pmatrix, 3)

        map = 0
        !$omp parallel do reduction(+:map) private(idetector, isample, ipixel)
        do idetector = 1, ndetectors
            do isample = 1, nsamples
                do ipixel = 1, npixels
                    map(pmatrix(ipixel,isample,idetector)%pixel) =     &
                        map(pmatrix(ipixel,isample,idetector)%pixel) + &
                        pmatrix(ipixel,isample,idetector)%weight * timeline(isample,idetector)
                end do
            end do
        end do
        !$omp end parallel do

    end subroutine pmatrix_transpose


    !---------------------------------------------------------------------------


    subroutine backprojection_weighted(pmatrix, timeline, mask, map, threshold)
        type(pointingelement), intent(in)     :: pmatrix(:,:,:)
        real(kind=p), intent(in)              :: timeline(:,:)
        real(kind=p), intent(out)             :: map(0:)
        logical(kind=1), intent(in), optional :: mask(:,:)
        real(kind=p), intent(in), optional    :: threshold
        real(kind=p)                          :: weight(0:ubound(map,1))
        integer                               :: npixels, nsamples, ndetectors
        integer                               :: ipixel, isample, idetector,imap
        real(kind=p)                          :: threshold_

        npixels    = size(pmatrix, 1)
        nsamples   = size(pmatrix, 2)
        ndetectors = size(pmatrix, 3)

        map = 0
        weight = 0
        !$omp parallel do reduction(+:map,weight) & 
        !$omp private(idetector,isample,ipixel,imap)
        do idetector = 1, ndetectors
            do isample = 1, nsamples
                if (present(mask)) then
                   if (mask(isample,idetector)) cycle
                end if
                do ipixel = 1, npixels
                    imap = pmatrix(ipixel,isample,idetector)%pixel
                    map(imap) = map(imap) + timeline(isample,idetector) * &
                        pmatrix(ipixel,isample,idetector)%weight
                    weight(imap) = weight(imap) + &
                        pmatrix(ipixel,isample,idetector)%weight
                end do
            end do
        end do
        !$omp end parallel do

        map = map / weight

        if (present(threshold)) then
            threshold_ = threshold
        else
            threshold_ = 0.d0
        endif

        where (weight <= threshold_)
            map = 0
        end where

    end subroutine backprojection_weighted


    !---------------------------------------------------------------------------


    ! backproject a frame into a minimap
    subroutine backprojection_weighted_roi(pmatrix, nx, timeline, mask, itime, &
                                           map, roi)
        type(pointingelement), intent(in) :: pmatrix(:,:,:)
        integer, intent(in)               :: nx
        real(kind=p), intent(in)          :: timeline(:,:)
        logical(kind=1), intent(in)       :: mask(:,:)
        integer, intent(in)               :: itime
        real(kind=p), intent(out)         :: map(0:)
        integer, intent(out)              :: roi(2,2)
        real(kind=sp)                     :: weight(0:size(map)-1)
        integer :: nxmap, xmap, ymap, imap
        integer :: ndetectors, npixels_per_sample, idetector, ipixel

        npixels_per_sample = size(pmatrix,1)
        ndetectors = size(pmatrix,3)

        roi(1,1) = minval(modulo(pmatrix(:,itime,:)%pixel,nx)) + 1 ! xmin
        roi(1,2) = maxval(modulo(pmatrix(:,itime,:)%pixel,nx)) + 1 ! xmax
        roi(2,1) = minval(pmatrix(:,itime,:)%pixel / nx) + 1       ! ymin
        roi(2,2) = maxval(pmatrix(:,itime,:)%pixel / nx) + 1       ! ymax

        nxmap = roi(1,2) - roi(1,1) + 1

        ! backprojection of the timeline and weights
        map = 0
        weight = 0
        do idetector = 1, ndetectors

            if (mask(itime,idetector)) cycle

            do ipixel = 1, npixels_per_sample

                if (pmatrix(ipixel,itime,idetector)%weight <= 0) exit

                xmap = mod(pmatrix(ipixel,itime,idetector)%pixel, nx)-roi(1,1)+1
                ymap = pmatrix(ipixel,itime,idetector)%pixel / nx - roi(2,1) + 1
                imap = xmap + ymap * nxmap
                map(imap) = map(imap) + timeline(itime,idetector) * &
                            pmatrix(ipixel,itime,idetector)%weight
                weight(imap) = weight(imap) + pmatrix(ipixel,itime,idetector)%weight
            end do

        end do

        map = map / weight
        where (weight <= 0)
            map = 0
        end where

    end subroutine backprojection_weighted_roi


end module module_pointingmatrix
