module module_pointingmatrix

    use module_math,       only : NaN, nint_down, nint_up
    use module_projection, only : intersection_polygon_unity_square
    use module_precision,  only : p, sp
    implicit none
    private

    public :: pointingelement
    public :: pmatrix_direct
    public :: pmatrix_transpose
    public :: pmatrix_ptp
    public :: xy2roi
    public :: roi2pmatrix
    public :: backprojection_weighted
    public :: backprojection_weighted_roi

    type pointingelement
        real*4    :: weight
        integer*4 :: pixel
    end type pointingelement

contains


    subroutine pmatrix_direct(pmatrix, map, timeline)

        type(pointingelement), intent(in) :: pmatrix(:,:,:)
        real(kind=p), intent(in)          :: map(0:)
        real(kind=p), intent(inout)       :: timeline(:,:)
        integer                           :: ipixel, isample, idetector, npixels_per_sample, nsamples, ndetectors

        npixels_per_sample = size(pmatrix,1)
        nsamples   = size(pmatrix, 2)
        ndetectors = size(pmatrix, 3)

        !$omp parallel do private(idetector, isample, ipixel)
        do idetector = 1, ndetectors
            do isample = 1, nsamples
!!$                timeline(isample,idetector) = sum(map(pmatrix(:,isample,idetector)%pixel) * pmatrix(:,isample,idetector)%weight)
                timeline(isample,idetector) = 0._p
                do ipixel = 1, npixels_per_sample
                    if (pmatrix(ipixel,isample,idetector)%pixel == -1) exit
                    timeline(isample,idetector) = timeline(isample,idetector) + map(pmatrix(ipixel,isample,idetector)%pixel) *     &
                        pmatrix(ipixel,isample,idetector)%weight
                end do
            end do
        end do
        !$omp end parallel do

    end subroutine pmatrix_direct


    !-------------------------------------------------------------------------------------------------------------------------------


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
                    if (pmatrix(ipixel,isample,idetector)%pixel == -1) exit
                    map(pmatrix(ipixel,isample,idetector)%pixel) = map(pmatrix(ipixel,isample,idetector)%pixel) +                  &
                        pmatrix(ipixel,isample,idetector)%weight * timeline(isample,idetector)
                end do
            end do
        end do
        !$omp end parallel do

    end subroutine pmatrix_transpose


    !-------------------------------------------------------------------------------------------------------------------------------
   
   
    subroutine pmatrix_ptp(pmatrix, ptp)
        type(pointingelement), intent(in) :: pmatrix(:,:,:)
        real(kind=p), intent(out)         :: ptp(0:,0:)
        integer                           :: idetector, isample
        integer                           :: ipixel, jpixel, i, j
        integer                           :: npixels, nsamples, ndetectors
        real(kind(pmatrix%weight))        :: pi, pj
       
        npixels    = size(pmatrix, 1)
        nsamples   = size(pmatrix, 2)
        ndetectors = size(pmatrix, 3)
       
        ptp = 0
        !$omp parallel do reduction(+:ptp) private(idetector, isample, ipixel, jpixel, i, j, pi, pj)
        do idetector = 1, ndetectors
            do isample = 1, nsamples
                do ipixel = 1, npixels
                    if (pmatrix(ipixel,isample,idetector)%pixel == -1) exit
                    i  = pmatrix(ipixel,isample,idetector)%pixel
                    pi = pmatrix(ipixel,isample,idetector)%weight
                    do jpixel = 1, npixels
                        if (pmatrix(jpixel,isample,idetector)%pixel == -1) exit
                        j  = pmatrix(jpixel,isample,idetector)%pixel
                        pj = pmatrix(jpixel,isample,idetector)%weight
                        ptp(i,j) = ptp(i,j) + pi * pj
                    end do
                end do
            end do
        end do
        !$omp end parallel do

    end subroutine pmatrix_ptp


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine backprojection_weighted(pmatrix, timeline, mask, map, weight, threshold)
        type(pointingelement), intent(in)     :: pmatrix(:,:,:)
        real(kind=p), intent(in)              :: timeline(:,:)
        real(kind=p), intent(out)             :: map(0:)
        real(kind=p), intent(out)             :: weight(0:)
        logical(kind=1), intent(in), optional :: mask(:,:)
        real(kind=p), intent(in), optional    :: threshold
        integer                               :: npixels, nsamples, ndetectors
        integer                               :: ipixel, isample, idetector,imap
        real(kind=p)                          :: threshold_
          
        logical :: domask

        npixels    = size(pmatrix, 1)
        nsamples   = size(pmatrix, 2)
        ndetectors = size(pmatrix, 3)
        domask     = present(mask)

        map    = 0.d0
        weight = 0.d0
        !$omp parallel do default(shared) reduction(+:map,weight) &
        !$omp private(idetector,isample,ipixel,imap)
        do idetector = 1, ndetectors
            do isample = 1, nsamples
                if (domask) then
                   if (mask(isample,idetector)) cycle
                end if
                do ipixel = 1, npixels
                    imap = pmatrix(ipixel,isample,idetector)%pixel
                    if (imap == -1) exit
                    map   (imap) = map   (imap) + pmatrix(ipixel,isample,idetector)%weight * timeline(isample,idetector)
                    weight(imap) = weight(imap) + pmatrix(ipixel,isample,idetector)%weight
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
            map = NaN
        end where

    end subroutine backprojection_weighted


    !-------------------------------------------------------------------------------------------------------------------------------


    ! backproject a frame into a minimap
    subroutine backprojection_weighted_roi(pmatrix, nx, timeline, mask, itime, map, roi)
        type(pointingelement), intent(in) :: pmatrix(:,:,:)
        integer, intent(in)               :: nx
        real(kind=p), intent(in)          :: timeline(:,:)
        logical(kind=1), intent(in)       :: mask(:,:)
        integer, intent(in)               :: itime
        real(kind=p), intent(out)         :: map(0:)
        integer, intent(out)              :: roi(2,2)
        real(kind=sp)                     :: weight(0:ubound(map,1))
        integer :: nxmap, xmap, ymap, imap
        integer :: ndetectors, npixels_per_sample, idetector, ipixel

        npixels_per_sample = size(pmatrix,1)
        ndetectors = size(pmatrix,3)

        roi(1,1) = minval(modulo(pmatrix(:,itime,:)%pixel,nx), pmatrix(:,itime,:)%pixel /= -1) + 1 ! xmin
        roi(1,2) = maxval(modulo(pmatrix(:,itime,:)%pixel,nx), pmatrix(:,itime,:)%pixel /= -1) + 1 ! xmax
        roi(2,1) = minval(pmatrix(:,itime,:)%pixel / nx, pmatrix(:,itime,:)%pixel /= -1) + 1       ! ymin
        roi(2,2) = maxval(pmatrix(:,itime,:)%pixel / nx, pmatrix(:,itime,:)%pixel /= -1) + 1       ! ymax

        nxmap = roi(1,2) - roi(1,1) + 1

        ! backprojection of the timeline and weights
        map = 0
        weight = 0
        do idetector = 1, ndetectors

            if (mask(itime,idetector)) cycle

            do ipixel = 1, npixels_per_sample

                if (pmatrix(ipixel,itime,idetector)%pixel == -1) exit

                xmap = mod(pmatrix(ipixel,itime,idetector)%pixel, nx) - roi(1,1) + 1
                ymap = pmatrix(ipixel,itime,idetector)%pixel / nx     - roi(2,1) + 1
                imap = xmap + ymap * nxmap
                map(imap) = map(imap) + timeline(itime,idetector) * pmatrix(ipixel,itime,idetector)%weight
                weight(imap) = weight(imap) + pmatrix(ipixel,itime,idetector)%weight

            end do

        end do

        map = map / weight
        where (weight <= 0)
            map = NaN
        end where

    end subroutine backprojection_weighted_roi


    !-------------------------------------------------------------------------------------------------------------------------------


    ! roi is a 3-dimensional array: [1=x|2=y,1=min|2=max,idetector]
    function xy2roi(xy, nvertices) result(roi)
        real*8, intent(in)  :: xy(:,:)
        integer, intent(in) :: nvertices
        real*8              :: roi(size(xy,1),2,size(xy,2)/nvertices)
        integer             :: idetector

        do idetector = 1, size(xy,2) / nvertices
            roi(1,1,idetector) = nint_up  (minval(xy(1,nvertices * (idetector-1)+1:nvertices*idetector)))
            roi(1,2,idetector) = nint_down(maxval(xy(1,nvertices * (idetector-1)+1:nvertices*idetector)))
            roi(2,1,idetector) = nint_up  (minval(xy(2,nvertices * (idetector-1)+1:nvertices*idetector)))
            roi(2,2,idetector) = nint_down(maxval(xy(2,nvertices * (idetector-1)+1:nvertices*idetector)))
        end do

    end function xy2roi


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine roi2pmatrix(roi, nvertices, coords, nx, ny, itime, nroi, out, pmatrix)
        integer, intent(in)                  :: roi(:,:,:)
        integer, intent(in)                  :: nvertices
        real*8, intent(in)                   :: coords(:,:)
        type(pointingelement), intent(inout) :: pmatrix(:,:,:)
        integer, intent(in)                  :: nx, ny, itime
        integer, intent(inout)               :: nroi
        logical, intent(inout)               :: out
        real*8                               :: polygon(size(roi,1),nvertices)

        integer                              :: npixels_per_sample, idetector, ix, iy, iroi
        real*4                               :: weight

        npixels_per_sample = size(pmatrix, 1)
        do idetector = 1, size(pmatrix,3)

            if (roi(2,1,idetector) < 1 .or. roi(2,2,idetector) > ny .or. roi(1,1,idetector) < 1 .or. roi(1,2,idetector) > nx) then
               out = .true.
            end if

            iroi = 1
            do iy = max(roi(2,1,idetector),1), min(roi(2,2,idetector),ny)

                do ix = max(roi(1,1,idetector),1), min(roi(1,2,idetector),nx)

                    polygon(1,:) = coords(1,(idetector-1)*nvertices+1:idetector*nvertices) - (ix-0.5d0)
                    polygon(2,:) = coords(2,(idetector-1)*nvertices+1:idetector*nvertices) - (iy-0.5d0)
                    weight = abs(intersection_polygon_unity_square(polygon, nvertices))
                    if (weight <= 0) cycle
                    if (iroi <= npixels_per_sample) then
                        pmatrix(iroi,itime,idetector)%pixel  = ix - 1 + (iy - 1) * nx
                        pmatrix(iroi,itime,idetector)%weight = weight
                    end if
                    iroi = iroi + 1

                end do

            end do

            ! fill the rest of the pointing matrix
            pmatrix(iroi:,itime,idetector)%pixel  = -1
            pmatrix(iroi:,itime,idetector)%weight = 0
            nroi = max(nroi, iroi-1)

        end do

    end subroutine roi2pmatrix


end module module_pointingmatrix
