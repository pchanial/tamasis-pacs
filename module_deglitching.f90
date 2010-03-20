module module_deglitching

    use precision, only : p, sp
    use module_math, only : moment
    use module_pointingelement, only : pointingelement
    implicit none
    private

    public :: deglitch_photproject


contains


    subroutine photproject_one_frame(pmatrix, nx, timeline, mask, itime,   &
                                     map, xmin, xmax, ymin, ymax)
        type(pointingelement), intent(in) :: pmatrix(:,:,:)
        integer, intent(in)               :: nx
        real(kind=p), intent(in)          :: timeline(:,:)
        logical, intent(in)               :: mask(:,:)
        integer, intent(in)               :: itime
        real(kind=p), intent(out)         :: map(0:)
        integer, intent(out)              :: xmin, xmax, ymin, ymax
        real(kind=sp)                     :: weight(0:size(map)-1)
        integer :: nxmap, xmap, ymap, imap
        integer :: ndetectors, npixels_per_sample, idetector, ipixel

        npixels_per_sample = size(pmatrix,1)
        ndetectors = size(pmatrix,3)

        xmin = minval(modulo(pmatrix(:,itime,:)%pixel,nx)) + 1
        xmax = maxval(modulo(pmatrix(:,itime,:)%pixel,nx)) + 1
        ymin = minval(pmatrix(:,itime,:)%pixel / nx) + 1
        ymax = maxval(pmatrix(:,itime,:)%pixel / nx) + 1

        nxmap = xmax - xmin + 1

        ! backprojection of the timeline and weights
        map = 0
        weight = 0
        do idetector = 1, ndetectors

            if (mask(itime,idetector)) cycle

            do ipixel = 1, npixels_per_sample
                xmap = modulo(pmatrix(ipixel,itime,idetector)%pixel, nx) - xmin
                ymap = pmatrix(ipixel,itime,idetector)%pixel / nx - ymin
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

    end subroutine photproject_one_frame


    !---------------------------------------------------------------------------


    subroutine deglitch_photproject(pmatrix, nx, ny, timeline, mask, &
                                    npixels_per_frame, nsigma)
        type(pointingelement), intent(in) :: pmatrix(:,:,:)
        integer, intent(in)               :: nx, ny
        real(kind=p), intent(in)          :: timeline(:,:)
        logical, intent(inout)            :: mask(:,:)
        integer, intent(in)               :: npixels_per_frame
        real(kind=p), intent(in)          :: nsigma
        integer, parameter                :: min_sample_size = 5
        integer         :: npixels_per_sample, ndetectors, ntimes
        real(kind=p)    :: map(0:npixels_per_frame-1,size(timeline,1))
        integer         :: hitmap(0:nx*ny-1)
        integer         :: xmin(size(timeline,1)), xmax(size(timeline,1))
        integer         :: ymin(size(timeline,1)), ymax(size(timeline,1))
        integer         :: i, j, ipixel, itime, idetector, isample, imap, iv
        integer         :: nv, nhits_max
        real(kind=p)    :: mean, stddev
        real(kind=p), allocatable  :: arrv(:)
        integer, allocatable :: arrt(:)
        logical, allocatable :: isglitch(:)
        
        npixels_per_sample = size(pmatrix,1)
        ntimes     = size(pmatrix,2)
        ndetectors = size(pmatrix,3)

        !$omp parallel do
        do itime=1, ntimes
            call photproject_one_frame(pmatrix, nx, timeline, mask, itime,     &
                                       map(:,itime), xmin(itime), xmax(itime), &
                                       ymin(itime), ymax(itime))
        end do
        !$omp end parallel do

        ! compute hit map, to determine the maximum number of hits
        hitmap = 0
        !$omp parallel do reduction(+:hitmap) private(idetector,isample,ipixel)
        do idetector = 1, ndetectors
            s: do isample = 1, ntimes
                do ipixel = 1, npixels_per_sample
                    if (pmatrix(ipixel,isample,idetector)%weight <= 0) cycle s
                    hitmap(pmatrix(ipixel,isample,idetector)%pixel) =     &
                        hitmap(pmatrix(ipixel,isample,idetector)%pixel) + 1
                end do
            end do s
        end do
        !$omp end parallel do

        nhits_max = maxval(hitmap)
        allocate(arrv(nhits_max))
        allocate(arrt(nhits_max))
        allocate(isglitch(nhits_max))

        !$omp parallel do private(i,j,itime,imap,nv,ipixel,arrv,arrt,isglitch)
        do j=1, ny
            do i=1, nx
                imap = (i-1) + (j-1)*nx
                nv = 0
                do itime=1, ntimes
                    if (i < xmin(itime) .or. i > xmax(itime) .or. &
                        j < ymin(itime) .or. j > ymax(itime)) cycle
                    nv = nv + 1
                    arrv(nv) = map(i-xmin(itime) + (xmax(itime)-xmin(itime)+1) * (j-ymin(itime)),itime)
                    arrt(nv) = itime
                end do

                ! check we have enough samples
                if (nv < min_sample_size) cycle

                ! sigma clip on arrv
                call moment(arrv(1:nv), mean, stddev=stddev)
                isglitch(1:nv) = abs(arrv(1:nv)-mean) > nsigma * stddev

                ! update mask
                do iv=1, nv
                    if (.not. isglitch(iv)) cycle
                    do idetector=1, ndetectors
                        do ipixel=1, npixels_per_sample
                            if (pmatrix(ipixel,arrt(iv),idetector)%pixel == imap) then
                                mask(arrt(iv),idetector) = .true.
                            end if
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        
    end subroutine deglitch_photproject

end module module_deglitching
