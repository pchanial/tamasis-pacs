module module_pointingmatrix

    use module_pointingelement

    implicit none

    public :: pmatrix_direct
    public :: pmatrix_transpose


contains


    subroutine pmatrix_direct(pmatrix, map, timeline)

        type(pointingelement), intent(in) :: pmatrix(:,:,:)
        real*8, intent(in)                :: map(0:)
        real*8, intent(inout)             :: timeline(:,:)
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
        real*8, intent(in)                :: timeline(:,:)
        real*8, intent(out)               :: map(0:)
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


end module module_pointingmatrix
