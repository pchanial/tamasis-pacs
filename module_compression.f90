module module_compression

    implicit none

    public :: compression_average_direct
    public :: compression_average_transpose


contains


    subroutine compression_average_direct(data, compressed, factor)

        real*8, intent(in)  :: data(:,:)
        real*8, intent(out) :: compressed(size(data,1)/factor,size(data,2))
        integer, intent(in) :: factor
        integer             :: nsamples, ndetectors, isample, idetector

        nsamples   = size(compressed,1)
        ndetectors = size(compressed,2)
        !$omp parallel do default(shared) private(idetector,isample)
        do idetector=1, ndetectors
            do isample = 1, nsamples
                compressed(isample, idetector) = sum(data((isample-1)*factor+1:isample*factor,idetector)) / factor
            end do
        end do
        !$omp end parallel do

    end subroutine compression_average_direct


    !-------------------------------------------------------------------------------


    subroutine compression_average_transpose(compressed, data, factor)

        real*8, intent(in)  :: compressed(:,:)
        real*8, intent(out) :: data(size(compressed,1)*factor,size(compressed,2))
        integer, intent(in) :: factor
        integer             :: nsamples, ndetectors, isample, idetector

        nsamples   = size(compressed,1)
        ndetectors = size(compressed,2)
        !$omp parallel do default(shared) private(idetector,isample)
        do idetector=1, ndetectors
            do isample = 1, nsamples
                data((isample-1)*factor+1:isample*factor,idetector) = compressed(isample, idetector) / factor
            end do
        end do
        !$omp end parallel do

    end subroutine compression_average_transpose


end module module_compression
