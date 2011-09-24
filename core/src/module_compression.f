! Copyright 2010-2011 Pierre Chanial
! All rights reserved
!
module module_compression

    use module_tamasis, only : p
    implicit none

    public :: compression_average_direct
    public :: compression_average_transpose
    public :: downsampling_direct
    public :: downsampling_transpose


contains


    subroutine compression_average_direct(data, compressed, factor)

        real(p), intent(in)  :: data(:)
        integer, intent(in)  :: factor
        real(p), intent(out) :: compressed(size(data)/factor)
        integer              :: isample

        do isample = 1, size(compressed)
            compressed(isample) = sum(data((isample-1)*factor+1:isample*factor)) / factor
        end do

    end subroutine compression_average_direct


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine compression_average_transpose(compressed, data, factor)

        real(p), intent(in)  :: compressed(:)
        integer, intent(in)  :: factor
        real(p), intent(out) :: data(size(compressed)*factor)
        integer              :: isample

        do isample = 1, size(compressed)
            data((isample-1)*factor+1:isample*factor) = compressed(isample) / factor
        end do

    end subroutine compression_average_transpose


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine downsampling_direct(data, compressed, factor)

        real(p), intent(in)  :: data(:)
        integer, intent(in)  :: factor
        real(p), intent(out) :: compressed(size(data)/factor)

        compressed = data(::factor)

    end subroutine downsampling_direct


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine downsampling_transpose(compressed, data, factor)

        real(p), intent(in)  :: compressed(:)
        integer, intent(in)  :: factor
        real(p), intent(out) :: data(size(compressed)*factor)
        integer              :: isample

        do isample = 1, size(compressed)
            data((isample-1)*factor+1) = compressed(isample)
            data((isample-1)*factor+2:isample*factor) = 0
        end do

    end subroutine downsampling_transpose


end module module_compression
