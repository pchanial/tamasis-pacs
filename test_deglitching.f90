program test_deglitching

    use :: precision, only : p
    use :: module_deglitching, only : deglitch_photproject
    implicit none

    integer                    :: npixels_per_sample
    real(kind=p), dimension(5) :: x, y
    
    npixels_per_sample = 2
    x = [(i)]
    y = 1


contains


    subroutine get_pmatrix(pmatrix, x, y, detectorsize)
        type(pointingelement), intent(in) :: pmatrix(:,:)
        real(kind=p), intent(in)          :: x, y


    end subroutine get_pmatrix


end program test_deglitching

