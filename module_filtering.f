module module_filtering

    use module_precision, only : p
    implicit none
    private

    public :: filterset

    type filterset
        integer                :: ncorrelations
        integer                :: ndetectors
        integer                :: nslices
        integer                :: bandwidth
        integer*8, allocatable :: first(:), last(:)
        real(p), allocatable   :: data(:,:,:)
    end type filterset

end module module_filtering
