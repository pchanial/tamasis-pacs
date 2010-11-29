program test_pointingmatrix

    use iso_fortran_env,      only : ERROR_UNIT
    use module_math,          only : sum_kahan, neq_real
    use module_pointingmatrix
    use module_tamasis,       only : p
    use module_projection,    only : surface_convex_polygon
    implicit none

    integer, parameter :: nvertices = 4
    real(p), allocatable, dimension(:)   :: x_vect, y_vect
    real(p), allocatable, dimension(:,:) :: xy
    integer, allocatable :: roi(:,:,:)
    integer i
    integer :: npixels_per_sample, ntimes, ndetectors, nroi, nx, ny, itime
    logical :: out
    type(pointingelement), allocatable :: pmatrix(:,:,:)

    ndetectors = 100
    ntimes = 1000
    npixels_per_sample = 4
    nx = 100
    ny = 200

    allocate(x_vect(ndetectors*nvertices))
    allocate(y_vect(ndetectors*nvertices))
    allocate(xy(2,ndetectors*nvertices))
    allocate(roi(2,2,ndetectors))
    allocate(pmatrix(npixels_per_sample,ntimes,ndetectors))
    x_vect = [([0.6_p,2.4_p,2.4_p,0.6_p], i=1,ndetectors)]
    y_vect = [([0.6_p,0.6_p,2.4_p,2.4_p], i=1,ndetectors)]
    xy(1,:) = x_vect
    xy(2,:) = y_vect
    roi = xy2roi(xy, 4)

    if (any(roi(:,1,:) /= 1 .or. roi(:,2,:) /= 2)) call failure('xy2roi')

    nroi = 0
    out = .false.
    do itime = 1, ntimes
       call roi2pmatrix(roi, nvertices, xy, nx, ny, nroi, out, pmatrix(:,itime,:))
    end do
    if (nroi /= npixels_per_sample) call failure('update npixels_per_sample')
    if (out) call failure('roi2pmatrix out')
    if (any(pmatrix%pixel /= 0 .and. pmatrix%pixel /= 1 .and. pmatrix%pixel /= 100 .and. pmatrix%pixel /= 101)) then
        call failure('roi2pmatrix1')
    end if

    if (neq_real(sum_kahan(real(pmatrix%weight,kind=p)), surface_convex_polygon(real(xy(:,1:4), p)) * ndetectors * ntimes,         &
        10._p * epsilon(1.0))) then
        call failure('roi2pmatrix2')
    end if

contains

    subroutine failure(errmsg)
        character(len=*), intent(in) :: errmsg
        write (ERROR_UNIT,'(a)'), 'FAILED: ' // errmsg
        stop 1
    end subroutine failure
 
end program test_pointingmatrix
