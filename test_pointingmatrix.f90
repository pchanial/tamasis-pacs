program test_pointingmatrix

    use module_math,          only : sum_kahan, neq_real
    use module_pointingmatrix
    use module_precision,     only : p
    use module_projection,    only : surface_convex_polygon
    implicit none

    integer, parameter :: nvertices = 4
    real(p), allocatable, dimension(:)  :: x_vect, y_vect
    real(p), allocatable, dimension(:,:):: xy
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

    if (any(roi(:,1,:) /= 1 .or. roi(:,2,:) /= 2)) stop 'FAILED: xy2roi'

    do itime = 1, ntimes
       call roi2pmatrix(roi, nvertices, xy, nx, ny, itime, nroi, out, pmatrix)
    end do
    if (nroi > npixels_per_sample) write (*,*) 'update npixels_per_sample:',nroi
    if (out) stop 'FAILED: roi2pmatrix out'
    if (any(pmatrix%pixel /= 0 .and. pmatrix%pixel /= 1 .and. pmatrix%pixel /= 100 .and. pmatrix%pixel /= 101)) then
        stop 'FAILED: roi2pmatrix1'
    end if

    if (neq_real(sum_kahan(real(pmatrix%weight,kind=p)), surface_convex_polygon(xy(:,1:4))*ndetectors*ntimes, 8)) then
        stop 'FAILED: roi2pmatrix2'
    end if

    stop "OK."
 
end program test_pointingmatrix
