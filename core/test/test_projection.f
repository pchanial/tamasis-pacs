program test_projection

    use module_projection
    implicit none

    real*8, parameter    :: square(2,4) = &
        reshape([0d0, 0d0, 1d0, 0d0, 1d0, 1d0, 0d0, 1d0], [2,4])
    real*8, parameter    :: poly(2,5)   =                                      &
        reshape([-0.1d0, -0.1d0, 1.01d0, -0.1d0, 2.5d0, 2.0d0, 0.9d0, 2.1d0,   &
        0.01d0, 3.5d0], [2,5])
    real*8, parameter    :: points(2,5) =                                      &
        reshape([1.d0, 1.d0, 1.d0, 0.d0, 2.d0, 2.d0, 0.3d0, 0.3d0, 0.d0, 1d0], &
        [2,5])
    integer              :: index(5), i, result
    integer, allocatable :: iconvex(:)

    if (.not. surface_parallelogram([0.1d0, 0.1d0], [2.0d0, 1.0d0], [1.0d0, 1.0d0]) >  0) stop 'FAILED: surface_parallelogram1'
    if (.not. surface_parallelogram([0.0d0, 1.0d0], [2.0d0, 1.0d0], [1.0d0, 1.0d0]) == 0) stop 'FAILED: surface_parallelogram2'
    if (.not. surface_parallelogram([0.0d0, 0.0d0], [0.0d0, 1.0d0], [1.0d0, 1.0d0]) <  0) stop 'FAILED: surface_parallelogram3'
    if (abs(surface_convex_polygon(reshape([-4d0, 3d0, 2d0, 5d0, 5d0, 1d0], [2,3]))) /= 15d0) stop 'FAILED: surface_convex_polygon'

    ! ordering points
    call set_pivot(points(:,2))
    call qsorti_point(points, index)
    if (any(index .ne. [2,3,1,5,4])) stop "FAILED: ordering points"

    ! ordering poly
    call set_pivot(poly(:,1))
    call qsorti_point(poly, index)
    if (any(index .ne. [1,2,3,4,5])) stop "FAILED: ordering poly"

    if (point_in_polygon([0.1d0, 0.2d0], square) /=  1) stop 'FAILED: square1'
    if (point_in_polygon([1.2d0, 0.5d0], square) /= -1) stop 'FAILED: square2'
    if (point_in_polygon([0.5d0, 1.2d0], square) /= -1) stop 'FAILED: square3'
    if (point_in_polygon([-1.d0,-3.0d0], square) /= -1) stop 'FAILED: square4'
    if (point_in_polygon([ 1.d0, 0.4d0], square) /=  0) stop 'FAILED: square5'
    if (point_in_polygon([0.3d0, 0.0d0], square) /=  0) stop 'FAILED: square6'

    if (point_in_polygon([1.0d0,-0.2d0], poly) /= -1) stop 'FAILED: conc1'
    if (point_in_polygon([1.1d0, 2.5d0], poly) /= -1) stop 'FAILED: conc2'
    if (point_in_polygon([3.0d0, 3.0d0], poly) /= -1) stop 'FAILED: conc3'
    if (point_in_polygon([2.5d0, 2.0d0], poly) /=  0) stop 'FAILED: conc4'

    call convex_hull(poly, iconvex)
    if (any(iconvex .ne. [5,1,2,3])) stop "FAILED: convex hull"

    do i = 1, size(poly,2)
        result = point_in_polygon(poly(:,i), poly(:,iconvex))
        if (result == 0 .and. i /= 4 .or. result == 1 .and. i == 4) cycle
        stop 'FAILED: pip'
    end do
    if (point_in_polygon([10.d0,1.d0], poly(:,iconvex)) /= -1) stop 'FAILED pip'
    stop 'OK.'

end program test_projection
