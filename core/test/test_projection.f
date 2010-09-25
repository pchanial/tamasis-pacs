program test_projection

    use module_projection
    use module_tamasis, only : p
    implicit none

    real(p), parameter   :: square(2,4) = reshape([0._p, 0.0_p, 1.0_p, 0.0_p, 1.0_p, 1.0_p, 0.0_p, 1.0_p], [2,4])
    real(p), parameter   :: poly(2,5)   = reshape([-0.1_p, -0.1_p, 1.01_p, -0.1_p, 2.5_p, 2.0_p, 0.9_p, 2.1_p, 0.01_p, 3.5_p],[2,5])
    real(p), parameter   :: points(2,5) = reshape([1._p, 1.0_p, 1.0_p, 0.0_p, 2._p, 2._p, 0.3_p, 0.3_p, 0.0_p, 1.0_p], [2,5])
    integer              :: index(5), i, result
    integer, allocatable :: iconvex(:)

    if (.not. surface_parallelogram([0.1_p, 0.1_p], [2.0_p, 1.0_p], [1.0_p, 1.0_p]) >  0) stop 'FAILED: surface_parallelogram1'
    if (.not. surface_parallelogram([0.0_p, 1.0_p], [2.0_p, 1.0_p], [1.0_p, 1.0_p]) == 0) stop 'FAILED: surface_parallelogram2'
    if (.not. surface_parallelogram([0.0_p, 0.0_p], [0.0_p, 1.0_p], [1.0_p, 1.0_p]) <  0) stop 'FAILED: surface_parallelogram3'
    if (abs(surface_convex_polygon(reshape([-4.0_p, 3.0_p, 2.0_p, 5.0_p, 5.0_p, 1.0_p], [2,3]))) /= 15._p)                         &
        stop 'FAILED: surface_convex_polygon'

    ! ordering points
    call set_pivot(points(:,2))
    call qsorti_point(points, index)
    if (any(index .ne. [2,3,1,5,4])) stop "FAILED: ordering points"

    ! ordering poly
    call set_pivot(poly(:,1))
    call qsorti_point(poly, index)
    if (any(index .ne. [1,2,3,4,5])) stop "FAILED: ordering poly"

    if (point_in_polygon([0.1_p, 0.2_p], square) /=  1) stop 'FAILED: square1'
    if (point_in_polygon([1.2_p, 0.5_p], square) /= -1) stop 'FAILED: square2'
    if (point_in_polygon([0.5_p, 1.2_p], square) /= -1) stop 'FAILED: square3'
    if (point_in_polygon([-1._p,-3.0_p], square) /= -1) stop 'FAILED: square4'
    if (point_in_polygon([1.0_p, 0.4_p], square) /=  0) stop 'FAILED: square5'
    if (point_in_polygon([0.3_p, 0.0_p], square) /=  0) stop 'FAILED: square6'

    if (point_in_polygon([1.0_p,-0.2_p], poly) /= -1) stop 'FAILED: conc1'
    if (point_in_polygon([1.1_p, 2.5_p], poly) /= -1) stop 'FAILED: conc2'
    if (point_in_polygon([3.0_p, 3.0_p], poly) /= -1) stop 'FAILED: conc3'
    if (point_in_polygon([2.5_p, 2.0_p], poly) /=  0) stop 'FAILED: conc4'

    call convex_hull(poly, iconvex)
    if (any(iconvex .ne. [5,1,2,3])) stop "FAILED: convex hull"

    do i = 1, size(poly,2)
        result = point_in_polygon(poly(:,i), poly(:,iconvex))
        if (result == 0 .and. i /= 4 .or. result == 1 .and. i == 4) cycle
        stop 'FAILED: pip'
    end do
    if (point_in_polygon([10._p, 1.0_p], poly(:,iconvex)) /= -1) stop 'FAILED pip'
    stop 'OK.'

end program test_projection
