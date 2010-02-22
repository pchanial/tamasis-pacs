program test_projection

    use module_projection
    implicit none

    real*8, parameter :: square(2,4) = reshape([0d0, 0d0, 1d0, 0d0, 1d0, 1d0, 0d0, 1d0], [2,4])
    real*8, parameter :: poly(2,5)   = reshape([-0.1d0, -0.1d0, 1.01d0, -0.1d0, 2.5d0, 2.0d0, 0.9d0, 2.1d0, 0.01d0, 3.5d0], [2,5])
    real*8, parameter :: points(2,5) = reshape([1.d0, 1.d0, 1.d0, 0.d0, 2.d0, 2.d0, 0.3d0, 0.3d0, 0.d0, 1d0], [2,5])
    integer           :: index(5), i, result
    integer, allocatable :: iconvex(:)

    write(*,*) .true. .eqv. surface_parallelogram([0.1d0, 0.1d0], [2.0d0, 1.0d0], [1.0d0, 1.0d0]) >  0
    write(*,*) .true. .eqv. surface_parallelogram([0.0d0, 1.0d0], [2.0d0, 1.0d0], [1.0d0, 1.0d0]) == 0
    write(*,*) .true. .eqv. surface_parallelogram([0.0d0, 0.0d0], [0.0d0, 1.0d0], [1.0d0, 1.0d0]) <  0
    write(*,*) abs(surface_convex_polygon(reshape([-4d0, 3d0, 2d0, 5d0, 5d0, 1d0], [2,3]))) == 15d0

    write (*,*) 'Ordering points:'
    call set_pivot(points(:,2))
    call qsorti_point(points, index)
    do i = 1, 5
       write(*,*) index(i), points(1, index(i)), points(2, index(i))
    end do

    write (*,*) 'Ordering poly:'
    call set_pivot(poly(:,1))
    call qsorti_point(poly, index)
    do i = 1, 5
       write(*,*) index(i), poly(1, index(i)), poly(2, index(i)), compare_point(1, index(i))
    end do

    write (*,*) 'Square:'
    write (*,*) point_in_polygon([0.1d0, 0.2d0], square) ==  1
    write (*,*) point_in_polygon([1.2d0, 0.5d0], square) == -1
    write (*,*) point_in_polygon([0.5d0, 1.2d0], square) == -1
    write (*,*) point_in_polygon([-1.d0,-3.0d0], square) == -1
    write (*,*) point_in_polygon([ 1.d0, 0.4d0], square) ==  0
    write (*,*) point_in_polygon([0.3d0, 0.0d0], square) ==  0

    write (*,*) 'Concave polygon:'
    write (*,*) point_in_polygon([1.0d0,-0.2d0], poly) == -1
    write (*,*) point_in_polygon([1.1d0, 2.5d0], poly) == -1
    write (*,*) point_in_polygon([3.0d0, 3.0d0], poly) == -1
    write (*,*) point_in_polygon([2.5d0, 2.0d0], poly) ==  0

    call convex_hull(poly, iconvex)
    write (*,*) 'iconvex: ', iconvex

    do i = 1, size(poly,2)
        result = point_in_polygon(poly(:,i), poly(:,iconvex))
        write (*,*) result
    end do

    stop "OK."

contains

    subroutine test1(array, n)
       integer :: array(2,n)
       integer :: n
       write(*,*) shape(array)
    end subroutine test1

end program test_projection
