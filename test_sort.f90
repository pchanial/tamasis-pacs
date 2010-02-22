program test_sort

    use module_sort

    implicit none

    real*8  :: a(10)
    integer :: b(10)
    integer :: index(10)

    a = [ 2.d0, 4.d0, 2.4d0, 1.3d0, -3.d0, 10.d0, 8.d0, 100.d0, 9.d0, 7.d0 ]
    b = [ 2, 4, 2, 1, 3, -8, -203, 990, 0, 3]
    call QSORTID(a, index)
    write (*,*) a(index)

    call QSORTI(a, index)
    write (*,*) a(index)

    call QSORTI(b, index)
    write (*,*) b(index)

end program test_sort
