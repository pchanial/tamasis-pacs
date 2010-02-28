program test_sort

    use module_sort

    implicit none

    real*8  :: a(10)
    integer :: b(10)
    integer :: index(10), i

    a = [ 2.d0, 4.d0, 2.4d0, 1.3d0, -3.d0, 10.d0, 8.d0, 100.d0, 9.d0, 7.d0 ]
    b = [ 2, 4, 2, 1, 3, -8, -203, 990, 0, 3]
    call QSORTID(a, index)
    do i = 2, 10
       if (a(index(i-1)) > a(index(i))) stop "Failed: qsortid."
    end do

    call QSORTI(a, index)
    do i = 2, 10
       if (a(index(i-1)) > a(index(i))) stop "Failed: qsorti_int."
    end do

    call QSORTI(b, index)
    do i = 2, 10
       if (b(index(i-1)) > b(index(i))) stop "Failed: qsorti_double."
    end do

    stop 'OK.'

end program test_sort
