program test_sort

    use module_fitstools, only : ft_read_image
    use module_math,      only : eq_real, neq_real
    use module_sort
    implicit none

    real*8  :: a(10)
    integer :: b(10)
    integer :: index(10), i, j, precision, nuniqs, status
    real*8, allocatable  :: table(:), timeline(:)
    integer, allocatable :: iuniq(:), hist(:)
    integer, allocatable :: order(:), isort(:)
    integer              :: count
    integer, allocatable :: i1(:), i2(:), i3(:)

    a = [ 2.d0, 4.d0, 2.d0, 1.3d0, -3.d0, 10.d0, 8.d0, 100.d0, 4.00000001d0, 7.d0 ]
    b = [ 2, 4, 2, 1, 3, -8, -203, 990, 0, 3]
    call qsortid(a, index)
    do i = 2, 10
       if (a(index(i-1)) > a(index(i))) stop "FAILED: qsortid"
    end do

    call qsorti(a, index)
    do i = 2, 10
       if (a(index(i-1)) > a(index(i))) stop "FAILED: qsorti_int"
    end do

    precision = 6
    call uniq(a, index, iuniq, precision)
    if (size(iuniq) /= 8) stop 'FAILED: uniq1'
    do i = 2, size(iuniq)
       if (a(iuniq(i-1)) > a(iuniq(i)) .or. eq_real(a(iuniq(i-1)),a(iuniq(i)),precision)) stop "FAILED: uniq2"
    end do

    precision = 10
    call uniq(a, index, iuniq, precision)
    if (size(iuniq) /= 9) stop 'FAILED: uniq3'
    do i = 2, size(iuniq)
       if (a(iuniq(i-1)) > a(iuniq(i)) .or. eq_real(a(iuniq(i-1)),a(iuniq(i)),precision)) stop "FAILED: uniq4"
    end do

    precision = 6
    call reorder(a, index, nuniqs, table, precision)
    if (nuniqs /= 8) stop 'FAILED: reorder1'
    do i = 1, nuniqs
        do j = 1, 10
            if (index(j) /= i) cycle
            if (neq_real(a(j), table(i), precision)) stop 'FAILED: reorder2'
        end do
    end do

    allocate (hist(nuniqs))
    hist = histogram(index, nuniqs)
    if (any(hist /= [1,1,2,2,1,1,1,1])) stop 'FAILED: histogram1'
    deallocate (hist)

    precision = 10
    call reorder(a, index, nuniqs, table, precision)
    if (nuniqs /= 9) stop 'FAILED: reorder3'
    do i = 1, nuniqs
        do j = 1, 10
            if (index(j) /= i) cycle
            if (neq_real(a(j), table(i), precision)) stop 'FAILED: reorder4'
        end do
    end do

    allocate (hist(nuniqs))
    hist = histogram(index, nuniqs)
    if (any(hist /= [1,1,2,1,1,1,1,1,1])) stop 'FAILED: histogram2'
    deallocate (hist)

    call qsorti(b, index)
    do i = 2, 10
       if (b(index(i-1)) > b(index(i))) stop "FAILED: qsorti_double."
    end do

    call ft_read_image('tests/timeline_transparent_mode.fits', timeline, status)
    if (status /= 0) stop 'FAILED ft_read_image'

    precision = 12
    allocate (order(size(timeline)), isort(size(timeline)))
    call reorder(timeline, order, nuniqs, table, precision)
    call qsorti(timeline, isort)
    do i = 2, size(timeline)
        if (order(isort(i)) < order(isort(i-1))) stop 'FAILED: reorder'
        if (order(isort(i)) == order(isort(i-1))) then
            if (neq_real(timeline(isort(i)), timeline(isort(i-1)), precision)) stop 'FAILED: reorder =='
        end if
    end do

    call uniq([2,2,3,4,4,5,8,8,9,9,9], isort)
    if (size(isort) /= 6 .or. any(isort /= [2,3,5,6,8,11])) stop 'FAILED: uniq_int'


    call where([.true., .false., .true., .true., .false., .false.], i1, count)
    if (count /= 3 .or. any(i1 /= [1, 3, 4])) stop 'FAILED: where 1d'

    call where(reshape([.true., .false., .true., .true., .false., .false.], [3,2]), i1, i2, count)
    if (count /= 3 .or. any(i1 /= [1, 3, 1]) .or. any(i2 /= [1,1,2])) stop 'FAILED: where 2d'

    call where(reshape([.true., .false., .true., .true., .false., .false.], [3,2]), i1, count)
    if (count /= 3 .or. any(i1 /= [1, 3, 4])) stop 'FAILED: where 2d-1d'

    call where(reshape([.true., .false., .true., .true., .false., .false.], [3,2,1]), i1, i2, i3, count)
    if (count /= 3 .or. any(i1 /= [1, 3, 1]) .or. any(i2 /= [1,1,2]) .or. any(i3 /= [1,1,1])) stop 'FAILED: where 3d'

    stop 'OK.'

end program test_sort
