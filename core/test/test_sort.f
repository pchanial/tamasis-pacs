program test_sort

    use module_fitstools, only : ft_read_image
    use module_math,      only : NaN, eq_real, neq_real
    use module_sort
    use module_tamasis,   only : p
    implicit none

    real(p) :: a(10), a_NaN(14), a_mask(14)
    integer :: b(10)
    real(p) :: rtol
    integer :: index(10), index_NaN(14), index_mask(14), i, j, nuniqs, status
    logical*1 :: mask(14)
    real(p), allocatable :: table(:), timeline(:)
    integer, allocatable :: iuniq(:), hist(:)
    integer, allocatable :: order(:), isort(:)
    integer              :: count
    integer, allocatable :: i1(:), i2(:), i3(:)

    a = [ 2._p, 4._p, 2._p, 1.3_p, -3._p, 10._p, 8._p, 100._p, 4.00000001_p, 7._p ]
    call qsortid(a, index)
    do i = 2, 10
       if (a(index(i-1)) > a(index(i))) stop "FAILED: qsortid"
    end do

    call qsorti(a, index)
    do i = 2, 10
       if (a(index(i-1)) > a(index(i))) stop "FAILED: qsorti_double"
    end do

    a_NaN = [ 2._p, 4._p, 2._p, NaN, 1.3_p, NaN, NaN, NaN, -3._p, 10._p, 8._p, 100._p, 4.00000001_p, 7._p ]
    call qsorti(a_NaN, index_NaN)
    do i = 2, 10
       if (a_NaN(index_NaN(i-1)) > a_NaN(index_NaN(i))) stop "FAILED: qsorti_double NaN 1"
    end do
    do i = 11, 14
        if (a_NaN(index_NaN(i)) == a_NaN(index_NaN(i))) stop "Failed: qsorti_double NaN 2"
    end do

    a_mask = [ 2._p, 4._p, 2._p, 0._p, 1.3_p, 0._p, NaN, 0._p, -3._p, 10._p, 8._p, 100._p, 4.00000001_p, 7._p ]
    mask =  logical(a_mask == 0,1)
    call qsorti(a_mask, mask, index_mask)
    do i = 2, 10
       if (a_mask(index_mask(i-1)) > a_mask(index_mask(i))) stop "FAILED: qsorti_double mask 1"
    end do
    do i = 11, 14
        if (.not. mask(index_mask(i)) .and. a_mask(index_mask(i)) == a_mask(index_mask(i))) stop "Failed: qsorti_double mask 2"
    end do

    rtol = 1.e-6_p
    call uniq(a, index, iuniq, rtol)
    if (size(iuniq) /= 8) stop 'FAILED: uniq1'
    do i = 2, size(iuniq)
       if (a(iuniq(i-1)) > a(iuniq(i)) .or. eq_real(a(iuniq(i-1)),a(iuniq(i)),rtol)) stop "FAILED: uniq2"
    end do
    call uniq(a_NaN, index_NaN, iuniq, rtol)
    if (size(iuniq) /= 8) stop 'FAILED: uniq1 NaN'
    do i = 2, size(iuniq)
       if (a_NaN(iuniq(i-1)) > a_NaN(iuniq(i)) .or. eq_real(a_NaN(iuniq(i-1)),a_NaN(iuniq(i)),rtol)) stop "FAILED: uniq2 NaN"
    end do

#if PRECISION_REAL != 4
    rtol = 1.e-10_p
    call uniq(a, index, iuniq, rtol)
    if (size(iuniq) /= 9) stop 'FAILED: uniq3'
    do i = 2, size(iuniq)
       if (a(iuniq(i-1)) > a(iuniq(i)) .or. eq_real(a(iuniq(i-1)),a(iuniq(i)),rtol)) stop "FAILED: uniq4"
    end do
    call uniq(a_NaN, index_NaN, iuniq, rtol)
    if (size(iuniq) /= 9) stop 'FAILED: uniq3 NaN'
    do i = 2, size(iuniq)
       if (a_NaN(iuniq(i-1)) > a_NaN(iuniq(i)) .or. eq_real(a_NaN(iuniq(i-1)),a_NaN(iuniq(i)),rtol)) stop "FAILED: uniq4 NaN"
    end do
#endif

    rtol = 1.e-6_p
    call reorder(a, index, nuniqs, table, rtol)
    if (nuniqs /= 8) stop 'FAILED: reorder1'
    do i = 1, nuniqs
        do j = 1, 10
            if (index(j) /= i) cycle
            if (neq_real(a(j), table(i), rtol)) stop 'FAILED: reorder2'
        end do
    end do
    call reorder(a_NaN, index_NaN, nuniqs, table, rtol)
    if (nuniqs /= 8) stop 'FAILED: reorder1 NaN'
    do i = 1, nuniqs
        do j = 1, 10
            if (index_NaN(j) /= i) cycle
            if (neq_real(a_NaN(j), table(i), rtol)) stop 'FAILED: reorder2 NaN 1'
        end do
    end do
    if (any(index_NaN([4,6,7,8]) /= huge(index_NaN))) stop 'FAILED: reorder 2 NaN 2'
    call reorder(a_mask, mask, index_mask, nuniqs, table, rtol)
    if (nuniqs /= 8) stop 'FAILED: reorder1 mask'
    do i = 1, nuniqs
        do j = 1, 10
            if (index_mask(j) /= i) cycle
            if (neq_real(a_mask(j), table(i), rtol)) stop 'FAILED: reorder2 mask 1'
        end do
    end do
    if (any(index_mask([4,6,7,8]) /= huge(index_mask))) stop 'FAILED: reorder 2 mask 2'

    allocate (hist(nuniqs))
    hist = histogram(index, nuniqs)
    if (any(hist /= [1,1,2,2,1,1,1,1])) stop 'FAILED: histogram1'
    hist = histogram(index_NaN, nuniqs)
    if (any(hist /= [1,1,2,2,1,1,1,1])) stop 'FAILED: histogram1 NaN'
    hist = histogram(index_mask, nuniqs)
    if (any(hist /= [1,1,2,2,1,1,1,1])) stop 'FAILED: histogram1 mask'
    deallocate (hist)

#if PRECISION_REAL != 4
    rtol = 1.e-10_p
    call reorder(a, index, nuniqs, table, rtol)
    if (nuniqs /= 9) stop 'FAILED: reorder3'
    do i = 1, nuniqs
        do j = 1, 10
            if (index(j) /= i) cycle
            if (neq_real(a(j), table(i), rtol)) stop 'FAILED: reorder4'
        end do
    end do
    call reorder(a_NaN, index_NaN, nuniqs, table, rtol)
    if (nuniqs /= 9) stop 'FAILED: reorder3 NaN'
    do i = 1, nuniqs
        do j = 1, 10
            if (index_NaN(j) /= i) cycle
            if (neq_real(a_NaN(j), table(i), rtol)) stop 'FAILED: reorder4 NaN 1'
        end do
    end do
    if (any(index_NaN([4,6,7,8]) /= huge(index_NaN))) stop 'FAILED: reorder 4 NaN 2'
    call reorder(a_mask, mask, index_mask, nuniqs, table, rtol)
    if (nuniqs /= 9) stop 'FAILED: reorder3 mask'
    do i = 1, nuniqs
        do j = 1, 10
            if (index_mask(j) /= i) cycle
            if (neq_real(a_mask(j), table(i), rtol)) stop 'FAILED: reorder4 mask 1'
        end do
    end do
    if (any(index_mask([4,6,7,8]) /= huge(index_mask))) stop 'FAILED: reorder 4 mask 2'

    allocate (hist(nuniqs))
    hist = histogram(index, nuniqs)
    if (any(hist /= [1,1,2,1,1,1,1,1,1])) stop 'FAILED: histogram2'
    deallocate (hist)

#endif

    b = [ 2, 4, 2, 1, 3, -8, -203, 990, 0, 3]
    call qsorti(b, index)
    do i = 2, 10
       if (b(index(i-1)) > b(index(i))) stop "FAILED: qsorti_double."
    end do

    call ft_read_image('core/test/data/timeline_transparent_mode.fits', timeline, status)
    if (status /= 0) stop 'FAILED ft_read_image'

#if PRECISION_REAL != 4
    rtol = 1.e-12_p
    allocate (order(size(timeline)), isort(size(timeline)))
    call reorder(timeline, order, nuniqs, table, rtol)
    call qsorti(timeline, isort)
    do i = 2, size(timeline)
        if (order(isort(i)) < order(isort(i-1))) stop 'FAILED: reorder'
        if (order(isort(i)) == order(isort(i-1))) then
            if (neq_real(timeline(isort(i)), timeline(isort(i-1)), rtol)) stop 'FAILED: reorder =='
        end if
    end do
#endif

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
