program test_preprocessor

    use iso_fortran_env,  only : ERROR_UNIT
    use module_fitstools, only : ft_read_image
    use module_math,      only : NaN, median_copy, neq_real
    use module_preprocessor
    use module_tamasis,   only : p
    implicit none

    real(p)                :: small(10)
    real(p), allocatable   :: timeline(:), timeline_NaN(:), filtered(:), reference(:), rnd(:)
    logical*1, allocatable :: mask(:)
    integer                :: count1, count2, count_rate, count_max, status, i

    small = [7.1_p,1.2_p,5._p, 3._p, 4._p, 5._p, 3._p,10._p,22._p,32._p]

    allocate (filtered(10))
    filtered = small
    call median_filtering(filtered, 3)
    !print(scipy.signalmedfilt(small, kernel_size=3))
    if (any(neq_real(small-filtered,[1.2_p,5._p,3._p,4._p,4._p,4._p,5._p,10._p,22._p,22._p]))) call failure('medfilt1')

    filtered = small
    call median_filtering(filtered, 5)
    !print(scipy.signalmedfilt(small, kernel_size=3))
    if (any(neq_real(small-filtered,[5._p,3._p,4._p,4._p,4._p,4._p,5._p,10._p,10._p,22._p]))) call failure('medfilt2')

    deallocate (filtered)

    small = [ NaN, NaN, 2._p, NaN, NaN, 5._p, NaN, 5._p, 4.5_p, NaN]
    call interpolate_linear(small)
    if (any(neq_real(small, [0._p, 1._p, 2._p, 3._p, 4._p, 5._p, 5._p, 5._p, 4.5_p, 4._p]))) call failure('interpolate linear')

    call ft_read_image('core/test/data/timeline_transparent_mode.fits', timeline, status)
    if (status /= 0) call failure('ft_read_image')

    allocate (filtered(size(timeline)), reference(size(timeline)), mask(size(timeline)))
    filtered = timeline
    
    call system_clock(count1, count_rate, count_max)
    call median_filtering(filtered, 99)
    call system_clock(count2, count_rate, count_max)
    write(*,'(f6.3,a)') real(count2-count1)/count_rate, 's'

    do i=1, size(timeline)
        reference(i) = median_copy(timeline(max(i-49,1):min(i+49,size(timeline))), .false.)
    end do

    do i=1, size(timeline)
        if (neq_real(reference(i)+filtered(i), timeline(i), 100._p * epsilon(1._p))) then
            print *, 'median_filtering:', i, reference(i)+filtered(i), timeline(i)
        end if
    end do
    if (any(neq_real(reference+filtered, timeline, 100._p * epsilon(1._p)))) call failure('median_filtering')

    allocate (rnd(size(timeline)/2))
    call random_number(rnd)
    allocate (timeline_NaN(size(timeline)))
    timeline_NaN = timeline
    timeline_NaN(ceiling(rnd * size(filtered))) = NaN
    filtered = timeline_NaN

    call system_clock(count1, count_rate, count_max)
    call median_filtering(filtered, 3)
    call system_clock(count2, count_rate, count_max)
    write(*,'(f6.3,a)') real(count2-count1)/count_rate, 's'

    
    do i=1, size(timeline_NaN)
        reference(i) = median_copy(timeline_NaN(max(i-1,1):min(i+1,size(timeline_NaN))), .true.)
    end do
    call interpolate_linear(reference)

    do i=1, size(timeline_NaN)
        if (neq_real(timeline_NaN(i)-reference(i), filtered(i), 100._p * epsilon(1._p))) then
            print *, 'median_filtering NaN:', i, filtered(i), reference(i)
        end if
    end do
    if (any(neq_real(timeline_NaN-reference, filtered, 100._p * epsilon(1._p)))) call failure('median_filtering mask')

    mask = .false.
    mask(ceiling(rnd * size(filtered))) = .true.
    filtered = timeline

    call system_clock(count1, count_rate, count_max)
    call median_filtering(filtered, mask, 3)
    call system_clock(count2, count_rate, count_max)
    write(*,'(f6.3,a)') real(count2-count1)/count_rate, 's'

    if (any(neq_real(timeline-reference, filtered, 100._p * epsilon(1._p)))) call failure('median_filtering mask')

contains

    subroutine failure(errmsg)
        character(len=*), intent(in) :: errmsg
        write (ERROR_UNIT,'(a)'), 'FAILED: ' // errmsg
        stop 1
    end subroutine failure

end program test_preprocessor
