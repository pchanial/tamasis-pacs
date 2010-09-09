program test_preprocessor

    use module_fitstools, only : ft_read_image
    use module_math, only      : median, neq_real
    use module_preprocessor
    implicit none

    real*8 :: small(10)
    real*8, allocatable :: timeline(:), filtered(:), reference(:)
    integer             :: count1, count2, count_rate, count_max, status, i

    small = [7.1d0,1.2d0,5.d0, 3.d0, 4.d0, 5.d0, 3.d0,10.d0,22.d0,32.d0]

    allocate (filtered(10))
    filtered = small
    call median_filtering(filtered, 3)
    !print(scipy.signalmedfilt(small, kernel_size=3))
    if (any(neq_real(small-filtered,[1.2d0,5.d0,3.d0,4.d0,4.d0,4.d0,5.d0,10d0,22.d0,22.d0],15))) stop 'FAILED: medfilt1'

    filtered = small
    call median_filtering(filtered, 5)
    !print(scipy.signalmedfilt(small, kernel_size=3))
    if (any(neq_real(small-filtered,[5.d0,3.d0,4.d0,4.d0,4.d0,4.d0,5.d0,10.d0,10.d0,22.d0],15))) stop 'FAILED: medfilt2'

    deallocate (filtered)

    call ft_read_image('core/test/data/timeline_transparent_mode.fits', timeline, status)
    if (status /= 0) stop 'FAILED ft_read_image'

    allocate (filtered(size(timeline)), reference(size(timeline)))
    filtered = timeline
    
    call system_clock(count1, count_rate, count_max)
    call median_filtering(filtered, 99)
    call system_clock(count2, count_rate, count_max)
    write(*,'(f6.3,a)') real(count2-count1)/count_rate, 's'

    do i=1, size(timeline)
        reference(i) = median(timeline(max(i-49,1):min(i+49,size(timeline))))
    end do

    do i=1, size(timeline)
        if (neq_real(reference(i), timeline(i)-filtered(i), 12)) then
            print *, 'FAILED:', i, reference(i), timeline(i)-filtered(i)
        end if
    end do

    stop 'OK.'

end program test_preprocessor
