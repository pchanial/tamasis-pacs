program test_math

    use module_math,      only : NaN, linspace, logspace, median, moment, nint_down, nint_up, neq_real
    use module_precision, only : p
    implicit none

    real(kind=p), parameter :: sample(5) = [0.204803265708662_p, 0.911376768702354_p,0.892714242396155_p,0.848272032735291_p,      &
                                   0.232127193083888_p]
    real(kind=p), parameter :: mresults(6) = [0.617858700525270_p, 0.133548079242052_p,-0.285146600004998_p,-2.246479349723560_p,  &
                                   0.365442306311204_p,0.319514776903196_p]
    real(kind=p) :: mean, variance, skewness, kurtosis, stddev, meandev
    real(kind=p), allocatable :: x(:)
    integer                   :: i

    call moment(sample, mean, variance, skewness, kurtosis, stddev=stddev, meandev=meandev)

    if (any(neq_real([mean,variance,skewness,kurtosis,stddev,meandev], mresults, 14))) then
        stop 'FAILED: moment 1'
    end if

    !XXX no [] constructor?
    call moment(sample(1:0), mean, variance, skewness, kurtosis, stddev=stddev, meandev=meandev)
    if (any(neq_real([mean,variance,skewness,kurtosis,stddev,meandev], NaN, 14))) then
        stop 'FAILED: moment 2'
    end if

    call moment([2.1_p], mean, variance, skewness, kurtosis, stddev=stddev, meandev=meandev)
    if (any(neq_real([mean,variance,skewness,kurtosis,stddev,meandev], [2.1_p,NaN,NaN,NaN,NaN,0._p], 14))) then
        stop 'FAILED: moment 3'
    end if

    call moment([2.1_p,2.1_p], mean, variance, skewness, kurtosis, stddev=stddev, meandev=meandev)
    if (any(neq_real([mean,variance,skewness,kurtosis,stddev,meandev], [2.1_p,0._p,NaN,NaN,0._p,0._p], 14))) then
        stop 'FAILED: moment 4'
    end if

    allocate(x(3))
    x = linspace(1._p, 2._p, 3)
    if (any(neq_real(x, [(1._p+i/2._p, i=0,2)], 14))) then
        stop 'FAILED: linspace'
    end if

    x = logspace(1._p, 2._p, 3)
    if (any(neq_real(x, [1._p, 1.41421356237309_p, 2._p], 14))) then
        stop 'FAILED: logspace'
    end if
    deallocate(x)

    allocate(x(8))
    x = [0._p,0.2_p,0.5_p,0.8_p,1._p,1.2_p,1.5_p,1.7_p]
    if (any(nint_down( x) /= [0,0,0,1,1,1,1,2])) stop 'FAILED: nint_down 1'
    if (any(nint_down(-x) /=-[0,0,1,1,1,1,2,2])) stop 'FAILED: nint_down 2'
    if (any(nint_up  ( x) /= [0,0,1,1,1,1,2,2])) stop 'FAILED: nint_up 1'
    if (any(nint_up  (-x) /=-[0,0,0,1,1,1,1,2])) stop 'FAILED: nint_up 2'
    deallocate(x)

    allocate(x(size(mresults)+1))
    x = [mresults, -999._p]
    if (median(x) /= mresults(2)) stop 'FAILED: median'
    deallocate(x)

    stop 'OK.'

end program test_math
