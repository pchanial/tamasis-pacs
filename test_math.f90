program test_math

    use precision, only : p
    use module_math, only : moment, linspace, logspace, test_real_eq
    implicit none

    real(kind=p), parameter :: sample(5) = [0.204803265708662_p,           &
            0.911376768702354_p,0.892714242396155_p,0.848272032735291_p,   &
            0.232127193083888_p]
    real(kind=p), parameter :: mresults(6) = [0.617858700525270_p,         &
            0.133548079242052_p,-0.285146600004998_p,-2.246479349723560_p, &
            0.365442306311204_p,0.319514776903196_p]
    real(kind=p) :: mean, variance, skewness, kurtosis, stddev, meandev, nan, zero
    real(kind=p), allocatable :: x(:)
    integer                   :: i

    !XXX ISO_IEEE...
    zero = 0._p
    NaN  = 0._p / zero

    call moment(sample, mean, variance, skewness, kurtosis, stddev=stddev, &
                meandev=meandev)

    if (any(.not.test_real_eq([mean,variance,skewness,kurtosis,stddev,meandev],&
                              mresults, 14))) then
        stop 'FAILED: moment 1'
    end if

    !XXX no [] constructor?
    call moment(sample(1:0), mean, variance, skewness, kurtosis, stddev=stddev,&
                meandev=meandev)
    if (any(.not.test_real_eq([mean,variance,skewness,kurtosis,stddev,meandev],&
                              NaN, 14))) then
        stop 'FAILED: moment 2'
    end if

    call moment([2.1_p], mean, variance, skewness, kurtosis, stddev=stddev, &
                meandev=meandev)
    if (any(.not.test_real_eq([mean,variance,skewness,kurtosis,stddev,meandev],&
                              [2.1_p,NaN,NaN,NaN,NaN,0._p], 14))) then
        stop 'FAILED: moment 3'
    end if

    call moment([2.1_p,2.1_p], mean, variance, skewness, kurtosis, &
                stddev=stddev, meandev=meandev)
    if (any(.not.test_real_eq([mean,variance,skewness,kurtosis,stddev,meandev],&
                              [2.1_p,0._p,NaN,NaN,0._p,0._p], 14))) then
        stop 'FAILED: moment 4'
    end if

    allocate(x(1))
    x = linspace(1._p, 2._p, 3)
    if (any(.not. test_real_eq(x, [(1._p+i/2._p, i=0,2)], 14))) then
        stop 'FAILED: linspace'
    end if

    x = logspace(1._p, 2._p, 3)
    if (any(.not. test_real_eq(x, [1._p, 1.41421356237309_p, 2._p], 14))) then
        stop 'FAILED: logspace'
    end if

    stop 'OK.'

end program test_math
