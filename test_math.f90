program test_math

    use module_math, only : moment, test_real_eq
    implicit none

    real*8, parameter :: sample(5) = [0.204803265708662d0,                 &
            0.911376768702354d0,0.892714242396155d0,0.848272032735291d0,   &
            0.232127193083888d0]
    real*8, parameter :: mresults(6) = [0.617858700525270d0,               &
            0.133548079242052d0,-0.285146600004998d0,-2.246479349723560d0, &
            0.365442306311204d0,0.319514776903196d0]
    real*8 :: mean, variance, skewness, kurtosis, stddev, meandev, nan, zero

    !XXX ISO_IEEE...
    zero = 0.d0
    NaN  = 0.d0 / zero

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

    call moment([2.1d0], mean, variance, skewness, kurtosis, stddev=stddev, &
                meandev=meandev)
    if (any(.not.test_real_eq([mean,variance,skewness,kurtosis,stddev,meandev],&
                              [2.1d0,NaN,NaN,NaN,NaN,0.d0], 14))) then
        stop 'FAILED: moment 3'
    end if

    call moment([2.1d0,2.1d0], mean, variance, skewness, kurtosis, &
                stddev=stddev, meandev=meandev)
    if (any(.not.test_real_eq([mean,variance,skewness,kurtosis,stddev,meandev],&
                              [2.1d0,0.d0,NaN,NaN,0.d0,0.d0], 14))) then
        stop 'FAILED: moment 4'
    end if

    stop 'OK.'

end program test_math
