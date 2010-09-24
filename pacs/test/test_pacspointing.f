program test_pacspointing

    use module_math,            only : NaN, neq_real
    use module_pacsobservation, only : PacsObservation, MaskPolicy
    use module_precision,       only : p
    implicit none

    class(pacsobservation), allocatable :: obs
    type(MaskPolicy)            :: policy
    character(len=*), parameter :: filename(1) = 'pacs/test/data/frames_blue.fits'
    integer                     :: i, index
    real*8                      :: time(12), ra(12), dec(12), pa(12), chop(12)
    real*8, parameter           :: timetest(6) = [1634799159.471412_p, 1634799159.5713959_p, 1634799159.6713789_p,                 &
                                       1634799159.771354_p, 1634799159.871338_p, 1634799159.9713209_p]
    real*8, parameter           :: ratest(6)  = [246.00222169646082_p, 246.00162449475283_p, 246.00103010675886_p,                 &
                                       246.00043574779528_p, 245.99983813619218_p, 245.99923984440272_p]
    real*8, parameter           :: dectest(6) = [61.500650548286465_p, 61.501139611012334_p, 61.501628617896479_p,                 &
                                       61.502117583112039_p, 61.502603243153743_p, 61.503088192854655_p]
    real*8, parameter           :: patest(6)  = [219.90452446277453_p, 219.90399442217031_p,219.90348064274335_p,                  &
                                       219.90296688571451_p, 219.90247630714885_p, 219.90199059816473_p] 
    real*8, parameter           :: choptest(6) = [0._p, 0._p, 0._p, 0._p, 0._p, -0.000325298332441215088_p]
    integer                     :: status

    allocate(obs)
    call obs%init(filename, policy, status)
    if (status /= 0) stop 'FAILED: init_pacsobservation'

    index = 0
    do i = 1, 6
        call obs%get_position_time(1, timetest(i), ra(i), dec(i), pa(i), chop(i), index)
    end do
    if (any(neq_real(ra  (1:6), ratest  ))) stop 'FAILED: ra'
    if (any(neq_real(dec (1:6), dectest ))) stop 'FAILED: dec'
    if (any(neq_real(pa  (1:6), patest  ))) stop 'FAILED: pa'
    if (any(neq_real(chop(1:6), choptest))) stop 'FAILED: chop'

    call obs%destructor()

    time = [(i * 0.5_p, i = -2, 9)]

    ! test interpolation (evenly spaced sampling)
    call obs%init_sim([1._p, 2._p], [0._p, 1._p], [3._p, 3._p], [2._p, 1._p], [0._p, 0._p], status)
    if (status /= 0) stop 'FAILED: init_sim 1'

    index = 0
    do i = 1, size(time)
       call obs%get_position_time(1, time(i), ra(i), dec(i), pa(i), chop(i), index)
    end do
    if (any(neq_real(ra,  [nan, nan,-1.0_p,-0.5_p, 0.0_p, 0.5_p, 1.0_p, 1.5_p, 2.0_p, nan, nan, nan]) .or.                         &
            neq_real(dec, [nan, nan, 3.0_p, 3.0_p, 3.0_p, 3.0_p, 3.0_p, 3.0_p, 3.0_p, nan, nan, nan]) .or.                         &
            neq_real(pa,  [nan, nan, 3.0_p, 2.5_p, 2.0_p, 1.5_p, 1.0_p, 0.5_p, 0.0_p, nan, nan, nan])))                            &
        stop 'FAILED: get_position_time 1'
    call obs%destructor()

    ! test slow interpolation (unevenly spaced sampling)
    call obs%init_sim([0.5_p, 1.5_p, 2.5_p, 3._p], [0._p, 1._p, 2._p, 2.5_p], [3._p, 3._p, 3._p, 3._p], [3._p, 2._p, 1._p, 0.5_p], &
                      [0._p, 0._p, 0._p, 1._p], status)
    if (status /= 0) stop 'FAILED: init_sim 2'

    index = 0
    do i = 1, size(time)
        call obs%get_position_time(1, time(i), ra(i), dec(i), pa(i), chop(i), index)
    end do

    ! time starts from -1 to 4.5
    if (any(neq_real(ra,  [nan, -1._p, -0.5_p, 0.0_p, 0.5_p, 1.0_p, 1.5_p, 2.0_p, 2.5_p, 3.0_p, nan, nan]) .or.                    &
            neq_real(dec, [nan,  3.0_p, 3.0_p, 3.0_p, 3.0_p, 3.0_p, 3.0_p, 3.0_p, 3.0_p, 3.0_p, nan, nan]) .or.                    &
            neq_real(pa,  [nan,  4.0_p, 3.5_p, 3.0_p, 2.5_p, 2.0_p, 1.5_p, 1.0_p, 0.5_p, 0.0_p, nan, nan]) .or.                    &
            neq_real(chop,[nan,  0.0_p, 0.0_p, 0.0_p, 0.0_p, 0.0_p, 0.0_p, 0.0_p, 1.0_p, 2.0_p, nan, nan]))) then
        stop 'FAILED: get_position_time 2'
    end if

    do i = 1, size(time)
        call obs%get_position_index(1, i, 3, ra(i), dec(i), pa(i), chop(i))
    end do

    if (any(neq_real(ra, [(0._p+i/3._p,i=0,2), (1._p+i/3._p,i=0,2), (2._p+0.5_p*i/3._p,i=0,2),(2.5_p+0.5_p*i/3._p,i=0,2)])         &
        .or. neq_real(dec, [(3.0_p,i=1,12)])                                                                                       &
        .or. neq_real(pa, [(3._p-i/3._p,i=0,2), (2._p-i/3._p,i=0,2), (1.-0.5_p*i/3._p,i=0,2), (0.5_p-0.5_p*i/3._p, i=0,2)])        &
        .or. neq_real(chop, [(0.0_p, i=1,6), (0.0_p+i/3.0_p, i=0,5)]))) stop 'FAILED: get_position_index'

    call obs%destructor()

    stop "OK."

end program test_pacspointing
