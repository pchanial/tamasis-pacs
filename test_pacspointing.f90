program test_pacspointing
    use module_math, only : NaN, test_real_eq, test_real_neq
    use module_pacsobservation, only : pacsobservation, init_pacsobservation
    use module_pacspointing
    implicit none

    type(pacspointing)          :: ptg
    type(pacsobservation)       :: obs(1)
    character(len=*), parameter :: filename(1) = 'tests/frames_blue.fits'
    integer                     :: i, index
    real*8                      :: time(12), ra(12), dec(12), pa(12), chop(12)
    real*8, parameter           :: timetest(6) = [1634799159.471412d0,         &
                                    1634799159.5713959d0, 1634799159.6713789d0,&
                                    1634799159.771354d0, 1634799159.871338d0,  &
                                    1634799159.9713209d0]
    real*8, parameter           :: ratest(6)  = [246.00222169646082d0,         &
                                    246.00162449475283d0, 246.00103010675886d0,&
                                    246.00043574779528d0, 245.99983813619218d0,&
                                    245.99923984440272d0]
    real*8, parameter           :: dectest(6) = [61.500650548286465d0,         &
                                    61.501139611012334d0, 61.501628617896479d0,&
                                    61.502117583112039d0, 61.502603243153743d0,&
                                    61.503088192854655d0]
    real*8, parameter           :: patest(6)  = [219.90452446277453d0,         &
                                    219.90399442217031d0,219.90348064274335d0, &
                                    219.90296688571451d0, 219.90247630714885d0,&
                                    219.90199059816473d0] 
    real*8, parameter           :: choptest(6) = [0.d0, 0.d0, 0.d0, 0.d0, 0.d0,&
                                    -0.000325298332441215088d0]
    integer                     :: status

    call init_pacsobservation(obs, filename, status)
    if (status /= 0) stop 'FAILED: init_pacsobservation'

    call ptg%init(obs, status)
    if (status /= 0) stop 'FAILED: ptg%init'

    index = 0
    do i = 1, 6
        call ptg%get_position(1, timetest(i), ra(i), dec(i), pa(i), chop(i), index)
    end do
    if (any(test_real_neq(ra  (1:6), ratest,   15))) stop 'FAILED: ra'
    if (any(test_real_neq(dec (1:6), dectest,  15))) stop 'FAILED: dec'
    if (any(test_real_neq(pa  (1:6), patest,   15))) stop 'FAILED: pa'
    if (any(test_real_neq(chop(1:6), choptest, 15))) stop 'FAILED: chop'

    call ptg%destructor()

    time = [(i * 0.5d0, i = -2, 9)]

    ! test interpolation (evenly spaced sampling)
    call ptg%init_sim([1.d0, 2.d0], [0.d0, 1.d0], [3.d0, 3.d0], [2.d0, 1.d0], [0.d0, 0.d0], status)
    if (status /= 0) stop 'FAILED: init_sim 1'

    index = 0
    do i = 1, size(time)
       call ptg%get_position(1, time(i), ra(i), dec(i), pa(i), chop(i), index)
    end do
    if (any(test_real_neq(ra,  [nan, nan,-1.0d0,-0.5d0, 0.0d0, 0.5d0, 1.0d0,   &
                                1.5d0, 2.0d0, nan, nan, nan], 10) .or.         &
            test_real_neq(dec, [nan, nan, 3.0d0, 3.0d0, 3.0d0, 3.0d0, 3.0d0,   &
                                3.0d0, 3.0d0, nan, nan, nan], 10) .or.         &
            test_real_neq(pa,  [nan, nan, 3.0d0, 2.5d0, 2.0d0, 1.5d0, 1.0d0,   &
                                0.5d0, 0.0d0, nan, nan, nan], 10)))            &
        stop 'FAILED: get_position 1'
    call ptg%destructor()

    ! test slow interpolation (unevenly spaced sampling)
    call ptg%init_sim([0.5d0, 1.5d0, 2.5d0, 3.d0], [0.d0, 1.d0, 2.d0, 2.5d0],  &
                      [3.d0, 3.d0, 3.d0, 3.d0], [3.d0, 2.d0, 1.d0, 0.5d0],     &
                      [0.d0, 0.d0, 0.d0, 1.d0], status)
    if (status /= 0) stop 'FAILED: init_sim 2'

    index = 0
    do i = 1, size(time)
        call ptg%get_position(1, time(i), ra(i), dec(i), pa(i), chop(i), index)
    end do

    ! time starts from -1 to 4.5
    if (any(test_real_neq(ra,  [nan, -1.d0, -0.5d0, 0.0d0, 0.5d0, 1.0d0, 1.5d0,&
                                2.0d0, 2.5d0, 3.0d0, 3.5d0, nan], 10) .or.     &
            test_real_neq(dec, [nan,  3.0d0, 3.0d0, 3.0d0, 3.0d0, 3.0d0, 3.0d0,&
                                3.0d0, 3.0d0, 3.0d0, 3.0d0, nan], 10) .or.     &
            test_real_neq(pa,  [nan, 4.0d0, 3.5d0, 3.0d0, 2.5d0, 2.0d0, 1.5d0, &
                                1.0d0, 0.5d0, 0.0d0, -0.5d0, nan], 10) .or.    &
            test_real_neq(chop,[nan, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
                                0.0d0, 1.0d0, 2.0d0, 3.0d0, nan], 10)))        &
       stop 'FAILED: get_position 2'

    call ptg%destructor()

    stop "OK."

end program test_pacspointing
