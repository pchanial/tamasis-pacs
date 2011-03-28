program test_observation

    use iso_fortran_env,        only : ERROR_UNIT
    use module_math,            only : NaN, neq_real
    use module_observation,     only : Observation
    use module_tamasis,         only : p
    implicit none

    class(Observation), allocatable :: obs
    logical*1, parameter            :: FALSE = .false.
    integer                         :: i, index
    real(p)                         :: time(12), ra(12), dec(12), pa(12), chop(12)
    integer                         :: status


    time = [(i * 0.5_p, i = -2, 9)]

    ! test interpolation (evenly spaced sampling)
    allocate (obs)
    call obs%init([1._p, 2._p], [0._p, 1._p], [3._p, 3._p], [2._p, 1._p], [0._p, 0._p], [FALSE, FALSE], [FALSE, FALSE], [2], [1],  &
                  [0._p], status)
    if (status /= 0) call failure('init 1')

    index = 0
    do i = 1, size(time)
       call obs%get_position_time(1, time(i), ra(i), dec(i), pa(i), chop(i), index)
    end do
    if (any(neq_real(ra,  [nan, nan,-1.0_p,-0.5_p, 0.0_p, 0.5_p, 1.0_p, 1.5_p, 2.0_p, nan, nan, nan]) .or.                         &
            neq_real(dec, [nan, nan, 3.0_p, 3.0_p, 3.0_p, 3.0_p, 3.0_p, 3.0_p, 3.0_p, nan, nan, nan]) .or.                         &
            neq_real(pa,  [nan, nan, 3.0_p, 2.5_p, 2.0_p, 1.5_p, 1.0_p, 0.5_p, 0.0_p, nan, nan, nan])))                            &
        call failure('get_position_time 1')

    ! test slow interpolation (unevenly spaced sampling)
    deallocate (obs%slice, obs%pointing)
    call obs%init([0.5_p, 1.5_p, 2.5_p, 3._p], [0._p, 1._p, 2._p, 2.5_p], [3._p, 3._p, 3._p, 3._p],                                &
                  [3._p, 2._p, 1._p, 0.5_p], [0._p, 0._p, 0._p, 1._p], [(FALSE, i=1, 4)], [(FALSE, i=1, 4)], [4],                  &
                  [1], [0._p], status)
    if (status /= 0) call failure('init 2')

    index = 0
    do i = 1, size(time)
        call obs%get_position_time(1, time(i), ra(i), dec(i), pa(i), chop(i), index)
    end do

    ! time starts from -1 to 4.5
    if (any(neq_real(ra,  [nan, -1.0_p, -0.5_p, 0.0_p, 0.5_p, 1.0_p, 1.5_p, 2.0_p, 2.5_p, 3.0_p, nan, nan]) .or.                   &
            neq_real(dec, [nan,  3.0_p,  3.0_p, 3.0_p, 3.0_p, 3.0_p, 3.0_p, 3.0_p, 3.0_p, 3.0_p, nan, nan]) .or.                   &
            neq_real(pa,  [nan,  4.0_p,  3.5_p, 3.0_p, 2.5_p, 2.0_p, 1.5_p, 1.0_p, 0.5_p, 0.0_p, nan, nan]) .or.                   &
            neq_real(chop,[nan,  0.0_p,  0.0_p, 0.0_p, 0.0_p, 0.0_p, 0.0_p, 0.0_p, 1.0_p, 2.0_p, nan, nan]))) then
        call failure('get_position_time 2')
    end if

    do i = 1, size(time)
        call obs%get_position_index(1, i, 3, 0._p, ra(i), dec(i), pa(i), chop(i))
    end do

    if (any( neq_real(ra, [(0+i/3._p,i=0,2), (1+i/3._p,i=0,2), (2+0.5_p*i/3._p,i=0,2),(2.5_p+0.5_p*i/3._p,i=0,2)])                 &
        .or. neq_real(dec, [(3.0_p,i=1,12)])                                                                                       &
        .or. neq_real(pa, [(3-i/3._p,i=0,2), (2-i/3._p,i=0,2), (1-0.5_p*i/3._p,i=0,2), (0.5_p-0.5_p*i/3._p, i=0,2)])               &
        .or. neq_real(chop, [(0._p, i=1,6), (i/3.0_p, i=0,5)]))) call failure('get_position_index')

contains

    subroutine failure(errmsg)
        character(len=*), intent(in) :: errmsg
        write (ERROR_UNIT,'(a)'), 'FAILED: ' // errmsg
        stop 1
    end subroutine failure

end program test_observation
