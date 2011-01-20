program test_pacspointing

    use iso_fortran_env,        only : ERROR_UNIT
    use module_math,            only : NaN, neq_real
    use module_pacsobservation, only : PacsObservation, MaskPolicy
    use module_tamasis,         only : p
    implicit none

    class(pacsobservation), allocatable :: obs
    type(MaskPolicy)            :: policy
    character(len=*), parameter :: filename(1) = 'pacs/test/data/frames_blue.fits'
    logical*1, parameter        :: FALSE = .false.
    real(p), parameter          :: rtol = 100._p * epsilon(1._p)
    integer                     :: i, index
    real(p)                     :: time(12), ra(12), dec(12), pa(12), chop(12)
    real(p), parameter          :: timetest(6) = [0._p, 9.9984e-2_p, 0.199967_p, 0.299942_p, 0.399926_p, 0.499909_p]
    real(p), parameter          :: ratest(6)   = [246.00222169646082_p, 246.00162449475283_p, 246.00103010675886_p,                &
                                                  246.00043574779528_p, 245.99983813619218_p, 245.99923984440272_p]
    real(p), parameter          :: dectest(6)  = [61.500650548286465_p, 61.501139611012334_p, 61.501628617896479_p,                &
                                                  61.502117583112039_p, 61.502603243153743_p, 61.503088192854655_p]
    real(p), parameter          :: patest(6)   = [219.90452446277453_p, 219.90399442217031_p,219.90348064274335_p,                 &
                                                  219.90296688571451_p, 219.90247630714885_p, 219.90199059816473_p] 
    real(p), parameter          :: choptest(6) = [0._p, 0._p, 0._p, 0._p, 0._p, -0.000325298332441215088_p]
    integer                     :: status

    allocate(obs)
    call obs%init(filename, policy, status)
    if (status /= 0) call failure('init_pacsobservation')

    index = 0
    do i = 1, 6
        call obs%get_position_time(1, timetest(i)-timetest(1), ra(i), dec(i), pa(i), chop(i), index)
    end do
    if (any(neq_real(ra  (1:6), ratest))) call failure('ra')
    if (any(neq_real(dec (1:6), dectest))) call failure('dec')
    if (any(neq_real(pa  (1:6), patest ))) call failure('pa')
    if (any(neq_real(chop(1:6), choptest))) call failure('chop')

    call obs%destructor()

    time = [(i * 0.5_p, i = -2, 9)]

    ! test interpolation (evenly spaced sampling)
    call obs%init_with_variables([1._p, 2._p], [0._p, 1._p], [3._p, 3._p], [2._p, 1._p], [0._p, 0._p], [FALSE, FALSE],             &
                                 [FALSE, FALSE], [2], [1], [0._p], status)
    if (status /= 0) call failure('init_with_var 1')

    index = 0
    do i = 1, size(time)
       call obs%get_position_time(1, time(i), ra(i), dec(i), pa(i), chop(i), index)
    end do
    if (any(neq_real(ra,  [nan, nan,-1.0_p,-0.5_p, 0.0_p, 0.5_p, 1.0_p, 1.5_p, 2.0_p, nan, nan, nan]) .or.                         &
            neq_real(dec, [nan, nan, 3.0_p, 3.0_p, 3.0_p, 3.0_p, 3.0_p, 3.0_p, 3.0_p, nan, nan, nan]) .or.                         &
            neq_real(pa,  [nan, nan, 3.0_p, 2.5_p, 2.0_p, 1.5_p, 1.0_p, 0.5_p, 0.0_p, nan, nan, nan])))                            &
        call failure('get_position_time 1')
    call obs%destructor()

    ! test slow interpolation (unevenly spaced sampling)
    call obs%init_with_variables([0.5_p, 1.5_p, 2.5_p, 3._p], [0._p, 1._p, 2._p, 2.5_p], [3._p, 3._p, 3._p, 3._p],                 &
                                 [3._p, 2._p, 1._p, 0.5_p], [0._p, 0._p, 0._p, 1._p], [(FALSE, i=1, 4)], [(FALSE, i=1, 4)], [4],   &
                                 [1], [0._p], status)
    if (status /= 0) call failure('init_with_var 2')

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

    call obs%destructor()

contains

    subroutine failure(errmsg)
        character(len=*), intent(in) :: errmsg
        write (ERROR_UNIT,'(a)'), 'FAILED: ' // errmsg
        stop 1
    end subroutine failure

end program test_pacspointing
