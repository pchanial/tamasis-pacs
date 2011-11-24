program test_pacsobservation

    use iso_fortran_env,        only : ERROR_UNIT, OUTPUT_UNIT
    use module_math,            only : neq_real
    use module_observation,     only : MaskPolicy, Observation
    use module_pacsobservation, only : PacsObservation
    use module_tamasis,         only : p, POLICY_REMOVE
    implicit none

    class(Observation), allocatable     :: obs
    class(PacsObservation), allocatable :: pacsobs
    type(MaskPolicy)                    :: policy
    character(len=*), parameter         :: filename = 'pacs/test/data/frames_blue.fits'
    character(len=255), allocatable     :: afilename(:)
    logical*1, parameter        :: FALSE = .false.
    real(p)                     :: ra(12), dec(12), pa(12), chop(12)
    real(p), parameter          :: timetest(6) = [0._p, 9.9984e-2_p, 0.199967_p, 0.299942_p, 0.399926_p, 0.499909_p]
    real(p), parameter          :: ratest(6)   = [246.00222169646082_p, 246.00162449475283_p, 246.00103010675886_p,                &
                                                  246.00043574779528_p, 245.99983813619218_p, 245.99923984440272_p]
    real(p), parameter          :: dectest(6)  = [61.500650548286465_p, 61.501139611012334_p, 61.501628617896479_p,                &
                                                  61.502117583112039_p, 61.502603243153743_p, 61.503088192854655_p]
    real(p), parameter          :: patest(6)   = [219.90452446277453_p, 219.90399442217031_p,219.90348064274335_p,                 &
                                                  219.90296688571451_p, 219.90247630714885_p, 219.90199059816473_p] 
    real(p), parameter          :: choptest(6) = [0._p, 0._p, 0._p, 0._p, 0._p, -0.000325298332441215088_p]
    integer                     :: status, itime, index
    integer                     :: first, last
 
    allocate (afilename(1))
    allocate(pacsobs)

    ! valid calls
    afilename(1) = filename
    call pacsobs%init(afilename, policy, status)
    if (status /= 0) call failure('init')
    call get_first_last(pacsobs%slice(1)%p%removed, first, last)
    if (first /= 1 .or. last /= 360) call failure('init1')

    afilename(1) = filename // '[23:23]'
    call pacsobs%init(afilename, policy, status)
    call get_first_last(pacsobs%slice(1)%p%removed, first, last)
    if (status /= 0 .or. first /= 23 .or. last /= 23) call failure('init2')

    afilename(1) = filename // '[200:]'
    call pacsobs%init(afilename, policy, status)
    call get_first_last(pacsobs%slice(1)%p%removed, first, last)
    if (status /= 0 .or. first /= 200 .or. last /= 360) call failure('init3')

    afilename(1) = filename // '[:130]'
    call pacsobs%init(afilename, policy, status)
    call get_first_last(pacsobs%slice(1)%p%removed, first, last)
    if (status /= 0 .or. first /= 1 .or. last /= 130) call failure('init4')

    afilename(1) = filename // '[:]'
    call pacsobs%init(afilename, policy, status)
    call get_first_last(pacsobs%slice(1)%p%removed, first, last)
    if (status /= 0 .or. first /= 1 .or. last /= 360) call failure('init5')

    afilename(1) = filename // '[]'
    call pacsobs%init(afilename, policy, status)
    call get_first_last(pacsobs%slice(1)%p%removed, first, last)
    if (status /= 0 .or. first /= 1 .or. last /= 360) call failure('init6')
   
    ! invalid calls
    write (OUTPUT_UNIT,'(a)') 'Testing error...'
    afilename(1) = filename // '[32:23]'
    call pacsobs%init(afilename, policy, status)
    if (status == 0) call failure('badinit1')
    write (OUTPUT_UNIT,'(a,/)') 'OK.'

    write (OUTPUT_UNIT,'(a)') 'Testing error...'
    afilename(1) = filename // '[32:361]'
    call pacsobs%init(afilename, policy, status)
    if (status == 0) call failure('badinit2')
    write (OUTPUT_UNIT,'(a,/)') 'OK.'

    write (OUTPUT_UNIT,'(a)') 'Testing error...'
    afilename(1) = filename // '[ljdf]'
    call pacsobs%init(afilename, policy, status)
    if (status == 0) call failure('badinit4')
    write (OUTPUT_UNIT,'(a,/)') 'OK.'

    write (OUTPUT_UNIT,'(a)') 'Testing error...'
    afilename(1) = filename // '[:l+]'
    call pacsobs%init(afilename, policy, status)
    if (status == 0) call failure('badinit5')
    write (OUTPUT_UNIT,'(a,/)') 'OK.'

    write (OUTPUT_UNIT,'(a)') 'Testing error...'
    afilename(1) = filename // '[370:]'
    call pacsobs%init(afilename, policy, status)
    if (status == 0) call failure('badinit6')
    write (OUTPUT_UNIT,'(a,/)') 'OK.'

    deallocate (afilename)

    write (OUTPUT_UNIT,'(a)') 'Testing error...'
    allocate(afilename(0))
    call pacsobs%init(afilename, policy, status)
    if (status == 0) call failure('ainit2')
    write (OUTPUT_UNIT,'(a,/)') 'OK.'
    deallocate (afilename)

    ! calls with array
    allocate(afilename(2))
    afilename(1) = filename // '[3:100]'
    afilename(2) = filename // '[201:300]'
    call pacsobs%init(afilename, policy, status)
    if (status /= 0 .or. sum(pacsobs%slice%nsamples) /= 360*2 .or. &
         count(.not. pacsobs%slice(1)%p%removed) + count(.not. pacsobs%slice(2)%p%removed) /= 198) stop 'FAILED:ainit1'
    call get_first_last(pacsobs%slice(1)%p%removed, first, last)
    if (first /= 3 .or. last /= 100) call failure('ainit1b')
    call get_first_last(pacsobs%slice(2)%p%removed, first, last)
    if (first /= 201 .or. last /= 300) call failure('ainit1c')
    deallocate(afilename)
 
    ! test get_position_time
    call pacsobs2obs(pacsobs, 1, obs, status)
    if (status /= 0) call failure('pacsobs2obs')

    index = 0
    do itime = 1, 6
        call obs%get_position_time(timetest(itime)-timetest(1), ra(itime), dec(itime), pa(itime), chop(itime), index)
    end do
    if (any(neq_real(ra  (1:6), ratest))) call failure('ra')
    if (any(neq_real(dec (1:6), dectest))) call failure('dec')
    if (any(neq_real(pa  (1:6), patest ))) call failure('pa')
    if (any(neq_real(chop(1:6), choptest))) call failure('chop')

    ! test reading observation
    allocate (afilename(1))
    afilename(1) = filename
    call pacsobs%init(afilename, policy, status)
    if (status /= 0) call failure('init_pacsobservation')
    if (pacsobs%band /= 'blue' .or. pacsobs%slice(1)%observing_mode /= 'prime') call failure('obs%init')

contains

    subroutine get_first_last(removed, first, last)

        logical*1, intent(in) :: removed(:)
        integer, intent(out)  :: first, last

        do first = 1, size(removed)
            if (.not. removed(first)) exit
        end do

        do last = size(removed), 1, -1
            if (.not. removed(last)) exit
        end do

    end subroutine get_first_last

    subroutine failure(errmsg)
        character(len=*), intent(in) :: errmsg
        write (ERROR_UNIT,'(a)'), 'FAILED: ' // errmsg
        stop 1
    end subroutine failure

    subroutine pacsobs2obs(pacsobs, islice, obs, status)
        class(PacsObservation), intent(in)           :: pacsobs
        integer, intent(in)                          :: islice
        class(Observation), allocatable, intent(out) :: obs
        integer, intent(out)                         :: status
        
        allocate (obs)
        call obs%init(pacsobs%slice(islice)%p%time, pacsobs%slice(islice)%p%ra, pacsobs%slice(islice)%p%dec,                       &
                      pacsobs%slice(islice)%p%pa, pacsobs%slice(islice)%p%chop, pacsobs%slice(islice)%p%masked,                    &
                      pacsobs%slice(islice)%p%removed, pacsobs%slice(islice)%compression_factor,                                   &
                      (pacsobs%slice(islice)%compression_factor - 1) / (2._p * pacsobs%slice(islice)%compression_factor), status)

    end subroutine pacsobs2obs

end program test_pacsobservation
