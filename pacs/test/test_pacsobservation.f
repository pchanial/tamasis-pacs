program test_pacsobservation

    use iso_fortran_env,        only : ERROR_UNIT
    use module_pacsinstrument,  only : PacsInstrument
    use module_pacsobservation, only : MaskPolicy, PacsObservation
    use module_observation,     only : Pointing
    use module_string,          only : strsection
    use module_tamasis,         only : p, POLICY_REMOVE
    implicit none

    class(PacsObservation), allocatable :: obs
    class(PacsInstrument), allocatable  :: pacs
    type(MaskPolicy)                    :: policy
    character(len=*), parameter         :: filename = 'pacs/test/data/frames_blue.fits'
    character(len=255), allocatable     :: afilename(:)
    integer                :: status, idetector
    logical*1, allocatable :: detector_mask(:,:)
    real(p), allocatable   :: signal(:,:)
    logical*1, allocatable :: mask(:,:)
    integer                :: first, last
    real(p)                :: signal_ref(360)
    logical*1              :: mask_ref(360)
 
    allocate (afilename(1))
    ! valid calls
    allocate(obs)

    afilename(1) = filename
    call obs%init(afilename, policy, status)
    if (status /= 0) call failure('init')
    call get_first_last(obs%slice(1)%p%removed, first, last)
    if (first /= 1 .or. last /= 360) call failure('init1')

    afilename(1) = filename // '[23:23]'
    call obs%init(afilename, policy, status)
    call get_first_last(obs%slice(1)%p%removed, first, last)
    if (status /= 0 .or. first /= 23 .or. last /= 23) call failure('init2')

    afilename(1) = filename // '[200:]'
    call obs%init(afilename, policy, status)
    call get_first_last(obs%slice(1)%p%removed, first, last)
    if (status /= 0 .or. first /= 200 .or. last /= 360) call failure('init3')

    afilename(1) = filename // '[:130]'
    call obs%init(afilename, policy, status)
    call get_first_last(obs%slice(1)%p%removed, first, last)
    if (status /= 0 .or. first /= 1 .or. last /= 130) call failure('init4')

    afilename(1) = filename // '[:]'
    call obs%init(afilename, policy, status)
    call get_first_last(obs%slice(1)%p%removed, first, last)
    if (status /= 0 .or. first /= 1 .or. last /= 360) call failure('init5')

    afilename(1) = filename // '[]'
    call obs%init(afilename, policy, status)
    call get_first_last(obs%slice(1)%p%removed, first, last)
    if (status /= 0 .or. first /= 1 .or. last /= 360) call failure('init6')
   
    ! invalid calls
    afilename(1) = filename // '[32:23]'
    call obs%init(afilename, policy, status)
    if (status == 0) call failure('badinit1')

    afilename(1) = filename // '[32:361]'
    call obs%init(afilename, policy, status)
    if (status == 0) call failure('badinit2')

    afilename(1) = filename // '[ljdf]'
    call obs%init(afilename, policy, status)
    if (status == 0) call failure('badinit4')

    afilename(1) = filename // '[:l+]'
    call obs%init(afilename, policy, status)
    if (status == 0) call failure('badinit5')
    deallocate (afilename)

    ! array
    allocate(afilename(2))
    afilename(1) = filename // '[3:100]'
    afilename(2) = filename // '[201:300]'
    call obs%init(afilename, policy, status)
    if (status /= 0 .or. sum(obs%slice%nsamples) /= 360*2 .or. &
         count(.not. obs%slice(1)%p%removed) + count(.not. obs%slice(2)%p%removed) /= 198) stop 'FAILED:ainit1'
    call get_first_last(obs%slice(1)%p%removed, first, last)
    if (first /= 3 .or. last /= 100) call failure('ainit1b')
    call get_first_last(obs%slice(2)%p%removed, first, last)
    if (first /= 201 .or. last /= 300) call failure('ainit1c')
    deallocate(afilename)
 
    allocate(afilename(0))
    call obs%init(afilename, policy, status)
    if (status == 0) call failure('ainit2')
    deallocate(afilename)

    ! test reading observation
    allocate (afilename(1))
    afilename(1) = filename
    call obs%init(afilename, policy, status)
    if (status /= 0) call failure('init_pacsobservation loop')
    if (obs%band /= 'blue' .or. obs%slice(1)%observing_mode /= 'prime') call failure('obs%init')
    deallocate (afilename)

    allocate (pacs)
    call pacs%read_detector_mask(obs%band, detector_mask, status,                                                                  &
         transparent_mode=obs%slice(1)%observing_mode=='transparent')
    if (status /= 0) call failure('pacs%read_detector_mask')
    call pacs%init_with_calfiles(obs%band, detector_mask, 1, status)
    if (status /= 0) call failure('pacs%init_with_calfiles')

    allocate (signal(obs%slice(1)%nsamples,pacs%ndetectors))
    allocate (mask  (obs%slice(1)%nsamples,pacs%ndetectors))
    call pacs%read(obs, signal, mask, status)

    ! get a detector that has been hit by a glitch
    do idetector = 1, pacs%ndetectors
        if (pacs%pq(1,idetector)+1 == 27.and.pacs%pq(2,idetector)+1 == 64) exit
    end do

    signal_ref = signal(:,idetector)
    mask_ref = .false.
    mask_ref([5,6,7,151]) = .true. ! (27,64)
    if (any(mask_ref .neqv. mask(:,idetector))) call failure('mask ref')
    deallocate (mask, signal)
    
    do first = 1, 360, 7
        do last = first, 360, 11
           obs%slice(1)%p%removed = .true.
           obs%slice(1)%p(first:last)%removed = .false.
           obs%slice(1)%nvalids = count(.not. obs%slice(1)%p%removed)
           obs%nvalids = obs%slice(1)%nvalids
           allocate(signal(obs%slice(1)%nvalids,pacs%ndetectors))
           allocate(mask  (obs%slice(1)%nvalids,pacs%ndetectors))
           call pacs%read(obs, signal, mask, status)
           if (status /= 0) call failure('read_pacsobservation loop')
           if (any(signal(:,idetector) /= signal_ref(first:last))) call failure('read signal')
           if (any(mask(:,idetector) .neqv. mask_ref(first:last))) then
               print *,first, last
               print *,shape(mask(:,idetector))
               print *,shape(mask_ref)
               print *,mask(:,idetector)
               print *,mask_ref(first:last)
               call failure('read mask')
           end if
           deallocate (signal)
           deallocate (mask)
        end do
    end do

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

end program test_pacsobservation
