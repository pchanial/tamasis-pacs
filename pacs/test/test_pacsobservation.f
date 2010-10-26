program test_pacsobservation

    use module_pacsinstrument,  only : PacsInstrument
    use module_pacsobservation, only : MaskPolicy, PacsObservation
    use module_observation,     only : Pointing
    use module_string,          only : strsection
    use module_tamasis,         only : init_tamasis, get_tamasis_path, p, POLICY_REMOVE
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
 
    ! initialise tamasis
    call init_tamasis

    allocate (afilename(1))
    ! valid calls
    allocate(obs)
    afilename(1) = get_tamasis_path() // filename

    call obs%init(afilename, policy, status)
    call get_first_last(obs%slice(1)%p%removed, first, last)
    if (status /= 0 .or. first /= 1 .or. last /= 360) stop 'FAILED: init1'

    afilename(1) = get_tamasis_path() // filename // '[23:23]'
    call obs%init(afilename, policy, status)
    call get_first_last(obs%slice(1)%p%removed, first, last)
    if (status /= 0 .or. first /= 23 .or. last /= 23) stop 'FAILED: init2'

    afilename(1) = get_tamasis_path() // filename // '[200:]'
    call obs%init(afilename, policy, status)
    call get_first_last(obs%slice(1)%p%removed, first, last)
    if (status /= 0 .or. first /= 200 .or. last /= 360) stop 'FAILED: init3'

    afilename(1) = get_tamasis_path() // filename // '[:130]'
    call obs%init(afilename, policy, status)
    call get_first_last(obs%slice(1)%p%removed, first, last)
    if (status /= 0 .or. first /= 1 .or. last /= 130) stop 'FAILED: init4'

    afilename(1) = get_tamasis_path() // filename // '[:]'
    call obs%init(afilename, policy, status)
    call get_first_last(obs%slice(1)%p%removed, first, last)
    if (status /= 0 .or. first /= 1 .or. last /= 360) stop 'FAILED: init5'

    afilename(1) = get_tamasis_path() // filename // '[]'
    call obs%init(afilename, policy, status)
    call get_first_last(obs%slice(1)%p%removed, first, last)
    if (status /= 0 .or. first /= 1 .or. last /= 360) stop 'FAILED: init6'
   
    ! invalid calls
    afilename(1) = get_tamasis_path() // filename // '[32:23]'
    call obs%init(afilename, policy, status)
    if (status == 0) stop 'FAILED: badinit1'

    afilename(1) = get_tamasis_path() // filename // '[32:361]'
    call obs%init(afilename, policy, status)
    if (status == 0) stop 'FAILED: badinit2'

    afilename(1) = get_tamasis_path() // filename // '[ljdf]'
    call obs%init(afilename, policy, status)
    if (status == 0) stop 'FAILED: badinit4'

    afilename(1) = get_tamasis_path() // filename // '[:l+]'
    call obs%init(afilename, policy, status)
    if (status == 0) stop 'FAILED: badinit5'
    deallocate (afilename)

    ! array
    allocate(afilename(2))
    afilename(1) = get_tamasis_path() // filename // '[3:100]'
    afilename(2) = get_tamasis_path() // filename // '[201:300]'
    call obs%init(afilename, policy, status)
    if (status /= 0 .or. sum(obs%slice%nsamples) /= 360*2 .or. &
         count(.not. obs%slice(1)%p%removed) + count(.not. obs%slice(2)%p%removed) /= 198) stop 'FAILED:ainit1'
    call get_first_last(obs%slice(1)%p%removed, first, last)
    if (first /= 3 .or. last /= 100) stop 'FAILED: ainit1b'
    call get_first_last(obs%slice(2)%p%removed, first, last)
    if (first /= 201 .or. last /= 300) stop 'FAILED: ainit1c'
    deallocate(afilename)
 
    allocate(afilename(0))
    call obs%init(afilename, policy, status)
    if (status == 0) stop 'FAILED: ainit2'
    deallocate(afilename)

    ! test reading observation
    allocate (afilename(1))
    afilename(1) = get_tamasis_path() // filename
    call obs%init(afilename, policy, status)
    if (status /= 0) stop 'FAILED: init_pacsobservation loop'
    if (obs%band /= 'blue' .or. obs%slice(1)%observing_mode /= 'prime') stop 'FAILED: obs%init'
    deallocate (afilename)

    allocate (pacs)
    call pacs%read_detector_mask(obs%band, detector_mask, status,                                                                  &
         transparent_mode=obs%slice(1)%observing_mode=='transparent')
    if (status /= 0) stop 'FAILED: pacs%read_detector_mask'
    call pacs%init_with_calfiles(obs%band, detector_mask, 1, status)
    if (status /= 0) stop 'FAILED: pacs%init_with_calfiles'

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
    if (any(mask_ref .neqv. mask(:,idetector))) stop 'FAILED: mask ref'
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
           if (status /= 0) stop 'FAILED: read_pacsobservation loop'
           if (any(signal(:,idetector) /= signal_ref(first:last))) stop 'FAILED: read signal'
           if (any(mask(:,idetector) .neqv. mask_ref(first:last))) then
               print *,first, last
               print *,shape(mask(:,idetector))
               print *,shape(mask_ref)
               print *,mask(:,idetector)
               print *,mask_ref(first:last)
               stop 'FAILED: read mask'
           end if
           deallocate (signal)
           deallocate (mask)
        end do
    end do

    stop 'OK.'

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

end program test_pacsobservation
