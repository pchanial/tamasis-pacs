program test_pacsobservation

    use module_pacsinstrument,  only : PacsInstrument
    use module_pacsobservation, only : MaskPolicy, PacsObservation
    use module_observation,     only : Pointing
    use module_string,          only : strsection
    use module_tamasis,         only : init_tamasis, get_tamasis_path, POLICY_REMOVE
    implicit none

    class(PacsObservation), allocatable :: obs
    class(PacsInstrument), allocatable  :: pacs
    type(Pointing), allocatable         :: p_(:)
    type(MaskPolicy)                    :: policy
    character(len=*), parameter         :: filename = 'pacs/test/data/frames_blue.fits'
    character(len=255), allocatable     :: afilename(:)
    integer                :: status, idetector
    real*8, allocatable    :: signal(:,:)
    logical*1, allocatable :: mask(:,:)
    integer                :: first, last
    real*8                 :: signal_ref(360)
    logical*1              :: mask_ref(360)

    ! initialise tamasis
    call init_tamasis

    allocate (afilename(1))
    ! valid calls
    allocate(obs)
    afilename(1) = get_tamasis_path() // filename

    call obs%init(afilename, policy, status)
    if (status /= 0 .or. obs%slice(1)%first /= 1 .or. obs%slice(1)%last /= 360) stop 'FAILED: init1'

    afilename(1) = get_tamasis_path() // filename // '[23:23]'
    call obs%init(afilename, policy, status)
    if (status /= 0 .or. obs%slice(1)%first /= 23 .or. obs%slice(1)%last /= 23) stop 'FAILED: init2'

    afilename(1) = get_tamasis_path() // filename // '[200:]'
    call obs%init(afilename, policy, status)
    if (status /= 0 .or. obs%slice(1)%first /= 200 .or. obs%slice(1)%last /= 360) stop 'FAILED: init3'

    afilename(1) = get_tamasis_path() // filename // '[:130]'
    call obs%init(afilename, policy, status)
    if (status /= 0 .or. obs%slice(1)%first /= 1 .or. obs%slice(1)%last /= 130) stop 'FAILED: init4'

    afilename(1) = get_tamasis_path() // filename // '[:]'
    call obs%init(afilename, policy, status)
    if (status /= 0 .or. obs%slice(1)%first /= 1 .or. obs%slice(1)%last /= 360) stop 'FAILED: init5'

    afilename(1) = get_tamasis_path() // filename // '[]'
    call obs%init(afilename, policy, status)
    if (status /= 0 .or. obs%slice(1)%first /= 1 .or. obs%slice(1)%last /= 360) stop 'FAILED: init6'
   
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
    if (status /= 0 .or. sum(obs%slice%nsamples) /= 198 .or.                     &
        obs%slice(1)%first /= 3   .or. obs%slice(1)%last /= 100 .or.             &
        obs%slice(2)%first /= 201 .or. obs%slice(2)%last /= 300) stop 'FAILED: ainit1'
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
    deallocate (afilename)

    allocate (pacs)
    call pacs%init(obs%channel, obs%observing_mode == 'Transparent', 1, status=status)
    if (status /= 0 .or. pacs%channel /= 'b' .or. pacs%transparent_mode) stop 'FAILED: pacs%init'

    allocate(signal(obs%slice(1)%nsamples,pacs%ndetectors))
    allocate(mask  (obs%slice(1)%nsamples,pacs%ndetectors))
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
    
    call move_alloc (obs%slice(1)%p, p_)
    
    do first = 1, 360, 7
        do last = first, 360, 11
           allocate (obs%slice(1)%p(last-first+1))
           obs%slice(1)%p = p_(first:last)
           obs%slice(1)%first = first
           obs%slice(1)%last = last
           obs%slice(1)%nsamples = last - first + 1
           obs%slice(1)%nvalids = count(.not. obs%slice(1)%p%invalid)
           obs%nsamples = obs%slice(1)%nsamples
           obs%nvalids = obs%slice(1)%nvalids
           if (allocated(obs%slice(1)%p)) deallocate (obs%slice(1)%p)
           allocate (obs%slice(1)%p(last-first+1))
           obs%slice(1)%p%invalid = .false.
           allocate(signal(obs%slice(1)%nsamples,pacs%ndetectors))
           allocate(mask  (obs%slice(1)%nsamples,pacs%ndetectors))
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
           deallocate (obs%slice(1)%p)
        end do
    end do
    

    stop 'OK.'
   
end program test_pacsobservation
