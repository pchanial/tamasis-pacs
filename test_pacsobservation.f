program test_pacsobservation

    use module_pacsinstrument,  only : pacsinstrument
    use module_pacsobservation, only : pacsobservation, maskarray
    use module_string,          only : strsection
    use module_tamasis,         only : init_tamasis, get_tamasis_path
    implicit none

    class(pacsobservation), allocatable :: obs
    class(pacsinstrument), allocatable  :: pacs
    type(maskarray)                     :: maskarray_policy
    character(len=*), parameter         :: filename = 'tests/frames_blue.fits'
    character(len=50), allocatable      :: afilename(:)
    integer                :: status, idetector
    real*8, allocatable    :: signal(:,:)
    logical*1, allocatable :: mask(:,:)
    integer*8              :: first, last
    logical*1              :: mask_ref(360), mask_frame(360)

    ! initialise tamasis
    call init_tamasis

    ! valid calls
    allocate(obs)
    call obs%init([filename], maskarray_policy, status)
    if (status /= 0 .or. obs%slice(1)%first /= 1 .or. obs%slice(1)%last /= 360) stop 'FAILED: init1'

    call obs%init([filename // '[23:23]'], maskarray_policy, status)
    if (status == 0 .or. obs%slice(1)%first /= 23 .or. obs%slice(1)%last /= 23) stop 'FAILED: init2'

    call obs%init([filename // '[200:]'], maskarray_policy, status)
    if (status /= 0 .or. obs%slice(1)%first /= 200 .or. obs%slice(1)%last /= 360) stop 'FAILED: init3'

    call obs%init([filename // '[:130]'], maskarray_policy, status)
    if (status /= 0 .or. obs%slice(1)%first /= 1 .or. obs%slice(1)%last /= 130) stop 'FAILED: init4'

    call obs%init([filename // '[:]'], maskarray_policy, status)
    if (status /= 0 .or. obs%slice(1)%first /= 1 .or. obs%slice(1)%last /= 360) stop 'FAILED: init5'

    call obs%init([filename // '[]'], maskarray_policy, status)
    if (status /= 0 .or. obs%slice(1)%first /= 1 .or. obs%slice(1)%last /= 360) stop 'FAILED: init6'
   
    ! invalid calls
    call obs%init([filename // '[32:23]'], maskarray_policy, status)
    if (status == 0) stop 'FAILED: badinit1'

    call obs%init([filename // '[32:361]'], maskarray_policy, status)
    if (status == 0) stop 'FAILED: badinit2'

    call obs%init([filename // '[ljdf]'], maskarray_policy, status)
    if (status == 0) stop 'FAILED: badinit4'

    call obs%init([filename // '[:l+]'], maskarray_policy, status)
    if (status == 0) stop 'FAILED: badinit5'

    ! array
    allocate(afilename(2))
    afilename(1) = trim(filename) // '[3:100]'
    afilename(2) = trim(filename) // '[201:300]'
    call obs%init(afilename, maskarray_policy, status)
    if (status /= 0 .or. sum(obs%slice%nsamples) /= 198 .or.                     &
        obs%slice(1)%first /= 3   .or. obs%slice(1)%last /= 100 .or.             &
        obs%slice(2)%first /= 201 .or. obs%slice(2)%last /= 300) stop 'FAILED: ainit1'
    deallocate(afilename)
 
    allocate(afilename(0))
    call obs%init(afilename, maskarray_policy, status)
    if (status == 0) stop 'FAILED: ainit2'
    deallocate(afilename)

    
    ! test reading observation
    call obs%init([filename], maskarray_policy, status)
    if (status /= 0) stop 'FAILED: init_pacsobservation loop'

    allocate(pacs)
    call pacs%init(obs%channel, obs%observing_mode == 'Transparent', 1, .false., .false., status)
    if (status /= 0 .or. pacs%channel /= 'b' .or. pacs%transparent_mode) stop 'FAILED: pacs%init'

    ! get a detector that has been hit by a glitch
    do idetector = 1, pacs%ndetectors
        if (pacs%pq(1,idetector)+1 == 27.and.pacs%pq(2,idetector)+1 == 64) exit
    end do

    mask_ref = .false.
    mask_ref([5,6,7,151]) = .true. ! 27 64
    mask_frame = .false.
    
    do first = 1, 360, 7
        do last = first, 360, 11
           obs%slice(1)%first = first
           obs%slice(1)%last = last
           obs%slice(1)%nsamples = last - first + 1
           obs%nsamples_tot = last - first + 1
           if (allocated(obs%slice(1)%p)) deallocate (obs%slice(1)%p)
           allocate (obs%slice(1)%p(last-first+1))
           obs%slice(1)%p%invalid = .false.
           allocate(signal(obs%slice(1)%nsamples,pacs%ndetectors))
           allocate(mask  (obs%slice(1)%nsamples,pacs%ndetectors))
           call pacs%read(obs, signal, mask, status)
           if (status /= 0) stop 'FAILED: read_pacsobservation loop'
           if (any(mask(1:last-first+1,idetector) .neqv. mask_ref(first:last))) then
               stop 'FAILED: read mask'
           end if
           deallocate(signal)
           deallocate(mask  )
        end do
    end do

    stop 'OK.'
   
end program test_pacsobservation
