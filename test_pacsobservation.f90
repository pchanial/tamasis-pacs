program test_pacsobservation

    use module_pacsinstrument, only : pacsinstrument
    use module_pacsobservation
    use string, only : strsection
    implicit none

    type(pacsobservation)              :: obs
    type(pacsobservation), allocatable :: aobs(:)
    type(pacsinstrument)        :: pacs
    character(len=*), parameter :: filename = 'tests/frames_blue.fits'
    character(len=50), allocatable :: afilename(:)
    integer                     :: status, idetector
    real*8, allocatable    :: signal(:,:)
    logical*1, allocatable :: mask(:,:)
    integer*8              :: first, last
    logical*1              :: mask_ref(360)

    ! valid calls
    call init_pacsobservation(obs, filename, status)
    if (status /= 0 .or. obs%first /= 1 .or. obs%last /= 360)                  &
        stop 'FAILED: init1'

    call init_pacsobservation(obs, filename // '[23:23]', status)
    if (status /= 0 .or. obs%first /= 23 .or. obs%last /= 23)                  &
        stop 'FAILED: init2'

    call init_pacsobservation(obs, filename // '[200:]', status)
    if (status /= 0 .or. obs%first /= 200 .or. obs%last /= 360)                &
        stop 'FAILED: init3'

    call init_pacsobservation(obs, filename // '[:130]', status)
    if (status /= 0 .or. obs%first /= 1 .or. obs%last /= 130)                  &
        stop 'FAILED: init4'

    call init_pacsobservation(obs, filename // '[:]', status)
    if (status /= 0 .or. obs%first /= 1 .or. obs%last /= 360)                  &
        stop 'FAILED: init5'

    call init_pacsobservation(obs, filename // '[]', status)
    if (status /= 0 .or. obs%first /= 1 .or. obs%last /= 360)                  &
        stop 'FAILED: init6'
   
    ! invalid calls
    call init_pacsobservation(obs, filename // '[32:23]', status)
    if (status == 0) stop 'FAILED: badinit1'

    call init_pacsobservation(obs, filename // '[32:361]', status)
    if (status == 0) stop 'FAILED: badinit2'

    call init_pacsobservation(obs, filename // '[ljdf]', status)
    if (status == 0) stop 'FAILED: badinit4'

    call init_pacsobservation(obs, filename // '[:l+]', status)
    if (status == 0) stop 'FAILED: badinit5'

    ! array
    allocate(aobs(2))
    allocate(afilename(2))
    afilename(1) = trim(filename) // '[3:100]'
    afilename(2) = trim(filename) // '[201:300]'
    call init_pacsobservation(aobs, afilename, status)
    if (status /= 0 .or. sum(aobs%nsamples) /= 198 .or.                        &
        aobs(1)%first /= 3   .or. aobs(1)%last /= 100 .or.                     &
        aobs(2)%first /= 201 .or. aobs(2)%last /= 300) stop 'FAILED: ainit1'
    deallocate(aobs)
    deallocate(afilename)
 
    allocate(aobs(0))
    allocate(afilename(0))
    call init_pacsobservation(aobs, afilename, status)
    if (status == 0) stop 'FAILED: ainit2'
    deallocate(aobs)
    deallocate(afilename)

    
    ! test reading observation
    call init_pacsobservation(obs, filename, status)

    call pacs%init_scalar(obs, 1, .false., status)
    if (status /= 0 .or. pacs%channel /= 'b' .or. pacs%transparent_mode)       &
        stop 'FAILED: pacs%init'

    call init_pacsobservation(obs, filename, status)
    if (status /= 0) stop 'FAILED: init_pacsobservation loop'

    ! get a detector that has been hit by a glitch
    do idetector = 1, pacs%ndetectors
        if (pacs%pq(1,idetector)+1 == 27.and.pacs%pq(2,idetector)+1 == 64) exit
    end do

    mask_ref = .false.
    mask_ref([5,6,7,151]) = .true. ! 27 64
    
    do first = 1, 360,7
        do last = first, 360, 11
           obs%first = first
           obs%last = last
           obs%nsamples = last - first + 1
           allocate(signal(obs%nsamples,pacs%ndetectors))
           allocate(mask  (obs%nsamples,pacs%ndetectors))
           call read_pacsobservation(obs, pacs%pq, signal, mask, status)
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
