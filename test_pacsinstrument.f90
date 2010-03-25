program test_pacsinstrument

    use module_pacsinstrument
    implicit none

    type(pacsinstrument)        :: pacs
    character(len=*), parameter :: filename = 'tests/frames_blue.fits'
    integer                :: status, idetector
    real*8, allocatable    :: signal(:,:)
    logical*1, allocatable :: mask(:,:)
    integer*8              :: first, last
    logical*1              :: mask_ref(360)

    first = 1
    last  = 360
    call pacs%init_filename(filename, .false., status)
    if (status /= 0 .or. pacs%channel /= 'b' .or. pacs%transparent_mode)       &
        stop 'FAILED: pacs%init'

    allocate(signal(last-first+1,pacs%ndetectors))
    allocate(mask  (last-first+1,pacs%ndetectors))

    call pacs%read_tod_file(filename, 0_8,   10_8,  signal, mask, status)
    if (status == 0) stop 'FAILED: read_tod_file 1'

    call pacs%read_tod_file(filename, 34_8,  12_8,  signal, mask, status)
    if (status == 0) stop 'FAILED: read_tod_file 2'

    call pacs%read_tod_file(filename, 129_8, 361_8, signal, mask, status)
    if (status == 0) stop 'FAILED: read_tod_file 3'

    ! get a detector that has been hit by a glitch
    do idetector = 1, pacs%ndetectors
        if (pacs%pq(1,idetector)+1 == 27.and.pacs%pq(2,idetector)+1 == 64) exit
    end do

    mask_ref = .false.
    mask_ref([5,6,7,151]) = .true. ! 27 64
    
    do first = 1, 360,7
        do last = first, 360, 11
           call pacs%read_tod_file(filename, first, last, signal, mask, status)
           if (status /= 0) stop 'FAILED: read_tod_file'
           if (any(mask(1:last-first+1,idetector) .neqv. mask_ref(first:last))) then
               stop 'FAILED: read mask'
           end if
        end do
    end do

    stop 'OK.'
    
end program test_pacsinstrument
