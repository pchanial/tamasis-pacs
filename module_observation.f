! Because of bugs in gfortran and ifort related to array of classes
! the architecture is (temporarily hopefully) a mess
! this module should define a derived type: observationslice and the
! module module_pacsobservation should extend it and define pacsobservationslice
! for now, many routines that should belong to module_pacsobservation are here...

module module_observation

    use iso_fortran_env,  only : ERROR_UNIT, OUTPUT_UNIT
    use module_fitstools, only : ft_open, ft_open_bintable, ft_read_column, ft_read_image, ft_close
    use module_math,      only : NaN, median, neq_real
    use module_precision, only : dp, p
    use module_string,    only : strinteger, strlowcase, strreal, strsection, strternary
    implicit none
    private

    public :: maskarray
    public :: observation
!!$    public :: observationslice
public :: pacsobservationslice
public :: pacsobservation    
    public :: pointing

    type maskarray
        logical :: remove_invalid = .true.
        logical :: off_scan = .true.
        logical :: off_target = .false.
        logical :: wrong_time = .true.
    end type maskarray

    type pointing
        real(dp) :: time
        real(dp) :: ra
        real(dp) :: dec
        real(dp) :: pa
        real(dp) :: chop
        logical  :: invalid
        logical  :: off_scan
        logical  :: off_target
        logical  :: wrong_time
    end type pointing

!!$    type observationslice
!!$        integer            :: first, last, nsamples
!!$        character(len=256) :: filename
!!$        character          :: channel
!!$        character(len=80)  :: observing_mode
!!$        integer            :: compression_factor
!!$        real(dp)           :: sampling_interval
!!$        type(pointing), allocatable :: p(:)
!!$
!!$    end type observationslice

!!$ type, extends(observationslice) :: pacsobservationslice
    type :: pacsobservationslice

        integer            :: first, last, nsamples, nvalids
        character(len=256) :: filename
        character          :: channel
        character(len=70)  :: observing_mode
        character(len=70)  :: unit
        integer            :: compression_factor
        real(dp)           :: sampling_interval
        type(pointing), allocatable :: p(:)

    contains
        
        procedure :: set_channel
        procedure :: set_compression_mode
        procedure :: set_filename
        procedure :: set_invalid
        procedure :: set_observing_mode
        procedure :: set_astrometry
        procedure :: set_unit
        procedure :: set_flags
        procedure :: set_sampling_interval
        procedure :: set_astrometry_oldstyle

    end type pacsobservationslice

    type :: observation
        integer           :: nslices
        integer           :: nsamples   ! number of pointings
        integer           :: nvalids    ! number of pointings for which the signal and mask should be read
        character         :: channel
        character(len=70) :: observing_mode
        character(len=70) :: unit
        type(maskarray)   :: policy
        type(pacsobservationslice), allocatable :: slice(:)
    contains
        procedure :: compute_center
        procedure :: get_position_time
        procedure :: get_position_index
        procedure :: destructor
    end type observation

    ! the following belongs to module_pacsobservation.f
    ! moved here because of gfortran bug #44065
    type, extends(observation) :: pacsobservation

    contains

        procedure :: init
        procedure :: init_sim
        procedure :: print

    end type pacsobservation


contains


    ! linear interpolation if time samples are not evenly spaced.
    ! for the first call, index should be set to zero
    subroutine get_position_time(this, islice, time, ra, dec, pa, chop, index)
        class(observation), intent(in) :: this
        integer, intent(in)            :: islice
        real*8, intent(in)             :: time
        real*8, intent(out)            :: ra, dec, pa, chop
        integer, intent(inout)         :: index
        integer                        :: i, first, last, nsamples
        real*8                         :: frac

        if (islice < 1 .or. islice > this%nslices) then
            write (ERROR_UNIT,'(a,i0,a)') "GET_POSITION: Invalid slice number '", islice, "'."
            stop
        end if

        first = this%slice(islice)%first
        last  = this%slice(islice)%last
        nsamples = last - first + 1

        if (time > 2.d0 * this%slice(islice)%p(nsamples)%time - this%slice(islice)%p(nsamples-1)%time .or.                         &
            time < 2.d0 * this%slice(islice)%p(1)%time - this%slice(islice)%p(2)%time) then
            ra   = NaN
            dec  = NaN
            pa   = NaN
            chop = NaN
            return
        end if

        if (index == 0 .or. index-1 > nsamples) then
            index = first + 1
        else if (time <= this%slice(islice)%p(index-1)%time) then
            index = first + 1
        end if

        do i = index, last
            if (time <= this%slice(islice)%p(i)%time) go to 100
        end do
        i = i - 1

    100 frac = (time - this%slice(islice)%p(i-1)%time) / (this%slice(islice)%p(i)%time - this%slice(islice)%p(i-1)%time)

        ra   = this%slice(islice)%p(i-1)%ra   * (1 - frac) + this%slice(islice)%p(i)%ra   * frac
        dec  = this%slice(islice)%p(i-1)%dec  * (1 - frac) + this%slice(islice)%p(i)%dec  * frac
        pa   = this%slice(islice)%p(i-1)%pa   * (1 - frac) + this%slice(islice)%p(i)%pa   * frac
        chop = this%slice(islice)%p(i-1)%chop * (1 - frac) + this%slice(islice)%p(i)%chop * frac

        index = i

    end subroutine get_position_time


    !-------------------------------------------------------------------------------------------------------------------------------


    ! The sampling factor is the compression factor times the fine sampling
    ! factor. itime is a fine sampling index
    ! if there is a gap in the timeline (which should not happen, and should be
    ! corrected beforehand), it will not be taken into account so that finer
    ! sampling will be interpolated using the start and end of the gap.
    subroutine get_position_index(this, islice, itime, sampling_factor, ra, dec, pa, chop)
        class(observation), intent(in) :: this
        integer, intent(in)            :: islice, itime, sampling_factor
        real*8, intent(out)            :: ra, dec, pa, chop
        integer                        :: i, itime_max
        real*8                         :: frac

        if (islice < 1 .or. islice > this%nslices) then
            write (ERROR_UNIT,'(a,i0,a)') "GET_POSITION: Invalid slice number '", islice, "'."
            stop
        end if

        itime_max = this%slice(islice)%nsamples * sampling_factor
        if (itime < 1 .or. itime > itime_max) then
            write (ERROR_UNIT,'(a,i0,a,i0,a)') "GET_POSITION: Invalid time index '", itime, "'. Valid range is [1:", itime_max,"]."
            stop
        end if

        i    = min((itime - 1) / sampling_factor, this%slice(islice)%nsamples-2)
        frac = real(itime - 1 - sampling_factor * i, p) / sampling_factor
        i    = i + 1

        ra   = this%slice(islice)%p(i)%ra   * (1 - frac) + this%slice(islice)%p(i+1)%ra   * frac
        dec  = this%slice(islice)%p(i)%dec  * (1 - frac) + this%slice(islice)%p(i+1)%dec  * frac
        pa   = this%slice(islice)%p(i)%pa   * (1 - frac) + this%slice(islice)%p(i+1)%pa   * frac
        chop = this%slice(islice)%p(i)%chop * (1 - frac) + this%slice(islice)%p(i+1)%chop * frac

    end subroutine get_position_index


    !-------------------------------------------------------------------------------------------------------------------------------


    ! compute the mean value of R.A. and dec, taking into account RA's singularity at 0
    subroutine compute_center(this, ra0, dec0)

        class(observation), intent(in) :: this
        real*8, intent(out)            :: ra0, dec0
        integer                        :: isample, islice, n180, nvalids
        real*8                         :: ra, dec
        logical                        :: zero_minus, zero_plus

        ra0  = 0.0d0
        dec0 = 0.0d0
        zero_minus = .false.
        zero_plus  = .false.
        n180 = 0
        nvalids = 0

        do islice = 1, this%nslices
            !$omp parallel do default(shared) reduction(+:n180, nvalids,ra0,dec0)           &
            !$omp reduction(.or.:zero_minus,zero_plus) private(isample, ra, dec)
            do isample=1, this%slice(islice)%nsamples
                if (this%slice(islice)%p(isample)%invalid) cycle
                ra  = this%slice(islice)%p(isample)%ra
                dec = this%slice(islice)%p(isample)%dec
                zero_minus = zero_minus .or. ra > 270.d0
                zero_plus  = zero_plus  .or. ra <= 90.d0
                if (ra >= 180.d0) n180 = n180 + 1
                ra0  = ra0  + ra
                dec0 = dec0 + dec
                nvalids = nvalids + 1
            end do
            !$omp end parallel do
        end do

        if (zero_minus .and. zero_plus) ra0 = ra0 - 360.d0 * n180

        ra0  = ra0  / nvalids
        dec0 = dec0 / nvalids
            
    end subroutine compute_center


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine print_pointing(this)
        type(pointing), intent(in) :: this(:)

        write (*,*) 'Time: ', minval(this%time), '...', maxval(this%time)
        write (*,*) 'RA  : ', minval(this%ra)  , '...', maxval(this%ra)
        write (*,*) 'Dec : ', minval(this%dec) , '...', maxval(this%dec)
        write (*,*) 'PA  : ', minval(this%pa)  , '...', maxval(this%pa)
        write (*,*) 'Chop: ', minval(this%chop), '...', maxval(this%chop)

    end subroutine print_pointing


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine destructor(this)

        class(observation), intent(inout) :: this

        integer                           :: islice

        if (.not. allocated(this%slice)) return
        do islice = 1, this%nslices
            if (allocated(this%slice(islice)%p)) deallocate (this%slice(islice)%p)
        end do
        deallocate (this%slice)

    end subroutine


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine set_filename(this, filename, first, last, status)

        class(pacsobservationslice), intent(inout) :: this
        character(len=*), intent(in)               :: filename
        integer, intent(out)                       :: first, last
        integer, intent(out)                       :: status

        integer :: pos, length, delim
        logical :: is_slice
 
        status = 0
        first = 0
        last  = 0

        length = len_trim(filename)
        if (length > len(this%filename)) then
            status = 1
            write (ERROR_UNIT,'(a,i0,a,i0,a)') 'ERROR: Input filename length is too long ', length, '. Maximum is ',               &
                  len(this%filename), '.'
            return
        end if

        this%filename = filename

        if (filename(length:length) /= ']') return

        do pos = length-1, 2, -1
            if (filename(pos:pos) == '[') exit
        end do

        if (pos == 1) then
            status = 1
            write (ERROR_UNIT,'(a)') "ERROR: Missing opening bracket in '" // trim(filename) // "'."
            return
        end if
      
        pos = pos + 1
        length = length - 1

        ! find delimiter ':'
        delim = index(filename(pos:length), ':')
        is_slice = delim > 0
        if (is_slice) then
            delim = pos + delim - 1
        else 
            delim = length + 1
        end if

        ! read first sample
        if (delim > pos) then
            read (filename(pos:delim-1), '(i20)', iostat=status) first
            if (status /= 0) then
                write (ERROR_UNIT,'(a)') "ERROR: Invalid first sample: '" // filename(pos-1:length+1) // "'."
                return
            end if
        end if

        ! read last sample
        if (is_slice) then
            if (delim < length) then
                read (filename(delim+1:length), '(i20)', iostat=status) last
                if (status /= 0) then
                    write (ERROR_UNIT,'(a)') "ERROR: Invalid last sample: '" // filename(pos-1:length+1) // "'."
                    return
                end if
            end if
        else
            last = first
        end if

        if (last /= 0 .and. last < first) then
            status = 1
            write (ERROR_UNIT,'(a,2(i0,a))')"ERROR: Last sample '", last,  "' is less than the first sample '", first, "'."
            return
        end if

        if (first < 0) then
            status = 1
            write (ERROR_UNIT,'(a,i0,a)')"ERROR: The first sample '", first, "' is less than 1."
            return
        end if

        ! range specifications could be read, remove it
        this%filename(pos-1:) = ' '

    end subroutine set_filename


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine set_invalid(this, policy, first, last, status)

        class(pacsobservationslice), intent(inout) :: this
        type(maskarray), intent(in)                :: policy
        integer, intent(in)                        :: first, last
        integer, intent(out)                       :: status

        type(pointing), allocatable :: p_(:)
        integer                     :: isample

        status = 0

        ! set invalid flag according to the mask policy
        if (policy%off_scan) then
            where (this%p%off_scan)
                this%p%invalid = .true.
            end where
        end if

        if (policy%off_target) then
            where (this%p%off_target)
                this%p%invalid = .true.
            end where
        end if

        if (policy%wrong_time) then
            where (this%p%wrong_time)
                this%p%invalid = .true.
            end where
        end if
        
        ! if a slice is specified, mask its complement otherwise, search for it
        this%first = first
        if (first /= 0) then
            this%p(1:first-1)%invalid = .true.
        else
            do isample = 1, size(this%p)
                if (.not. this%p(isample)%invalid) exit
            end do
            this%first = isample
        end if

        this%last = last
        if (last /= 0) then
            if (last > this%nsamples) then
                status = 1
                write (ERROR_UNIT,'(a,2(i0,a))') "ERROR: The last sample '", last, "' exceeds the size of the observation '",      &
                     this%nsamples, "'."
                return
            end if
            this%p(last+1:)%invalid = .true.
        else
            do isample = size(this%p), 1, -1
                if (.not. this%p(isample)%invalid) exit
            end do
            this%last = isample
        end if

        if (this%first > this%last) then
            this%first = 1
            this%last  = this%nsamples
        end if
        
        ! retain pointings bracketed by [first,last]
        this%nsamples = this%last - this%first + 1
        allocate (p_(this%nsamples))
        p_ = this%p(this%first:this%last)
        call move_alloc (p_, this%p)

        if (policy%remove_invalid) then
            this%nvalids = count(.not. this%p%invalid)
        else
            this%nvalids = this%nsamples
        end if

    end subroutine set_invalid


    !-------------------------------------------------------------------------------------------------------------------------------


    ! sets 'b', 'g' or 'r' for a blue, green or red channel observation
    subroutine set_channel(this, status)

        class(pacsobservationslice), intent(inout) :: this
        integer, intent(out)                       :: status

        integer                       :: length, unit, ivalid, nsamples
        character(len=5), allocatable :: channels(:)

        this%channel = ' '

        ! old style file format
        length = len_trim(this%filename)
        if (this%filename(length-4:length) /= '.fits') then
            status = 0
            if (strlowcase(this%filename(length-3:length)) == 'blue') then
               this%channel = 'b'
            else if (strlowcase(this%filename(length-4:length)) == 'green') then
               this%channel = 'g'
            else if (strlowcase(this%filename(length-2:length)) == 'red') then
               this%channel = 'r'
            else
               status = 1
               write (ERROR_UNIT,'(a)') 'File name does not contain the array channel identifier (blue, green, red).'
            end if
            return
        endif

        call ft_open_bintable(trim(this%filename) // '[Status]', unit, nsamples, status)
        if (status /= 0) return

        allocate(channels(nsamples))
        call ft_read_column(unit, 'BAND', 1, nsamples, channels, status)
        if (status /= 0) return

        call ft_close(unit, status)
        if (status /= 0) return

        ! get first defined BAND
        do ivalid=1, nsamples
            if (channels(ivalid) /= 'UNDEF') exit
        end do

        ! check that there is at least one defind BAND
        if (ivalid > nsamples) then
            status = 1
            write (ERROR_UNIT, '(a)') 'All observation samples have UNDEF band.'
            return
        end if

        ! check the observation only has one BAND
        if (any(channels /= channels(ivalid) .and. channels /= 'UNDEF')) then
            status = 1
            write (ERROR_UNIT,'(a)') 'Observation is BANDSWITCH.'
            return
        end if

        ! mark UNDEF BAND as invalid
        this%p%invalid = this%p%invalid .or. channels == 'UNDEF'

        select case (channels(ivalid))
        case ('BS')
            this%channel = 'b'
        case ('BL')
            this%channel = 'g'
        case ('R')
            this%channel = 'r'
        case default
            status = 1
            write (ERROR_UNIT,'(a)') 'Invalid array BAND value: ' // channels(ivalid)
        end select

    end subroutine set_channel


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine set_compression_mode(this, status)

        class(pacsobservationslice), intent(inout) :: this
        integer, intent(out)                       :: status

        integer           :: unit, ikey, length, status_close
        character(len=70) :: algorithm, keyword, comment

        this%compression_factor = 0

        ! old style file format
        length = len_trim(this%filename)
        if (this%filename(length-4:length) /= '.fits') then
            status = 0
            this%compression_factor = 1
            write (*,'(a)') 'Warning: Assuming compression factor of one for ob&
                &solete file format.'
            return
        end if

        call ft_open(trim(this%filename), unit, status)
        if (status /= 0) return

        ikey = 1
        do
            call ftgkys(unit, 'key.META_' // strinteger(ikey), keyword, comment, status)
            if (status /= 0) then
                write (ERROR_UNIT,'(a)') "ERROR: FITS keyword 'algorithm' is not found."
                go to 999
            end if
            if (keyword == 'algorithm') exit
            ikey = ikey + 1
        end do

        call ftgkys(unit, 'META_' // strinteger(ikey), algorithm, comment, status)

        if (algorithm(1:19) == 'Floating Average  :') then
            read (algorithm(20:),'(i3)', iostat=status) this%compression_factor
            if (status /= 0) then
                write (ERROR_UNIT, '(a)') "ERROR: The compression algorithm '" // trim(algorithm) // "' is not understood."
                go to 999
            end if
        else if (algorithm == 'None') then
            this%compression_factor = 1
        end if
        
    999 call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine set_compression_mode


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine set_observing_mode(this, status)
        class(pacsobservationslice), intent(inout) :: this
        integer, intent(out)                       :: status
        integer           :: unit, ikey, length, status_close
        character(len=70) :: compression, keyword, comment

        this%observing_mode = 'Unknown'

        ! old style file format
        length = len_trim(this%filename)
        if (this%filename(length-4:length) /= '.fits') then
            status = 0
            this%observing_mode = 'Transparent'
            write (*,'(a)') 'Warning: Transparent mode assumed for obsolete file format.'
            return
        end if

        call ft_open(trim(this%filename), unit, status)
        if (status /= 0) return

        ikey = 1
        do
            call ftgkys(unit, 'key.META_'//strinteger(ikey), keyword, comment, status)
            if (status /= 0) then
                status = 0
                write (ERROR_UNIT,'(a)') "Warning: FITS keyword 'compMode' is n&
                    &ot found. Assuming normal mode."
                go to 999
            end if
            if (keyword == 'compMode') exit
            ikey = ikey + 1
        end do

        call ftgkys(unit, 'META_'//strinteger(ikey), compression, comment, status)
        if (compression == 'Photometry Lossless Compression Mode') then
            this%observing_mode = 'Transparent'
        else if (compression == 'Photometry Default Mode') then
            this%observing_mode = 'Prime'
        else if (compression == 'Photometry Double Compression Mode') then
            this%observing_mode = 'Parallel'
        else
            write (OUTPUT_UNIT,'(a)') "Warning: Unknown compression mode: '" // trim(compression) // "'."
        end if

    999 call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine set_observing_mode


    !-------------------------------------------------------------------------------------------------------------------------------


    ! if the first sample to consider is not specified, then:
    !     - if the bbid is set, we search starting from the beginning the first valid sample
    !     - otherwise, we search starting from the end for the first sample that satisfies abs(chopfpuangle) > 0.01.
    ! if invalid samples have been found, we then discard the first 10, that might be affected by relaxation.
    subroutine set_flags(this, status)

        class(pacsobservationslice), intent(inout) :: this
        integer, intent(out)                       :: status

        logical(1), allocatable :: off_target(:)
        logical(1), allocatable :: off_scan(:)
        integer                 :: length
        integer                 :: unit
        integer                 :: isample, itarget
        integer*8, allocatable  :: bbid(:)
        real*8, allocatable     :: chop(:)

        length = len_trim(this%filename)
        if (this%filename(length-4:length) /= '.fits') then
            write (OUTPUT_UNIT,'(a)') 'Warning: Obsolete file format.'
            call ft_read_image(trim(this%filename) // '_ChopFpuAngle.fits', chop, status)
            if (status /= 0) return
            this%nsamples = size(chop)
            allocate (bbid(this%nsamples))
            bbid = 0
        else

            call ft_open_bintable(trim(this%filename) // '[Status]', unit, this%nsamples, status)
            if (status /= 0) return

            allocate (bbid(this%nsamples))
            call ft_read_column(unit, 'BBID', 1, this%nsamples, bbid, status)
            if (status /= 0) return
            bbid = ishft(bbid, -16)

            allocate (chop(this%nsamples))
            call ft_read_column(unit, 'CHOPFPUANGLE', 1, this%nsamples, chop, status)
            if (status /= 0) return

            call ft_close(unit, status)
            if (status /= 0) return

        end if

        if (this%nsamples == 0) then
            status = 1
            write (ERROR_UNIT, '(a)') 'ERROR: Status extension is empty.'
            return
        end if

        if (maxval(bbid) == 0) then
            where (abs(chop) > 0.01_p)
                bbid = z'4004'
            elsewhere
                bbid = z'cd2'
            end where
            do isample = 1, this%nsamples
                if (bbid(isample) == z'4004') exit
            end do
            if (isample <= this%nsamples) then
                bbid(1:isample-1) = z'4000'
            end if
        end if

        allocate (off_target(this%nsamples))
        allocate (off_scan(this%nsamples))

        off_target = bbid /= z'cd2'
        off_scan = bbid /= z'cd2' .and. bbid /= z'4000'

        ! search for the first on-target sample
        do itarget = 1, this%nsamples
            if (.not. off_target(itarget)) exit
        end do
        if (itarget <= this%nsamples) then
            ! backward search first off-scan sample
            do isample = itarget-1, 1, -1
                if (off_scan(isample)) exit
            end do
            ! mask as off-scan all samples before the last 0x4000 block
            off_scan(1:isample) = .true.
        end if
        
        if (allocated(this%p)) deallocate (this%p)
        allocate (this%p(this%nsamples))
        this%p%invalid = .false.
        this%p%off_scan = off_scan
        this%p%off_target = off_target
        this%p%wrong_time = .false.

    end subroutine set_flags


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine set_astrometry(this, status)

        class(pacsobservationslice), intent(inout) :: this
        integer, intent(out)                       :: status

        integer*8, allocatable :: timeus(:)
        integer                :: nsamples
        integer                :: length, unit
        
        length = len_trim(this%filename)
        if (this%filename(length-4:length) /= '.fits') then
            call this%set_astrometry_oldstyle(status)
            return
        end if
        
        call ft_open_bintable(trim(this%filename) // '[Status]', unit, nsamples, status)
        if (status /= 0) return
        
        if (nsamples <= 1) then
            write (ERROR_UNIT,'(a)') 'Input time has less than 2 samples.'
            call ft_close(unit, status)
            status = 1
            return
        end if

        allocate(timeus(nsamples))
        
        call ft_read_column(unit, 'FINETIME', 1, nsamples, timeus, status)
        if (status /= 0) return
        this%p%time = timeus * 1.d-6
        
        call ft_read_column(unit, 'RaArray', 1, nsamples, this%p%ra, status)
        if (status /= 0) return
        
        call ft_read_column(unit, 'DecArray', 1, nsamples, this%p%dec, status)
        if (status /= 0) return
        
        call ft_read_column(unit, 'PaArray', 1, nsamples, this%p%pa, status)
        if (status /= 0) return
                
        call ft_read_column(unit, 'CHOPFPUANGLE', 1, nsamples, this%p%chop, status)
        if (status /= 0) return

        call ft_close(unit, status)

    end subroutine set_astrometry


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine set_astrometry_oldstyle(this, status)

        class(pacsobservationslice), intent(inout) :: this
        integer, intent(out)                       :: status

        real(dp), allocatable  :: buffer(:)
        integer*8, allocatable :: timeus(:)
        
        call ft_read_image(trim(this%filename)//'_Time.fits', timeus, status)
        if (status /= 0) return
        allocate(buffer(size(timeus)))
        buffer = timeus * 1.0d-6
        this%p%time = buffer

        call ft_read_image(trim(this%filename) // '_RaArray.fits', buffer, status)
        if (status /= 0) return
        this%p%ra = buffer

        call ft_read_image(trim(this%filename) // '_DecArray.fits', buffer, status)
        if (status /= 0) return
        this%p%dec = buffer

        call ft_read_image(trim(this%filename) // '_PaArray.fits', buffer, status)
        if (status /= 0) return
        this%p%pa = buffer

        call ft_read_image(trim(this%filename) // '_ChopFpuAngle.fits', buffer, status)
        if (status /= 0) return
        this%p%chop = buffer

    end subroutine set_astrometry_oldstyle


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine set_unit(this, status)

        class(pacsobservationslice), intent(inout) :: this
        integer, intent(out)                       :: status

        character(len=70) :: comment
        integer           :: length, unit

        this%unit = ''
        length = len_trim(this%filename)
        if (this%filename(length-4:length) /= '.fits') then
            return
        end if

        call ft_open(trim(this%filename) // '[Signal]', unit, status)
        if (status /= 0) return

        call ftgkys(unit, 'QTTY____', this%unit, comment, status)
        status = 0
        this%unit = trim(this%unit)

        call ft_close(unit, status)

    end subroutine set_unit

    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine set_sampling_interval(this, status, verbose)

        class(pacsobservationslice), intent(inout) :: this
        integer, intent(out)                       :: status
        logical, intent(in)                        :: verbose

        integer                :: isample, njumps
        integer                :: compression_factor
        real(dp), allocatable  :: delta(:)
        real(dp)               :: delta_max

        ! special case if there is only one sample
        if (this%nsamples == 1) then
            this%sampling_interval = this%compression_factor * 0.024996_dp
            return
        end if

        status = 1

        ! check that the input time is monotonous (and increasing)
        allocate(delta(this%nsamples-1))
        delta = this%p(2:this%nsamples)%time - this%p(1:this%nsamples-1)%time
        this%sampling_interval = median(delta)
        njumps = 0
        do isample = 2, this%nsamples
            if (delta(isample-1) <= 0) then
                if (isample < this%nsamples) then
                    if (delta(isample) <= 0) then
                        write (ERROR_UNIT,'(a)') "ERROR: The pointing time is not strictly increasing."
                        return
                    end if
                end if
                njumps = njumps + 1
                this%p(isample)%wrong_time = .true.
                this%p(isample)%time = this%p(isample-1)%time + this%sampling_interval
                this%p(isample)%ra   = this%p(isample-1)%ra  
                this%p(isample)%dec  = this%p(isample-1)%dec 
                this%p(isample)%pa   = this%p(isample-1)%pa  
                this%p(isample)%chop = this%p(isample-1)%chop
            end if
        end do
        
        if (njumps > 0) then
            if (verbose) then
                write (OUTPUT_UNIT,'(a,i0,a)') 'Warning: The pointing fine time has ', njumps, ' negative jump(s). The affected fra&
                     &mes have been masked.'
            end if
            delta = this%p(2:this%nsamples)%time - this%p(1:this%nsamples-1)%time
        end if
        
        ! check if there are gaps
        delta_max = maxval(abs(delta))
        if (verbose .and. any(neq_real(delta, this%sampling_interval, 3))) then
            write (*,'(a,$)') "Warning: In observation '" // trim(this%filename) //"', the pointing time is not evenly spaced."
            if (delta_max > 1.5_p * this%sampling_interval) then
                write (*,'(a)') ' Largest gap is ' // strreal(delta_max*1000._p,1) // 'ms.'
            else
                write (*,*)
            end if
        end if
        
        ! check the compression factor from the data themselves
        compression_factor = nint(this%sampling_interval / 0.024996_dp)
        if (neq_real(compression_factor * 0.024996_dp, this%sampling_interval, 2)) then
            write (*,'(a)') 'Error: The sampling time is not an integer number of PACS sampling time (40Hz).'
            return
        end if
        
        if (compression_factor /= this%compression_factor) then
            write (*,'(a,2(i0,a))') "Error: The compression_factor determined from the observation header '",                      &
                  this%compression_factor, "' is different from the one inferred from the time in the Status table '",             &
                  compression_factor, "'."
            return
        end if

        status = 0

    end subroutine set_sampling_interval


    ! the following belongs to module_pacsobservation.f
    ! moved here because of gfortran bug #44065

    subroutine init(this, filename, policy, status, verbose)

        class(pacsobservation), intent(inout) :: this
        character(len=*), intent(in)          :: filename(:)
        type(maskarray), intent(in)           :: policy
        integer, intent(out)                  :: status
        logical, intent(in), optional         :: verbose
        integer                               :: first, last
        integer                               :: islice
        logical                               :: verbose_

        ! parameter checking
        if (size(filename) == 0) then 
            status = 1
            write (ERROR_UNIT,'(a)') 'INIT_PACSOBSERVATION: Array has zero size.'
            return
        end if

        if (present(verbose)) then
            verbose_ = verbose
        else
            verbose_ = .false.
        end if

        ! set number of observations
        this%nslices = size(filename)
        if (allocated(this%slice)) deallocate(this%slice)
        allocate(this%slice(this%nslices))

        ! set the mask policy '.true.' means rejected or masked out
        this%policy = policy        
        
        do islice = 1, this%nslices

            call this%slice(islice)%set_filename(filename(islice), first, last, status)
            if (status /= 0) return

            call this%slice(islice)%set_flags(status)
            if (status /= 0) return

            call this%slice(islice)%set_channel(status)
            if (status /= 0) return
                
            call this%slice(islice)%set_compression_mode(status)
            if (status /= 0) return
        
            call this%slice(islice)%set_observing_mode(status)
            if (status /= 0) return

            call this%slice(islice)%set_unit(status)
            if (status /= 0) return

            call this%slice(islice)%set_astrometry(status)
            if (status /= 0) return
            
            call this%slice(islice)%set_sampling_interval(status, verbose_)
            if (status /= 0) return

            call this%slice(islice)%set_invalid(this%policy, first, last, status)
            if (status /= 0) return

        end do

        ! make sure the channel is the same for all observations
        this%channel = this%slice(1)%channel
        if (any(this%slice%channel /= this%channel)) then
            status = 1
            write (ERROR_UNIT,'(a)') 'ERROR: Observations do not have the same channel: ', this%slice%channel, '.'
            return
        end if

        ! make sure the observing mode is the same for all observations
        ! that could be relaxed, but we would need to properly handle the bad detector mask for transparent observations
        this%observing_mode = this%slice(1)%observing_mode
        if (any(this%slice%observing_mode /= this%observing_mode)) then
            status = 1
            write (ERROR_UNIT,'(a)') 'ERROR: Observations do not have the same observing mode.'
            return
        end if

        ! make sure the unit is the same for all observations
        this%unit = ''
        do islice = 1, this%nslices
            if (this%slice(islice)%unit /= '') then
                this%unit = this%slice(islice)%unit
                exit
            end if
        end do

        if (this%unit == '') then
            if (verbose_) then
                write (OUTPUT_UNIT,'(a)') "Warning: Observation has no unit defined. Assuming 'Jy'."
            end if
            this%unit = 'Jy'
            this%slice%unit = 'Jy'
        else
            do islice = 1, this%nslices
                if (this%slice(islice)%unit == '' .and. status == 0) then
                    if (verbose_) then
                        write (OUTPUT_UNIT,'(a,i0,a)') 'Warning: Observation ', islice, " has no units. Assuming '" //             &
                              trim(this%unit) // "'."
                    end if
                    this%slice(islice)%unit = this%unit
                else if (this%slice(islice)%unit /= this%unit) then
                    write (ERROR_UNIT,'(a,i0,a)') 'Error: Observation ', islice, " has a unit '" // this%slice(islice)%unit // "' i&
                          &ncompatible with '" // this%unit // "'."
                    status = 1
                end if
            end do
            if (status /= 0) return
        end if

        this%nsamples = sum(this%slice%nsamples)
        this%nvalids  = sum(this%slice%nvalids)

    end subroutine init


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine init_sim(this, time, ra, dec, pa, chop, status, verbose)
        class(pacsobservation), intent(inout) :: this
        real*8, intent(in)                    :: time(:),ra(:),dec(:),pa(:),chop(:)
        integer, intent(out)                  :: status
        logical, intent(in), optional         :: verbose
        integer                               :: nsamples
        logical                               :: verbose_

        if (present(verbose)) then
            verbose_ = verbose
        else
            verbose_ = .false.
        end if

        status = 1

        ! check conformity of time, ra, dec, pa
        nsamples = size(time)

        if (size(ra) /= nsamples) then
            write (ERROR_UNIT,'(a)') "Input R.A. has an invalid number of samples."
            return
        endif
        if (size(dec) /= nsamples) then
            write (ERROR_UNIT,'(a)') "Input declination has an invalid number of samples."
            return
        endif
        if (size(pa) /= nsamples) then
            write (ERROR_UNIT,'(a)') "Input P.A. has an invalid number of samples."
            return
        endif
        if (size(chop) /= nsamples) then
            write (ERROR_UNIT,'(a)') "Input chop angle has an invalid number of samples."
            return
        endif

        this%nslices      = 1
        this%nsamples = nsamples
        this%channel      = ' '
        this%observing_mode = 'Unknown'
        this%policy = maskarray(wrong_time=.true., off_target=.true.)

        allocate (this%slice(1))
        this%slice(1)%first    = 1
        this%slice(1)%last     = nsamples
        this%slice(1)%nsamples = nsamples
        this%slice(1)%filename = 'simulation'
        this%slice(1)%channel  = ' '
        this%slice(1)%observing_mode = 'Unknown'
        this%slice(1)%compression_factor = nint((time(2)-time(1)) / 0.024996_dp)
        this%slice(1)%sampling_interval = time(2)-time(1)
        allocate (this%slice(1)%p(nsamples))

        this%slice(1)%p%time       = time
        this%slice(1)%p%ra         = ra
        this%slice(1)%p%dec        = dec
        this%slice(1)%p%pa         = pa
        this%slice(1)%p%chop       = chop
        this%slice(1)%p%invalid    = .false.
        this%slice(1)%p%wrong_time = .false.
        this%slice(1)%p%off_target = .false.

        call this%slice(1)%set_sampling_interval(status, verbose_)
        if (status /= 0) return

        call this%slice(1)%set_invalid(this%policy, 1, nsamples, status)

    end subroutine init_sim


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine print(this)

        class(pacsobservation), intent(in) :: this

        integer                :: islice, nrejected, nsamples
        logical*1, allocatable :: invalid(:)
        character(len=8)       :: mode

        write (*,*)

        do islice = 1, this%nslices

            nsamples = this%slice(islice)%nsamples

            ! observation number & file name
            write (OUTPUT_UNIT,'(a)') 'Info: Observation' // strternary(this%nslices>1, ' ' // strinteger(islice), '') // ': ' //  &
                  trim(this%slice(islice)%filename)
            write (OUTPUT_UNIT,'(a)') '      Section: [' // strsection(this%slice(islice)%first,this%slice(islice)%last) // ']'
            
            ! channel
            write (OUTPUT_UNIT,'(a,$)') "      Channel: "
            select case (this%slice(islice)%channel)
                case ('b')
                    write (OUTPUT_UNIT,'(a)') 'Blue'
                case ('g')
                    write (OUTPUT_UNIT,'(a)') 'Green'
                case ('r')
                    write (OUTPUT_UNIT,'(a)') 'Red'
                case default
                    write (OUTPUT_UNIT,'(a)') 'Unknown'
            end select

            ! observing mode
            write (OUTPUT_UNIT,'(a)') '      Observing mode: ' // trim(this%observing_mode)

            ! compression factor
            write (OUTPUT_UNIT,'(a,i0)') '      Compression factor: ', this%slice(islice)%compression_factor

            ! unit
            write (OUTPUT_UNIT,'(a)') '      Unit: ' // trim(this%slice(islice)%unit)

            ! recompute invalid samples to 
            allocate (invalid(nsamples))
            invalid = .false.

            ! print maskarray information

            mode = strternary(this%policy%remove_invalid, 'Removed ', 'Masked ')
            if (this%policy%off_scan) then
                invalid = invalid .or. this%slice(islice)%p%off_scan
                nrejected = count(this%slice(islice)%p%off_scan)
                if (nrejected > 0) then
                    write (OUTPUT_UNIT,'(3a,2(i0,a))') "      Mask 'off scan': ", trim(mode), ' ', nrejected, ' / ', nsamples, '.'
                end if
            end if

            if (this%policy%off_target) then
                invalid = invalid .or. this%slice(islice)%p%off_target
                nrejected = count(this%slice(islice)%p%off_target)
                if (nrejected > 0) then
                    write (OUTPUT_UNIT,'(3a,2(i0,a))') "      Mask 'off target': ", trim(mode), ' ', nrejected, ' / ', nsamples, '.'
                end if
            end if

            if (this%policy%wrong_time) then
                invalid = invalid .or. this%slice(islice)%p%wrong_time
                nrejected = count(this%slice(islice)%p%wrong_time)
                if (nrejected > 0) then
                    write (OUTPUT_UNIT,'(3a,2(i0,a))') "      Mask 'wrong time': ", trim(mode), ' ', nrejected, ' / ', nsamples, '.'
                end if
            end if

            nrejected = count(invalid) - count(this%slice(islice)%p%invalid)
            if (nrejected > 0) then
                write (OUTPUT_UNIT,'(3a,2(i0,a))') "      Other: ", trim(mode), ' ', nrejected, ' / ', nsamples, '.'
            end if

            deallocate (invalid)

            write (*,*)

        end do

    end subroutine print

end module module_observation
