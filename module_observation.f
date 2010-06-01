! Because of bugs in gfortran and ifort related to array of classes
! the architecture is (temporarily hopefully) a mess
! this module should define a derived type: observationslice and the
! module module_pacsobservation should extend it and define pacsobservationslice
! for now, many routines that should belong to module_pacsobservation are here...

module module_observation

    use iso_fortran_env,  only : ERROR_UNIT, OUTPUT_UNIT
    use module_fitstools, only : ft_open, ft_open_bintable, ft_read_column, ft_read_extension, ft_close
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
        logical :: wrong_time = .true.
        logical :: off_target = .true.
    end type maskarray

    type pointing
        real(dp) :: time
        real(dp) :: ra
        real(dp) :: dec
        real(dp) :: pa
        real(dp) :: chop
        logical  :: invalid
        logical  :: wrong_time
        logical  :: off_target
    end type pointing

!!$    type observationslice
!!$        integer*8          :: first, last, nsamples
!!$        character(len=256) :: filename
!!$        character          :: channel
!!$        character(len=80)  :: observing_mode
!!$        integer            :: compression_factor
!!$        real(dp)           :: delta
!!$        type(pointing), allocatable :: p(:)
!!$
!!$    end type observationslice

!!$ type, extends(observationslice) :: pacsobservationslice
    type :: pacsobservationslice

        integer*8          :: first, last, nsamples
        character(len=256) :: filename
        character          :: channel
        character(len=70)  :: observing_mode
        character(len=70)  :: unit
        integer            :: compression_factor
        real(dp)           :: delta
        type(pointing), allocatable :: p(:)

    contains
        
        procedure :: set_filename
        procedure :: set_valid_slice
        procedure :: set_channel
        procedure :: set_compression_mode
        procedure :: set_observing_mode
        procedure :: set_pointing
        procedure :: set_unit
        procedure :: validate_pointing
        procedure :: set_pointing_oldstyle

    end type pacsobservationslice

    type :: observation
        integer           :: nslices
        integer           :: nsamples_tot
        character         :: channel
        character(len=70) :: observing_mode
        character(len=70) :: unit
        type(maskarray)   :: maskarray_policy
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
        integer                        :: isample, islice, n180
        real*8                         :: ra, dec
        logical                        :: zero_minus, zero_plus

        ra0  = 0.0d0
        dec0 = 0.0d0
        zero_minus = .false.
        zero_plus  = .false.
        n180 = 0

        do islice = 1, this%nslices
            !$omp parallel do default(shared) reduction(+:n180,ra0,dec0)           &
            !$omp reduction(.or.:zero_minus,zero_plus) private(isample, ra, dec)
            do isample=1, this%slice(islice)%nsamples
                ra  = this%slice(islice)%p(isample)%ra
                dec = this%slice(islice)%p(isample)%dec
                zero_minus = zero_minus .or. ra > 270.d0
                zero_plus  = zero_plus  .or. ra <= 90.d0
                if (ra >= 180.d0) n180 = n180 + 1
                ra0  = ra0  + ra
                dec0 = dec0 + dec
            end do
            !$omp end parallel do
        end do

        if (zero_minus .and. zero_plus) ra0 = ra0 - 360.d0 * n180

        ra0  = ra0  / this%nsamples_tot
        dec0 = dec0 / this%nsamples_tot
            
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
        integer*8, intent(out)                     :: first, last
        integer, intent(out)                       :: status

        integer :: pos, length, delim
 
        status = 0
        first = 0
        last  = 0

        length = len_trim(filename)
        if (length > len(this%filename)) then
            status = 1
            write (ERROR_UNIT,'(a,i0,a,i0,a)') 'ERROR: Input filename length is&
                 & too long ', length, '. Maximum is ', len(this%filename), '.'
            return
        end if

        this%filename = filename

        if (filename(length:length) /= ']') return

        do pos = length-1, 2, -1
            if (filename(pos:pos) == '[') exit
        end do

        if (pos == 1) then
            status = 1
            write (ERROR_UNIT,'(a)') "ERROR: Missing opening bracket in '" //  &
                 trim(filename) // "'."
            return
        end if
      
        pos = pos + 1
        length = length - 1

        ! find delimiter ':'
        delim = index(filename(pos:length), ':')
        if (delim == 0) then
            delim = length + 1
        else 
            delim = pos + delim - 1
        end if

        ! read last sample
        if (delim < length) then
            read (filename(delim+1:length), '(i20)', iostat=status) last
            if (status /= 0) then
                write (ERROR_UNIT,'(a)') "ERROR: Invalid last sample: '" //    &
                      filename(pos-1:length+1) // "'."
                return
            end if
        end if

        ! read first sample
        if (delim > pos) then
            read (filename(pos:delim-1), '(i20)', iostat=status) first
            if (status /= 0) then
                write (ERROR_UNIT,'(a)') "ERROR: Invalid first sample: '" //    &
                      filename(pos-1:length+1) // "'."
                return
            end if
        end if

        if (last /= 0 .and. last < first) then
            status = 1
            write (ERROR_UNIT,'(a,2(i0,a))')"ERROR: Last sample '", last,  &
                  "' is less than the first sample '", first, "'."
            return
        end if

        if (first < 0) then
            status = 1
            write (ERROR_UNIT,'(a,i0,a)')"ERROR: The first sample '", first,&
                  "' is less than 1."
            return
        end if

        ! range specifications could be read, remove it
        this%filename(pos-1:) = ' '

    end subroutine set_filename


    !-------------------------------------------------------------------------------------------------------------------------------


    ! sets 'b', 'g' or 'r' for a blue, green or red channel observation
    subroutine set_channel(this, status)

        class(pacsobservationslice), intent(inout) :: this
        integer, intent(out)                       :: status

        integer                       :: length, unit, nsamples
        integer                       :: status_close
        character(len=2), allocatable :: channels(:)

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

        if (any(channels /= channels(1))) then
            status = 1
            write (ERROR_UNIT,'(a)') 'Observation is BANDSWITCH.'
        else
           select case (channels(1))
               case ('BS')
                   this%channel = 'b'
               case ('BL')
                   this%channel = 'g'
               case ('R ')
                   this%channel = 'r'
               case default
                   status = 1
                   write (ERROR_UNIT,'(a)') 'Invalid array BAND value: ' // channels(1)
           end select
        end if

        call ft_close(unit, status_close)
        if (status == 0) status = status_close

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
                write (ERROR_UNIT,'(a)') "ERROR: FITS keyword 'algorithm' is &
                    &not found."
                go to 999
            end if
            if (keyword == 'algorithm') exit
            ikey = ikey + 1
        end do

        call ftgkys(unit, 'META_' // strinteger(ikey), algorithm, comment, status)

        if (algorithm(1:19) == 'Floating Average  :') then
            read (algorithm(20:),'(i3)', iostat=status) this%compression_factor
            if (status /= 0) then
                write (ERROR_UNIT, '(a)') "ERROR: The compression algorithm '" &
                    // trim(algorithm) // "' is not understood."
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
    !     - if the bbtype is set, we search starting from the beginning the first valid sample
    !     - otherwise, we search starting from the end for the first sample that satisfies abs(chopfpuangle) > 0.01.
    ! if invalid samples have been found, we then discard the first 10, that might be affected by relaxation.
    subroutine set_valid_slice(this, first, last, maskarray_policy, status)

        class(pacsobservationslice), intent(inout) :: this
        integer*8, intent(in)                      :: first, last
        type(maskarray), intent(in)                :: maskarray_policy
        integer, intent(out)                       :: status

        logical(1), allocatable :: off_target(:)
        integer*8               :: nsamples
        integer                 :: length
        integer*8               :: isample
        integer*8, allocatable  :: bbtype(:)
        integer*8, allocatable  :: chop(:)

        length = len_trim(this%filename)
        if (this%filename(length-4:length) /= '.fits') then
            write (OUTPUT_UNIT,'(a)') 'Warning: Obsolete file format.'
            call ft_read_extension(trim(this%filename) // '_ChopFpuAngle.fits', chop, status)
            if (status /= 0) return
            allocate (bbtype(size(chop)))
            bbtype = 0
        else
            call ft_read_column(trim(this%filename)//'[Status]', 'BBTYPE', bbtype, status)
            if (status /= 0) return
            call ft_read_column(trim(this%filename)//'[Status]', 'CHOPFPUANGLE', chop, status)
            if (status /= 0) return
        end if

        nsamples = size(bbtype)

        if (nsamples == 0) then
            status = 1
            write (ERROR_UNIT) 'ERROR: Status extension is empty.'
            return
        end if

        allocate (off_target(nsamples))

        if (maxval(bbtype) == 0) then
            off_target = abs(chop) > 1._p
        else
            off_target = bbtype /= 3282
        end if

        if (last > nsamples) then
           status = 1
           write (ERROR_UNIT,'(a,2(i0,a))') "ERROR: The last sample '", last, "' exceeds the size of the observation '",       &
                nsamples, "'."
           return
        end if

        ! set last valid sample
        if (last == 0) then
           do isample = nsamples, 1, -1
               if (.not. off_target(isample)) exit
           end do
           if (isample == 0) then
               write (OUTPUT_UNIT,'(a)') 'Warning: No in-scan sample. Automatic search for valid slice is disabled.'
               this%first = 1
               this%last  = nsamples
               go to 999
           end if
           this%last = isample
        else
           this%last = last
        end if

        ! set first valid sample. 
        if (first == 0) then
            if (maxval(bbtype) == 0) then
                do isample = this%last, 1, -1
                    if (off_target(isample)) exit
                end do
                this%first = isample + 1
            else
                do isample = 1, this%last
                    if (.not. off_target(isample)) exit
                end do
                this%first = isample
            end if

            ! discard first 10 samples if invalid samples have been detected
            if (this%first /= 1) then
                if (this%first + 10 > this%last) then
                    write (OUTPUT_UNIT,'(a)') 'Warning: It is not possible to discard the first 10 scan samples. The observation is&
                          & too short or the specified last sample is not large enough.'
                else
                    this%first = this%first + 10
                end if
            end if
        else
            this%first = first
        end if
        
    999 this%nsamples = this%last - this%first + 1
        if (allocated(this%p)) deallocate (this%p)
        allocate (this%p(this%nsamples))
        this%p%invalid    = .false.
        this%p%wrong_time = .false.
        this%p%off_target = off_target(this%first:this%last)
        if (maskarray_policy%off_target) this%p%invalid = this%p%off_target

    end subroutine set_valid_slice


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine set_pointing(this, status)

        class(pacsobservationslice), intent(inout) :: this
        integer, intent(out)                       :: status

        real(dp), allocatable  :: buffer(:)
        integer*8, allocatable :: timeus(:)
        integer                :: nsamples, first, last, status_close
        integer                :: length, unit
        
        length = len_trim(this%filename)
        if (this%filename(length-4:length) /= '.fits') then
            call this%set_pointing_oldstyle(status)
            return
        end if
        
        if (this%nsamples <= 1) then
            status = 1
            write (ERROR_UNIT,'(a)') 'Input time has less than 2 samples.'
            return
        end if

        allocate(timeus(this%nsamples))
        allocate(buffer(this%nsamples))
        
        call ft_open_bintable(trim(this%filename) // '[Status]', unit, nsamples, status)
        if (status /= 0) return
        if (nsamples < this%nsamples) then
            status = 1
            write (ERROR_UNIT,'(a)') 'ERROR: There is not pointing information for every frame.'
            go to 999
        end if
        
        ! get first and last sample in file
        first = this%first
        last  = this%last
        
        call ft_read_column(unit, 'FINETIME', first, last, timeus, status)
        if (status /= 0) go to 999
        buffer = timeus * 1.d-6
        this%p%time = buffer
        
        call ft_read_column(unit, 'RaArray', first, last, buffer, status)
        if (status /= 0) go to 999
        this%p%ra = buffer
        
        call ft_read_column(unit, 'DecArray', first, last, buffer, status)
        if (status /= 0) go to 999
        this%p%dec = buffer
        
        call ft_read_column(unit, 'PaArray', first, last, buffer, status)
        if (status /= 0) go to 999
        this%p%pa = buffer
        
        call ft_read_column(unit, 'CHOPFPUANGLE', first, last,buffer,status)
        if (status /= 0) go to 999
        this%p%chop = buffer

    999 call ft_close(unit, status_close)
        if (status == 0) status = status_close
        if (status /= 0) return

    end subroutine set_pointing


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine set_pointing_oldstyle(this, status)

        class(pacsobservationslice), intent(inout) :: this
        integer, intent(out)                       :: status

        real(dp), allocatable  :: buffer(:)
        integer*8, allocatable :: timeus(:)
        integer                :: first, last
        
        if (this%nsamples <= 1) then
            status = 1
            write (ERROR_UNIT,'(a)') 'Input time has less than 2 samples.'
            return
        end if

        ! get first and last sample in file
        first = this%first
        last  = this%last
        
        call ft_read_extension(trim(this%filename)//'_Time.fits', timeus, status)
        if (status /= 0) return
        allocate(buffer(size(timeus)))
        buffer = timeus * 1.0d-6
        this%p%time = buffer(first:last)

        call ft_read_extension(trim(this%filename) // '_RaArray.fits', buffer, status)
        if (status /= 0) return
        this%p%ra = buffer(first:last)

        call ft_read_extension(trim(this%filename) // '_DecArray.fits', buffer, status)
        if (status /= 0) return
        this%p%dec = buffer(first:last)

        call ft_read_extension(trim(this%filename) // '_PaArray.fits', buffer, status)
        if (status /= 0) return
        this%p%pa = buffer(first:last)

        call ft_read_extension(trim(this%filename) // '_ChopFpuAngle.fits', buffer, status)
        if (status /= 0) return
        this%p%chop = buffer(first:last)

    end subroutine set_pointing_oldstyle


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

    end subroutine set_unit

    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine validate_pointing(this, maskarray_policy, status, verbose)

        class(pacsobservationslice), intent(inout) :: this
        type(maskarray), intent(in)                :: maskarray_policy
        integer, intent(out)                       :: status
        logical, intent(in)                        :: verbose

        integer                :: isample, njumps
        integer                :: compression_factor
        real(dp), allocatable  :: delta(:)
        real(dp)               :: delta_max

        status = 1

        ! check that the input time is monotonous (and increasing)
        allocate(delta(this%nsamples-1))
        delta = this%p(2:this%nsamples)%time - this%p(1:this%nsamples-1)%time
        this%delta = median(delta)
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
                if (maskarray_policy%wrong_time) this%p(isample)%invalid = .true.
                this%p(isample)%time = this%p(isample-1)%time + this%delta
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
        if (verbose .and. any(neq_real(delta, this%delta, 3))) then
            write (*,'(a,$)') "Warning: In observation '" // trim(this%filename) //"', the pointing time is not evenly spaced."
            if (delta_max > 1.5_p * this%delta) then
                write (*,'(a)') ' Largest gap is ' // strreal(delta_max*1000._p,1) // 'ms.'
            else
                write (*,*)
            end if
        end if
        
        ! check the compression factor from the data themselves
        compression_factor = nint(this%delta / 0.024996_dp)
        if (neq_real(compression_factor * 0.024996_dp, this%delta, 2)) then
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

    end subroutine validate_pointing


    ! the following belongs to module_pacsobservation.f
    ! moved here because of gfortran bug #44065

    subroutine init(this, filename, maskarray_policy, status, verbose)

        class(pacsobservation), intent(inout) :: this
        character(len=*), intent(in)          :: filename(:)
        type(maskarray), intent(in)           :: maskarray_policy
        integer, intent(out)                  :: status
        logical, intent(in), optional         :: verbose
        integer*8                             :: first, last
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

        ! set the mask policy '.true.' means 'taken into account'
        this%maskarray_policy = maskarray_policy        
        
        if (allocated(this%slice)) then
            do islice = 1, this%nslices
                if (allocated(this%slice(islice)%p)) deallocate (this%slice(islice)%p)
            end do
            deallocate(this%slice)
        end if
        allocate(this%slice(this%nslices))

        do islice = 1, this%nslices

            call this%slice(islice)%set_filename(filename(islice), first, last, status)
            if (status /= 0) return

            call this%slice(islice)%set_valid_slice(first, last, this%maskarray_policy, status)
            if (status /= 0) return

            call this%slice(islice)%set_channel(status)
            if (status /= 0) return
                
            call this%slice(islice)%set_compression_mode(status)
            if (status /= 0) return
        
            call this%slice(islice)%set_observing_mode(status)
            if (status /= 0) return

            call this%slice(islice)%set_unit(status)
            if (status /= 0) return

            call this%slice(islice)%set_pointing(status)
            if (status /= 0) return
                
            call this%slice(islice)%validate_pointing(this%maskarray_policy, status, verbose_)
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

        this%nsamples_tot = sum(this%slice%nsamples)

    end subroutine init


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine init_sim(this, time, ra, dec, pa, chop, status, verbose)
        class(pacsobservation), intent(inout) :: this
        real*8, intent(in)                    :: time(:),ra(:),dec(:),pa(:),chop(:)
        integer, intent(out)                  :: status
        logical, intent(in), optional         :: verbose
        integer*8                             :: nsamples
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
        this%nsamples_tot = nsamples
        this%channel      = ' '
        this%observing_mode = 'Unknown'
        this%maskarray_policy = maskarray(wrong_time=.true., off_target=.true.)

        allocate (this%slice(1))
        this%slice(1)%first    = 1
        this%slice(1)%last     = nsamples
        this%slice(1)%nsamples = nsamples
        this%slice(1)%filename = 'simulation'
        this%slice(1)%channel  = ' '
        this%slice(1)%observing_mode = 'Unknown'
        this%slice(1)%compression_factor = nint((time(2)-time(1)) / 0.024996_dp)
        this%slice(1)%delta = time(2)-time(1)
        allocate (this%slice(1)%p(nsamples))

        this%slice(1)%p%time       = time
        this%slice(1)%p%ra         = ra
        this%slice(1)%p%dec        = dec
        this%slice(1)%p%pa         = pa
        this%slice(1)%p%chop       = chop
        this%slice(1)%p%invalid    = .false.
        this%slice(1)%p%wrong_time = .false.
        this%slice(1)%p%off_target = .false.

        call this%slice(1)%validate_pointing(this%maskarray_policy, status, verbose_)

    end subroutine init_sim


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine print(this)

        class(pacsobservation), intent(in) :: this

        integer :: islice, nrejected, nsamples

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

            ! compression factor
            if (this%unit /= 'Jy') then
                write (OUTPUT_UNIT,'(a)') '      Unit: ' // this%slice(islice)%unit
            end if

            ! print maskarray information
            if (this%maskarray_policy%off_target) then
                nrejected = count(this%slice(islice)%p%off_target)
                if (nrejected > 0) then
                    write (OUTPUT_UNIT,'(a,2(i0,a))') "      Mask 'off target': Rejecting ", nrejected, ' / ', nsamples, '.'
                end if
            end if

            if (this%maskarray_policy%wrong_time) then
                nrejected = count(this%slice(islice)%p%wrong_time)
                if (nrejected > 0) then
                    write (OUTPUT_UNIT,'(a,2(i0,a))') "      Mask 'wrong time': Rejecting ", nrejected, ' / ', nsamples, '.'
                end if
            end if

            write (*,*)

        end do

    end subroutine print

end module module_observation
