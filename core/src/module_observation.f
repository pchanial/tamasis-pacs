! Because of bugs in gfortran and ifort related to array of classes
! the architecture is (temporarily hopefully) a mess
! this module should define a derived type: observationslice and the
! module module_PacsObservation should extend it and define PacsObservationSlice
! for now, many routines that should belong to module_PacsObservation are here...

module module_observation

    use iso_fortran_env,  only : ERROR_UNIT, OUTPUT_UNIT
    use module_fitstools, only : FLEN_VALUE, ft_close, ft_open, ft_open_bintable, ft_read_column, ft_read_image, ft_read_keyword,  &
                                 ft_read_keyword_hcss, ft_check_error_cfitsio
    use module_math,      only : NaN, barycenter_lonlat, mean, mean_degrees, median, neq_real
    use module_string,    only : strinteger, strlowcase, strreal, strsection, strternary
    use module_tamasis,   only : p, POLICY_KEEP, POLICY_MASK, POLICY_REMOVE
    implicit none
    private

    public :: MaskPolicy
    public :: Observation
!!$    public :: ObservationSlice
public :: PacsObservation    
public :: PacsObservationSlice
    public :: Pointing
    public :: POINTING_INSCAN, POINTING_TURNAROUND, POINTING_OTHER

    integer, parameter :: POINTING_INSCAN = 1
    integer, parameter :: POINTING_TURNAROUND = 2
    integer, parameter :: POINTING_OTHER = 3

    type MaskPolicy
        integer :: inscan = POLICY_KEEP
        integer :: turnaround = POLICY_KEEP
        integer :: other = POLICY_REMOVE
        integer :: invalid = POLICY_MASK
    end type MaskPolicy

    type Pointing
        real(p)   :: time
        real(p)   :: ra
        real(p)   :: dec
        real(p)   :: pa
        real(p)   :: chop

        logical*1 :: inscan
        logical*1 :: turnaround
        logical*1 :: other      ! not inscan and not turnaround
        logical*1 :: invalid    ! wrong finetime, undef. band

        logical*1 :: removed
        logical*1 :: masked

        logical*1 :: notused1
        logical*1 :: notused2
    end type Pointing

!!$    type observationslice
!!$        integer                     :: nsamples, nvalids
!!$        character(len=256)          :: filename
!!$        character                   :: band
!!$        character(len=FLEN_VALUE)   :: observing_mode
!!$        integer                     :: compression_factor
!!$        real(dp)                    :: sampling_interval
!!$        type(pointing), allocatable :: p(:)
!!$    end type observationslice

    real(p), parameter :: PACS_SAMPLING = 0.024996_p

!!$ type, extends(observationslice) :: PacsObservationSlice
    type PacsObservationSlice

        integer                   :: nsamples, nvalids
        character(len=256)        :: filename
        character(len=FLEN_VALUE) :: band
        character(len=FLEN_VALUE) :: observing_mode
        integer                   :: obsid
        real(p)                   :: ra, dec
        real(p)                   :: cam_angle, scan_angle, scan_length, scan_step
        integer                   :: scan_nlegs
        character(len=FLEN_VALUE) :: unit
        integer                   :: compression_factor
        real(p)                   :: offset
        real(p)                   :: sampling_interval
        integer                   :: nmasks
        character(len=FLEN_VALUE), allocatable :: mask_name(:)
        logical, allocatable                   :: mask_activated(:)
        type(pointing), allocatable            :: p(:)

    contains
        
        procedure :: set_band
        procedure :: set_filename
        procedure :: set_policy
        procedure :: set_observing_mode
        procedure :: set_astrometry
        procedure :: set_mask
        procedure :: set_unit
        procedure :: set_flags
        procedure :: set_sampling_interval
        procedure :: set_offset

    end type PacsObservationSlice

    type Observation
        integer                   :: nslices
        integer                   :: nsamples  ! number of pointings
        integer                   :: nvalids   ! number of non-removed pointings, to be included in tod & pmatrix
        character(len=FLEN_VALUE) :: band
        character(len=FLEN_VALUE) :: unit
        type(MaskPolicy)          :: policy
        type(PacsObservationSlice), allocatable :: slice(:)
    contains
        procedure :: compute_center
        procedure :: get_position_time
        procedure :: get_position_index
        procedure :: destructor
    end type Observation

    ! the following belongs to module_PacsObservation.f
    ! moved here because of gfortran bug #44065
    type, extends(Observation) :: PacsObservation

    contains

        procedure :: init
        procedure :: init_with_variables
        procedure :: print

    end type PacsObservation


contains


    ! linear interpolation if time samples are not evenly spaced.
    ! for the first call, index should be set to zero
    subroutine get_position_time(this, islice, time, ra, dec, pa, chop, index)
        class(observation), intent(in) :: this
        integer, intent(in)            :: islice
        real(p), intent(in)            :: time
        real(p), intent(out)           :: ra, dec, pa, chop
        integer, intent(inout)         :: index
        integer                        :: i, nsamples
        real(p)                        :: frac

        if (islice < 1 .or. islice > this%nslices) then
            write (ERROR_UNIT,'(a,i0,a)') "GET_POSITION: Invalid slice number '", islice, "'."
            stop
        end if

        nsamples = this%slice(islice)%nsamples

        if (time > 2._p * this%slice(islice)%p(nsamples)%time - this%slice(islice)%p(nsamples-1)%time .or.                         &
            time < 2._p * this%slice(islice)%p(1)%time - this%slice(islice)%p(2)%time) then
            ra   = NaN
            dec  = NaN
            pa   = NaN
            chop = NaN
            return
        end if

        if (index == 0 .or. index-1 > nsamples) then
            index = 2
        else if (time <= this%slice(islice)%p(index-1)%time) then
            index = 2
        end if

        do i = index, nsamples
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
    subroutine get_position_index(this, islice, itime, sampling_factor, offset, ra, dec, pa, chop)

        class(observation), intent(in) :: this
        integer, intent(in)            :: islice, itime, sampling_factor
        real(p), intent(in)            :: offset
        real(p), intent(out)           :: ra, dec, pa, chop
        integer                        :: i, itime_max
        real(p)                        :: frac

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
        frac = real(itime - 1 - sampling_factor * i, p) / sampling_factor - offset
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
        real(p), intent(out)           :: ra0, dec0

        real(p), allocatable           :: ra(:), dec(:)
        integer                        :: islice

        allocate (ra(this%nslices), dec(this%nslices))
        do islice = 1, this%nslices
            call barycenter_lonlat(pack(this%slice(islice)%p%ra,  .not. this%slice(islice)%p%removed),                             &
                                   pack(this%slice(islice)%p%dec, .not. this%slice(islice)%p%removed),                             &
                                   ra(islice), dec(islice))
        end do

        call barycenter_lonlat(ra, dec, ra0, dec0)
        
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

        class(PacsObservationSlice), intent(inout) :: this
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
            write (ERROR_UNIT,'(a,i0,a,i0,a)') 'Error: Input filename length is too long ', length, '. Maximum is ',               &
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
            write (ERROR_UNIT,'(a)') "Error: Missing opening bracket in '" // trim(filename) // "'."
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
                write (ERROR_UNIT,'(a)') "Error: Invalid first sample: '" // filename(pos-1:length+1) // "'."
                return
            end if
        end if

        ! read last sample
        if (is_slice) then
            if (delim < length) then
                read (filename(delim+1:length), '(i20)', iostat=status) last
                if (status /= 0) then
                    write (ERROR_UNIT,'(a)') "Error: Invalid last sample: '" // filename(pos-1:length+1) // "'."
                    return
                end if
            end if
        else
            last = first
        end if

        if (last /= 0 .and. last < first) then
            status = 1
            write (ERROR_UNIT,'(a,2(i0,a))')"Error: Last sample '", last,  "' is less than the first sample '", first, "'."
            return
        end if

        if (first < 0) then
            status = 1
            write (ERROR_UNIT,'(a,i0,a)')"Error: The first sample '", first, "' is less than 1."
            return
        end if

        ! range specifications could be read, remove it
        this%filename(pos-1:) = ' '

    end subroutine set_filename


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine set_policy(this, policy, first, last, status)

        class(PacsObservationSlice), intent(inout) :: this
        type(MaskPolicy), intent(in)               :: policy
        integer, intent(in)                        :: first, last
        integer, intent(out)                       :: status

        status = 0

        this%p%masked  = (policy%inscan     == POLICY_MASK)   .and. this%p%inscan     .or.                                         &
                         (policy%turnaround == POLICY_MASK)   .and. this%p%turnaround .or.                                         &
                         (policy%other      == POLICY_MASK)   .and. this%p%other      .or.                                         &
                         (policy%invalid    == POLICY_MASK)   .and. this%p%invalid

        this%p%removed = (policy%inscan     == POLICY_REMOVE) .and. this%p%inscan     .or.                                         &
                         (policy%turnaround == POLICY_REMOVE) .and. this%p%turnaround .or.                                         &
                         (policy%other      == POLICY_REMOVE) .and. this%p%other      .or.                                         &
                         (policy%invalid    == POLICY_REMOVE) .and. this%p%invalid
        
        if (first > this%nsamples) then
            write (ERROR_UNIT,'(a,2(i0,a))') "Error: In file '" // trim(this%filename) // "', the lower bound '", first,           &
                 "' exceeds the size of the observation '", this%nsamples, "'."
            status = 1
            return
        end if

        ! if a slice is specified, remove its complement
        this%p(1:first-1)%removed = .true.
        if (last > this%nsamples) then
            status = 1
            write (ERROR_UNIT,'(a,2(i0,a))') "Error: In file '" // trim(this%filename) // "', the upper bound '", last,            &
                 "' exceeds the size of the observation '", this%nsamples, "'."
            return
        end if
        if (last /= 0) then
            this%p(last+1:)%removed = .true.
        end if
        
        ! get number of valid pointings
        this%nvalids  = count(.not. this%p%removed)

    end subroutine set_policy


    !-------------------------------------------------------------------------------------------------------------------------------


    ! sets 'b', 'g' or 'r' for a blue, green or red band observation
    subroutine set_band(this, status)

        class(PacsObservationSlice), intent(inout) :: this
        integer, intent(out)                       :: status

        integer                       :: nsamples, unit
        character(len=5), allocatable :: bands(:)
        character(len=FLEN_VALUE)     :: camname, blue
        logical                       :: found


        call ft_open(trim(this%filename), unit, status)
        if (status /= 0) return

        call ft_read_keyword_hcss(unit, 'camName', camname, status=status)
        if (status /= 0) return

        call ft_read_keyword_hcss(unit, 'blue', blue, found, status)
        if (status /= 0) return
        if (.not. found) then
            call ft_read_keyword(unit, 'FILTER', blue, status=status)
            if (status /= 0) return
        end if

        call ft_close(unit, status)
        if (status /= 0) return

        if (camname == 'Blue Photometer') then
            if (blue == 'blue1' .or. blue == 'blue70um') then
                this%band = 'blue'
            else
                this%band = 'green'
            end if
        else
            this%band = 'red'
        end if

        call ft_open_bintable(trim(this%filename) // '[Status]', unit, nsamples, status)
        if (status /= 0) return

        allocate(bands(nsamples))
        call ft_read_column(unit, 'BAND', 1, nsamples, bands, status)
        if (status /= 0) return

        call ft_close(unit, status)
        if (status /= 0) return

        select case (this%band)
            case ('blue')
                this%p%invalid = this%p%invalid .or. bands /= 'BS'
            case ('green')
                this%p%invalid = this%p%invalid .or. bands /= 'BL'
            case ('red')
                this%p%invalid = this%p%invalid .or. bands /= 'R'
        end select

    end subroutine set_band


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine set_observing_mode(this, status)

        class(PacsObservationSlice), intent(inout) :: this
        integer, intent(out)                       :: status

        logical                   :: found
        integer                   :: unit
        character(len=FLEN_VALUE) :: algorithm, mode

        call ft_open(trim(this%filename), unit, status)
        if (status /= 0) return

        call ft_read_keyword(unit, 'obs_id', this%obsid, found, status)
        if (status /= 0) return

        call ft_read_keyword_hcss(unit, 'algorithm', algorithm, status=status)
        if (status /= 0) return

        if (algorithm(1:19) == 'Floating Average  :') then
            read (algorithm(20:),'(i3)', iostat=status) this%compression_factor
            if (status /= 0) then
                write (ERROR_UNIT, '(a)') "Error: The compression algorithm '" // trim(algorithm) // "' is not understood."
                return
            end if
        else if (algorithm == 'None') then
            this%compression_factor = 1
        else
            status = 1
            write (ERROR_UNIT, '(a)') "Error: The compression algorithm '" // trim(algorithm) // "' is not understood."
        end if

        call ft_read_keyword_hcss(unit, 'mapScanAngle',     this%scan_angle, found, status)
        if (status /= 0) return

        call ft_read_keyword_hcss(unit, 'mapScanLegLength', this%scan_length, found, status)
        if (status /= 0) return        

        call ft_read_keyword_hcss(unit, 'mapScanNumLegs',   this%scan_nlegs, found, status)
        if (status /= 0) return

        call ft_read_keyword(unit, 'ra', this%ra, found, status)
        if (status /= 0) return

        call ft_read_keyword(unit, 'dec', this%dec, found, status)
        if (status /= 0) return

        call ft_read_keyword(unit, 'cusmode', mode, found, status)
        if (status /= 0) return
        if (.not. found) then

            call ft_read_keyword_hcss(unit, 'compMode', mode, status=status)
            if (status /= 0) return

            select case (mode)
                case ('Photometry Lossless Compression Mode')
                    this%observing_mode = 'transparent'
                case ('Photometry Default Mode')
                    if (this%band == 'red') then
                        write (OUTPUT_UNIT,'(a)') 'Warning: Prime or parallel mode cannot be inferred.'
                    end if
                    this%observing_mode = 'prime'
                case ('Photometry Double Compression Mode')
                    this%observing_mode = 'parallel'
                case default
                    status = 1
                    write (ERROR_UNIT,'(a)') "SET_OBSERVING_MODE: In file '" // trim(this%filename) // "', unhandled value '" //   &
                          trim(mode) // "' for keyword 'compMode'."
                    return
            end select

        else
            
            select case (mode)
                case ('PacsPhoto')
                    this%observing_mode = 'prime'
                case ('SpirePacsParallel')
                    this%observing_mode = 'parallel'
                case ('__PacsTranspScan')
                    this%observing_mode = 'transparent'
                case ('__Calibration')
                    this%observing_mode = 'unknown'
                case default
                    status = 1
                    write (ERROR_UNIT,'(a)') "SET_OBSERVING_MODE: In file '" // this%filename // "', unhandled value '" //         &
                          trim(mode) // "' for keyword 'CUSMODE'."
                    return
            end select

        end if

        this%cam_angle = 0.
        this%scan_step = 0.

        call ft_close(unit, status)
        if (status /= 0) return

    end subroutine set_observing_mode


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine set_flags(this, status)

        class(PacsObservationSlice), intent(inout) :: this
        integer, intent(out)                       :: status

        logical(1), allocatable :: inscan(:)
        logical(1), allocatable :: turnaround(:)
        integer                 :: unit
        integer                 :: isample, itarget
        integer*8, allocatable  :: bbid(:)
        real(p), allocatable    :: chop(:)

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

        if (this%nsamples == 0) then
            status = 1
            write (ERROR_UNIT, '(a)') 'Error: Status extension is empty.'
            return
        end if

        ! handle observations with missing BBIDs
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

        allocate (inscan(this%nsamples))
        allocate (turnaround(this%nsamples))

        inscan = bbid == z'cd2'
        turnaround = bbid == z'4000'

        ! search for the first inscan sample
        do itarget = 1, this%nsamples
            if (inscan(itarget)) exit
        end do
        if (itarget <= this%nsamples) then
            ! backward search of the first non turnaround sample (but don't go back too much, we would include too much of the slew)
            do isample = itarget-1, max(itarget-201,1), -1
                if (.not. turnaround(isample)) exit
            end do
            ! mask as off-scan all samples before the 0x4000 block that precedes the first inscan block
            turnaround(1:isample) = .false.
        end if
        
        if (allocated(this%p)) deallocate (this%p)
        allocate (this%p(this%nsamples))
        this%p%invalid = .false.
        this%p%inscan = inscan
        this%p%turnaround = turnaround
        this%p%other = .not. inscan .and. .not. turnaround

    end subroutine set_flags


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine set_astrometry(this, status)

        class(PacsObservationSlice), intent(inout) :: this
        integer, intent(out)                       :: status

        integer*8, allocatable :: timeus(:)
        integer                :: nsamples
        integer                :: unit
        
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
        this%p%time = (timeus - timeus(1)) * 1.e-6_p
        
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


    subroutine set_mask(this, status)

        class(PacsObservationSlice), intent(inout) :: this
        integer, intent(out)                       :: status

        integer :: unit, hdutype, i
        logical :: found
        character(len=8*FLEN_VALUE) :: info
        integer, allocatable        :: mask_extension(:)

        call ft_open(trim(this%filename) // '[MASK]', unit, found, status)
        if (status /= 0) return
        if (.not. found) then
            this%nmasks = 0
            allocate (this%mask_name(0), this%mask_activated(0))
            return
        end if

        call ft_read_keyword(unit, 'DSETS___', this%nmasks, status=status)
        if (status /= 0) return

        allocate (mask_extension(this%nmasks), this%mask_name(this%nmasks), this%mask_activated(this%nmasks))
        do i=1, this%nmasks
            call ft_read_keyword(unit, 'DS_' // strinteger(i-1), mask_extension(i), status=status)
            if (status /= 0) return
        end do
        
        do i=1, this%nmasks
            call FTMAHD(unit, mask_extension(i)+1, hdutype, status)
            if (ft_check_error_cfitsio(status, unit, trim(this%filename))) return

            call ft_read_keyword(unit, 'EXTNAME', this%mask_name(i), status=status)
            if (status /= 0) return
            call ft_read_keyword(unit, 'INFO____', info, status=status)
            if (status /= 0) return
            this%mask_activated(i) = .not. index(info, 'deactivated') > 0
            
        end do

        call ft_close(unit, status)
        
    end subroutine set_mask


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine set_unit(this, status)

        class(PacsObservationSlice), intent(inout) :: this
        integer, intent(out)                       :: status

        integer :: unit
        logical :: found

        this%unit = ''

        call ft_open(trim(this%filename) // '[Signal]', unit, found, status)
        if (status /= 0 .or. .not. found) return

        call ft_read_keyword(unit, 'QTTY____', this%unit, found, status)
        if (status /= 0) return

        call ft_close(unit, status)

        select case (this%unit)
            case ('Jy')
                this%unit = 'Jy / detector_reference'
            case ('V')
                this%unit = 'V / detector_reference'
        end select

    end subroutine set_unit


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine set_sampling_interval(this, status, verbose)

        class(PacsObservationSlice), intent(inout) :: this
        integer, intent(out)                       :: status
        logical, intent(in)                        :: verbose

        integer              :: isample, njumps
        integer              :: compression_factor
        real(p), allocatable :: delta(:)
        real(p)              :: delta_max

        ! special case if there is only one sample
        if (this%nsamples == 1) then
            this%sampling_interval = this%compression_factor * PACS_SAMPLING
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
                        write (ERROR_UNIT,'(a)') "Error: The pointing time is not strictly increasing."
                        return
                    end if
                end if
                njumps = njumps + 1
                this%p(isample)%invalid = .true.
                this%p(isample)%time = this%p(isample-1)%time + this%sampling_interval
                this%p(isample)%ra   = this%p(isample-1)%ra  
                this%p(isample)%dec  = this%p(isample-1)%dec 
                this%p(isample)%pa   = this%p(isample-1)%pa  
                this%p(isample)%chop = this%p(isample-1)%chop
            end if
        end do
        
        if (njumps > 0) then
            if (verbose) then
                write (OUTPUT_UNIT,'(a,i0,a)') "Warning: In file '" // trim(this%filename) //                                      &
                     "', the pointing fine time has ", njumps, ' negative jump(s). The affected fra&
                     &mes have been marked as invalid.'
            end if
            delta = this%p(2:this%nsamples)%time - this%p(1:this%nsamples-1)%time
        end if
        
        ! check if there are gaps
        delta_max = maxval(abs(delta))
        if (verbose .and. any(neq_real(delta, this%sampling_interval, 1.e-3_p))) then
            write (OUTPUT_UNIT,'(a,$)') "Warning: In file '" // trim(this%filename) //"', the pointing time is not evenly spaced."
            if (delta_max > 1.5_p * this%sampling_interval) then
                write (OUTPUT_UNIT,'(a)') ' Largest gap is ' // strreal(delta_max*1000._p,1) // 'ms.'
            else
                write (OUTPUT_UNIT,*)
            end if
        end if
        
        ! check the compression factor from the data themselves
        compression_factor = nint(this%sampling_interval / PACS_SAMPLING)
        if (neq_real(compression_factor * PACS_SAMPLING, this%sampling_interval, 1e-2_p)) then
            write (ERROR_UNIT,'(a)') 'Error: The sampling time is not an integer number of PACS sampling time (40Hz).'
            return
        end if
        
        if (compression_factor /= this%compression_factor) then
            write (OUTPUT_UNIT,'(a,2(i0,a))') "Warning: The compression factor determined from the observation header '",          &
                  this%compression_factor, "' is different from the one inferred from the fine time in the Status table '",        &
                  compression_factor, "'. The latter is assumed to be true."
            this%compression_factor = compression_factor
        end if

        status = 0

    end subroutine set_sampling_interval


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine set_offset(this)
    ! offset of the fine time of the compressed frame wrt the fine time of the first frame
    ! the fine time of a compressed frame is the average of the fine times of the contributing frames
    ! the offset is in fraction of sampling interval

        class(PacsObservationSlice), intent(inout) :: this

        this%offset = (this%compression_factor - 1) / (2._p * this%compression_factor)

    end subroutine set_offset


    !-------------------------------------------------------------------------------------------------------------------------------


    ! the following belongs to module_PacsObservation.f
    ! moved here because of gfortran bug #44065
    subroutine init(this, filename, policy, status, verbose)

        class(PacsObservation), intent(inout) :: this
        character(len=*), intent(in)          :: filename(:)
        type(MaskPolicy), intent(in)          :: policy
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

            call this%slice(islice)%set_band(status)
            if (status /= 0) return
                
            call this%slice(islice)%set_observing_mode(status)
            if (status /= 0) return

            call this%slice(islice)%set_mask(status)
            if (status /= 0) return

            call this%slice(islice)%set_unit(status)
            if (status /= 0) return

            call this%slice(islice)%set_astrometry(status)
            if (status /= 0) return
            
            call this%slice(islice)%set_sampling_interval(status, verbose_)
            if (status /= 0) return

            call this%slice(islice)%set_offset()

            call this%slice(islice)%set_policy(this%policy, first, last, status)
            if (status /= 0) return

        end do

        ! make sure the band is the same for all observations
        this%band = this%slice(1)%band
        if (any(this%slice%band /= this%band)) then
            status = 1
            write (ERROR_UNIT,'(a)') 'Error: Observations do not have the same band: '
            do islice = 1, this%nslices
                write (ERROR_UNIT,'(a)') '- ' // trim(this%slice(islice)%band)
            end do
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
                write (OUTPUT_UNIT,'(a)') "Warning: Observation has no units. Assuming 'Jy / detector_reference'."
            end if
            this%unit = 'Jy / detector_reference'
            this%slice%unit = 'Jy / detector_reference'
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


    subroutine init_with_variables(this, time, ra, dec, pa, chop, masked, removed, nsamples, compression_factor, delay, status)
        ! delay is the instrument lag wrt the telescope

        class(PacsObservation), intent(inout) :: this
        real(p), intent(in)                   :: time(:),ra(:),dec(:),pa(:),chop(:)
        logical*1, intent(in)                 :: masked(:), removed(:)
        integer, intent(in)                   :: nsamples(:)
        integer, intent(in)                   :: compression_factor(:)
        real(p), intent(in)                   :: delay(:)
        integer, intent(out)                  :: status

        integer :: islice, dest, nsamples_tot

        status = 1

        ! check conformity of time, ra, dec, pa, etc.
        nsamples_tot = size(time)

        if (size(ra) /= nsamples_tot) then
            write (ERROR_UNIT,'(a)') "INIT_WITH_VAR: Input argument 'ra' has an invalid number of samples."
            return
        end if
        if (size(dec) /= nsamples_tot) then
            write (ERROR_UNIT,'(a)') "INIT_WITH_VAR: Input argument 'dec' has an invalid number of samples."
            return
        end if
        if (size(pa) /= nsamples_tot) then
            write (ERROR_UNIT,'(a)') "INIT_WITH_VAR: Input argument 'pa' has an invalid number of samples."
            return
        end if
        if (size(chop) /= nsamples_tot) then
            write (ERROR_UNIT,'(a)') "INIT_WITH_VAR: Input argument 'chop' has an invalid number of samples."
            return
        end if
        if (size(masked) /= nsamples_tot) then
            write (ERROR_UNIT,'(a)') "INIT_WITH_VAR: Input argument 'masked' has an invalid number of samples."
            return
        end if
        if (size(removed) /= nsamples_tot) then
            write (ERROR_UNIT,'(a)') "INIT_WITH_VAR: Input argument 'removed' has an invalid number of samples."
            return
        end if
        if (sum(nsamples) /= nsamples_tot) then
            write (ERROR_UNIT,'(a)') "INIT_WITH_VAR: Input argument 'nsamples' has an invalid number of samples."
            return
        end if
        if (size(compression_factor) /= size(nsamples)) then
            write (ERROR_UNIT,'(a)') "INIT_WITH_VAR: Input argument 'compression_factor' has an invalid number of elements."
            return
        end if

        allocate (this%slice(size(nsamples)))
        this%slice%nsamples = nsamples
        this%slice%compression_factor = compression_factor
        this%slice%offset = (compression_factor-1) / (2._p * compression_factor) +                                                 &
                            delay / (1000 * PACS_SAMPLING * compression_factor)
        this%slice%filename = ''
        this%slice%band  = ''
        this%slice%observing_mode = ''
        this%slice%sampling_interval = 0._p

        dest = 1
        do islice = 1, size(nsamples)

            this%slice(islice)%nvalids  = count(.not. removed(dest:dest+nsamples(islice)-1))
            allocate (this%slice(islice)%p(nsamples(islice)))
            this%slice(islice)%p%time    = time   (dest:dest+nsamples(islice)-1)
            this%slice(islice)%p%ra      = ra     (dest:dest+nsamples(islice)-1)
            this%slice(islice)%p%dec     = dec    (dest:dest+nsamples(islice)-1)
            this%slice(islice)%p%pa      = pa     (dest:dest+nsamples(islice)-1)
            this%slice(islice)%p%chop    = chop   (dest:dest+nsamples(islice)-1)
            this%slice(islice)%p%masked  = masked (dest:dest+nsamples(islice)-1)
            this%slice(islice)%p%removed = removed(dest:dest+nsamples(islice)-1)
            dest = dest + nsamples(islice)

        end do

        this%nslices  = size(nsamples)
        this%nsamples = nsamples_tot
        this%nvalids  = sum(this%slice%nvalids)
        this%band     = ''
        this%unit     = ''

        status = 0

    end subroutine init_with_variables


    !-------------------------------------------------------------------------------------------------------------------------------


    pure function strpolicy_len(policy)
        
        integer, intent(in) :: policy
        integer             :: strpolicy_len

        select case (policy)
            case (POLICY_KEEP)
                strpolicy_len = 4
            case (POLICY_MASK)
                strpolicy_len = 6
            case (POLICY_REMOVE)
                strpolicy_len = 7
            case default
                strpolicy_len = 0
        end select

    end function strpolicy_len


    !-------------------------------------------------------------------------------------------------------------------------------


    function strpolicy(policy)

        integer, intent(in)                  :: policy
        character(len=strpolicy_len(policy)) :: strpolicy

        select case (policy)
            case (POLICY_KEEP)
                strpolicy = 'kept'
            case (POLICY_MASK)
                strpolicy = 'masked'
            case (POLICY_REMOVE)
                strpolicy = 'removed'
        end select

    end function strpolicy


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine print(this)

        class(PacsObservation), intent(in) :: this

        integer :: islice, nsamples

        do islice = 1, this%nslices

            nsamples = this%slice(islice)%nsamples

            ! observation number & file name
            write (OUTPUT_UNIT,'(a6,a)') strternary(this%nslices>1, '  #' // strinteger(islice,3), ''), 'Observation: ' //         &
                  trim(this%slice(islice)%filename)
            
            ! band
            write (OUTPUT_UNIT,'(a,$)') "      Band: "
            select case (this%slice(islice)%band)
                case ('blue')
                    write (OUTPUT_UNIT,'(a)') 'Blue'
                case ('green')
                    write (OUTPUT_UNIT,'(a)') 'Green'
                case ('red')
                    write (OUTPUT_UNIT,'(a)') 'Red'
                case default
                    write (OUTPUT_UNIT,'(a)') 'Unknown'
            end select

            ! observing mode
            write (OUTPUT_UNIT,'(a)') '      Observing mode: ' // trim(this%slice(islice)%observing_mode)

            ! compression factor
            write (OUTPUT_UNIT,'(a,i0)') '      Compression factor: ', this%slice(islice)%compression_factor

            ! unit
            write (OUTPUT_UNIT,'(a)') '      Unit: ' // trim(this%slice(islice)%unit)

            ! print mask information
            write (OUTPUT_UNIT,'(a,i0,a)') '      In-scan:    ', count(this%slice(islice)%p%inscan),     ' ' //                    &
                  strpolicy(this%policy%inscan)
            write (OUTPUT_UNIT,'(a,i0,a)') '      Turnaround: ', count(this%slice(islice)%p%turnaround), ' ' //                    &
                  strpolicy(this%policy%turnaround)
            write (OUTPUT_UNIT,'(a,i0,a)') '      Other:      ', count(this%slice(islice)%p%other),      ' ' //                    &
                  strpolicy(this%policy%other)
            write (OUTPUT_UNIT,'(a,i0,a)') '      Invalid:    ', count(this%slice(islice)%p%invalid),    ' ' //                    &
                  strpolicy(this%policy%invalid)

        end do

    end subroutine print

end module module_observation
