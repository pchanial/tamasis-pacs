!-------------------------------------------------------------------------------
!
! This module defines the pacsobservation derived type.
! It contains all the information of a PACS observation relevant to its data
! processing. In particular:
!    - file name of the observation
!    - channel ('b', 'g' or 'r')
!    - transparent_mode (false or true)
!    - section to be processed
!
! Author: Pierre Chanial
!
!-------------------------------------------------------------------------------

module module_pacsobservation

    use ISO_FORTRAN_ENV,  only : ERROR_UNIT, OUTPUT_UNIT
    use module_fitstools, only : ft_open, ft_open_image, ft_open_bintable,  &
                                 ft_read_column, ft_readextension,          &
                                 ft_readslice, ft_test_extension,ft_close
    use precision,        only : p
    use string,           only : strinteger, strlowcase, strsection, strternary
    implicit none
    private

    public :: pacsobservation
    public :: init_pacsobservation
    public :: read_pacsobservation

    type pacsobservation

        private
        character(len=256), public :: filename
        character,          public :: channel
        logical,            public :: transparent_mode
        integer,            public :: compression_factor
        integer*8,          public :: nsamples
        integer*8,          public :: first
        integer*8,          public :: last

    contains

        private
        procedure :: init
        procedure :: read
        procedure :: read_oldstyle
        procedure :: set_channel
        procedure :: set_compression_mode
        procedure :: set_filename
        procedure :: set_transparent_mode
        procedure :: set_valid_slice

    end type pacsobservation

    interface init_pacsobservation
        module procedure init_scalar, init_array
    end interface init_pacsobservation

    interface read_pacsobservation
        module procedure read_scalar, read_array
    end interface read_pacsobservation


contains


    subroutine init(this, filename, status)
        class(pacsobservation), intent(inout) :: this
        character(len=*), intent(in)          :: filename
        integer, intent(out)                  :: status
        integer*8                             :: first, last

        call this%set_filename(filename, first, last, status)
        if (status /= 0) return

        call this%set_valid_slice(first, last, status)
        if (status /= 0) return

        this%nsamples = this%last - this%first + 1

        call this%set_channel(status)
        if (status /= 0) return

        call this%set_compression_mode(status)
        if (status /= 0) return
        
        call this%set_transparent_mode(status)
        if (status /= 0) return
        
    end subroutine init


    !---------------------------------------------------------------------------


    subroutine init_scalar(obs, filename, status)
        class(pacsobservation), intent(inout) :: obs
        character(len=*), intent(in)          :: filename
        integer, intent(out)                  :: status

        ! initialise observation
        call obs%init(filename, status)
        if (status /= 0) return

        ! print some info
        write (OUTPUT_UNIT,'(a)') 'Info: Section: ' //           &
             strsection(obs%first, obs%last)

    end subroutine init_scalar


    !---------------------------------------------------------------------------


    subroutine init_array(obs, filename, status)
        type(pacsobservation), intent(inout) :: obs(:)
        character(len=*), intent(in)         :: filename(:)
        integer, intent(out)                 :: status
        integer   :: iobs

        ! parameter checking
        if (size(obs) == 0) then 
            status = 1
            write (ERROR_UNIT,'(a)')'INIT_PACSOBSERVATION: Array has zero size.'
            return
        end if
        if (size(obs) /= size(filename)) then
            status = 1
            write (ERROR_UNIT,'(a,2(i0,a))') "INIT_PACSOBSERVATION: Invalid num&
                &ber of input files '", size(filename), "' instead of '",      &
                size(obs), "'."
            return
        end if

        ! loop over the observations
        do iobs = 1, size(obs)
            call obs(iobs)%init(filename(iobs), status)
            if (status /= 0) return
        end do

        ! print some info
        do iobs = 1, size(obs)
            write (OUTPUT_UNIT,'(a)') 'Info: Section' // strternary(iobs>1, ' '&
                // strinteger(iobs), ''), strsection(obs(iobs)%first,          &
                obs(iobs)%last)
        end do

    end subroutine init_array


    !---------------------------------------------------------------------------


    ! returns 'b', 'g' or 'r' for a blue, green or red channel observation
    subroutine set_channel(this, status)
        class(pacsobservation), intent(inout) :: this
        integer, intent(out)                  :: status
        integer                               :: length, unit, nsamples
        integer                               :: status_close
        character(len=2), allocatable         :: channels(:)

        this%channel = ' '

        ! old style file format
        length = len_trim(this%filename)
        if (this%filename(length-4:length) /= '.fits') then
            status = 0
            write (ERROR_UNIT,'(a)') 'SET_CHANNEL: obsolete file format.'
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
            go to 999
        endif

        call ft_open_bintable(trim(this%filename) // '[Status]', unit,         &
                              nsamples, status)
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

    999 if (status == 0) then
            write (OUTPUT_UNIT,'(a,$)') "Info: Channel: "
            select case (this%channel)
                case ('b')
                    write (OUTPUT_UNIT,'(a)') "'Blue'"
                case ('g')
                    write (OUTPUT_UNIT,'(a)') "'Green'"
                case ('r')
                    write (OUTPUT_UNIT,'(a)') "'Red'"
            end select
        end if

    end subroutine set_channel


    !---------------------------------------------------------------------------


    subroutine set_compression_mode(this, status)
        class(pacsobservation), intent(inout) :: this
        integer, intent(out)                  :: status
        integer           :: unit, ikey, length, status_close
        character(len=72) :: algorithm, keyword, comment

        this%compression_factor = 0

        ! old style file format
        length = len_trim(this%filename)
        if (this%filename(length-4:length) /= '.fits') then
            status = 0
            this%compression_factor = 1
            write (*,'(a)') 'Info: Compression factor is one for obsolete file &
                &format.'
            return
        end if

        call ft_open(trim(this%filename), unit, status)
        if (status /= 0) return

        ikey = 1
        do
            call ftgkys(unit, 'key.META_'//strinteger(ikey), keyword, comment, status)
            if (status /= 0) then
                write (ERROR_UNIT,'(a)') "ERROR: FITS keyword 'algorithm' is &
                    &not found."
                go to 999
            end if
            if (keyword == 'algorithm') exit
            ikey = ikey + 1
        end do

        call ftgkys(unit, 'META_'//strinteger(ikey), algorithm, comment, status)

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
        
        write (OUTPUT_UNIT,'(a,i0)') "Info: Compression factor: ",             &
            this%compression_factor

    999 call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine set_compression_mode


    !---------------------------------------------------------------------------


    subroutine set_transparent_mode(this, status)
        class(pacsobservation), intent(inout) :: this
        integer, intent(out)                  :: status
        integer                               :: unit, ikey, length, status_close
        character(len=72)                     :: compression, keyword, comment

        this%transparent_mode = .false.

        ! old style file format
        length = len_trim(this%filename)
        if (this%filename(length-4:length) /= '.fits') then
            status = 0
            this%transparent_mode = .true.
            write (*,'(a)') 'Info: Transparent mode assumed for obsolete file f&
                &ormat.'
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
        if (compression == 'Photometry Default Mode') then
            this%transparent_mode = .false.
        else if (compression == 'Photometry Lossless Compression Mode') then
            this%transparent_mode = .true.
        else
            status = 1
            write (ERROR_UNIT,'(a)') "ERROR: Unknown compression mode: '" // trim(compression) // "'."
            go to 999
        end if

        write (OUTPUT_UNIT,'(a)') "Info: Compression mode: '" // trim(compression) // "'"

    999 call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine set_transparent_mode


    !---------------------------------------------------------------------------


    subroutine set_valid_slice(this, first, last, status)
        class(pacsobservation), intent(inout) :: this
        integer*8, intent(in)                 :: first, last
        integer, intent(out)                  :: status
        integer*8              :: nsamples
        integer                :: length
        integer*8              :: isample
        integer*8, allocatable :: chop(:)
        integer*4, allocatable :: lbl(:)

        length = len_trim(this%filename)
        if (this%filename(length-4:length) /= '.fits') then
            write(*,'(a)') 'Info: Obsolete file format.'

            call ft_readextension(trim(this%filename) // '_ChopFpuAngle.fits', &
                                  chop, status)
            if (status /= 0) return

            nsamples = size(chop)

            if (nsamples == 0) then
               status = 1
               write (ERROR_UNIT) 'ERROR: Status extension is empty.'
               return
            end if

            if (last > nsamples) then
               status = 1
               write (ERROR_UNIT,'(a,2(i0,a))') "ERROR: The last sample '",    &
                    last, "' exceeds the size of the observation '",        &
                    nsamples, "'."
               return
            end if

            ! set last valid sample
            if (last == 0) then
               this%last = nsamples
            else
               this%last = last
            end if

            ! set first valid sample
            if (first == 0) then
                if (abs(chop(this%last)) > 1.d0) then
                    write (OUTPUT_UNIT,'(a)') 'Info: last sample is invalid. Au&
                        &tomatic search for valid slice is disabled.'
                    this%first = 1
                else
                    do isample = this%last-1, 1, -1
                        if (abs(chop(isample)) > 1.d0) exit
                    end do
                  this%first = isample + 1
                end if
            else
                this%first = first
            end if
            
            return

        end if
  
        ! read status column 'LBL', which discriminates valid samples
        call ft_read_column(trim(this%filename)//'[Status]', 'LBL', lbl, status)
        if (status /= 0) return

        nsamples = size(lbl)

        ! make sure that the slice is inside the observation range
        if (first < 0 .or. last > nsamples) then
            status = 1
            write (ERROR_UNIT,'(a)') "ERROR: Slice '[" //                      &
                 strsection(first, last) // "]' is not inside observation range&
                 & '[" // strsection(1_8, nsamples) // "]'."
            return
        end if
        
        ! if first or last are specified, use them
        if (last /= 0) then
            this%last = last
        else
            this%last = nsamples
        end if
        if (first /= 0) then
            this%first = first
            return
        endif

        ! first = 0, let's determine first valid sample
        ! if we only have one valid sample, keep it
        if (nsamples == 1) then
            this%first = this%last
            return
        end if

        ! if last is already invalid, keep everything
        if (lbl(nsamples) /= 0) then
            write (OUTPUT_UNIT,'(a)') 'Info: last sample is invalid. Automatic &
                &search for valid slice is disabled.'
            this%first = 1
            return
        end if

        ! find last valid sample and exclude it
        do isample = nsamples-1, 1, -1
            if (lbl(isample) /= 0) exit
        end do
        this%first = isample + 1

    end subroutine set_valid_slice


    !---------------------------------------------------------------------------


    subroutine set_filename(this, filename, first, last, status)
        class(pacsobservation), intent(inout) :: this
        character(len=*), intent(in)          :: filename
        integer*8, intent(out)                :: first, last
        integer, intent(out)                  :: status
        integer                               :: pos, length, delim
 
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
                write (ERROR_UNIT,'(a)') 'ERROR: Invalid last sample: ' //    &
                    filename(pos-1:length+1)
                return
            end if
        end if

        ! read first sample
        if (delim > pos) then
            read (filename(pos:delim-1), '(i20)', iostat=status) first
            if (status /= 0) then
                write (ERROR_UNIT,'(a)') 'ERROR: Invalid first sample: ' //    &
                    filename(pos-1:length+1)
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


    !---------------------------------------------------------------------------


    subroutine read_scalar(obs, pq, signal, mask, status)
        type(pacsobservation), intent(inout) :: obs
        integer, intent(in)                  :: pq(:,:)
        real(p), intent(out)                 :: signal(:,:)
        logical(1), intent(out)              :: mask(:,:)
        integer, intent(out)                 :: status


        ! check that the number of samples in the observation is equal
        ! to the number of samples in the signal and mask arrays
        if (obs%nsamples /= size(signal,1) .or.                                &
            obs%nsamples /= size(mask,1)) then
            status = 1
            write (ERROR_UNIT,'(a)') 'READ_TOD: invalid dimensions.'
            return
        end if

        call obs%read(pq, 1_8, signal, mask, status)

    end subroutine read_scalar


    !---------------------------------------------------------------------------


    subroutine read_array(obs, pq, signal, mask, status)
        type(pacsobservation), intent(in) :: obs(:)
        integer, intent(in)               :: pq(:,:)
        real(p), intent(out)              :: signal(:,:)
        logical(1), intent(out)           :: mask(:,:)
        integer, intent(out)              :: status
        integer*8 :: nsamples, destination
        integer   :: nobs, iobs

        nobs     = size(obs)
        nsamples = sum(obs%nsamples)

        ! check that the total number of samples in the observations is equal
        ! to the number of samples in the signal and mask arrays
        if (nsamples /= size(signal,1) .or. nsamples /= size(mask,1)) then
            status = 1
            write (ERROR_UNIT,'(a)') 'READ_TOD: invalid dimensions.'
            return
        end if

        ! loop over the PACS observations
        destination = 1
        do iobs = 1, nobs
            call obs(iobs)%read(pq, destination, signal, mask, status)
            if (status /= 0) return
            destination = destination + obs(iobs)%nsamples
        end do
   
    end subroutine read_array


    !---------------------------------------------------------------------------


    subroutine read(this, pq, destination, signal, mask, status)
        class(pacsobservation), intent(in) :: this
        integer, intent(in)                :: pq(:,:)
        integer*8, intent(in)              :: destination
        real*8, intent(inout)              :: signal(:,:)
        logical*1, intent(inout)           :: mask  (:,:)
        integer, intent(out)               :: status
        integer*8                          :: p, q
        integer                            :: idetector, unit, length
        integer                            :: status_close
        integer, allocatable               :: imageshape(:)
        logical                            :: mask_found
        integer*4, allocatable             :: maskcompressed(:)
        integer*4                          :: maskval
        integer*8                          :: ncompressed
        integer*8                          :: isample, icompressed, ibit
        integer*8                          :: firstcompressed, lastcompressed

        ! old style file format
        length = len_trim(this%filename)
        if (this%filename(length-4:length) /= '.fits') then
            write (ERROR_UNIT,'(a)') 'READ: obsolete file format.'
            call read_oldstyle(this, pq, destination, signal, mask, status)
            return
        end if

        ! read signal HDU
        call ft_open_image(trim(this%filename) // '[Signal]', unit, 3,         &
                           imageshape, status)
        if (status /= 0) return

        do idetector = 1, size(signal,2)

            p = pq(1,idetector)
            q = pq(2,idetector)

            call ft_readslice(unit, this%first, this%last, q+1, p+1,           &
                 imageshape, signal(destination:destination+this%nsamples-1,   &
                 idetector),status)
            if (status /= 0) go to 999

        end do

        call ft_close(unit, status)
        if (status /= 0) return


        ! read Mask HDU
        mask(destination:destination+this%nsamples-1,:) = .false.
        mask_found = ft_test_extension(trim(this%filename)//'[Master]', status)
        if (status /= 0) return

        if (.not. mask_found) then
            write (*,'(a)') 'Info: mask Master is not found.'
            return
        end if

        call ft_open_image(trim(this%filename) // '[Master]', unit, 3,         &
                           imageshape, status)
        if (status /= 0) return

        allocate(maskcompressed(imageshape(1)))

        do idetector = 1, size(mask,2)

            p = pq(1,idetector)
            q = pq(2,idetector)


            firstcompressed = (this%first - 1) / 32 + 1
            lastcompressed  = (this%last  - 1) / 32 + 1
            ncompressed = lastcompressed - firstcompressed + 1

            call ft_readslice(unit, firstcompressed, lastcompressed,       &
                 q+1, p+1, imageshape, maskcompressed(1:ncompressed),status)
            if (status /= 0) go to 999

            ! loop over the bytes of the compressed mask
            do icompressed = firstcompressed, lastcompressed

                maskval = maskcompressed(icompressed-firstcompressed+1)
                if (maskval == 0) cycle

                isample = (icompressed-1)*32 - this%first + destination + 1

                ! loop over the bits of a compressed mask byte
                do ibit = max(0, this%first - (icompressed-1)*32-1),           &
                          min(31, this%last - (icompressed-1)*32-1)
                    mask(isample+ibit,idetector) =                             &
                         mask(isample+ibit,idetector) .or. btest(maskval,ibit)
                end do

            end do

        end do

    999 call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine read


    !---------------------------------------------------------------------------


    subroutine read_oldstyle(this, pq, dest, signal, mask, status)
        class(pacsobservation), intent(in) :: this
        integer, intent(in)                :: pq(:,:)
        integer*8, intent(in)              :: dest
        real*8, intent(inout)              :: signal(:,:)
        logical*1, intent(inout)           :: mask(:,:)
        integer, intent(out)               :: status
        integer*8                          :: p, q
        integer                            :: idetector, unit, status_close
        integer, allocatable               :: imageshape(:)

        ! handle signal
        call ft_open_image(trim(this%filename) // '_Signal.fits', unit, 3,     &
                           imageshape, status)
        if (status /= 0) return

        do idetector = 1, size(signal,2)

            p = pq(1,idetector)
            q = pq(2,idetector)
            call ft_readslice(unit, this%first, this%last, q+1, p+1,           &
                 imageshape, signal(dest:dest+this%nsamples-1,idetector),status)
            if (status /= 0) go to 999

        end do

        call ft_close(unit, status)
        if (status /= 0) return

        ! handle mask
        call ft_open_image(trim(this%filename) // '_Mask.fits', unit, 3,       &
                           imageshape, status)
        if (status /= 0) return

        do idetector = 1, size(mask,2)

            p = pq(1,idetector)
            q = pq(2,idetector)
            call ft_readslice(unit, this%first, this%last, q+1, p+1,           &
                 imageshape, mask(dest:dest+this%nsamples-1,idetector), status)
            if (status /= 0) go to 999

        end do

    999 call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine read_oldstyle


 end module module_pacsobservation
