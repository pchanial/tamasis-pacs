module module_pacsinstrument

    use, intrinsic :: ISO_FORTRAN_ENV, only : ERROR_UNIT, OUTPUT_UNIT
    use module_fitstools
    use module_math, only : nint_down, nint_up
    use module_pacspointing
    use module_pointingmatrix
    use module_projection
    use string, only : strinteger, strlowcase
    implicit none
    private

    public :: ndims
    public :: nvertices
    public :: pacsinstrument
    public :: xy2roi
    public :: roi2pmatrix
    public :: multiplexing_direct
    public :: multiplexing_transpose

    integer, parameter :: ndims = 2
    integer, parameter :: nvertices = 4
    integer, parameter :: shape_blue(2) = [32, 64]
    integer, parameter :: shape_red (2) = [16, 32]
    integer, parameter :: distortion_degree = 3
    real*8, parameter  :: sampling = 0.025d0

    integer :: tamasis_dir_len ! length of the path to the tamasis installation directory

    type pacsinstrument
 
        logical*1            :: mask_blue(shape_blue(1), shape_blue(2))
        logical*1            :: mask_red (shape_red (1), shape_red (2))

        real*8               :: flatfield_blue (shape_blue(1), shape_blue(2))
        real*8               :: flatfield_green(shape_blue(1), shape_blue(2))
        real*8               :: flatfield_red  (shape_red (1), shape_red (2))

        integer              :: ndetectors
        integer              :: nrows
        integer              :: ncolumns
        integer              :: fine_sampling_factor
        character            :: channel
        logical              :: transparent_mode

        logical*1, pointer   :: mask(:,:)
        integer, allocatable :: ij(:,:)
        integer, allocatable :: pq(:,:)

        real*8, allocatable  :: flatfield(:)

        real*8               :: corners_uv_blue(ndims, nvertices, shape_blue(1), shape_blue(2))
        real*8               :: corners_uv_red (ndims, nvertices, shape_red (1), shape_red (2))
        real*8, allocatable  :: corners_uv(:,:)

        real*8               :: corners_yz_blue(ndims, nvertices, shape_blue(1), shape_blue(2))
        real*8               :: corners_yz_red (ndims, nvertices, shape_red (1), shape_red (2))
        real*8, allocatable  :: corners_yz(:,:)

        real*8               :: distortion_yz_blue(ndims, distortion_degree, distortion_degree, distortion_degree)
        real*8               :: distortion_yz_red (ndims, distortion_degree, distortion_degree, distortion_degree)
        real*8               :: distortion_yz     (ndims, distortion_degree, distortion_degree, distortion_degree)

    contains

        procedure          :: init
        procedure          :: init_filename
        procedure          :: read_calibration_files
        procedure          :: filter_detectors
        procedure          :: read_tod_file
        procedure          :: compute_mapheader
        procedure          :: find_minmax
        procedure          :: compute_projection_sharp_edges
        procedure, nopass  :: uv2yz
        procedure, nopass  :: yz2ad
        procedure, nopass  :: xy2roi
        procedure, nopass  :: roi2pmatrix

        procedure, private :: filter_detectors_array
        procedure, private :: read_signal_file_oldstyle
        procedure, private :: read_mask_file_oldstyle

    end type pacsinstrument


contains

    
    subroutine init(this, channel, transparent_mode, keep_bad_detectors,       &
                    status, bad_detector_mask)
        class(pacsinstrument), intent(inout) :: this
        character, intent(in)                :: channel
        logical, intent(in)                  :: transparent_mode
        logical, intent(in)                  :: keep_bad_detectors
        integer, intent(out)                 :: status
        logical*1, intent(in), optional      :: bad_detector_mask(:,:)

        this%channel = channel
        this%transparent_mode = transparent_mode

        call this%read_calibration_files(status)
        if (status /= 0) return

        if (present(bad_detector_mask)) then
            if (channel == 'r') then
                !XXX gfortran bug test_bug_shape
                ! should be any(shape(bad_detector_mask) /= shape(this%mask_red))
                if (size(bad_detector_mask,1) /= size(this%mask_red,1) .or.    &
                    size(bad_detector_mask,2) /= size(this%mask_red,2)) then
                   status = 1
                   write (ERROR_UNIT,'(a)') 'INIT: input red mask has invalid size.'
                   return
                end if
                this%mask_red = bad_detector_mask
            else
                if (channel == 'g') write (*,*) 'pacs%init: check me.'
                if (size(bad_detector_mask,1) /= size(this%mask_blue,1) .or.   &
                    size(bad_detector_mask,2) /= size(this%mask_blue,2)) then
                   status = 1
                   write (ERROR_UNIT,'(a)') 'INIT: input blue mask has invalid size.'
                   return
                end if
                this%mask_blue = bad_detector_mask
            end if

        end if
        
        call this%filter_detectors(channel, transparent_mode,                  &
                                   keep_bad_detectors, status)

    end subroutine init


    !---------------------------------------------------------------------------


    subroutine init_filename(this, filename, keep_bad_detectors, status,       &
                            bad_detector_mask)
        class(pacsinstrument), intent(inout) :: this
        character(len=*), intent(in)         :: filename
        logical, intent(in)                  :: keep_bad_detectors
        integer, intent(out)                 :: status
        logical*1, intent(in), optional      :: bad_detector_mask(:,:)
        character :: channel
        logical   :: transparent_mode

        channel = get_channel(filename, status)
        if (status /= 0) return

        transparent_mode = get_transparent_mode(filename, status)
        if (status /= 0) return

        call this%init(channel, transparent_mode, keep_bad_detectors, status,  &
                       bad_detector_mask)

    end subroutine init_filename


    !---------------------------------------------------------------------------


    subroutine read_calibration_files(this, status)

        class(pacsinstrument), intent(inout) :: this
        integer, intent(out)                 :: status

        character(len=*), parameter :: filename_saa = 'PCalPhotometer_SubArrayArray_FM_v5.fits'
        character(len=*), parameter :: filename_ai  = 'PCalPhotometer_ArrayInstrument_FM_v4.fits'
        character(len=*), parameter :: filename_bpm = 'PCalPhotometer_BadPixelMask_FM_v3.fits'
        character(len=*), parameter :: filename_ff  = 'PCalPhotometer_FlatField_FM_v1.fits'

        integer,          parameter :: hdu_blue(4) = [8, 12, 16, 20]
        integer,          parameter :: hdu_red (4) = [6, 10, 14, 18]

        integer                     :: ivertex
        logical*1, allocatable      :: tmplogical(:,:)
        real*8, allocatable         :: tmp2(:,:)
        real*8, allocatable         :: tmp3(:,:,:)
        character(len=100)          :: value

        ! might be moved elsewhere
        call get_environment_variable('TAMASIS_DIR', value, length=tamasis_dir_len)
        if (tamasis_dir_len == 0) tamasis_dir_len = len('/home/pchanial/work/tamasis')


        ! read bad pixel mask

        call ft_readextension(get_calfile(filename_bpm) // '[blue]', tmplogical, status)
        if (status /= 0) return
        this%mask_blue = transpose(tmplogical)

        call ft_readextension(get_calfile(filename_bpm) // '[red]', tmplogical, status)
        if (status /= 0) return
        this%mask_red = transpose(tmplogical)


        ! read flat fields

        call ft_readextension(get_calfile(filename_ff) // '+5',  tmp2, status)
        if (status /= 0) return
        this%flatfield_blue = transpose(tmp2)
        call ft_readextension(get_calfile(filename_ff) // '+8', tmp2, status)
        if (status /= 0) return
        this%flatfield_green = transpose(tmp2)
        call ft_readextension(get_calfile(filename_ff) // '+2',   tmp2, status)
        if (status /= 0) return
        this%flatfield_red = transpose(tmp2)


        ! read detector corners in the (u,v) plane

        do ivertex=1, nvertices

            call ft_readextension(get_calfile(filename_saa), tmp2, status, hdu_blue(ivertex)  )
            if (status /= 0) return
            this%corners_uv_blue(1,ivertex,:,:) = transpose(tmp2)
            call ft_readextension(get_calfile(filename_saa), tmp2, status, hdu_blue(ivertex)+1)
            if (status /= 0) return
            this%corners_uv_blue(2,ivertex,:,:) = transpose(tmp2)
            call ft_readextension(get_calfile(filename_saa), tmp2, status, hdu_red (ivertex)  )
            if (status /= 0) return
            this%corners_uv_red (1,ivertex,:,:) = transpose(tmp2)
            call ft_readextension(get_calfile(filename_saa), tmp2, status, hdu_red (ivertex)+1)
            if (status /= 0) return
            this%corners_uv_red (2,ivertex,:,:) = transpose(tmp2)

        end do


        ! read the distortion coefficients in the (y,z) plane

        call ft_readextension(get_calfile(filename_ai) // '[ycoeffblue]', tmp3, status)
        if (status /= 0) return
        this%distortion_yz_blue(1,:,:,:) = tmp3
        call ft_readextension(get_calfile(filename_ai) // '[zcoeffblue]', tmp3, status)
        if (status /= 0) return
        this%distortion_yz_blue(2,:,:,:) = tmp3
        call ft_readextension(get_calfile(filename_ai) // '[ycoeffred]', tmp3, status)
        if (status /= 0) return
        this%distortion_yz_red (1,:,:,:) = tmp3
        call ft_readextension(get_calfile(filename_ai) // '[zcoeffred]', tmp3, status)
        if (status /= 0) return
        this%distortion_yz_red (2,:,:,:) = tmp3

    end subroutine read_calibration_files


    !---------------------------------------------------------------------------


    function get_calfile(filename)
        character(len=*), intent(in)                   :: filename
        character(len=tamasis_dir_len+6+len(filename)) :: get_calfile
        integer                                        :: length

        ! Environment variables are not visible in photran...
        call get_environment_variable('TAMASIS_DIR', get_calfile, length)
        if (length == 0) get_calfile(1:tamasis_dir_len) = '/home/pchanial/work/tamasis'
        get_calfile(tamasis_dir_len+1:) = '/data/' // filename

    end function get_calfile


    !---------------------------------------------------------------------------


    subroutine read_signal_file_oldstyle(this, filename, first, last, signal, status)
        class(pacsinstrument), intent(in) :: this
        character(len=*), intent(in)      :: filename
        integer*8, intent(in)             :: first, last
        real*8, intent(out)               :: signal(:,:)
        integer, intent(out)              :: status
        integer*8                         :: p, q
        integer                           :: idetector, unit
        integer, allocatable              :: imageshape(:)
        logical                           :: keep_bad_detectors

        keep_bad_detectors = size(this%mask) > this%ndetectors
        status = 0
        call ft_open_image(filename // '_Signal.fits', unit, 3, imageshape, status)
        if (status /= 0) return

        do idetector = 1, this%ndetectors

            p = this%pq(1,idetector)
            q = this%pq(2,idetector)
            call ft_readslice(unit, first, last, q+1, p+1, imageshape, signal(:,idetector), status)
            if (status /= 0) return

        end do

        call ft_close(unit, status)

    end subroutine read_signal_file_oldstyle


    !---------------------------------------------------------------------------


    subroutine read_mask_file_oldstyle(this, filename, first, last, mask, status)

        class(pacsinstrument), intent(in) :: this
        character(len=*), intent(in)      :: filename
        integer*8, intent(in)             :: first, last
        logical*1, intent(out)            :: mask(:,:)
        integer, intent(out)              :: status
        integer*8                         :: p, q
        integer                           :: idetector, unit
        integer, allocatable              :: imageshape(:)

        call ft_open_image(filename // '_Mask.fits', unit, 3, imageshape, status)
        if (status /= 0) return

        do idetector = 1, this%ndetectors

            p = this%pq(1,idetector)
            q = this%pq(2,idetector)
            call ft_readslice(unit, first, last, q+1, p+1, imageshape, mask(:,idetector), status)
            if (status /= 0) return

        end do

        call ft_close(unit, status)

    end subroutine read_mask_file_oldstyle


    !---------------------------------------------------------------------------


    subroutine read_tod_file(this, filename, first, last, signal, mask, status)
        class(pacsinstrument), intent(in) :: this
        character(len=*), intent(in)      :: filename
        integer*8, intent(in)             :: first, last
        real*8, intent(out)               :: signal(:,:)
        logical*1, intent(out)            :: mask  (:,:)
        integer, intent(out)              :: status
        integer*8                         :: p, q
        integer                           :: idetector, unit, pos
        integer, allocatable              :: imageshape(:)
        logical                           :: keep_bad_detectors, mask_found
        integer*4, allocatable            :: maskcompressed(:)
        integer*4                         :: maskval
        integer                           :: ntimes, ncompressed
        integer                           :: itime, icompressed, ibit
        integer*8                         :: firstcompressed, lastcompressed

        if (last < first) then
            status = 1
            write (ERROR_UNIT,'(a,2(i0,a))')"READ_TOD_FILE: The last sample '",&
                  last, "' is less than the first sample '", first, "'."
            return
        end if

        if (first < 1) then
            status = 1
            write (ERROR_UNIT,'(a,i0,a)') "READ_TOD_FILE: The first sample '", &
                  first, "' is less than 1."
            return
        end if
 
        pos = len_trim(filename)
        if (filename(pos-4:pos) /= '.fits') then
            write (ERROR_UNIT,'(a)') 'READ_TOD_FILE: obsolete file format.'
            call this%read_signal_file_oldstyle(filename, first, last, signal, status)
            if (status /= 0) return
            call this%read_mask_file_oldstyle(filename, first, last, mask, status)
            return
        endif

        keep_bad_detectors = size(this%mask) > this%ndetectors

        ! read signal HDU
        call ft_open_image(filename // '[Signal]', unit, 3, imageshape, status)
        if (status /= 0) return

        if (last > imageshape(1)) then
            status = 1
            write (ERROR_UNIT,'(a,2(i0,a))')"READ_TOD_FILE: The last sample '",&
                  last, "' exceeds the size of the observation '",             &
                  imageshape(1), "'."
            return
        end if
        ntimes = last - first + 1

        do idetector = 1, this%ndetectors

            p = this%pq(1,idetector)
            q = this%pq(2,idetector)
            call ft_readslice(unit, first, last, q+1, p+1, imageshape, signal(:,idetector), status)
            if (status /= 0) return

        end do

        call ft_close(unit, status)
        if (status /= 0) return

        ! read MMT_Glitchmask HDU
        mask = .false.
        mask_found = ft_test_extension(filename // '[MMT_Glitchmask]', status)
        if (status /= 0) return

        if (.not. mask_found) then
            write (*,'(a)') 'Info: mask MMT_Glitchmask is not found.'
            return
        end if

        call ft_open_image(filename // '[MMT_Glitchmask]', unit, 3, imageshape, status)
        if (status /= 0) return

        firstcompressed = (first - 1) / 32 + 1
        lastcompressed  = (last  - 1) / 32 + 1
        ncompressed = lastcompressed - firstcompressed + 1
        allocate(maskcompressed(ncompressed))

        do idetector = 1, this%ndetectors

            p = this%pq(1,idetector)
            q = this%pq(2,idetector)
            call ft_readslice(unit, firstcompressed, lastcompressed, q+1, p+1, &
                              imageshape, maskcompressed, status)
            if (status /= 0) go to 999

            ! loop over the bytes of the compressed mask
            do icompressed = firstcompressed, lastcompressed

               maskval = maskcompressed(icompressed-firstcompressed+1)
               if (maskval == 0) cycle

               itime = (icompressed-1)*32 - first + 2

               ! loop over the bits of a compressed mask byte
               do ibit= max(0, first - (icompressed-1)*32-1),                  &
                        min(31, last - (icompressed-1)*32-1)
                    mask(itime+ibit,idetector) = mask(itime+ibit,idetector)    &
                                                 .or. btest(maskval, ibit)
               end do

            end do

        end do

    999 call ft_close(unit, status)

    end subroutine read_tod_file


    !---------------------------------------------------------------------------


    subroutine filter_detectors(this, channel, transparent_mode, keep_bad_detectors, status)
        class(pacsinstrument), intent(inout)   :: this
        character, intent(in)                  :: channel
        logical, intent(in), optional          :: transparent_mode
        logical, intent(in), optional          :: keep_bad_detectors
        integer, intent(out)                   :: status

        select case (channel)
           case ('r')
              call this%filter_detectors_array(this%mask_red,                  &
                   this%corners_uv_red, this%distortion_yz_red,                &
                   this%flatfield_red, transparent_mode=transparent_mode,      &
                   keep_bad_detectors=keep_bad_detectors)
           case ('g')
              write(*,*) 'FILTER_DETECTORS: Check calibration files for the green band...'
              call this%filter_detectors_array(this%mask_blue,                 &
                   this%corners_uv_blue, this%distortion_yz_blue,              &
                   this%flatfield_green, transparent_mode=transparent_mode,    &
                   keep_bad_detectors=keep_bad_detectors)
           case ('b')
              call filter_detectors_array(this, this%mask_blue,                &
                   this%corners_uv_blue, this%distortion_yz_blue,              &
                   this%flatfield_blue, transparent_mode=transparent_mode,     &
                   keep_bad_detectors=keep_bad_detectors)

           case default
                status = 1
                write (ERROR_UNIT,'(a)') "FILTER_DETECTORS: invalid array &
                    &channel ('blue', 'green' or 'red'): " // channel
                return
        end select

        status = 0

    end subroutine filter_detectors


    !---------------------------------------------------------------------------


    subroutine filter_detectors_array(this, mask, uv, distortion, flatfield, transparent_mode, keep_bad_detectors)

        class(pacsinstrument), intent(inout)    :: this
        logical*1, intent(in), target :: mask(:,:)
        real*8, intent(in)            :: uv(:,:,:,:)
        real*8, intent(in)            :: distortion(:,:,:,:)
        real*8, intent(in)            :: flatfield(:,:)
        logical, intent(in), optional :: transparent_mode
        logical, intent(in), optional :: keep_bad_detectors

        integer                       :: idetector, p, q
        logical                       :: transmode, keepbaddetectors

        if (present(transparent_mode)) then
            transmode = transparent_mode
        else
            transmode = .false.
        end if

        if (present(keep_bad_detectors)) then
            keepbaddetectors = keep_bad_detectors
        else
            keepbaddetectors = .false.
        end if

        this%nrows    = size(mask, 1)
        this%ncolumns = size(mask, 2)

        allocate(this%mask(this%nrows, this%ncolumns))
        this%mask => mask

        if (transmode) then
            this%mask(1:16,1:16) = .true.
            this%mask(1:16,33:)  = .true.
            this%mask(17:,:)     = .true.
        end if
        this%ndetectors = size(this%mask)
        if (.not. keepbaddetectors) then
            this%ndetectors = this%ndetectors - count(this%mask)
        end if

        allocate(this%ij(ndims, this%ndetectors))
        allocate(this%pq(ndims, this%ndetectors))
        allocate(this%corners_uv(ndims, nvertices * this%ndetectors))
        allocate(this%corners_yz(ndims, nvertices * this%ndetectors))
        allocate(this%flatfield(this%ndetectors))
        this%distortion_yz = distortion

        idetector = 1

        do p = 0, this%nrows - 1
            do q = 0, this%ncolumns - 1
                if (this%mask(p+1,q+1) .and. .not. keepbaddetectors) cycle
                this%pq(1, idetector) = p
                this%pq(2, idetector) = q
                this%ij(1, idetector) = mod(p, 16)
                this%ij(2, idetector) = mod(q, 16)
                this%corners_uv(:,nvertices * (idetector-1)+1:nvertices*idetector) = uv(:,:,p+1,q+1)
                this%flatfield(idetector) = flatfield(p+1,q+1)
                idetector = idetector + 1
            end do
        end do

    end subroutine filter_detectors_array


    !---------------------------------------------------------------------------


    ! returns 'b', 'g' or 'r' for a blue, green or red channel observation
    function get_channel(filename, status) result(channel)
        character                     :: channel
        character(len=*), intent(in)  :: filename
        integer, intent(out)          :: status
        integer                       :: length, unit, ncolumns, nrecords
        integer                       :: status_close
        character(len=2), allocatable :: channels(:)

        channel = ' '
        length = len_trim(filename)
        if (filename(length-4:length) /= '.fits') then
            status = 0
            write (ERROR_UNIT,'(a)') 'GET_CHANNEL: obsolete file format.'
            if (strlowcase(filename(length-3:length)) == 'blue') then
               channel = 'b'
            else if (strlowcase(filename(length-4:length)) == 'green') then
               channel = 'g'
            else if (strlowcase(filename(length-2:length)) == 'red') then
               channel = 'r'
            else
               status = 1
               write (ERROR_UNIT,'(a)') 'File name does not contain the array channel identifier (blue, green, red).'
            end if
            go to 999
        endif

        call ft_open_bintable(filename // '[Status]', unit, ncolumns, nrecords, status)
        if (status /= 0) return

        allocate(channels(nrecords))
        call ft_read_column(unit, 'BAND', channels, status)
        if (status /= 0) return

        if (any(channels /= channels(1))) then
            status = 1
            write (ERROR_UNIT,'(a)') 'Observation is BANDSWITCH.'
        else
           select case (channels(1))
               case ('BS')
                   channel = 'b'
               case ('BL')
                   channel = 'g'
               case ('R ')
                   channel = 'r'
               case default
                   status = 1
                   write (ERROR_UNIT,'(a)') 'Invalid array BAND value: ' // channels(1)
           end select
        end if

        call ft_close(unit, status_close)
        if (status == 0) status = status_close

    999 if (status == 0) then
            write (OUTPUT_UNIT,'(a,$)') "Info: Channel: '"
            select case (channel)
                case ('b')
                    write (OUTPUT_UNIT,'(a,$)') 'Blue'
                case ('g')
                    write (OUTPUT_UNIT,'(a,$)') 'Green'
                case ('r')
                    write (OUTPUT_UNIT,'(a,$)') 'Red'
            end select
            write (OUTPUT_UNIT,'(a)') "'"
        end if

    end function get_channel


    !---------------------------------------------------------------------------


    function get_transparent_mode(filename, status) result(tmode)
        logical                      :: tmode
        character(len=*), intent(in) :: filename
        integer, intent(out)         :: status
        integer           :: unit, ikey, length
        character(len=72) :: compression, keyword, comment

        tmode = .true.

        length = len_trim(filename)
        if (filename(length-4:length) /= '.fits') then
            write(*,'(a)') 'GET_TRANSPARENT_MODE: Obsolete file format.'
            status = 0
            return
        end if

        call ft_open(filename, unit, status)
        if (status /= 0) return

        ikey = 1
        do
            call ftgkys(unit, 'key.META_'//strinteger(ikey), keyword, comment, status)
            if (status /= 0) go to 999
            if (keyword == 'compMode') exit
            ikey = ikey + 1
        end do

        call ftgkys(unit, 'META_'//strinteger(ikey), compression, comment, status)
        if (compression == 'Photometry Default Mode') then
            tmode = .false.
        else if (compression == 'Photometry Lossless Compression Mode') then
            tmode = .true.
        else
            status = 1
            write (ERROR_UNIT,'(a)') "ERROR: Unknown compression mode: '" // trim(compression) // "'."
            return
        end if

        write (OUTPUT_UNIT,'(a)') "Info: Compression mode: '" // trim(compression) // "'"

    999 call ft_close(unit, status)

    end function get_transparent_mode


    !---------------------------------------------------------------------------


    function uv2yz(uv, distortion_yz, chop) result(yz)

        real*8, intent(in) :: uv(:,:)
        real*8, intent(in) :: distortion_yz(ndims, distortion_degree, distortion_degree, distortion_degree)
        real*8, intent(in) :: chop ! chop angle in degree
        real*8             :: yz(ndims, size(uv,2))

        real*8             :: a_pow, ratio
        integer            :: n, i, j, k, l
        real*8             :: u_pow(size(uv,2)), v_pow(size(uv,2))

        yz = 0
        n = size(uv,2)
        do k = 1, distortion_degree
            if (k /= 1) then
                u_pow = u_pow * uv(1,:)
            else
                u_pow = 1.d0
            endif
            do j = 1, distortion_degree
              if (j /= 1) then
                  v_pow = v_pow * uv(2,:)
              else
                  v_pow = 1.d0
              end if
              do i = 1, distortion_degree
                  if (i /= 1) then
                      a_pow = a_pow * chop
                  else
                      a_pow = 1.d0
                  end if
                  do l = 1, n
                      ratio = - u_pow(l) * v_pow(l) * a_pow
                      yz(1,l) = yz(1,l) + distortion_yz(1, i, j, k) * ratio
                      yz(2,l) = yz(2,l) + distortion_yz(2, i, j, k) * ratio
                  end do
              end do
           end do
        end do
    
        ! convert from arcsecond to degrees
        yz = - yz / 3600.d0

    end function uv2yz
 

    !---------------------------------------------------------------------------


    function yz2ad(yz, ra0, dec0, pa0) result (ad)

        real*8, intent(in) :: yz(:,:) ! in degrees
        real*8, intent(in) :: ra0, dec0, pa0
        real*8             :: ad(ndims,size(yz,2))

        integer            :: i
        real*8             :: cospa, sinpa
        real*8, parameter  :: pi = 4.0d0 * atan(1.0d0)
        real*8, parameter  :: radeg = pi / 180d0


        cospa =  cos(pa0*radeg)
        sinpa = -sin(pa0*radeg)

        do i=1, size(yz, 2)

            ad(2,i) = dec0 + (yz(1,i) * sinpa + yz(2,i) * cospa)
            ad(1,i) = ra0  + (yz(1,i) * cospa - yz(2,i) * sinpa) / cos(ad(2,i) * radeg)

        end do

    end function yz2ad


    !---------------------------------------------------------------------------


    function xy2roi(xy) result(roi)
        real*8, intent(in) :: xy(:,:)
        real*8             :: roi(ndims,2,size(xy,2)/nvertices)
        integer            :: idetector

        do idetector = 1, size(xy,2) / nvertices
            roi(:,1,idetector) = nint_up  (minval(xy(:,nvertices * (idetector-1)+1:nvertices*idetector),2))
            roi(:,2,idetector) = nint_down(maxval(xy(:,nvertices * (idetector-1)+1:nvertices*idetector),2))
        end do

    end function xy2roi


    !---------------------------------------------------------------------------


    subroutine roi2pmatrix(roi, coords, nx, ny, isample, nroi, pmatrix)
        integer, intent(in)                  :: roi(:,:,:)
        real*8, intent(in)                   :: coords(:,:)
        type(pointingelement), intent(inout) :: pmatrix(:,:,:)
        integer, intent(in)                  :: nx, ny, isample
        integer, intent(out)                 :: nroi
        real*8                               :: polygon(2,nvertices)

        integer                              :: npixels_per_sample, idetector, ix, iy, iroi, ipixel
        real*4                               :: weight

        ipixel = 0
        nroi = 0
        npixels_per_sample = size(pmatrix, 1)
        do idetector = 1, size(pmatrix,3)
            iroi = 1
if (roi(2,1,idetector) < 1 .or. roi(2,2,idetector) > ny) then
    write(*,*) 'roi2pmatrix: map y too small', roi(2,:,idetector), ny
end if
            do iy = max(roi(2,1,idetector),1), min(roi(2,2,idetector),ny)
if (roi(1,1,idetector) < 1 .or. roi(1,2,idetector) > nx) then
    write(*,*) 'roi2pmatrix: map x too small', roi(1,:,idetector), nx
end if
                do ix = max(roi(1,1,idetector),1), min(roi(1,2,idetector),nx)
                    ipixel = ix - 1 + (iy - 1) * nx
                    polygon(1,:) = coords(1,(idetector-1)*nvertices+1:idetector*nvertices) - (ix-0.5d0)
                    polygon(2,:) = coords(2,(idetector-1)*nvertices+1:idetector*nvertices) - (iy-0.5d0)
                    weight = abs(intersection_polygon_unity_square(polygon, nvertices))
                    if (weight <= 0) cycle
                    if (iroi <= npixels_per_sample) then
                        pmatrix(iroi,isample,idetector)%pixel  = ipixel
                        pmatrix(iroi,isample,idetector)%weight = weight
                    end if
                    iroi = iroi + 1
                end do
            end do
            ! fill the rest of the pointing matrix
            pmatrix(iroi:,isample,idetector)%pixel  = ipixel
            pmatrix(iroi:,isample,idetector)%weight = 0
            nroi = max(nroi, iroi-1)
        end do

    end subroutine roi2pmatrix


    !---------------------------------------------------------------------------


    subroutine compute_mapheader(this, pointing, times, resolution, header, status)

        use module_fitstools, only : ft_create_header
        use module_wcs, only : init_astrometry

        class(pacsinstrument), intent(in) :: this
        class(pacspointing), intent(in)   :: pointing
        real*8, intent(in)                :: times(:)  ! times at which the PACS array should be considered as inside the map
        real*8, intent(in)                :: resolution
        character(len=2880), intent(out)  :: header
        integer, intent(out)              :: status

        integer                           :: nx, ny
        integer                           :: ixmin, ixmax, iymin, iymax
        real*8                            :: ra0, dec0, xmin, xmax, ymin, ymax

        call pointing%compute_center(times, ra0, dec0)

        call ft_create_header(0, 0, -resolution/3600.d0, resolution/3600.d0, 0.d0, ra0, dec0, 1.d0, 1.d0, header)

        call init_astrometry(header, status=status)
        if (status /= 0) return

        call this%find_minmax(pointing, times, xmin, xmax, ymin, ymax)

        ixmin = nint(xmin)
        ixmax = nint(xmax)
        iymin = nint(ymin)
        iymax = nint(ymax)

        nx = ixmax-ixmin+1
        ny = iymax-iymin+1

        ! move the reference pixels (not the reference values!)
        call ft_create_header(nx, ny, -resolution/3600.d0, resolution/3600.d0, 0.d0, &
                              ra0, dec0, -ixmin+2.d0, -iymin+2.d0, header)

    end subroutine compute_mapheader


    !---------------------------------------------------------------------------


    ! find minimum and maximum pixel coordinates in maps
    subroutine find_minmax(this, pointing, times, xmin, xmax, ymin, ymax)

        use module_projection, only: convex_hull
        use module_wcs, only : ad2xy_gnomonic
        class(pacsinstrument), intent(in) :: this
        class(pacspointing), intent(in)   :: pointing
        real*8, intent(in)                :: times(:)  ! times at which the PACS array should be considered as inside the map
        real*8, intent(out)               :: xmin, xmax, ymin, ymax

        real*8                            :: ra, dec, pa, chop
        real*8                            :: infinity, zero
        real*8, allocatable               :: hull_uv(:,:), hull(:,:)
        integer, allocatable              :: ihull(:)
        integer                           :: nsamples, isample, index

        zero = 0.d0
        infinity = 1.d0 / zero

        xmin =  infinity
        xmax = -infinity
        ymin =  infinity
        ymax = -infinity

        nsamples = size(times)
        call convex_hull(this%corners_uv, ihull)

        allocate(hull_uv(ndims,size(ihull)))
        allocate(hull(ndims, size(ihull)))
        hull_uv = this%corners_uv(:,ihull)
        index = 2

        ! with gfortran 4.5, the following openmp construct corrupts the memory...
        !!$omp parallel do default(shared) reduction(min:xmin,ymin) reduction(max:xmax,ymax) &
        !!$omp private(isample,ra,dec,pa,chop,hull) firstprivate(index)
        do isample = 1, nsamples
           call pointing%get_position(times(isample), ra, dec, pa, chop, index)
           hull = uv2yz(hull_uv, this%distortion_yz_blue, chop)
           hull = yz2ad(hull, ra, dec, pa)
           hull = ad2xy_gnomonic(hull)
           xmin = min(xmin, minval(hull(1,:)))
           xmax = max(xmax, maxval(hull(1,:)))
           ymin = min(ymin, minval(hull(2,:)))
           ymax = max(ymax, maxval(hull(2,:)))
        end do
        !!$omp end parallel do

    end subroutine find_minmax


    !---------------------------------------------------------------------------


    subroutine compute_projection_sharp_edges(this, pointing, time, header, nx, ny, pmatrix, status)

        use module_wcs, only : init_astrometry, ad2xy_gnomonic
        class(pacsinstrument), intent(in)    :: this
        class(pacspointing), intent(in)      :: pointing
        character(len=*), intent(in)         :: header
        real*8, intent(in)                   :: time(:)
        integer, intent(in)                  :: nx, ny
        type(pointingelement), intent(out)   :: pmatrix(:,:,:)
        integer, intent(out)                 :: status

        real*8  :: coords(ndims,this%ndetectors*nvertices), coords_yz(ndims,this%ndetectors*nvertices)
        real*8  :: ra, dec, pa, chop, chop_old
        integer :: roi(ndims,2,this%ndetectors), isample, nsamples, npixels_per_sample, nroi, index

        call init_astrometry(header, status=status)
        if (status /= 0) return

        nsamples = size(time)
        npixels_per_sample = -1

        chop_old = -9999999999.0d0 ! XXX should be NaN
        index = 2
        !$omp parallel do default(shared) firstprivate(index, chop_old)   &
        !$omp private(isample, ra, dec, pa, chop, coords, coords_yz, roi) &
        !$omp reduction(max : npixels_per_sample)
        do isample = 1, nsamples
            call pointing%get_position(time(isample), ra, dec, pa, chop, index)
            if (abs(chop-chop_old) > 1.d-2) then
                coords_yz = this%uv2yz(this%corners_uv, this%distortion_yz, chop)
                chop_old = chop
            end if
            coords = this%yz2ad(coords_yz, ra, dec, pa)
            coords = ad2xy_gnomonic(coords)
            roi    = xy2roi(coords) ! [1=x|2=y,1=min|2=max,idetector]
            call roi2pmatrix(roi, coords, nx, ny, isample, nroi, pmatrix)
            npixels_per_sample = max(npixels_per_sample, nroi)
        end do
        !$omp end parallel do
        if (npixels_per_sample /= size(pmatrix,1)) then
            write(*,'(a,i0,a)') 'Warning: to compute the Pointing Matrix, npixels_per_sample may be updated to ', npixels_per_sample, '.'
        end if

    end subroutine compute_projection_sharp_edges


    !---------------------------------------------------------------------------


    subroutine multiplexing_direct(signal, sampled_signal, sampling, ij)
        real*8, intent(in)  :: signal(:,:)
        integer, intent(in) :: sampling
        real*8, intent(out) :: sampled_signal(size(signal,1)/sampling,size(signal,2))
        integer, intent(in) :: ij(2,size(signal,2))
        integer             :: ndetectors, idetector, isample, j
        real*8              :: frac

        ndetectors = size(signal, 2)
        if (sampling == 16) then
            !$omp parallel do default(shared) private(idetector)
            do idetector = 1, ndetectors
                sampled_signal(:,idetector) = signal(ij(2,idetector)+1::16,idetector)
            end do
            !$omp end parallel do
        else
            ! for a fine sampling of 4
            ! sampling 16*40Hz: 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16
            !                j:            |          7|           |
            !          isample: 1          |2          |3          |4
            !$omp parallel do default(shared) private(idetector,isample,frac,j)
            do idetector = 1, ndetectors
                j = ij(2,idetector)
                isample = j * sampling / 16 + 1
                frac = (j * sampling / 16.d0) - floor(j * sampling / 16.d0)
                sampled_signal(:,idetector) = (1-frac) * signal(isample  ::sampling,idetector) + &
                                                 frac  * signal(isample+1::sampling,idetector)
            end do
            !$omp end parallel do
        end if

    end subroutine multiplexing_direct


    !---------------------------------------------------------------------------


    subroutine multiplexing_transpose(sampled_signal, signal, sampling, ij)
        real*8, intent(in)  :: sampled_signal(:,:)
        integer, intent(in) :: sampling
        real*8, intent(out) :: signal(size(sampled_signal,1)*sampling,size(sampled_signal,2))
        integer, intent(in) :: ij(2,size(sampled_signal,2))
        integer             :: ndetectors, isample, idetector, j
        real*8              :: frac

        ndetectors = size(sampled_signal, 2)
        if (sampling == 16) then
            signal = 0.d0
            !$omp parallel do default(shared) private(idetector)
            do idetector = 1, ndetectors
                signal(ij(2,idetector)+1::16,idetector) = sampled_signal(:,idetector)
            end do
            !$omp end parallel do
        else
            !$omp parallel default(shared) private(idetector,isample,frac,j)
            !$omp workshare
            signal = 0.d0
            !$omp end workshare
            !$omp do
            do idetector = 1, ndetectors
                j = ij(2,idetector)
                isample = j * sampling / 16 + 1 ! isample+1 is guaranteed to exist
                frac = (j * sampling / 16.d0) - floor(j * sampling / 16.d0)
                signal(isample  ::sampling,idetector) = (1-frac) * sampled_signal(isample::sampling,idetector)
                signal(isample+1::sampling,idetector) =    frac  * sampled_signal(isample::sampling,idetector)
            end do
            !$omp end do
            !$omp end parallel
        end if

    end subroutine multiplexing_transpose


end module module_pacsinstrument
