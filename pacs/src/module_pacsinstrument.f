module module_pacsinstrument

    use iso_fortran_env,        only : ERROR_UNIT, OUTPUT_UNIT
    use module_filtering,       only : FilterUncorrelated
    use module_fitstools,       only : ft_check_error_cfitsio, ft_close, ft_create_header, ft_open, ft_open_image,                 &
                                       ft_read_image, ft_read_keyword_hcss, ft_read_slice, ft_test_extension
    use module_math,            only : DEG2RAD, mInf, pInf, mean, nint_down, nint_up
    use module_pacsobservation, only : PacsObservationSlice, PacsObservation
    use module_pointingmatrix,  only : PointingElement, xy2pmatrix, xy2roi, roi2pmatrix
    use module_projection,      only : convex_hull, surface_convex_polygon
    use module_string,          only : strinteger
    use module_tamasis,         only : get_tamasis_path, tamasis_path_len, p, POLICY_KEEP, POLICY_MASK, POLICY_REMOVE
    use module_wcs,             only : init_astrometry, ad2xy_gnomonic, ad2xy_gnomonic_vect, ad2xys_gnomonic, refpix_area
    use omp_lib
    implicit none
    private

    public :: NDIMS
    public :: NVERTICES
    public :: NEAREST_NEIGHBOUR, SHARP_EDGES

    public :: PacsInstrument

    public :: get_calfile
    public :: multiplexing_direct
    public :: multiplexing_transpose
    public :: read_calibration_files_new
    public :: read_detector_mask
    public :: read_filter_calibration_ncorrelations
    public :: read_filter_calibration
    public :: read_filter_filename

    integer, parameter :: NDIMS = 2
    integer, parameter :: NVERTICES = 4
    integer, parameter :: SHAPE_BLUE(2) = [32, 64]
    integer, parameter :: SHAPE_RED (2) = [16, 32]
    integer, parameter :: NEAREST_NEIGHBOUR = 0
    integer, parameter :: SHARP_EDGES = 1
    integer, parameter :: DISTORTION_DEGREE = 3
    real(p), parameter :: SAMPLING = 0.024996_p

    ! these calibration files are taken from HCSS 4.0.70
    character(len=*), parameter :: FILENAME_SAA = 'PCalPhotometer_SubArrayArray_FM_v5.fits'
    character(len=*), parameter :: FILENAME_AI  = 'PCalPhotometer_ArrayInstrument_FM_v5.fits'
    character(len=*), parameter :: FILENAME_BPM = 'PCalPhotometer_BadPixelMask_FM_v5.fits'
    character(len=*), parameter :: FILENAME_FF  = 'PCalPhotometer_FlatField_FM_v3.fits'
    character(len=*), parameter :: FILENAME_IB  = 'PCalPhotometer_InvnttBS_FM_v1.fits[Contents]'
    character(len=*), parameter :: FILENAME_IG  = 'PCalPhotometer_InvnttBL_FM_v1.fits[Contents]'
    character(len=*), parameter :: FILENAME_IR  = 'PCalPhotometer_InvnttRed_FM_v1.fits[Contents]'
    character(len=*), parameter :: FILENAME_RES = 'PCalPhotometer_Responsivity_FM_v5.fits'

    type PacsInstrument
 
        character(len=5)       :: band
        integer                :: nrows
        integer                :: ncolumns
        integer                :: ndetectors
        integer                :: fine_sampling_factor
        logical                :: transparent_mode

        logical*1, allocatable :: bad(:,:)  ! .true. for bad detectors
        logical*1, allocatable :: mask(:,:) ! .true. if the bad detector has been filtered out in the flattened list of detectors
        integer, allocatable   :: ij(:,:)
        integer, allocatable   :: pq(:,:)

        real(p)                :: responsivity ! Jy/V
        real(p), allocatable   :: detector_center_all(:,:,:)   ! remove me
        real(p), allocatable   :: detector_corner_all(:,:,:,:) ! remove me
        real(p), allocatable   :: detector_area_all(:,:)       ! remove me
        real(p), allocatable   :: flatfield_optical_all(:,:)   ! remove me
        real(p), allocatable   :: flatfield_detector_all(:,:)  ! remove me
        real(p), allocatable   :: flatfield_total_all(:,:)     ! remove me

        real(p), allocatable   :: detector_center(:,:)
        real(p), allocatable   :: detector_corner(:,:)
        real(p), allocatable   :: detector_area(:)
        real(p), allocatable   :: flatfield_optical(:)
        real(p), allocatable   :: flatfield_detector(:)
        real(p), allocatable   :: flatfield_total(:)

        real(p)                :: distortion_yz(NDIMS, DISTORTION_DEGREE, DISTORTION_DEGREE, DISTORTION_DEGREE)

    contains

        private
        procedure, public :: init
        procedure, public :: init_with_calfiles
        procedure, public :: init_with_variables
        procedure, public :: compute_map_header
        procedure, public :: find_minmax
        procedure, public :: compute_projection
        procedure, public :: destroy
        procedure, public :: read
        procedure, public, nopass :: read_detector_mask

        procedure, public, nopass :: uv2yz
        procedure, public, nopass :: yz2ad

        procedure :: compute_projection_nearest_neighbour
        procedure :: compute_projection_sharp_edges
        procedure :: filter_detectors
        procedure :: read_one
        procedure :: read_calibration_files
        procedure, public, nopass :: read_calibration_files_new

    end type PacsInstrument


contains


    subroutine init(this, band, transparent_mode, fine_sampling_factor, detector_mask, status)

        class(PacsInstrument), intent(inout) :: this
        character(len=*), intent(in)         :: band
        logical, intent(in)                  :: transparent_mode
        integer, intent(in)                  :: fine_sampling_factor
        logical*1, intent(in), optional      :: detector_mask(:,:) ! if all bad, use the calibration mask
        integer, intent(out)                 :: status

print *, 'INIT CALIBRATION PACS'
        ! check band
        if (band /= 'blue' .and. band /= 'green' .and. band /= 'red') then
            status = 1
            write (ERROR_UNIT,'(a)') "ERROR: Invalid band '" // band // "'. Valid values are 'red', 'green' or 'blue'."
            return
        end if

        ! check fine sampling factor
        if (all(fine_sampling_factor /= [1,2,4,8,16,32])) then
            status = 1
            write (ERROR_UNIT,'(a)') "ERROR: Invalid sampling factor '" //     &
                strinteger(fine_sampling_factor) // "'. Valid values are '1', '2', '4', '8', '16' or '32'."
            return
        end if

        this%band                 = band
        this%fine_sampling_factor = fine_sampling_factor
        this%transparent_mode     = transparent_mode

        if (band /= 'red') then
            this%nrows = SHAPE_BLUE(1)
            this%ncolumns = SHAPE_BLUE(2)
        else
            this%nrows = SHAPE_RED(1)
            this%ncolumns = SHAPE_RED(2)
        end if

        call this%read_calibration_files(status)
        if (status /= 0) return

        if (present(detector_mask)) then
            if (any(shape(detector_mask) /= [this%nrows,this%ncolumns])) then
                status = 1
                write (ERROR_UNIT,'(a)') 'INIT: the input bad detector mask has an invalid size.'
                return
            end if
            this%mask = detector_mask
        end if

        call this%filter_detectors()

    end subroutine init


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine init_with_calfiles(this, band, detector_mask, fine_sampling_factor, status)

        class(PacsInstrument), intent(inout) :: this
        character(len=*), intent(in)         :: band
        logical*1, intent(in)                :: detector_mask(:,:)
        integer, intent(in)                  :: fine_sampling_factor
        integer, intent(out)                 :: status

        real(p), allocatable :: detector_center_all(:,:,:)
        real(p), allocatable :: detector_corner_all(:,:,:,:)
        real(p), allocatable :: detector_area_all(:,:)
        real(p), allocatable :: flatfield_optical_all(:,:)
        real(p), allocatable :: flatfield_detector_all(:,:)
        real(p)              :: distortion_yz(NDIMS,DISTORTION_DEGREE,DISTORTION_DEGREE,DISTORTION_DEGREE)
        real(p)              :: responsivity

        call this%read_calibration_files_new(band, detector_mask, detector_center_all, detector_corner_all, detector_area_all,     &
                                             flatfield_optical_all, flatfield_detector_all, distortion_yz, responsivity, status)
        if (status /= 0) return

        call this%init_with_variables(band, detector_mask, fine_sampling_factor, status, detector_center_all, detector_corner_all, &
                                      detector_area_all, flatfield_optical_all, flatfield_detector_all, distortion_yz,             &
                                      responsivity)
        if (status /= 0) return

    end subroutine init_with_calfiles


    !-------------------------------------------------------------------------------------------------------------------------------


    ! a flattened array of detector values is travelled through columns and then through rows
    subroutine init_with_variables(this, band, detector_mask, fine_sampling_factor, status, detector_center_all,                   &
                                   detector_corner_all, detector_area_all, flatfield_optical_all, flatfield_detector_all,          &
                                   distortion_yz, responsivity)

        class(PacsInstrument), intent(inout) :: this
        character(len=*), intent(in)         :: band
        logical*1, intent(in)                :: detector_mask(:,:)
        integer, intent(in)                  :: fine_sampling_factor
        integer, intent(out)                 :: status
        real(p), intent(in), optional        :: detector_center_all(:,:,:)
        real(p), intent(in), optional        :: detector_corner_all(:,:,:,:)
        real(p), intent(in), optional        :: detector_area_all(:,:)
        real(p), intent(in), optional        :: flatfield_optical_all(:,:)
        real(p), intent(in), optional        :: flatfield_detector_all(:,:)
        real(p), intent(in), optional        :: distortion_yz(:,:,:,:)
        real(p), intent(in), optional        :: responsivity

        integer :: nrows, ncolumns, idetector, ip, iq

        status = 1

        ! check band
        if (band /= 'blue' .and. band /= 'green' .and. band /= 'red') then
            write (ERROR_UNIT,'(a)') "INIT: Invalid band '" // band // "'. Valid values are 'blue', 'green' and 'red'."
            return
        end if

        ! check fine sampling factor
        if (all(fine_sampling_factor /= [1,2,4,8,16,32])) then
            write (ERROR_UNIT,'(a)') "INIT: Invalid sampling factor '" //     &
                strinteger(fine_sampling_factor) // "'. Valid values are '1', '2', '4', '8', '16' or '32'."
            return
        end if

        if (band == 'red') then
            nrows    = SHAPE_RED(1)
            ncolumns = SHAPE_RED(2)
        else
            nrows    = SHAPE_BLUE(1)
            ncolumns = SHAPE_BLUE(2)
        end if

        if (size(detector_mask,1) /= nrows .or. size(detector_mask,2) /= ncolumns) then
            write (ERROR_UNIT, '(a,4(i0,a))') 'Invalid shape of the input detector mask (', size(detector_mask,1), ',',            &
                  size(detector_mask,2), ') instead of (', nrows, ',', ncolumns, ')'
            return
        end if

        ! get basic information
        this%band       = band
        this%nrows      = nrows
        this%ncolumns   = ncolumns
        this%ndetectors = count(.not. detector_mask)
        this%fine_sampling_factor = fine_sampling_factor

        ! allocate space for the arrays
        allocate (this%mask(nrows,ncolumns))
        this%mask = detector_mask

        allocate (this%ij(NDIMS, this%ndetectors))
        allocate (this%pq(NDIMS, this%ndetectors))

        if (present(detector_center_all)) then
            if (any(shape(detector_center_all) /= [NDIMS,shape(detector_mask)])) then
                write (ERROR_UNIT,'(a)') 'Invalid shape for the detector center array.'
                return
            end if
            allocate (this%detector_center(NDIMS,this%ndetectors))
        end if

        if (present(detector_corner_all)) then
            if (any(shape(detector_corner_all) /= [NDIMS,NVERTICES,shape(detector_mask)])) then
                write (ERROR_UNIT,'(a)') 'Invalid shape for the detector corner array.'
                return
            end if
            allocate (this%detector_corner(NDIMS,NVERTICES*this%ndetectors))
        end if

        if (present(detector_area_all)) then
            if (any(shape(detector_area_all) /= shape(detector_mask))) then
                write (ERROR_UNIT,'(a)') 'Invalid shape for the detector area.'
                return
            end if
            allocate (this%detector_area(this%ndetectors))
        end if

        if (present(flatfield_optical_all)) then
            if (any(shape(flatfield_optical_all) /= shape(detector_mask))) then
                write (ERROR_UNIT,'(a)') 'Invalid shape for the optical flat field.'
                return
            end if
            allocate (this%flatfield_optical(this%ndetectors))
        end if

        if (present(flatfield_detector_all)) then
            if (any(shape(flatfield_detector_all) /= shape(detector_mask))) then
                write (ERROR_UNIT,'(a)') 'Invalid shape for the detector flat field.'
                return
            end if
            allocate (this%flatfield_detector(this%ndetectors))
        end if

        if (present(flatfield_optical_all) .and. present(flatfield_detector_all)) then
            allocate (this%flatfield_total(this%ndetectors))
        end if

        if (present(distortion_yz)) then
            if (any(shape(distortion_yz) /= [NDIMS,DISTORTION_DEGREE,DISTORTION_DEGREE,DISTORTION_DEGREE])) then
                write (ERROR_UNIT,'(a)') 'Invalid shape for the distortion array.'
                return
            end if
            this%distortion_yz = distortion_yz
        end if

        if (present(responsivity)) then
            this%responsivity = responsivity
        else
            this%responsivity = 1._p
        end if

        idetector = 1

        do ip = 1, nrows

            do iq = 1, ncolumns
                if (this%mask(ip,iq)) cycle
                this%pq(1,idetector) = ip-1
                this%pq(2,idetector) = iq-1
                this%ij(1,idetector) = mod(ip-1, 16)
                this%ij(2,idetector) = mod(iq-1, 16)
                if (present(detector_center_all)) then
                    this%detector_center(:,idetector) = detector_center_all(:,ip,iq)
                end if
                if (present(detector_corner_all)) then
                    this%detector_corner(:,NVERTICES * (idetector-1)+1:NVERTICES*idetector) = detector_corner_all(:,:,ip,iq)
                end if
                if (present(detector_area_all)) then
                    this%detector_area(idetector) = detector_area_all(ip,iq)
                end if
                if (present(flatfield_optical_all)) then
                    this%flatfield_optical(idetector) = flatfield_optical_all(ip,iq)
                end if
                if (present(flatfield_detector_all)) then
                    this%flatfield_detector(idetector) = flatfield_detector_all(ip,iq)
                end if
                idetector = idetector + 1
            end do
        end do
        
        if (present(flatfield_optical_all) .and. present(flatfield_detector_all)) then
            this%flatfield_total = this%flatfield_optical * this%flatfield_detector
        end if

        status = 0

    end subroutine init_with_variables


    !-------------------------------------------------------------------------------------------------------------------------------
 

    subroutine read_detector_mask(band, detector_mask, status, transparent_mode, reject_bad_line)

        character(len=*), intent(in)        :: band
        logical*1, intent(out), allocatable :: detector_mask(:,:)
        integer, intent(out)                :: status
        logical, intent(in), optional       :: transparent_mode
        logical, intent(in), optional       :: reject_bad_line

        logical*1, allocatable :: tmp(:,:)
        integer                :: nrows, ncolumns
        
        status = 1

        if (band /= 'blue' .and. band /= 'green' .and. band /= 'red') then
            write (ERROR_UNIT,'(a)') "READ_DETECTOR_MASK: Invalid band '" // band // "'. Valid values are 'blue', 'green' or 'red'."
            return
        end if

        if (band == 'red') then
            nrows    = SHAPE_RED(1)
            ncolumns = SHAPE_RED(2)
        else
            nrows    = SHAPE_BLUE(1)
            ncolumns = SHAPE_BLUE(2)
        end if

        call ft_read_image(get_calfile(FILENAME_BPM) // '[' // trim(band) // ']', tmp, status)
        if (status /= 0) return

        if (size(tmp,2) /= nrows .or. size(tmp,1) /= ncolumns) then
            write (ERROR_UNIT, '(a,4(i0,a))') 'Invalid shape of the calibration detector mask (', size(tmp,2), ',', size(tmp,1),   &
                  ') instead of (', nrows, ',', ncolumns, ')'
            return
        end if

        allocate (detector_mask(nrows,ncolumns))
        detector_mask = transpose(tmp)

        if (present(transparent_mode)) then
            if (transparent_mode) then
                if (band == 'red') then
                    detector_mask(1:8,1:8) = .true.
                    detector_mask(1:8,17:) = .true.
                    detector_mask(9:,:)    = .true.
                else
                    detector_mask(1:16,1:16) = .true.
                    detector_mask(1:16,33:)  = .true.
                    detector_mask(17:,:)     = .true.
                end if
            end if
        end if

        if (present(reject_bad_line)) then
            if (reject_bad_line .and. band /= 'red') then
                detector_mask(12,17:32) = .true.
            end if
        end if
        
        status = 0

    end subroutine read_detector_mask


    !-------------------------------------------------------------------------------------------------------------------------------
 

    subroutine read_calibration_files_new(band, detector_mask, detector_center_all, detector_corner_all, detector_area_all,  &
                                          flatfield_optical_all, flatfield_detector_all, distortion_yz, responsivity, status)

        character(len=*), intent(in)      :: band
        logical*1, intent(in)             :: detector_mask(:,:)
        real(p), intent(out), allocatable :: detector_center_all(:,:,:)
        real(p), intent(out), allocatable :: detector_corner_all(:,:,:,:)
        real(p), intent(out), allocatable :: detector_area_all(:,:)
        real(p), intent(out), allocatable :: flatfield_optical_all(:,:)
        real(p), intent(out), allocatable :: flatfield_detector_all(:,:)
        real(p), intent(out)              :: distortion_yz(NDIMS,DISTORTION_DEGREE,DISTORTION_DEGREE,DISTORTION_DEGREE)
        real(p), intent(out)              :: responsivity
        integer, intent(out)              :: status

        integer, parameter   :: HDU_CORNER_BLUE(4) = [8, 12, 16, 20]
        integer, parameter   :: HDU_CORNER_RED (4) = [6, 10, 14, 18]

        integer              :: ip, iq, ivertex, unit, ncolumns, nrows
        real(p), allocatable :: tmp2(:,:)
        real(p), allocatable :: tmp3(:,:,:)
        real(p), allocatable :: flatfield_total_all(:,:)
        real(p)              :: center(2,2)

        status = 1

        if (band /= 'blue' .and. band /= 'green' .and. band /= 'red') then
            write (ERROR_UNIT,'(a)') "READ_DETECTOR_MASK: Invalid band '" // band // "'. Valid values are 'blue', 'green' or 'red'."
            return
        end if

        if (band == 'red') then
            nrows    = SHAPE_RED(1)
            ncolumns = SHAPE_RED(2)
        else
            nrows    = SHAPE_BLUE(1)
            ncolumns = SHAPE_BLUE(2)
        end if

        if (size(detector_mask,1) /= nrows .or. size(detector_mask,2) /= ncolumns) then
            write (ERROR_UNIT, '(a,4(i0,a))') 'Invalid shape of the input detector mask (', size(detector_mask,1), ',',            &
                  size(detector_mask,2), ') instead of (', nrows, ',', ncolumns, ')'
            return
        end if

        allocate (detector_center_all   (NDIMS,nrows,ncolumns))
        allocate (detector_corner_all   (NDIMS,NVERTICES,nrows,ncolumns))
        allocate (detector_area_all     (nrows,ncolumns))
        allocate (flatfield_optical_all (nrows,ncolumns))
        allocate (flatfield_detector_all(nrows,ncolumns))
        allocate (flatfield_total_all   (nrows,ncolumns))

        select case (band)

            case ('blue')
                call ft_read_image(get_calfile(FILENAME_FF) // '+12',  tmp2, status)
                if (status /= 0) return
                flatfield_total_all = transpose(tmp2)

            case ('green')
                call ft_read_image(get_calfile(FILENAME_FF) // '+7',  tmp2, status)
                if (status /= 0) return
                flatfield_total_all = transpose(tmp2)

            case ('red')
                call ft_read_image(get_calfile(FILENAME_FF) // '+2',  tmp2, status)
                if (status /= 0) return
                flatfield_total_all = transpose(tmp2)

        end select

        if (band /= 'red') then

            ! UV centers
            call ft_read_image(get_calfile(FILENAME_SAA) // '[ublue]', tmp2, status)
            if (status /= 0) return
            detector_center_all(1,:,:) = transpose(tmp2)
            call ft_read_image(get_calfile(FILENAME_SAA) // '[vblue]', tmp2, status)
            if (status /= 0) return
            detector_center_all(2,:,:) = transpose(tmp2)

            ! UV corners
            do ivertex=1, NVERTICES
                call ft_read_image(get_calfile(FILENAME_SAA), tmp2, status, HDU_CORNER_BLUE(ivertex)  )
                if (status /= 0) return
                detector_corner_all(1,ivertex,:,:) = transpose(tmp2)
                call ft_read_image(get_calfile(FILENAME_SAA), tmp2, status, HDU_CORNER_BLUE(ivertex)+1)
                if (status /= 0) return
                detector_corner_all(2,ivertex,:,:) = transpose(tmp2)
            end do

            ! Distortion coefficients in the (y,z) plane
            call ft_read_image(get_calfile(FILENAME_AI) // '[ycoeffblue]', tmp3, status)
            if (status /= 0) return
            distortion_yz(1,:,:,:) = tmp3
            call ft_read_image(get_calfile(FILENAME_AI) // '[zcoeffblue]', tmp3, status)
            if (status /= 0) return
            distortion_yz(2,:,:,:) = tmp3

        else

            ! UV centers
            call ft_read_image(get_calfile(FILENAME_SAA) // '[ured]', tmp2, status)
            if (status /= 0) return
            detector_center_all(1,:,:) = transpose(tmp2)
            call ft_read_image(get_calfile(FILENAME_SAA) // '[vred]', tmp2, status)
            if (status /= 0) return
            detector_center_all(2,:,:) = transpose(tmp2)

            ! UV corners
            do ivertex=1, NVERTICES
                call ft_read_image(get_calfile(FILENAME_SAA), tmp2, status, HDU_CORNER_RED(ivertex)  )
                if (status /= 0) return
                detector_corner_all(1,ivertex,:,:) = transpose(tmp2)
                call ft_read_image(get_calfile(FILENAME_SAA), tmp2, status, HDU_CORNER_RED(ivertex)+1)
                if (status /= 0) return
                detector_corner_all(2,ivertex,:,:) = transpose(tmp2)
            end do

            ! Distortion coefficients in the (y,z) plane
            call ft_read_image(get_calfile(FILENAME_AI) // '[ycoeffred]', tmp3, status)
            if (status /= 0) return
            distortion_yz(1,:,:,:) = tmp3
            call ft_read_image(get_calfile(FILENAME_AI) // '[zcoeffred]', tmp3, status)
            if (status /= 0) return
            distortion_yz(2,:,:,:) = tmp3

        end if

        ! Responsivity
        call ft_open(get_calfile(FILENAME_RES) // '[' // trim(band) // ']', unit, status)
        if (status /= 0) return
        call ft_read_keyword_hcss(unit, 'Responsivity', responsivity, status=status)
        if (status /= 0) return
        call ft_close(unit, status)
        if (status /= 0) return

        ! detector area
        do ip = 1, nrows
            do iq = 1, ncolumns
                detector_area_all(ip,iq) = abs(surface_convex_polygon(real(uv2yz(detector_corner_all(:,:,ip,iq),         &
                    distortion_yz, 0._p), p))) * 3600._p**2
            end do
        end do
        
        ! optical and detector flat fields
        center = detector_area_all(nrows/2:nrows/2+1,ncolumns/2:ncolumns/2+1)
        flatfield_optical_all = detector_area_all / mean(reshape(center,[4]))
        flatfield_detector_all = flatfield_total_all / flatfield_optical_all
        where (detector_mask)
            flatfield_optical_all = 1._p
            flatfield_detector_all = 1._p
        end where

    end subroutine read_calibration_files_new


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine read_calibration_files(this, status)

        class(PacsInstrument), intent(inout) :: this
        integer, intent(out)                 :: status

        integer, parameter     :: HDU_CORNER_BLUE(4) = [8, 12, 16, 20]
        integer, parameter     :: HDU_CORNER_RED (4) = [6, 10, 14, 18]

        character(len=5)       :: band_name
        integer                :: ivertex, unit
        logical*1, allocatable :: tmplogical(:,:)
        real(p), allocatable   :: tmp2(:,:)
        real(p), allocatable   :: tmp3(:,:,:)

        allocate (this%flatfield_total_all(this%nrows,this%ncolumns))
        select case (this%band)

            case ('blue')
                band_name = 'blue'
                call ft_read_image(get_calfile(FILENAME_FF) // '+12',  tmp2, status)
                if (status /= 0) return
                this%flatfield_total_all = transpose(tmp2)

            case ('green')
                band_name = 'green'
                call ft_read_image(get_calfile(FILENAME_FF) // '+7',  tmp2, status)
                if (status /= 0) return
                this%flatfield_total_all = transpose(tmp2)

            case ('red')
                band_name = 'red'
                call ft_read_image(get_calfile(FILENAME_FF) // '+2',  tmp2, status)
                if (status /= 0) return
                this%flatfield_total_all = transpose(tmp2)

        end select

        allocate (this%detector_center_all(NDIMS,this%nrows,this%ncolumns))
        allocate (this%detector_corner_all(NDIMS,NVERTICES,this%nrows,this%ncolumns))
        if (this%band /= 'red') then

            ! UV centers
            call ft_read_image(get_calfile(FILENAME_SAA) // '[ublue]', tmp2, status)
            if (status /= 0) return
            this%detector_center_all(1,:,:) = transpose(tmp2)
            call ft_read_image(get_calfile(FILENAME_SAA) // '[vblue]', tmp2, status)
            if (status /= 0) return
            this%detector_center_all(2,:,:) = transpose(tmp2)

            ! UV corners
            do ivertex=1, NVERTICES
                call ft_read_image(get_calfile(FILENAME_SAA), tmp2, status, HDU_CORNER_BLUE(ivertex)  )
                if (status /= 0) return
                this%detector_corner_all(1,ivertex,:,:) = transpose(tmp2)
                call ft_read_image(get_calfile(FILENAME_SAA), tmp2, status, HDU_CORNER_BLUE(ivertex)+1)
                if (status /= 0) return
                this%detector_corner_all(2,ivertex,:,:) = transpose(tmp2)
            end do

            ! Distortion coefficients in the (y,z) plane
            call ft_read_image(get_calfile(FILENAME_AI) // '[ycoeffblue]', tmp3, status)
            if (status /= 0) return
            this%distortion_yz(1,:,:,:) = tmp3
            call ft_read_image(get_calfile(FILENAME_AI) // '[zcoeffblue]', tmp3, status)
            if (status /= 0) return
            this%distortion_yz(2,:,:,:) = tmp3

        else

            ! UV centers
            call ft_read_image(get_calfile(FILENAME_SAA) // '[ured]', tmp2, status)
            if (status /= 0) return
            this%detector_center_all(1,:,:) = transpose(tmp2)
            call ft_read_image(get_calfile(FILENAME_SAA) // '[vred]', tmp2, status)
            if (status /= 0) return
            this%detector_center_all(2,:,:) = transpose(tmp2)

            ! UV corners
            do ivertex=1, NVERTICES
                call ft_read_image(get_calfile(FILENAME_SAA), tmp2, status, HDU_CORNER_RED(ivertex)  )
                if (status /= 0) return
                this%detector_corner_all(1,ivertex,:,:) = transpose(tmp2)
                call ft_read_image(get_calfile(FILENAME_SAA), tmp2, status, HDU_CORNER_RED(ivertex)+1)
                if (status /= 0) return
                this%detector_corner_all(2,ivertex,:,:) = transpose(tmp2)
            end do

            ! Distortion coefficients in the (y,z) plane
            call ft_read_image(get_calfile(FILENAME_AI) // '[ycoeffred]', tmp3, status)
            if (status /= 0) return
            this%distortion_yz(1,:,:,:) = tmp3
            call ft_read_image(get_calfile(FILENAME_AI) // '[zcoeffred]', tmp3, status)
            if (status /= 0) return
            this%distortion_yz(2,:,:,:) = tmp3

        end if

        ! Bad pixel mask
        allocate (this%mask(this%nrows,this%ncolumns))
        call ft_read_image(get_calfile(FILENAME_BPM) // '[' // trim(band_name) // ']', tmplogical, status)
        if (status /= 0) return
        this%mask = transpose(tmplogical)

        ! mask detectors rejected in transparent mode
        if (this%transparent_mode) then
            if (this%band /= 'red') then
                this%mask(1:16,1:16) = .true.
                this%mask(1:16,33:)  = .true.
                this%mask(17:,:)     = .true.
            else
                this%mask(1:8,1:8) = .true.
                this%mask(1:8,17:) = .true.
                this%mask(9:,:)    = .true.
            end if
        end if

        ! Responsivity
        call ft_open(get_calfile(FILENAME_RES) // '[' // trim(band_name) // ']', unit, status)
        if (status /= 0) return
        call ft_read_keyword_hcss(unit, 'Responsivity', this%responsivity, status=status)
        if (status /= 0) return
        call ft_close(unit, status)
        if (status /= 0) return

    end subroutine read_calibration_files


    !-------------------------------------------------------------------------------------------------------------------------------


    function get_calfile(filename)

        character(len=*), intent(in)                     :: filename
        character(len=tamasis_path_len+10+len(filename)) :: get_calfile

        get_calfile = get_tamasis_path() // 'pacs/data/' // filename

    end function get_calfile


    !-------------------------------------------------------------------------------------------------------------------------------


    ! a flattened array of detector values is travelled through columns and then through rows
    subroutine filter_detectors(this)

        class(PacsInstrument), intent(inout) :: this

        integer :: idetector, ip, iq
        real(p) :: center(2,2)

        ! get the number of detectors
        this%ndetectors = count(.not. this%mask)

        allocate (this%ij(NDIMS, this%ndetectors))
        allocate (this%pq(NDIMS, this%ndetectors))
        allocate (this%detector_center(NDIMS, this%ndetectors))
        allocate (this%detector_corner(NDIMS, NVERTICES * this%ndetectors))
        allocate (this%detector_area  (this%ndetectors))

        allocate (this%flatfield_detector(this%ndetectors))
        allocate (this%flatfield_optical(this%ndetectors))
        allocate (this%flatfield_total(this%ndetectors))

        allocate (this%detector_area_all(this%nrows, this%ncolumns))
        allocate (this%flatfield_detector_all(this%nrows,this%ncolumns))
        allocate (this%flatfield_optical_all (this%nrows,this%ncolumns))

        do ip = 1, this%nrows
            do iq = 1, this%ncolumns
                this%detector_area_all(ip,iq) = abs(surface_convex_polygon(real(this%uv2yz(this%detector_corner_all(:,:,ip,iq),    &
                    this%distortion_yz, 0._p), p))) * 3600._p**2
            end do
        end do

        center = this%detector_area_all(this%nrows/2:this%nrows/2+1,this%ncolumns/2:this%ncolumns/2+1)
        this%flatfield_optical_all = this%detector_area_all / mean(reshape(center,[4]))
        this%flatfield_detector_all = this%flatfield_total_all / this%flatfield_optical_all
        where (this%mask)
            this%flatfield_total_all    = 1
            this%flatfield_detector_all = 1
            this%flatfield_optical_all  = 1
        end where
        
        idetector = 1
        do ip = 1, this%nrows
            do iq = 1, this%ncolumns
                if (this%mask(ip,iq)) cycle
                this%pq(1,idetector) = ip-1
                this%pq(2,idetector) = iq-1
                this%ij(1,idetector) = mod(ip-1, 16)
                this%ij(2,idetector) = mod(iq-1, 16)
                this%detector_center(:,idetector) = this%detector_center_all(:,ip,iq)
                this%detector_corner(:,NVERTICES * (idetector-1)+1:NVERTICES*idetector) = this%detector_corner_all(:,:,ip,iq)
                this%detector_area(idetector) = this%detector_area_all(ip,iq)
                this%flatfield_optical (idetector) = this%flatfield_optical_all (ip,iq)
                this%flatfield_detector(idetector) = this%flatfield_detector_all(ip,iq)
                this%flatfield_total   (idetector) = this%flatfield_total_all   (ip,iq)
                idetector = idetector + 1
            end do
        end do
        
    end subroutine filter_detectors


    !-------------------------------------------------------------------------------------------------------------------------------


    function uv2yz(uv, distortion_yz, chop) result(yz)

        real(p), intent(in) :: uv(:,:)
        real(p), intent(in) :: distortion_yz(NDIMS, DISTORTION_DEGREE, DISTORTION_DEGREE, DISTORTION_DEGREE)
        real(p), intent(in) :: chop ! chop angle in degree
        real(p)             :: yz(NDIMS, size(uv,2))

        real(p) :: a_pow, ratio
        real(p) :: u_pow(size(uv,2)), v_pow(size(uv,2))
        integer :: n, i, j, k, l

        yz = 0
        n = size(uv,2)
        do k = 1, DISTORTION_DEGREE
            if (k /= 1) then
                u_pow = u_pow * uv(1,:)
            else
                u_pow = 1
            endif
            do j = 1, DISTORTION_DEGREE
              if (j /= 1) then
                  v_pow = v_pow * uv(2,:)
              else
                  v_pow = 1
              end if
              do i = 1, DISTORTION_DEGREE
                  if (i /= 1) then
                      a_pow = a_pow * chop
                  else
                      a_pow = 1
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
        yz = - yz / 3600._p

    end function uv2yz
 

    !-------------------------------------------------------------------------------------------------------------------------------


    function yz2ad(yz, ra0, dec0, pa0) result (ad)

        real(p), intent(in) :: yz(:,:) ! in degrees
        real(p), intent(in) :: ra0, dec0, pa0
        real(p)             :: ad(NDIMS,size(yz,2))

        real(p) :: cospa, sinpa
        integer :: i

        cospa =  cos(pa0*DEG2RAD)
        sinpa = -sin(pa0*DEG2RAD)

        do i=1, size(yz, 2)

            ad(2,i) = dec0 + (yz(1,i) * sinpa + yz(2,i) * cospa)
            ad(1,i) = ra0  + (yz(1,i) * cospa - yz(2,i) * sinpa) / cos(ad(2,i) * DEG2RAD)

        end do

    end function yz2ad


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine compute_map_header(this, obs, oversampling, resolution, header, status)

        class(PacsInstrument), intent(in)  :: this
        class(PacsObservation), intent(in) :: obs
        logical, intent(in)                :: oversampling
        real(p), intent(in)                :: resolution
        character(len=2880), intent(out)   :: header
        integer, intent(out)               :: status

        integer :: nx, ny
        integer :: ixmin, ixmax, iymin, iymax
        real(p) :: ra0, dec0, xmin, xmax, ymin, ymax


        call obs%compute_center(ra0, dec0)

        call ft_create_header(0, 0, [[-resolution/3600._p, 0._p], [0._p, resolution/3600._p]], ra0, dec0, 1._p, 1._p,  &
             header)

        call init_astrometry(header, status=status)
        if (status /= 0) return

        call this%find_minmax(obs, oversampling, xmin, xmax, ymin, ymax)

        ixmin = nint(xmin)
        ixmax = nint(xmax)
        iymin = nint(ymin)
        iymax = nint(ymax)

        nx = ixmax-ixmin+1
        ny = iymax-iymin+1

        ! move the reference pixel (not the reference value!)
        call ft_create_header(nx, ny, [[-resolution/3600._p, 0._p], [0._p, resolution/3600._p]], ra0, dec0, -ixmin+2._p,           &
             -iymin+2._p, header)

    end subroutine compute_map_header


    !-------------------------------------------------------------------------------------------------------------------------------


    ! find minimum and maximum pixel coordinates in maps
    subroutine find_minmax(this, obs, oversampling, xmin, xmax, ymin, ymax)

        class(PacsInstrument), intent(in)  :: this
        class(PacsObservation), intent(in) :: obs
        logical, intent(in)                :: oversampling
        real(p), intent(out)               :: xmin, xmax, ymin, ymax

        real(p)              :: ra, dec, pa, chop
        real(p), allocatable :: hull_uv(:,:), hull(:,:)
        integer, allocatable :: ihull(:)
        integer              :: sampling_factor, isample, islice, itime

        xmin = pInf
        xmax = mInf
        ymin = pInf
        ymax = mInf

        call convex_hull(real(this%detector_corner, p), ihull)

        allocate(hull_uv(NDIMS,size(ihull)))
        allocate(hull(NDIMS, size(ihull)))
        hull_uv = this%detector_corner(:,ihull)

        do islice = 1, obs%nslices

            ! check if it is required to interpolate pointing positions
            if (oversampling) then
               sampling_factor = this%fine_sampling_factor * obs%slice(islice)%compression_factor
            else
               sampling_factor = 1
            end if

#ifndef IFORT
            !$omp parallel do default(shared) reduction(min:xmin,ymin) reduction(max:xmax,ymax) &
            !$omp private(itime,isample,ra,dec,pa,chop,hull)
#endif
            do itime = 1, obs%slice(islice)%nsamples * sampling_factor

                isample = (itime - 1) / sampling_factor + 1
                if (obs%slice(islice)%p(isample)%removed) cycle

                call obs%get_position_index(islice, itime, sampling_factor, ra, dec, pa, chop)
                hull = this%uv2yz(hull_uv, this%distortion_yz, chop)
                hull = this%yz2ad(hull, ra, dec, pa)
                hull = ad2xy_gnomonic(hull)
                xmin = min(xmin, minval(hull(1,:)))
                xmax = max(xmax, maxval(hull(1,:)))
                ymin = min(ymin, minval(hull(2,:)))
                ymax = max(ymax, maxval(hull(2,:)))

             end do
#ifndef IFORT
             !$omp end parallel do
#endif

        end do

    end subroutine find_minmax


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine compute_projection(this, method, obs, oversampling, header, nx, ny, pmatrix, status)

        class(PacsInstrument), intent(in)  :: this
        integer, intent(in)                :: method
        class(PacsObservation), intent(in) :: obs
        logical, intent(in)                :: oversampling
        character(len=*), intent(in)       :: header
        integer, intent(in)                :: nx, ny
        type(PointingElement), intent(out) :: pmatrix(:,:,:)
        integer, intent(out)               :: status

        integer :: npixels_per_sample
        integer :: count1, count2, count_rate, count_max
        logical :: out

        write(*,'(a)', advance='no') 'Info: Computing the projector... '
        call system_clock(count1, count_rate, count_max)

        call init_astrometry(header, status=status)
        if (status /= 0) return

        select case (method)

            case (NEAREST_NEIGHBOUR)
                call this%compute_projection_nearest_neighbour(obs, oversampling, nx, ny, pmatrix, npixels_per_sample, out)
            
            case (SHARP_EDGES)
               call this%compute_projection_sharp_edges(obs, oversampling, nx, ny, pmatrix, npixels_per_sample, out)

        end select

        call system_clock(count2, count_rate, count_max)
        write(*,'(f7.2,a)') real(count2-count1)/count_rate, 's'

        if (npixels_per_sample > size(pmatrix,1)) then
            status = 1
            write(ERROR_UNIT,'(a,i0,a)') 'Error: Please update npixels_per_sample to ', npixels_per_sample, '.'
        else if (npixels_per_sample < size(pmatrix,1)) then
            write(OUTPUT_UNIT,'(a,i0,a)') 'Warning: You may update npixels_per_sample to ', npixels_per_sample, '.'
        end if

        if (out) then
            write (OUTPUT_UNIT,'(a)') 'Warning: Some detectors fall outside the map.'
        end if

    end subroutine compute_projection


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine compute_projection_nearest_neighbour(this, obs, oversampling, nx, ny, pmatrix, npixels_per_sample, out)

        class(PacsInstrument), intent(in)  :: this
        class(PacsObservation), intent(in) :: obs
        logical, intent(in)                :: oversampling
        integer, intent(in)                :: nx, ny
        type(PointingElement), intent(out) :: pmatrix(:,:,:)
        integer, intent(out)               :: npixels_per_sample
        logical, intent(out)               :: out

        real(p) :: coords(NDIMS,this%ndetectors), coords_yz(NDIMS,this%ndetectors)
        real(p) :: x(this%ndetectors), y(this%ndetectors), s(this%ndetectors)
        real(p) :: ra, dec, pa, chop, chop_old, reference_area
        integer   :: ifine, isample, islice, itime, ivalid, nsamples, nvalids, sampling_factor, dest
        integer, allocatable :: valids(:)

        npixels_per_sample = 1
        reference_area = refpix_area()

        dest = 0
        out = .false.

        ! loop over the observations
        do islice = 1, obs%nslices

            nsamples = obs%slice(islice)%nsamples
            nvalids  = obs%slice(islice)%nvalids
            
            chop_old = pInf

            ! check if it is required to interpolate pointing positions
            if (oversampling) then
               sampling_factor = this%fine_sampling_factor * obs%slice(islice)%compression_factor
            else
               sampling_factor = 1
            end if

            allocate (valids(nvalids))
            ivalid = 1
            do isample = 1, nsamples
                if (obs%slice(islice)%p(isample)%removed) cycle
                valids(ivalid) = isample
                ivalid = ivalid + 1
            end do

            !$omp parallel do default(shared) firstprivate(chop_old)                            &
            !$omp private(ifine, isample, itime, ra, dec, pa, chop, coords, coords_yz, x, y, s) &
            !$omp reduction(max : npixels_per_sample) reduction(.or. : out)
            
            ! loop over the samples which have not been removed
            do ivalid = 1, nvalids
                isample = valids(ivalid)
                do ifine = 1, sampling_factor
                     itime = (isample - 1) * sampling_factor + ifine
                     call obs%get_position_index(islice, itime, sampling_factor, ra, dec, pa, chop)
                     if (abs(chop-chop_old) > 1.e-2_p) then
                         coords_yz = this%uv2yz(this%detector_center, this%distortion_yz, chop)
                         chop_old = chop
                     end if
                     coords = this%yz2ad(coords_yz, ra, dec, pa)
                     call ad2xys_gnomonic(coords, x, y, s)
                     itime  = (ivalid-1) * sampling_factor + ifine + dest
                     ! the input map has flux densities, not surface brightness
                     ! f_idetector = f_imap * weight
                     ! with weight = detector_area / pixel_area
                     ! and pixel_area = reference_area / s
                     call xy2pmatrix(x, y, nx, ny, out, pmatrix(1,itime,:))
                     pmatrix(1,itime,:)%weight = s * this%detector_area / reference_area
                 end do
             end do

             !$omp end parallel do
             dest = dest + nvalids * sampling_factor
             deallocate (valids)

        end do

        pmatrix(2:,:,:)%pixel  = -1
        pmatrix(2:,:,:)%weight = 0

    end subroutine compute_projection_nearest_neighbour


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine compute_projection_sharp_edges(this, obs, oversampling, nx, ny, pmatrix, npixels_per_sample, out)

        class(PacsInstrument), intent(in)  :: this
        class(PacsObservation), intent(in) :: obs
        logical, intent(in)                :: oversampling
        integer, intent(in)                :: nx, ny
        type(PointingElement), intent(out) :: pmatrix(:,:,:)
        integer, intent(out)               :: npixels_per_sample
        logical, intent(out)               :: out

        real(p) :: coords(NDIMS,this%ndetectors*NVERTICES), coords_yz(NDIMS,this%ndetectors*NVERTICES)
        real(p) :: ra, dec, pa, chop, chop_old
        integer :: roi(NDIMS,2,this%ndetectors)
        integer :: ifine, isample, islice, itime, ivalid, nsamples, nvalids, sampling_factor, dest
        integer, allocatable :: valids(:)

        npixels_per_sample = 0
        dest = 0
        out = .false.

        ! loop over the observations
        do islice = 1, obs%nslices

            nsamples = obs%slice(islice)%nsamples
            nvalids  = obs%slice(islice)%nvalids
            
            chop_old = pInf

            ! check if it is required to interpolate pointing positions
            if (oversampling) then
               sampling_factor = this%fine_sampling_factor * obs%slice(islice)%compression_factor
            else
               sampling_factor = 1
            end if

            allocate (valids(nvalids))
            ivalid = 1
            do isample = 1, nsamples
                if (obs%slice(islice)%p(isample)%removed) cycle
                valids(ivalid) = isample
                ivalid = ivalid + 1
            end do

            !$omp parallel do default(shared) firstprivate(chop_old)                        &
            !$omp private(ifine, isample, itime, ra, dec, pa, chop, coords, coords_yz, roi) &
            !$omp reduction(max : npixels_per_sample) reduction(.or. : out)
            
            ! loop over the samples which have not been removed
            do ivalid = 1, nvalids
                isample = valids(ivalid)
                do ifine = 1, sampling_factor
                     itime = (isample - 1) * sampling_factor + ifine
                     call obs%get_position_index(islice, itime, sampling_factor, ra, dec, pa, chop)
                     if (abs(chop-chop_old) > 1.e-2_p) then
                         coords_yz = this%uv2yz(this%detector_corner, this%distortion_yz, chop)
                         chop_old = chop
                     end if
                     itime  = (ivalid-1) * sampling_factor + ifine + dest
                     coords = this%yz2ad(coords_yz, ra, dec, pa)
                     coords = ad2xy_gnomonic(coords)
                     roi    = xy2roi(coords, NVERTICES)
                     call roi2pmatrix(roi, NVERTICES, coords, nx, ny, npixels_per_sample, out, pmatrix(:,itime,:))
                 end do
             end do

             !$omp end parallel do
             dest = dest + nvalids * sampling_factor
             deallocate (valids)

        end do

    end subroutine compute_projection_sharp_edges


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine read(this, obs, signal, mask, status, verbose)

        class(PacsInstrument), intent(in)  :: this
        class(PacsObservation), intent(in) :: obs
        real(p), intent(out)               :: signal(:,:)
        logical(1), intent(out)            :: mask(:,:)
        integer, intent(out)               :: status
        logical, intent(in), optional      :: verbose

        integer :: dest
        integer :: iobs, nobs
        integer :: count1, count2, count_rate, count_max
        logical :: verbose_

        verbose_ = .false.
        if (present(verbose)) verbose_ = verbose

        if (verbose_) then
            write(*,'(a)', advance='no') 'Info: Reading timeline... '
        end if
        call system_clock(count1, count_rate, count_max)

        status = 1
        nobs   = obs%nslices

        ! check that the total number of samples in the observations is equal
        ! to the number of samples in the signal and mask arrays
        if (obs%nvalids /= size(signal,1) .or. obs%nvalids /= size(mask,1)) then
            write (ERROR_UNIT,'(a)') 'Error: read: invalid dimensions.'
            return
        end if

        ! set mask to ok
        mask = .false.

        ! loop over the PACS observations
        dest = 1
        do iobs = 1, nobs
            call this%read_one(obs%slice(iobs), signal(dest:dest+obs%slice(iobs)%nvalids-1,:),                                     &
                 mask(dest:dest+obs%slice(iobs)%nvalids-1,:), status)
            if (status /= 0) return
            dest = dest + obs%slice(iobs)%nvalids
        end do

        call system_clock(count2, count_rate, count_max)
        if (verbose_) then
            write (*,'(f6.2,a)') real(count2-count1) / count_rate, 's'
        end if

    end subroutine read


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine read_one(this, obs, signal, mask, status)

        class(PacsInstrument), intent(in)       :: this
        class(PacsObservationSlice), intent(in) :: obs
        real(p), intent(out)                    :: signal(:,:)
        logical*1, intent(out)                  :: mask  (:,:)
        integer, intent(out)                    :: status

        integer                :: first, last
        integer                :: ip, iq
        integer                :: idetector, ndetectors, unit
        integer, allocatable   :: imageshape(:)
        logical                :: mask_found
        integer*4, allocatable :: maskcompressed(:)
        integer*4              :: maskval
        integer                :: ncompressed
        integer                :: isample, icompressed, ibit
        integer                :: firstcompressed, lastcompressed
        real(p), allocatable   :: signal_(:)
        logical*1, allocatable :: mask_(:)

        ndetectors = size(signal, 2)

        ! get first and last valid pointing, to decrease I/O
        do first = 1, obs%nsamples
            if (.not. obs%p(first)%removed) exit
        end do
        do last = obs%nsamples, 1, -1
            if (.not. obs%p(last)%removed) exit
        end do

        ! set mask from policy
        mask = spread(pack(obs%p(first:last)%masked, .not. obs%p(first:last)%removed), 2, ndetectors)

        ! read signal HDU
        call ft_open_image(trim(obs%filename) // '[Signal]', unit, 3, imageshape, status)
        if (status /= 0) return

        allocate (signal_(last - first + 1), mask_(last - first + 1))
        do idetector = 1, ndetectors

            ip = this%pq(1,idetector)
            iq = this%pq(2,idetector)

            call ft_read_slice(unit, first, last, iq+1, ip+1, imageshape, signal_, status)
            if (status /= 0) return

            signal(:,idetector) = pack(signal_, .not. obs%p(first:last)%removed)

        end do

        call ft_close(unit, status)
        if (status /= 0) return

        ! read Mask HDU
        mask_found = ft_test_extension(trim(obs%filename)//'[Master]', status)
        if (status /= 0) return

        if (.not. mask_found) then
!!$            write (*,'(a)') 'Info: mask Master is not found.'
            return
        end if

        call ft_open_image(trim(obs%filename) // '[Master]', unit, 3, imageshape, status)
        if (status /= 0) return

        firstcompressed = (first - 1) / 32 + 1
        lastcompressed  = (last  - 1) / 32 + 1
        ncompressed = lastcompressed - firstcompressed + 1

        if (lastcompressed > imageshape(1)) then
            status = 1
            write (ERROR_UNIT, '(a)') 'There is not enough samples in ' // trim(obs%filename) // '[Master] for this pointing.'
            return
        end if

        allocate (maskcompressed(ncompressed))

        do idetector = 1, size(mask,2)

            ip = this%pq(1,idetector)
            iq = this%pq(2,idetector)

            mask_ = .false.

            call ft_read_slice(unit, firstcompressed, lastcompressed, iq+1, ip+1, imageshape, maskcompressed, status)
            if (status /= 0) return

            ! loop over the bytes of the compressed mask
            do icompressed = firstcompressed, lastcompressed

                maskval = maskcompressed(icompressed-firstcompressed+1)
                if (maskval == 0) cycle

                isample = (icompressed-1)*32 - first + 2

                ! loop over the bits of a compressed mask byte
                do ibit = max(0, first - (icompressed-1)*32 - 1), min(31, last - (icompressed-1)*32 - 1)
                    mask_(isample+ibit) = btest(maskval,ibit)
                end do

            end do
            
            mask(:,idetector) = mask(:,idetector) .or. pack(mask_, .not. obs%p(first:last)%removed)

        end do

        call ft_close(unit, status)

    end subroutine read_one


    !-------------------------------------------------------------------------------------------------------------------------------


    function read_filter_calibration_ncorrelations(band, status) result (ncorrelations)

        character(len=*), intent(in) :: band
        integer, intent(out)         :: status
        integer                      :: ncorrelations

        integer           :: unit
        character(len=70) :: comment

        select case (band)
            case ('blue')
                call ft_open(get_calfile(FILENAME_IB), unit, status)
            case ('green')
                call ft_open(get_calfile(FILENAME_IG), unit, status)
            case ('red')
                call ft_open(get_calfile(FILENAME_IR), unit, status)
            case default
                status = 1
                write (ERROR_UNIT,'(a)') "READ_FILTER_CALIBRATION_NCORRELATIONS: invalid band: '" // band // "'."
        end select
        if (status /= 0) return

        call ftgkyj(unit, 'NAXIS1', ncorrelations, comment, status)
        if (ft_check_error_cfitsio(status, unit)) return

        ncorrelations = ncorrelations - 1

        call ft_close(unit, status)

    end function read_filter_calibration_ncorrelations


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine read_filter_calibration(band, mask, filter, status)

        character(len=*), intent(in)          :: band
        logical*1, intent(in)                 :: mask(:,:)
        type(FilterUncorrelated), intent(out) :: filter
        integer, intent(out)                  :: status

        select case (band)
            case ('blue')
                call read_filter_filename(get_calfile(FILENAME_IB), mask, filter, status)
            case ('green')
                call read_filter_filename(get_calfile(FILENAME_IG), mask, filter, status)
            case ('red')
                call read_filter_filename(get_calfile(FILENAME_IR), mask, filter, status)
            case default
                status = 1
                write (ERROR_UNIT,'(a)') "READ_FILTER_CALIBRATION: invalid band: '" // band // "'."
        end select

    end subroutine read_filter_calibration


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine read_filter_filename(filename, mask, filter, status)

        character(len=*), intent(in)          :: filename
        logical*1, intent(in)                 :: mask(:,:)
        type(FilterUncorrelated), intent(out) :: filter
        integer, intent(out)                  :: status

        real(p), allocatable :: data(:,:)
        integer              :: idetector, ndetectors, nrows, ncolumns, ip, iq

        ndetectors = size(mask) - count(mask)
        nrows = size(mask, 1)
        ncolumns = size(mask, 2)

        call ft_read_image(filename, data, status)
        if (status /= 0) return

        filter%ncorrelations = size(data, 1) - 1
        if (size(data, 2) /= size(mask)) then
            status = 1
            write (ERROR_UNIT,'(a,2(i0,a))') "Calibration file '" // filename // "' contains an incorrect number of detectors: '", &
                  size(data, 2), "' instead of '", size(mask), "'."
            return
        end if
        filter%ndetectors = ndetectors
        filter%bandwidth  = 2 * filter%ncorrelations + 1
        
        allocate (filter%data(filter%ncorrelations+1,filter%ndetectors))

        !XXX CHECK P, Q layout in calibration file
        idetector = 1
        do ip = 1, nrows
            do iq = 1, ncolumns
                if (mask(ip,iq)) cycle
                filter%data(:,idetector) = data(:, (ip-1) * ncolumns + iq)
                idetector = idetector + 1
            end do
        end do

    end subroutine read_filter_filename


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine multiplexing_direct(signal, sampled_signal, sampling, ij)

        real(p), intent(in)  :: signal(:,:)
        integer, intent(in)  :: sampling
        real(p), intent(out) :: sampled_signal(size(signal,1)/sampling,size(signal,2))
        integer, intent(in)  :: ij(2,size(signal,2))

        integer :: ndetectors, idetector, isample, j
        real(p) :: frac

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
                frac = (j * sampling / 16._p) - floor(j * sampling / 16._p)
                sampled_signal(:,idetector) = (1-frac) * signal(isample  ::sampling,idetector) + &
                                                 frac  * signal(isample+1::sampling,idetector)
            end do
            !$omp end parallel do
        end if

    end subroutine multiplexing_direct


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine multiplexing_transpose(sampled_signal, signal, sampling, ij)

        real(p), intent(in) :: sampled_signal(:,:)
        integer, intent(in) :: sampling
        real(p), intent(out):: signal(size(sampled_signal,1)*sampling,size(sampled_signal,2))
        integer, intent(in) :: ij(2,size(sampled_signal,2))

        integer :: ndetectors, isample, idetector, j
        real(p) :: frac

        ndetectors = size(sampled_signal, 2)
        if (sampling == 16) then
            signal = 0
            !$omp parallel do default(shared) private(idetector)
            do idetector = 1, ndetectors
                signal(ij(2,idetector)+1::16,idetector) = sampled_signal(:,idetector)
            end do
            !$omp end parallel do
        else
            !$omp parallel default(shared) private(idetector,isample,frac,j)
            !$omp workshare
            signal = 0
            !$omp end workshare
            !$omp do
            do idetector = 1, ndetectors
                j = ij(2,idetector)
                isample = j * sampling / 16 + 1 ! isample+1 is guaranteed to exist
                frac = (j * sampling / 16._p) - floor(j * sampling / 16._p)
                signal(isample  ::sampling,idetector) = (1-frac) * sampled_signal(isample::sampling,idetector)
                signal(isample+1::sampling,idetector) =    frac  * sampled_signal(isample::sampling,idetector)
            end do
            !$omp end do
            !$omp end parallel
        end if

    end subroutine multiplexing_transpose


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine destroy(this)

        class(PacsInstrument), intent(inout) :: this

        deallocate (this%mask)
        deallocate (this%ij)
        deallocate (this%pq)
        deallocate (this%flatfield_optical)
        deallocate (this%flatfield_detector)
        deallocate (this%flatfield_total)
        deallocate (this%detector_center)
        deallocate (this%detector_corner)
        deallocate (this%detector_area)

        if (allocated(this%flatfield_total_all))    deallocate (this%flatfield_total_all)
        if (allocated(this%flatfield_detector_all)) deallocate (this%flatfield_detector_all)
        if (allocated(this%flatfield_optical_all))  deallocate (this%flatfield_optical_all)
        if (allocated(this%detector_center_all))    deallocate (this%detector_center_all)
        if (allocated(this%detector_corner_all))    deallocate (this%detector_corner_all)
        if (allocated(this%detector_area_all))      deallocate (this%detector_area_all)

    end subroutine destroy


end module module_pacsinstrument
