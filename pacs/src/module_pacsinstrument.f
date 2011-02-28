module module_pacsinstrument

    use iso_fortran_env,        only : ERROR_UNIT, OUTPUT_UNIT
    use module_filtering,       only : FilterUncorrelated
    use module_fitstools,       only : FLEN_VALUE, ft_check_error_cfitsio, ft_close, ft_create_header, ft_open, ft_open_image,     &
                                       ft_read_image, ft_read_keyword, ft_read_keyword_hcss, ft_read_slice, ft_test_extension
    use module_math,            only : DEG2RAD, mInf, pInf, NaN, mean, nint_down, nint_up
    use module_pacsobservation, only : PacsObservationSlice, PacsObservation
    use module_pointingmatrix,  only : PointingElement, xy2pmatrix, xy2roi, roi2pmatrix
    use module_projection,      only : convex_hull, surface_convex_polygon
    use module_string,          only : strinteger, strlowcase
    use module_tamasis,         only : tamasis_dir, p, POLICY_KEEP, POLICY_MASK, POLICY_REMOVE, info_time
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
    public :: read_calibration_files
    public :: read_filter_calibration_ncorrelations
    public :: read_filter_calibration
    public :: read_filter_filename
    public :: uv2yz
    public :: yz2ad

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
    character(len=*), parameter :: FILENAME_AI  = 'PCalPhotometer_ArrayInstrument_FM_v6.fits'
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

        logical*1, allocatable :: mask(:,:) ! true if the detector has been filtered out in the flattened list of detectors
        integer, allocatable   :: ij(:,:)
        integer, allocatable   :: pq(:,:)

        logical*1, allocatable :: detector_bad(:)  ! true for bad detectors (returned timeline values are null)
        real(p), allocatable   :: detector_center(:,:)
        real(p), allocatable   :: detector_corner(:,:)
        real(p), allocatable   :: detector_area(:)
        real(p), allocatable   :: flatfield_optical(:)
        real(p), allocatable   :: flatfield_detector(:)
        real(p), allocatable   :: flatfield_total(:)

        real(p)                :: distortion_yz(NDIMS, DISTORTION_DEGREE, DISTORTION_DEGREE, DISTORTION_DEGREE)

    contains

        private
        procedure, public :: init_with_calfiles
        procedure, public :: init_with_variables
        procedure, public :: compute_map_header
        procedure, public :: find_minmax
        procedure, public :: compute_projection
        procedure, public :: destroy
        procedure, public :: read

        procedure, public, nopass :: uv2yz
        procedure, public, nopass :: yz2ad

        procedure :: compute_projection_nearest_neighbour
        procedure :: compute_projection_sharp_edges
        procedure :: read_one
        procedure, public, nopass :: read_calibration_files

    end type PacsInstrument


contains


    subroutine init_with_calfiles(this, band, detector_mask, fine_sampling_factor, status)

        class(PacsInstrument), intent(inout) :: this
        character(len=*), intent(in)         :: band
        logical*1, intent(in)                :: detector_mask(:,:)
        integer, intent(in)                  :: fine_sampling_factor
        integer, intent(out)                 :: status

        logical*1, allocatable :: detector_bad_all(:,:)
        real(p), allocatable   :: detector_center_all(:,:,:)
        real(p), allocatable   :: detector_corner_all(:,:,:,:)
        real(p), allocatable   :: detector_area_all(:,:)
        real(p), allocatable   :: flatfield_optical_all(:,:)
        real(p), allocatable   :: flatfield_detector_all(:,:)
        real(p)                :: distortion_yz(NDIMS,DISTORTION_DEGREE,DISTORTION_DEGREE,DISTORTION_DEGREE)
        real(p)                :: responsivity, active_fraction

        active_fraction = 0

        call this%read_calibration_files(band, detector_bad_all, detector_center_all, detector_corner_all, detector_area_all,      &
             flatfield_optical_all, flatfield_detector_all, distortion_yz, responsivity,  active_fraction, status)
        if (status /= 0) return

        call this%init_with_variables(band, detector_mask, fine_sampling_factor, status, detector_bad_all, detector_center_all,    &
             detector_corner_all, detector_area_all, flatfield_optical_all, flatfield_detector_all, distortion_yz)
        if (status /= 0) return

    end subroutine init_with_calfiles


    !-------------------------------------------------------------------------------------------------------------------------------


    ! a flattened array of detector values is travelled through columns and then through rows
    subroutine init_with_variables(this, band, detector_mask, fine_sampling_factor, status, detector_bad_all, detector_center_all, &
                                   detector_corner_all, detector_area_all, flatfield_optical_all, flatfield_detector_all,          &
                                   distortion_yz)

        class(PacsInstrument), intent(inout) :: this
        character(len=*), intent(in)         :: band
        logical*1, intent(in)                :: detector_mask(:,:)
        integer, intent(in)                  :: fine_sampling_factor
        integer, intent(out)                 :: status
        logical*1, intent(in), optional      :: detector_bad_all(:,:)
        real(p), intent(in), optional        :: detector_center_all(:,:,:)
        real(p), intent(in), optional        :: detector_corner_all(:,:,:,:)
        real(p), intent(in), optional        :: detector_area_all(:,:)
        real(p), intent(in), optional        :: flatfield_optical_all(:,:)
        real(p), intent(in), optional        :: flatfield_detector_all(:,:)
        real(p), intent(in), optional        :: distortion_yz(:,:,:,:)

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

        allocate (this%detector_bad(this%ndetectors))
        if (present(detector_bad_all)) then
            if (any(shape(detector_bad_all) /= shape(detector_mask))) then
                write (ERROR_UNIT,'(a)') 'Error: Invalid shape for the bad detector mask.'
                return
            end if
        else
            this%detector_bad = .false.
        end if

        if (present(detector_center_all)) then
            if (any(shape(detector_center_all) /= [NDIMS,shape(detector_mask)])) then
                write (ERROR_UNIT,'(a)') 'Error: Invalid shape for the detector center array.'
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

        idetector = 1

        do ip = 1, nrows

            do iq = 1, ncolumns
                if (this%mask(ip,iq)) cycle
                this%pq(1,idetector) = ip-1
                this%pq(2,idetector) = iq-1
                this%ij(1,idetector) = mod(ip-1, 16)
                this%ij(2,idetector) = mod(iq-1, 16)
                if (present(detector_bad_all)) then
                    this%detector_bad(idetector) = detector_bad_all(ip,iq)
                end if
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
 

    subroutine read_calibration_files(band, detector_bad_all, detector_center_all, detector_corner_all, detector_area_all,         &
                                      flatfield_optical_all, flatfield_detector_all, distortion_yz, responsivity, active_fraction, &
                                      status)

        character(len=*), intent(in)        :: band
        logical*1, intent(out), allocatable :: detector_bad_all(:,:)
        real(p), intent(out), allocatable   :: detector_center_all(:,:,:)
        real(p), intent(out), allocatable   :: detector_corner_all(:,:,:,:)
        real(p), intent(out), allocatable   :: detector_area_all(:,:)
        real(p), intent(out), allocatable   :: flatfield_optical_all(:,:)
        real(p), intent(out), allocatable   :: flatfield_detector_all(:,:)
        real(p), intent(out)                :: distortion_yz(NDIMS,DISTORTION_DEGREE,DISTORTION_DEGREE,DISTORTION_DEGREE)
        real(p), intent(out)                :: responsivity
        real(p), intent(inout)              :: active_fraction
        integer, intent(out)                :: status

        
        ! hdu extension for the blue, green and red band, in this order
        integer, parameter   :: HDU_FLATFIELD(3)    = [12, 7, 2]                                    ! v3
        integer, parameter   :: HDU_RESPONSIVITY(3) = [6, 4, 2]                                     ! v5
        integer, parameter   :: HDU_CORNER(4,2)     = reshape([7, 11, 15, 19, 5, 9, 13, 17], [4,2]) ! v5

        integer              :: iband, ichannel, ip, iq, iv, ncolumns, nrows, unit
        character(len=4)     :: channel
        real(p), allocatable :: tmp1(:)
        real(p), allocatable :: tmp2(:,:)
        real(p), allocatable :: tmp3(:,:,:)
        real(p), allocatable :: flatfield_total_all(:,:)
        real(p)              :: center(2,2)
        real(p)              :: calib_active_fraction, coef
        logical*1, allocatable :: tmpl(:,:)

        status = 1

        if (band /= 'blue' .and. band /= 'green' .and. band /= 'red') then
            write (ERROR_UNIT,'(a)') "Error: Invalid band '" // band // "'. Expected values are 'blue', 'green' or 'red'."
            return
        end if

        ! 1, 2, 3 for blue, green and red band
        iband = index('bgr', band(1:1))

        if (band == 'red') then
            nrows    = SHAPE_RED(1)
            ncolumns = SHAPE_RED(2)
            ichannel = 2
            channel  = 'red'
        else
            nrows    = SHAPE_BLUE(1)
            ncolumns = SHAPE_BLUE(2)
            ichannel = 1
            channel  = 'blue'
        end if

        allocate (detector_bad_all      (nrows,ncolumns))
        allocate (detector_center_all   (NDIMS,nrows,ncolumns))
        allocate (detector_corner_all   (NDIMS,NVERTICES,nrows,ncolumns))
        allocate (detector_area_all     (nrows,ncolumns))
        allocate (flatfield_optical_all (nrows,ncolumns))
        allocate (flatfield_detector_all(nrows,ncolumns))
        allocate (flatfield_total_all   (nrows,ncolumns))
        
        ! bad detectors
        call ft_read_image(get_calfile(FILENAME_BPM) // '[' // band // ']', tmpl, status)
        if (status /= 0) return
        detector_bad_all = transpose(tmpl)

        ! UV centers
        call ft_read_image(get_calfile(FILENAME_SAA) // '[u' // trim(channel) // ']', tmp2, status)
        if (status /= 0) return
        detector_center_all(1,:,:) = transpose(tmp2)
        call ft_read_image(get_calfile(FILENAME_SAA) // '[v' // trim(channel) // ']', tmp2, status)
        if (status /= 0) return
        detector_center_all(2,:,:) = transpose(tmp2)

        ! UV corners
        do iv=1, NVERTICES
            call ft_read_image(get_calfile(FILENAME_SAA) // '+' // strinteger(HDU_CORNER(iv,ichannel)), tmp2, status)
            if (status /= 0) return
            detector_corner_all(1,iv,:,:) = transpose(tmp2)
            call ft_read_image(get_calfile(FILENAME_SAA) // '+' // strinteger(HDU_CORNER(iv,ichannel)+1), tmp2, status)
            if (status /= 0) return
            detector_corner_all(2,iv,:,:) = transpose(tmp2)
            
        end do

        ! active fraction
        call ft_open(get_calfile(FILENAME_SAA), unit, status)
        if (status /= 0) return

        call ft_read_keyword_hcss(unit, 'activeFraction', calib_active_fraction, status=status)
        if (status /= 0) return

        call ft_close(unit, status)
        if (status /= 0) return

        ! distortion coefficients in the (y,z) plane
        call ft_read_image(get_calfile(FILENAME_AI) // '[ycoeff' // trim(channel) // ']', tmp3, status)
        if (status /= 0) return
        distortion_yz(1,:,:,:) = tmp3
        call ft_read_image(get_calfile(FILENAME_AI) // '[zcoeff' // trim(channel) // ']', tmp3, status)
        if (status /= 0) return
        distortion_yz(2,:,:,:) = tmp3

        ! detector flat field
        call ft_read_image(get_calfile(FILENAME_FF) // '+' // strinteger(HDU_FLATFIELD(iband)),  tmp2, status)
        if (status /= 0) return
        flatfield_total_all = transpose(tmp2)

        ! responsivity
        call ft_read_image(get_calfile(FILENAME_RES) // '+' // strinteger(HDU_RESPONSIVITY(iband)), tmp1, status)
        if (status /= 0) return
        responsivity = tmp1(1)

        ! modify detector size if active_fraction is not set to zero
        if (active_fraction /= 0) then
            coef = sqrt(active_fraction/calib_active_fraction)
            do ip = 1, nrows
                do iq = 1, ncolumns
                    do iv = 1, NVERTICES
                        detector_corner_all(:,iv,ip,iq) = (detector_corner_all(:,iv,ip,iq) - detector_center_all(:,ip,iq)) *       &
                                                           coef + detector_center_all(:,ip,iq)
                    end do
                end do
            end do
        else
            active_fraction = calib_active_fraction
        end if

        ! detector area
        do ip = 1, nrows
            do iq = 1, ncolumns
                detector_area_all(ip,iq) = abs(surface_convex_polygon(real(uv2yz(detector_corner_all(:,:,ip,iq), distortion_yz,    &
                     0._p), p))) * 3600._p**2
            end do
        end do
        
        ! optical and detector flat fields
        center = detector_area_all(nrows/2:nrows/2+1,ncolumns/2:ncolumns/2+1)
        flatfield_optical_all  = detector_area_all / mean(reshape(center, [4]))
        flatfield_detector_all = flatfield_total_all / flatfield_optical_all
        where (detector_bad_all)
            flatfield_detector_all = NaN
        end where

    end subroutine read_calibration_files


    !-------------------------------------------------------------------------------------------------------------------------------


    function get_calfile(filename)

        character(len=*), intent(in)                    :: filename
        character(len=len(tamasis_dir)+6+len(filename)) :: get_calfile

        get_calfile = tamasis_dir // '/pacs/' // filename

    end function get_calfile


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

        call ft_create_header(0, 0, [[-resolution/3600._p, 0._p], [0._p, resolution/3600._p]], ra0, dec0, 1._p, 1._p, header)

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

        real(p)              :: ra, dec, pa, chop, offset
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
               offset = obs%slice(islice)%offset
            else
               sampling_factor = 1
               offset = 0
            end if

#ifndef IFORT
            !$omp parallel do default(shared) reduction(min:xmin,ymin) reduction(max:xmax,ymax) &
            !$omp private(itime,isample,ra,dec,pa,chop,hull)
#endif
            do itime = 1, obs%slice(islice)%nsamples * sampling_factor

                isample = (itime - 1) / sampling_factor + 1
                if (obs%slice(islice)%p(isample)%removed) cycle

                call obs%get_position_index(islice, itime, sampling_factor, offset, ra, dec, pa, chop)
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


    subroutine compute_projection(this, method, obs, oversampling, header, nx, ny, pmatrix, npixels_per_sample, status)

        class(PacsInstrument), intent(in)  :: this
        integer, intent(in)                :: method
        class(PacsObservation), intent(in) :: obs
        logical, intent(in)                :: oversampling
        character(len=*), intent(in)       :: header
        integer, intent(in)                :: nx, ny
        type(PointingElement), intent(out) :: pmatrix(:,:,:)
        integer, intent(out)               :: npixels_per_sample
        integer, intent(out)               :: status

        integer :: count_start
        logical :: out

        call system_clock(count_start)

        call init_astrometry(header, status=status)
        if (status /= 0) return

        select case (method)

            case (NEAREST_NEIGHBOUR)
                call this%compute_projection_nearest_neighbour(obs, oversampling, nx, ny, pmatrix, npixels_per_sample, out)
            
            case (SHARP_EDGES)
               call this%compute_projection_sharp_edges(obs, oversampling, nx, ny, pmatrix, npixels_per_sample, out)

        end select

        if (size(pmatrix,1) > 0) then
            call info_time('Computing the projector', count_start)
        else
            call info_time('Computing the projector size', count_start)
        end if

        if (npixels_per_sample > size(pmatrix,1) .and. size(pmatrix,1) > 0) then
            status = 1
            write (ERROR_UNIT,'(a,i0,a)') 'Error: Please update npixels_per_sample to ', npixels_per_sample, '.'
        else if (npixels_per_sample < size(pmatrix,1) .or. size(pmatrix,1) == 0) then
            write (OUTPUT_UNIT,'(a,i0,a)') "Warning: For this observation, you can set the keyword 'npixels_per_sample' to ",      &
                  npixels_per_sample, ' for better performances.'
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
        real(p) :: ra, dec, pa, chop, chop_old, reference_area, offset
        integer :: ifine, isample, islice, itime, ivalid, nsamples, nvalids, sampling_factor, dest
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
               offset = obs%slice(islice)%offset
            else
               sampling_factor = 1
               offset = 0
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
                     call obs%get_position_index(islice, itime, sampling_factor, offset, ra, dec, pa, chop)
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
        real(p) :: ra, dec, pa, chop, chop_old, offset
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
               offset = obs%slice(islice)%offset
            else
               sampling_factor = 1
               offset = 0
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
                     call obs%get_position_index(islice, itime, sampling_factor, offset, ra, dec, pa, chop)
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


    subroutine read(this, obs, selected_masks, signal, mask, status, verbose)

        class(PacsInstrument), intent(in)  :: this
        class(PacsObservation), intent(in) :: obs
        character(len=*), intent(in)       :: selected_masks(:)
        real(p), intent(out)               :: signal(:,:)
        logical(1), intent(out)            :: mask(:,:)
        integer, intent(out)               :: status
        logical, intent(in), optional      :: verbose

        integer :: dest
        integer :: iobs, nobs
        integer :: count_start
        logical :: verbose_

        verbose_ = .false.
        if (present(verbose)) verbose_ = verbose

        call system_clock(count_start)

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
            call this%read_one(obs%slice(iobs), selected_masks, signal(dest:dest+obs%slice(iobs)%nvalids-1,:),                     &
                               mask(dest:dest+obs%slice(iobs)%nvalids-1,:), status)
            if (status /= 0) return
            dest = dest + obs%slice(iobs)%nvalids
        end do

        if (verbose_) then
            call info_time('Reading timeline', count_start)
        end if

    end subroutine read


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine read_one(this, obs, selected_mask, signal, mask, status)

        class(PacsInstrument), intent(in)       :: this
        class(PacsObservationSlice), intent(in) :: obs
        character(len=*), intent(in)            :: selected_mask(:)
        real(p), intent(out)                    :: signal(:,:)
        logical*1, intent(out)                  :: mask  (:,:)
        integer, intent(out)                    :: status

        integer                   :: first, last
        integer                   :: ip, iq
        integer                   :: idetector, ndetectors, unit
        integer, allocatable      :: signal_shape(:), mask_shape(:)
        integer                   :: imask, isample, icompressed, ibit
        integer                   :: firstcompressed, lastcompressed
        real(p), allocatable      :: signal_(:)
        logical*1, allocatable    :: mask_(:)
        integer, allocatable      :: mask_extension(:)
        character(len=FLEN_VALUE) :: mask_name
        integer                   :: nmasks
        integer*4, allocatable    :: maskcompressed(:)
        integer*4                 :: maskval
        logical                   :: mask_found
        integer                   :: hdutype, naxis1
        integer                   :: ncompressed

        ndetectors = size(signal, 2)

        ! get first and last valid pointing, to decrease I/O
        do first = 1, obs%nsamples
            if (.not. obs%p(first)%removed) exit
        end do
        do last = obs%nsamples, 1, -1
            if (.not. obs%p(last)%removed) exit
        end do
        
        ! test if all samples are removed in this slice
        if (first > last) then
            status = 0
            return
        end if

        ! set mask from policy
        mask = spread(pack(obs%p(first:last)%masked, .not. obs%p(first:last)%removed), 2, ndetectors)

        ! read signal HDU
        call ft_open_image(trim(obs%filename) // '[Signal]', unit, 3, signal_shape, status)
        if (status /= 0) return

        allocate (mask_shape(size(signal_shape)))
        mask_shape(1)  = (signal_shape(1) - 1) / 32 + 1
        mask_shape(2:) = signal_shape(2:)

        allocate (signal_(last - first + 1), mask_(last - first + 1))
        do idetector = 1, ndetectors

            ip = this%pq(1,idetector)
            iq = this%pq(2,idetector)

            call ft_read_slice(unit, first, last, iq+1, ip+1, signal_shape, signal_, status)
            if (status /= 0) return

            signal(:,idetector) = pack(signal_, .not. obs%p(first:last)%removed)

        end do
        
        call ft_close(unit, status)
        if (status /= 0) return

        ! get the mask tree
        call ft_open(trim(obs%filename) // '[MASK]', unit, mask_found, status)
        if (status /= 0 .or. .not. mask_found) return

        call ft_read_keyword(unit, 'DSETS___', nmasks, status=status)
        if (status /= 0) return

        allocate (mask_extension(nmasks))
        do imask=1, nmasks
            call ft_read_keyword(unit, 'DS_' // strinteger(imask-1), mask_extension(imask), status=status)
            if (status /= 0) return
        end do
        
        firstcompressed = (first - 1) / 32 + 1
        lastcompressed  = (last  - 1) / 32 + 1
        ncompressed     = lastcompressed - firstcompressed + 1
        allocate (maskcompressed(ncompressed))

        ! loop over all the masks in the file
        do imask=1, nmasks

            call FTMAHD(unit, mask_extension(imask)+1, hdutype, status)
            if (ft_check_error_cfitsio(status, unit, trim(obs%filename))) return

            call ft_read_keyword(unit, 'EXTNAME', mask_name, status=status)
            if (status /= 0) return

            ! check if the mask is selected
            if (all(strlowcase(mask_name) /= selected_mask)) cycle

            ! check mask dimension
            call ft_read_keyword(unit, 'NAXIS1', naxis1, status=status)
            if (status /= 0) return            
            if (mask_shape(1) /= naxis1) then
                status = 1
                write (ERROR_UNIT, '(a,2(i0,a))') "The observation '" // trim(obs%filename) // "', has an incompatible dimension in&
                       & mask '" // trim(mask_name) // "' (", naxis1, ' instead of ', mask_shape(1), ').'
                return
            end if

            ! read mask
            do idetector = 1, size(mask,2)

                ip = this%pq(1,idetector)
                iq = this%pq(2,idetector)

                mask_ = .false.

                call ft_read_slice(unit, firstcompressed, lastcompressed, iq+1, ip+1, mask_shape, maskcompressed, status)
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

        end do

        call ft_close(unit, status)

        ! apply bad detector mask
        do idetector = 1, size(mask,2)
            if (this%detector_bad(idetector)) mask(:,idetector) = .true.
        end do

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
        deallocate (this%detector_bad)
        deallocate (this%detector_center)
        deallocate (this%detector_corner)
        deallocate (this%detector_area)

    end subroutine destroy


end module module_pacsinstrument
