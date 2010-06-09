module module_pacsinstrument

    use iso_fortran_env,        only : ERROR_UNIT, OUTPUT_UNIT
    use module_fitstools,       only : ft_close, ft_create_header, ft_open_image, ft_read_extension, ft_read_slice,ft_test_extension
    use module_math,            only : DEG2RAD, pInf, mInf, NaN, mean, nint_down, nint_up
    use module_pacsobservation, only : pacsobservationslice, pacsobservation
    use module_pointingmatrix,  only : pointingelement, xy2roi, roi2pmatrix
    use module_precision,       only : p
    use module_projection,      only : convex_hull, surface_convex_polygon
    use module_string,          only : strinteger
    use module_tamasis,         only : get_tamasis_path, tamasis_path_len
    use module_wcs,             only : init_astrometry, ad2xy_gnomonic, ad2xy_gnomonic_vect
    implicit none
    private

    public :: ndims
    public :: nvertices
    public :: pacsinstrument
    public :: multiplexing_direct
    public :: multiplexing_transpose

    integer, parameter :: ndims = 2
    integer, parameter :: nvertices = 4
    integer, parameter :: shape_blue(2) = [32, 64]
    integer, parameter :: shape_red (2) = [16, 32]
    integer, parameter :: distortion_degree = 3
    real*8,  parameter :: sampling = 0.025d0

    type pacsinstrument
 
        logical*1            :: mask_blue (shape_blue(1), shape_blue(2))
        logical*1            :: mask_green(shape_blue(1), shape_blue(2))
        logical*1            :: mask_red  (shape_red (1), shape_red (2))

        real*8               :: flatfield_blue (shape_blue(1), shape_blue(2))
        real*8               :: flatfield_green(shape_blue(1), shape_blue(2))
        real*8               :: flatfield_red  (shape_red (1), shape_red (2))

        integer              :: ndetectors
        integer              :: nrows
        integer              :: ncolumns
        integer              :: fine_sampling_factor
        character            :: channel
        logical              :: transparent_mode
        logical              :: keep_bad_detectors
        logical              :: mask_bad_line

        logical*1, allocatable :: mask(:,:) ! .true. if the detector has been filtered out in the flattened list of detectors
        logical*1, allocatable :: bad(:,:)  ! .true. if the detector is flagged as bad (same as mask if keep_bad_detector is false)
        integer, allocatable   :: ij(:,:)
        integer, allocatable   :: pq(:,:)

        real*8, allocatable  :: flatfield_total(:,:)
        real*8, allocatable  :: flatfield_detector(:,:)
        real*8, allocatable  :: flatfield_optical(:,:)
        real*8, allocatable  :: detector_area(:,:)

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

        private
        procedure, public :: init
        procedure, public :: compute_map_header
        procedure, public :: find_minmax
        procedure, public :: compute_projection_sharp_edges
        procedure, public :: destroy
        procedure, public :: read

        procedure, nopass, public :: uv2yz
        procedure, nopass, public :: yz2ad

        procedure :: read_one
        procedure :: read_oldstyle
        procedure :: read_calibration_files
        procedure :: filter_detectors
        procedure :: filter_detectors_array

    end type pacsinstrument


contains


    subroutine init(this, channel, transparent_mode, fine_sampling_factor, keep_bad_detectors, mask_bad_line, status,              &
                    bad_detector_mask)
        class(pacsinstrument), intent(inout) :: this
        character, intent(in)                :: channel
        logical, intent(in)                  :: transparent_mode
        integer, intent(in)                  :: fine_sampling_factor
        logical, intent(in)                  :: keep_bad_detectors
        logical, intent(in)                  :: mask_bad_line
        integer, intent(out)                 :: status
        logical*1, intent(in), optional      :: bad_detector_mask(:,:)

        ! check channel
        if (index('rgb', channel) == 0) then
            status = 1
            write (ERROR_UNIT,'(a)') "ERROR: Invalid channel '" // channel // "'. Valid values are 'r', 'g' and 'b'."
            return
        end if

        ! check fine sampling factor
        if (all(fine_sampling_factor /= [1,2,4,8,16,32])) then
            status = 1
            write (ERROR_UNIT,'(a)') "ERROR: Invalid sampling factor '" //     &
                strinteger(fine_sampling_factor) // "'. Valid values are '1', '2', '4', '8', '16' or '32'."
            return
        end if

        this%channel              = channel
        this%fine_sampling_factor = fine_sampling_factor
        this%transparent_mode     = transparent_mode
        this%keep_bad_detectors   = keep_bad_detectors
        this%mask_bad_line        = mask_bad_line

        call this%read_calibration_files(status)
        if (status /= 0) return

        if (present(bad_detector_mask)) then

            select case (this%channel)
                case ('b')
                    if (any(shape(bad_detector_mask) /= shape_blue)) then
                        status = 1
                        write (ERROR_UNIT,'(a)') 'INIT: input blue bad detector mask has invalid size.'
                        return
                    end if
                    this%mask_blue = bad_detector_mask
                case ('g')
                    if (any(shape(bad_detector_mask) /= shape_blue)) then
                        status = 1
                        write (ERROR_UNIT,'(a)') 'INIT: input green bad detector mask has invalid size.'
                        return
                    end if
                    this%mask_green = bad_detector_mask
                case ('r')
                    if (any(shape(bad_detector_mask) /= shape_red)) then
                        status = 1
                        write (ERROR_UNIT,'(a)') 'INIT: input red bad detector mask has invalid size.'
                        return
                    end if
                    this%mask_red = bad_detector_mask
                case default
                    status = 1
                    write (ERROR_UNIT,'(a)') "INIT: invalid channel: '" // this%channel // "'."
            end select

        end if
        
        call this%filter_detectors()

    end subroutine init


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine read_calibration_files(this, status)

        class(pacsinstrument), intent(inout) :: this
        integer, intent(out)                 :: status

        ! these calibration files are taken from HCSS 4.0.70
        character(len=*), parameter :: filename_saa = 'PCalPhotometer_SubArrayArray_FM_v5.fits'
        character(len=*), parameter :: filename_ai  = 'PCalPhotometer_ArrayInstrument_FM_v5.fits'
        character(len=*), parameter :: filename_bpm = 'PCalPhotometer_BadPixelMask_FM_v5.fits'
        character(len=*), parameter :: filename_ff  = 'PCalPhotometer_FlatField_FM_v3.fits'

        integer,          parameter :: hdu_blue(4) = [8, 12, 16, 20]
        integer,          parameter :: hdu_red (4) = [6, 10, 14, 18]

        integer                     :: ivertex
        logical*1, allocatable      :: tmplogical(:,:)
        real*8, allocatable         :: tmp2(:,:)
        real*8, allocatable         :: tmp3(:,:,:)

        ! read bad pixel mask

        call ft_read_extension(get_calfile(filename_bpm) // '[blue]', tmplogical, status)
        if (status /= 0) return
        this%mask_blue = transpose(tmplogical)

        call ft_read_extension(get_calfile(filename_bpm) // '[green]', tmplogical, status)
        if (status /= 0) return
        this%mask_green = transpose(tmplogical)

        call ft_read_extension(get_calfile(filename_bpm) // '[red]', tmplogical, status)
        if (status /= 0) return
        this%mask_red = transpose(tmplogical)


        ! read flat fields

        call ft_read_extension(get_calfile(filename_ff) // '+12',  tmp2, status)
        if (status /= 0) return
        this%flatfield_blue = transpose(tmp2)
        call ft_read_extension(get_calfile(filename_ff) // '+7', tmp2, status)
        if (status /= 0) return
        this%flatfield_green = transpose(tmp2)
        call ft_read_extension(get_calfile(filename_ff) // '+2',   tmp2, status)
        if (status /= 0) return
        this%flatfield_red = transpose(tmp2)


        ! read detector corners in the (u,v) plane

        do ivertex=1, nvertices

            call ft_read_extension(get_calfile(filename_saa), tmp2, status, hdu_blue(ivertex)  )
            if (status /= 0) return
            this%corners_uv_blue(1,ivertex,:,:) = transpose(tmp2)
            call ft_read_extension(get_calfile(filename_saa), tmp2, status, hdu_blue(ivertex)+1)
            if (status /= 0) return
            this%corners_uv_blue(2,ivertex,:,:) = transpose(tmp2)
            call ft_read_extension(get_calfile(filename_saa), tmp2, status, hdu_red (ivertex)  )
            if (status /= 0) return
            this%corners_uv_red (1,ivertex,:,:) = transpose(tmp2)
            call ft_read_extension(get_calfile(filename_saa), tmp2, status, hdu_red (ivertex)+1)
            if (status /= 0) return
            this%corners_uv_red (2,ivertex,:,:) = transpose(tmp2)

        end do


        ! read the distortion coefficients in the (y,z) plane

        call ft_read_extension(get_calfile(filename_ai) // '[ycoeffblue]', tmp3, status)
        if (status /= 0) return
        this%distortion_yz_blue(1,:,:,:) = tmp3
        call ft_read_extension(get_calfile(filename_ai) // '[zcoeffblue]', tmp3, status)
        if (status /= 0) return
        this%distortion_yz_blue(2,:,:,:) = tmp3
        call ft_read_extension(get_calfile(filename_ai) // '[ycoeffred]', tmp3, status)
        if (status /= 0) return
        this%distortion_yz_red (1,:,:,:) = tmp3
        call ft_read_extension(get_calfile(filename_ai) // '[zcoeffred]', tmp3, status)
        if (status /= 0) return
        this%distortion_yz_red (2,:,:,:) = tmp3

    end subroutine read_calibration_files


    !-------------------------------------------------------------------------------------------------------------------------------


    function get_calfile(filename)
        character(len=*), intent(in)                     :: filename
        character(len=tamasis_path_len+10+len(filename)) :: get_calfile

        get_calfile = get_tamasis_path() // 'pacs/data/' // filename

    end function get_calfile


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine filter_detectors(this)

        class(pacsinstrument), intent(inout)   :: this

        select case (this%channel)
           case ('r')
              call this%filter_detectors_array(this%mask_red,   this%corners_uv_red,  this%distortion_yz_red,  this%flatfield_red)
           case ('g')
              call this%filter_detectors_array(this%mask_green, this%corners_uv_blue, this%distortion_yz_blue, this%flatfield_green)
           case ('b')
              call this%filter_detectors_array(this%mask_blue,  this%corners_uv_blue, this%distortion_yz_blue, this%flatfield_blue)
        end select

    end subroutine filter_detectors


    !-------------------------------------------------------------------------------------------------------------------------------


    ! a flatten array of detector values is travelled through columns and then through rows
    subroutine filter_detectors_array(this, mask, uv, distortion, flatfield)

        class(pacsinstrument), intent(inout) :: this
        logical*1, intent(in), target        :: mask(:,:)
        real*8, intent(in)                   :: uv(:,:,:,:)
        real*8, intent(in)                   :: distortion(:,:,:,:)
        real*8, intent(in)                   :: flatfield(:,:)

        integer   :: idetector, p, q
        real*8    :: center(2,2)

        this%nrows    = size(mask, 1)
        this%ncolumns = size(mask, 2)

        allocate(this%mask(this%nrows, this%ncolumns))
        allocate(this%bad (this%nrows, this%ncolumns))
        this%bad = mask

        ! mask detectors rejected in transparent mode
        if (this%transparent_mode) then
            if (this%channel /= 'r') then
                this%bad(1:16,1:16) = .true.
                this%bad(1:16,33:)  = .true.
                this%bad(17:,:)     = .true.
            else
                this%bad(1:8,1:8) = .true.
                this%bad(1:8,17:) = .true.
                this%bad(9:,:)    = .true.
            end if
        end if

        ! mask erratic line
        if (this%mask_bad_line .and. this%channel /= 'r') then
            this%bad(12,17:32) = .true.
        end if
        
        if (this%keep_bad_detectors) then
            this%mask = .false.
        else
            this%mask = this%bad
        end if

        ! get the number of detectors
        this%ndetectors = count(.not. this%mask)

        allocate (this%ij(ndims, this%ndetectors))
        allocate (this%pq(ndims, this%ndetectors))
        allocate (this%corners_uv(ndims, nvertices * this%ndetectors))
        allocate (this%corners_yz(ndims, nvertices * this%ndetectors))
        allocate (this%flatfield_total(this%nrows, this%ncolumns))
        allocate (this%flatfield_detector(this%nrows, this%ncolumns))
        allocate (this%flatfield_optical(this%nrows, this%ncolumns))
        allocate (this%detector_area(this%nrows, this%ncolumns))
        this%distortion_yz = distortion

        idetector = 1

        do p = 1, this%nrows
            do q = 1, this%ncolumns
                this%detector_area(p,q) = abs(surface_convex_polygon(this%uv2yz(uv(:,:,p,q), distortion, 0.d0))) * 3600.d0**2
                if (this%mask(p,q)) cycle
                this%pq(1, idetector) = p-1
                this%pq(2, idetector) = q-1
                this%ij(1, idetector) = mod(p-1, 16)
                this%ij(2, idetector) = mod(q-1, 16)
                this%corners_uv(:,nvertices * (idetector-1)+1:nvertices*idetector) = uv(:,:,p,q)
                idetector = idetector + 1
            end do
        end do
        
        center = this%detector_area(this%nrows/2:this%nrows/2+1,this%ncolumns/2:this%ncolumns/2)
        this%flatfield_optical = this%detector_area / mean(reshape(center,[4]))
        this%flatfield_total = flatfield
        this%flatfield_detector = this%flatfield_total / this%flatfield_optical

        where (mask)
            this%flatfield_total    = NaN
            this%flatfield_detector = NaN
            this%flatfield_optical  = NaN
        end where

    end subroutine filter_detectors_array


    !-------------------------------------------------------------------------------------------------------------------------------


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
 

    !-------------------------------------------------------------------------------------------------------------------------------


    function yz2ad(yz, ra0, dec0, pa0) result (ad)

        real*8, intent(in) :: yz(:,:) ! in degrees
        real*8, intent(in) :: ra0, dec0, pa0
        real*8             :: ad(ndims,size(yz,2))

        integer            :: i
        real*8             :: cospa, sinpa


        cospa =  cos(pa0*DEG2RAD)
        sinpa = -sin(pa0*DEG2RAD)

        do i=1, size(yz, 2)

            ad(2,i) = dec0 + (yz(1,i) * sinpa + yz(2,i) * cospa)
            ad(1,i) = ra0  + (yz(1,i) * cospa - yz(2,i) * sinpa) / cos(ad(2,i) * DEG2RAD)

        end do

    end function yz2ad


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine compute_map_header(this, obs, oversampling, resolution, header, status)

        class(pacsinstrument), intent(in)  :: this
        class(pacsobservation), intent(in) :: obs
        logical, intent(in)                :: oversampling
        real*8, intent(in)                 :: resolution
        character(len=2880), intent(out)   :: header
        integer, intent(out)               :: status

        integer :: nx, ny
        integer :: ixmin, ixmax, iymin, iymax
        real*8  :: ra0, dec0, xmin, xmax, ymin, ymax


        call obs%compute_center(ra0, dec0)

        call ft_create_header(0, 0, [[-resolution/3600.d0, 0.d0], [0.d0, resolution/3600.d0]], ra0, dec0, 1.d0, 1.d0, header)

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
        call ft_create_header(nx, ny, [[-resolution/3600.d0, 0.d0], [0.d0, resolution/3600.d0]], ra0, dec0, -ixmin+2.d0,           &
             -iymin+2.d0, header)

    end subroutine compute_map_header


    !-------------------------------------------------------------------------------------------------------------------------------


    ! find minimum and maximum pixel coordinates in maps
    subroutine find_minmax(this, obs, oversampling, xmin, xmax, ymin, ymax)
        class(pacsinstrument), intent(in)  :: this
        class(pacsobservation), intent(in) :: obs
        logical, intent(in)                :: oversampling
        real*8, intent(out)                :: xmin, xmax, ymin, ymax
        real*8               :: ra, dec, pa, chop
        real*8, allocatable  :: hull_uv(:,:), hull(:,:)
        integer, allocatable :: ihull(:)
        integer              :: sampling_factor, islice, itime

        xmin = pInf
        xmax = mInf
        ymin = pInf
        ymax = mInf

        call convex_hull(this%corners_uv, ihull)

        allocate(hull_uv(ndims,size(ihull)))
        allocate(hull(ndims, size(ihull)))
        hull_uv = this%corners_uv(:,ihull)

        do islice = 1, obs%nslices

            ! check if it is required to interpolate pointing positions
            if (oversampling) then
               sampling_factor = this%fine_sampling_factor * obs%slice(islice)%compression_factor
            else
               sampling_factor = 1
            end if

            !XXX bug IFORT
            !!$omp parallel do default(shared) reduction(min:xmin,ymin) &
            !!$omp reduction(max:xmax,ymax) &
            !!$omp private(itime,ra,dec,pa,chop,hull) firstprivate(index)
            do itime = 1, obs%slice(islice)%nsamples * sampling_factor

                call obs%get_position_index(islice, itime, sampling_factor, ra, dec, pa, chop)
                hull = this%uv2yz(hull_uv, this%distortion_yz_blue, chop)
                hull = this%yz2ad(hull, ra, dec, pa)
                hull = ad2xy_gnomonic(hull)
                xmin = min(xmin, minval(hull(1,:)))
                xmax = max(xmax, maxval(hull(1,:)))
                ymin = min(ymin, minval(hull(2,:)))
                ymax = max(ymax, maxval(hull(2,:)))

             end do
             !!$omp end parallel do

        end do

    end subroutine find_minmax


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine compute_projection_sharp_edges(this, obs, oversampling, header, nx, ny, pmatrix, status)
        class(pacsinstrument), intent(in)  :: this
        class(pacsobservation), intent(in) :: obs
        logical, intent(in)                :: oversampling
        character(len=*), intent(in)       :: header
        integer, intent(in)                :: nx, ny
        type(pointingelement), intent(out) :: pmatrix(:,:,:)
        integer, intent(out)               :: status

        real*8  :: coords(ndims,this%ndetectors*nvertices), coords_yz(ndims,this%ndetectors*nvertices)
        real*8  :: ra, dec, pa, chop, chop_old
        integer :: roi(ndims,2,this%ndetectors)
        integer :: nsamples, npixels_per_sample, islice, itime, sampling_factor, dest
        integer :: count1, count2, count_rate, count_max
        logical :: out

        write(*,'(a)', advance='no') 'Info: Computing the projector... '
        call system_clock(count1, count_rate, count_max)

        call init_astrometry(header, status=status)
        if (status /= 0) return

        npixels_per_sample = 0
        dest = 0
        out = .false.

        ! loop over the observations
        do islice = 1, obs%nslices

            nsamples = int(obs%slice(islice)%nsamples, kind=kind(0))
            chop_old = pInf

            ! check if it is required to interpolate pointing positions
            if (oversampling) then
               sampling_factor = this%fine_sampling_factor * obs%slice(islice)%compression_factor
            else
               sampling_factor = 1
            end if

            !$omp parallel do default(shared) firstprivate(chop_old)   &
            !$omp private(itime, ra, dec, pa, chop, coords, coords_yz, roi) &
            !$omp reduction(max : npixels_per_sample) reduction(.or. : out)
            do itime = 1, nsamples * sampling_factor
                call obs%get_position_index(islice, itime, sampling_factor, ra, dec, pa, chop)
                if (abs(chop-chop_old) > 1.d-2) then
                    coords_yz = this%uv2yz(this%corners_uv, this%distortion_yz, chop)
                    chop_old = chop
                 end if
                 coords = this%yz2ad(coords_yz, ra, dec, pa)
                 coords = ad2xy_gnomonic(coords)
                 roi    = xy2roi(coords, nvertices)
                 call roi2pmatrix(roi, nvertices, coords, nx, ny, itime + dest, npixels_per_sample, out, pmatrix)
             end do
             !$omp end parallel do
             dest = dest + nsamples * sampling_factor
        end do

        call system_clock(count2, count_rate, count_max)
        write(*,'(f6.2,a)') real(count2-count1)/count_rate, 's'

        if (npixels_per_sample > size(pmatrix,1)) then
            status = 1
            write(ERROR_UNIT,'(a,i0,a)') 'Error: Please update npixels_per_sample to ', npixels_per_sample, '.'
        else if (npixels_per_sample < size(pmatrix,1)) then
            write(OUTPUT_UNIT,'(a,i0,a)') 'Warning: You may update npixels_per_sample to ', npixels_per_sample, '.'
        end if

        if (out) then
            write (OUTPUT_UNIT,'(a)') 'Warning: Some detectors fall outside the map.'
        end if

    end subroutine compute_projection_sharp_edges


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine read(this, obs, signal, mask, status, verbose)

        class(pacsinstrument), intent(in)  :: this
        class(pacsobservation), intent(in) :: obs
        real(p), intent(out)               :: signal(:,:)
        logical(1), intent(out)            :: mask(:,:)
        integer, intent(out)               :: status
        logical, intent(in), optional      :: verbose

        integer*8 :: nsamples, destination
        integer   :: iobs, idetector, nobs
        integer   :: count1, count2, count_rate, count_max
        logical   :: verbose_

        verbose_ = .false.
        if (present(verbose)) verbose_ = verbose

        if (verbose_) then
            write(*,'(a)', advance='no') 'Info: Reading timeline... '
        end if
        call system_clock(count1, count_rate, count_max)

        status   = 1
        nobs     = obs%nslices
        nsamples = obs%nsamples_tot

        ! check that the total number of samples in the observations is equal
        ! to the number of samples in the signal and mask arrays
        if (nsamples /= size(signal,1) .or. nsamples /= size(mask,1)) then
            write (ERROR_UNIT,'(a)') 'Error: read: invalid dimensions.'
            return
        end if

        ! loop over the PACS observations
        destination = 1
        do iobs = 1, nobs
            call this%read_one(obs%slice(iobs), destination, signal, mask, status)
            if (status /= 0) return
            destination = destination + obs%slice(iobs)%nsamples
        end do

        ! mask bad detectors if they haven't been filtered out
        if (this%keep_bad_detectors) then
            do idetector = 1, this%ndetectors
                if (this%bad(this%pq(1,idetector)+1,this%pq(2,idetector)+1)) then
                    signal(:,idetector) = 0.d0
                    mask  (:,idetector) = .true.
                end if
            end do
        end if

        call system_clock(count2, count_rate, count_max)
        if (verbose_) then
            write(*,'(f6.2,a)') real(count2-count1)/count_rate, 's'
        end if

    end subroutine read


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine read_one(this, obs, destination, signal, mask, status)

        class(pacsinstrument), intent(in)   :: this
!!$        class(observationslice), intent(in) :: obs
        class(pacsobservationslice), intent(in) :: obs
        integer*8, intent(in)               :: destination
        real*8, intent(inout)               :: signal(:,:)
        logical*1, intent(inout)            :: mask  (:,:)
        integer, intent(out)                :: status

        integer*8              :: p, q
        integer                :: idetector, unit, length
        integer                :: status_close
        integer, allocatable   :: imageshape(:)
        logical                :: mask_found
        integer*4, allocatable :: maskcompressed(:)
        integer*4              :: maskval
        integer*8              :: ncompressed
        integer*8              :: isample, icompressed, ibit
        integer*8              :: firstcompressed, lastcompressed

        ! old style file format
        length = len_trim(obs%filename)
        if (obs%filename(length-4:length) /= '.fits') then
            call this%read_oldstyle(obs, destination, signal, mask, status)
            return
        end if

        ! read signal HDU
        call ft_open_image(trim(obs%filename) // '[Signal]', unit, 3, imageshape, status)
        if (status /= 0) return

        do idetector = 1, size(signal,2)

            p = this%pq(1,idetector)
            q = this%pq(2,idetector)

            call ft_read_slice(unit, obs%first, obs%last, q+1, p+1, imageshape, signal(destination:destination+obs%nsamples-1,     &
                 idetector),status)
            if (status /= 0) go to 999

        end do

        call ft_close(unit, status)
        if (status /= 0) return


        ! read Mask HDU
        mask(destination:destination+obs%nsamples-1,:) = .false.
        mask_found = ft_test_extension(trim(obs%filename)//'[Master]', status)
        if (status /= 0) return

        if (.not. mask_found) then
            write (*,'(a)') 'Info: mask Master is not found.'
            go to 100
        end if

        call ft_open_image(trim(obs%filename) // '[Master]', unit, 3, imageshape, status)
        if (status /= 0) return

        allocate(maskcompressed(imageshape(1)))

        do idetector = 1, size(mask,2)

            p = this%pq(1,idetector)
            q = this%pq(2,idetector)

            firstcompressed = (obs%first - 1) / 32 + 1
            lastcompressed  = (obs%last  - 1) / 32 + 1
            ncompressed = lastcompressed - firstcompressed + 1

            call ft_read_slice(unit, firstcompressed, lastcompressed, q+1, p+1, imageshape, maskcompressed(1:ncompressed),status)
            if (status /= 0) go to 999

            ! loop over the bytes of the compressed mask
            do icompressed = firstcompressed, lastcompressed

                maskval = maskcompressed(icompressed-firstcompressed+1)
                if (maskval == 0) cycle

                isample = (icompressed-1)*32 - obs%first + destination + 1

                ! loop over the bits of a compressed mask byte
                do ibit = max(0, obs%first - (icompressed-1)*32-1), min(31, obs%last - (icompressed-1)*32-1)
                    mask(isample+ibit,idetector) = mask(isample+ibit,idetector) .or. btest(maskval,ibit)
                end do

            end do

        end do

        call ft_close(unit, status)
        if (status /= 0) return

    100 do isample = 1, obs%nsamples
            if (obs%p(isample)%invalid) mask(destination+isample-1,:) = .true.
        end do
        return

    999 call ft_close(unit, status_close)

    end subroutine read_one


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine read_oldstyle(this, obs, dest, signal, mask, status)

        class(pacsinstrument), intent(in)   :: this
!!$        class(observationslice), intent(in) :: obs
        class(pacsobservationslice), intent(in) :: obs
        integer*8, intent(in)               :: dest
        real*8, intent(inout)               :: signal(:,:)
        logical*1, intent(inout)            :: mask(:,:)
        integer, intent(out)                :: status

        integer*8            :: p, q
        integer              :: idetector, unit, status_close
        integer, allocatable :: imageshape(:)

        ! handle signal
        call ft_open_image(trim(obs%filename) // '_Signal.fits', unit, 3, imageshape, status)
        if (status /= 0) return

        do idetector = 1, size(signal,2)

            p = this%pq(1,idetector)
            q = this%pq(2,idetector)
            call ft_read_slice(unit, obs%first, obs%last, q+1, p+1, imageshape, signal(dest:dest+obs%nsamples-1,idetector),status)
            if (status /= 0) go to 999

        end do

        call ft_close(unit, status)
        if (status /= 0) return

        ! handle mask
        call ft_open_image(trim(obs%filename) // '_Mask.fits', unit, 3, imageshape, status)
        if (status /= 0) return

        do idetector = 1, size(mask,2)

            p = this%pq(1,idetector)
            q = this%pq(2,idetector)

            call ft_read_slice(unit, obs%first, obs%last, q+1, p+1, imageshape, mask(dest:dest+obs%nsamples-1,idetector), status)
            if (status /= 0) go to 999

        end do

    999 call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine read_oldstyle


    !-------------------------------------------------------------------------------------------------------------------------------


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


    !-------------------------------------------------------------------------------------------------------------------------------


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


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine destroy(this)

        class(pacsinstrument), intent(inout) :: this

        deallocate (this%bad)
        deallocate (this%mask)
        deallocate (this%ij)
        deallocate (this%pq)
        deallocate (this%flatfield_total)
        deallocate (this%flatfield_detector)
        deallocate (this%flatfield_optical)
        deallocate (this%detector_area)
        deallocate (this%corners_uv)
        deallocate (this%corners_yz)

    end subroutine destroy


end module module_pacsinstrument
