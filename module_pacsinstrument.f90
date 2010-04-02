module module_pacsinstrument

    use ISO_FORTRAN_ENV,        only : ERROR_UNIT, OUTPUT_UNIT
    use module_fitstools,       only : ft_close, ft_create_header, ft_open_image, ft_readextension, ft_readslice, ft_test_extension
    use module_math,            only : DEG2RAD, pInf, mInf, NaN, nint_down, nint_up
    use module_pacsobservation, only : pacsobservation, pacsobsinfo
    use module_pacspointing,    only : pacspointing
    use module_pointingmatrix,  only : pointingelement, xy2roi, roi2pmatrix
    use module_projection,      only : convex_hull
    use precision,              only : p
    use string,                 only : strinteger
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
        logical              :: keep_bad_detectors

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

        private
        procedure, public :: init
        procedure, public :: compute_mapheader
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


    subroutine init(this, channel, transparent_mode, fine_sampling_factor, keep_bad_detectors, status, bad_detector_mask)
        class(pacsinstrument), intent(inout) :: this
        character, intent(in)                :: channel
        logical, intent(in)                  :: transparent_mode
        integer, intent(in)                  :: fine_sampling_factor
        logical, intent(in)                  :: keep_bad_detectors
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

        call this%read_calibration_files(status)
        if (status /= 0) return

        if (present(bad_detector_mask)) then
            if (this%channel == 'r') then
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
                if (this%channel == 'g') write (*,*) 'pacs%init: check me.'
                if (size(bad_detector_mask,1) /= size(this%mask_blue,1) .or.   &
                    size(bad_detector_mask,2) /= size(this%mask_blue,2)) then
                   status = 1
                   write (ERROR_UNIT,'(a)') 'INIT: input blue mask has invalid size.'
                   return
                end if
                this%mask_blue = bad_detector_mask
            end if

        end if
        
        call this%filter_detectors()

    end subroutine init


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


    subroutine filter_detectors(this)
        class(pacsinstrument), intent(inout)   :: this

        select case (this%channel)
           case ('r')
              call this%filter_detectors_array(this%mask_red,                  &
                   this%corners_uv_red, this%distortion_yz_red,                &
                   this%flatfield_red)
           case ('g')
              write(*,*) 'XXX FILTER_DETECTORS: Check calibration files for the green band...'
              call this%filter_detectors_array(this%mask_blue,                 &
                   this%corners_uv_blue, this%distortion_yz_blue,              &
                   this%flatfield_green)
           case ('b')
              call this%filter_detectors_array(this%mask_blue,                 &
                   this%corners_uv_blue, this%distortion_yz_blue,              &
                   this%flatfield_blue)

        end select

    end subroutine filter_detectors


    !---------------------------------------------------------------------------


    subroutine filter_detectors_array(this, mask, uv, distortion, flatfield)

        class(pacsinstrument), intent(inout)    :: this
        logical*1, intent(in), target :: mask(:,:)
        real*8, intent(in)            :: uv(:,:,:,:)
        real*8, intent(in)            :: distortion(:,:,:,:)
        real*8, intent(in)            :: flatfield(:,:)

        integer                       :: idetector, p, q

        this%nrows    = size(mask, 1)
        this%ncolumns = size(mask, 2)

        allocate(this%mask(this%nrows, this%ncolumns))
        this%mask => mask

        if (this%transparent_mode) then
            this%mask(1:16,1:16) = .true.
            this%mask(1:16,33:)  = .true.
            this%mask(17:,:)     = .true.
        end if
        this%ndetectors = size(this%mask)
        if (.not. this%keep_bad_detectors) then
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
                if (this%mask(p+1,q+1) .and. .not.this%keep_bad_detectors) cycle
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


        cospa =  cos(pa0*DEG2RAD)
        sinpa = -sin(pa0*DEG2RAD)

        do i=1, size(yz, 2)

            ad(2,i) = dec0 + (yz(1,i) * sinpa + yz(2,i) * cospa)
            ad(1,i) = ra0  + (yz(1,i) * cospa - yz(2,i) * sinpa) / cos(ad(2,i) * DEG2RAD)

        end do

    end function yz2ad


    !---------------------------------------------------------------------------


    subroutine compute_mapheader(this, pointing, finer_sampling, resolution,   &
                                 header, status)
        class(pacsinstrument), intent(in) :: this
        class(pacspointing), intent(in)   :: pointing
        logical, intent(in)               :: finer_sampling
        real*8, intent(in)                :: resolution
        character(len=2880), intent(out)  :: header
        integer, intent(out)              :: status
        integer                           :: nx, ny
        integer                           :: ixmin, ixmax, iymin, iymax
        real*8                            :: ra0, dec0, xmin, xmax, ymin, ymax

        call pointing%compute_center(ra0, dec0)

        call ft_create_header(0, 0, -resolution/3600.d0, resolution/3600.d0, 0.d0, ra0, dec0, 1.d0, 1.d0, header)

        call init_astrometry(header, status=status)
        if (status /= 0) return

        call this%find_minmax(pointing, finer_sampling, xmin, xmax, ymin, ymax)

        ixmin = nint(xmin)
        ixmax = nint(xmax)
        iymin = nint(ymin)
        iymax = nint(ymax)

        nx = ixmax-ixmin+1
        ny = iymax-iymin+1

        ! move the reference pixel (not the reference value!)
        call ft_create_header(nx, ny, -resolution/3600.d0, resolution/3600.d0, &
                              0.d0, ra0, dec0, -ixmin+2.d0, -iymin+2.d0, header)

    end subroutine compute_mapheader


    !---------------------------------------------------------------------------


    ! find minimum and maximum pixel coordinates in maps
    subroutine find_minmax(this, pointing, finer_sampling, xmin, xmax, ymin, ymax)
        class(pacsinstrument), intent(in) :: this
        class(pacspointing), intent(in)   :: pointing
        logical, intent(in)               :: finer_sampling
        real*8, intent(out)               :: xmin, xmax, ymin, ymax
        real*8                            :: ra, dec, pa, chop
        real*8, allocatable               :: hull_uv(:,:), hull(:,:)
        integer, allocatable              :: ihull(:)
        integer                           :: sampling_factor, islice, itime

        xmin = pInf
        xmax = mInf
        ymin = pInf
        ymax = mInf

        call convex_hull(this%corners_uv, ihull)

        allocate(hull_uv(ndims,size(ihull)))
        allocate(hull(ndims, size(ihull)))
        hull_uv = this%corners_uv(:,ihull)

        do islice = 1, pointing%nslices

            ! check if it is required to interpolate pointing positions
            if (finer_sampling) then
               sampling_factor = this%fine_sampling_factor * pointing%compression_factor(islice)
            else
               sampling_factor = 1
            end if

            !XXX bug IFORT
            !!$omp parallel do default(shared) reduction(min:xmin,ymin) &
            !!$omp reduction(max:xmax,ymax) &
            !!$omp private(itime,ra,dec,pa,chop,hull) firstprivate(index)
            do itime = 1, pointing%nsamples(islice) * sampling_factor

                call pointing%get_position_index(islice, itime, sampling_factor, ra, dec, pa, chop)
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


    !---------------------------------------------------------------------------


    subroutine compute_projection_sharp_edges(this, pointing, finer_sampling, header, nx, ny, pmatrix, status)
        class(pacsinstrument), intent(in)  :: this
        class(pacspointing), intent(in)    :: pointing
        logical, intent(in)                :: finer_sampling
        character(len=*), intent(in)       :: header
        integer, intent(in)                :: nx, ny
        type(pointingelement), intent(out) :: pmatrix(:,:,:)
        integer, intent(out)               :: status

        real*8  :: coords(ndims,this%ndetectors*nvertices), coords_yz(ndims,this%ndetectors*nvertices)
        real*8  :: ra, dec, pa, chop, chop_old
        integer :: roi(ndims,2,this%ndetectors)
        integer :: nsamples,  npixels_per_sample, islice, itime, sampling_factor, nroi, dest

        call init_astrometry(header, status=status)
        if (status /= 0) return

        npixels_per_sample = -1
        dest = 0

        ! loop over the observations
        do islice = 1, pointing%nslices

            nsamples = int(pointing%nsamples(islice), kind=kind(0))
            chop_old = pInf

            ! check if it is required to interpolate pointing positions
            if (finer_sampling) then
               sampling_factor = this%fine_sampling_factor * pointing%compression_factor(islice)
            else
               sampling_factor = 1
            end if

            !$omp parallel do default(shared) firstprivate(chop_old)   &
            !$omp private(itime, ra, dec, pa, chop, coords, coords_yz, roi) &
            !$omp reduction(max : npixels_per_sample)
            do itime = 1, nsamples * sampling_factor
                call pointing%get_position_index(islice, itime, sampling_factor, ra, dec, pa, chop)
                if (abs(chop-chop_old) > 1.d-2) then
                    coords_yz = this%uv2yz(this%corners_uv, this%distortion_yz, chop)
                    chop_old = chop
                 end if
                 coords = this%yz2ad(coords_yz, ra, dec, pa)
                 coords = ad2xy_gnomonic(coords)
                 roi    = xy2roi(coords, nvertices)
                 call roi2pmatrix(roi, nvertices, coords, nx, ny, itime + dest, nroi, pmatrix)
                 npixels_per_sample = max(npixels_per_sample, nroi)
             end do
             !$omp end parallel do
             dest = dest + nsamples * sampling_factor
        end do

        if (npixels_per_sample > size(pmatrix,1)) then
            status = 1
            write(*,'(a,i0,a)') 'Error: Please update npixels_per_sample to ', npixels_per_sample, '.'
        else if (npixels_per_sample < size(pmatrix,1)) then
            write(*,'(a,i0,a)') 'Warning: You may update npixels_per_sample to ', npixels_per_sample, '.'
        end if

    end subroutine compute_projection_sharp_edges


    !---------------------------------------------------------------------------


    subroutine read(this, obs, signal, mask, status)

        class(pacsinstrument), intent(in)  :: this
        class(pacsobservation), intent(in) :: obs
        real(p), intent(out)               :: signal(:,:)
        logical(1), intent(out)            :: mask(:,:)
        integer, intent(out)               :: status
        integer*8 :: nsamples, destination
        integer   :: nobs, iobs

        nobs     = size(obs%info)
        nsamples = sum(obs%info%nsamples)

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
            call this%read_one(obs%info(iobs), destination, signal, mask, status)
            if (status /= 0) return
            destination = destination + obs%info(iobs)%nsamples
        end do
   
    end subroutine read


    !---------------------------------------------------------------------------


    subroutine read_one(this, obs, destination, signal, mask, status)
        class(pacsinstrument), intent(in) :: this
        class(pacsobsinfo), intent(in)    :: obs
        integer*8, intent(in)         :: destination
        real*8, intent(inout)         :: signal(:,:)
        logical*1, intent(inout)      :: mask  (:,:)
        integer, intent(out)          :: status
        integer*8                     :: p, q
        integer                       :: idetector, unit, length
        integer                       :: status_close
        integer, allocatable          :: imageshape(:)
        logical                       :: mask_found
        integer*4, allocatable        :: maskcompressed(:)
        integer*4                     :: maskval
        integer*8                     :: ncompressed
        integer*8                     :: isample, icompressed, ibit
        integer*8                     :: firstcompressed, lastcompressed

        ! old style file format
        length = len_trim(obs%filename)
        if (obs%filename(length-4:length) /= '.fits') then
            call this%read_oldstyle(obs, destination, signal, mask, status)
            return
        end if

        ! read signal HDU
        call ft_open_image(trim(obs%filename) // '[Signal]', unit, 3,         &
                           imageshape, status)
        if (status /= 0) return

        do idetector = 1, size(signal,2)

            p = this%pq(1,idetector)
            q = this%pq(2,idetector)

            call ft_readslice(unit, obs%first, obs%last, q+1, p+1,           &
                 imageshape, signal(destination:destination+obs%nsamples-1,   &
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
            return
        end if

        call ft_open_image(trim(obs%filename) // '[Master]', unit, 3,         &
                           imageshape, status)
        if (status /= 0) return

        allocate(maskcompressed(imageshape(1)))

        do idetector = 1, size(mask,2)

            p = this%pq(1,idetector)
            q = this%pq(2,idetector)

            if (this%mask(p+1,q+1)) then
                mask(destination:destination+obs%nsamples-1,idetector) = .true.
                cycle
            end if

            firstcompressed = (obs%first - 1) / 32 + 1
            lastcompressed  = (obs%last  - 1) / 32 + 1
            ncompressed = lastcompressed - firstcompressed + 1

            call ft_readslice(unit, firstcompressed, lastcompressed,       &
                 q+1, p+1, imageshape, maskcompressed(1:ncompressed),status)
            if (status /= 0) go to 999

            ! loop over the bytes of the compressed mask
            do icompressed = firstcompressed, lastcompressed

                maskval = maskcompressed(icompressed-firstcompressed+1)
                if (maskval == 0) cycle

                isample = (icompressed-1)*32 - obs%first + destination + 1

                ! loop over the bits of a compressed mask byte
                do ibit = max(0, obs%first - (icompressed-1)*32-1),           &
                          min(31, obs%last - (icompressed-1)*32-1)
                    mask(isample+ibit,idetector) =                             &
                         mask(isample+ibit,idetector) .or. btest(maskval,ibit)
                end do

            end do

        end do

    999 call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine read_one


    !---------------------------------------------------------------------------


    subroutine read_oldstyle(this, obs, dest, signal, mask, status)
        class(pacsinstrument), intent(in) :: this
        class(pacsobsinfo), intent(in)    :: obs
        integer*8, intent(in)          :: dest
        real*8, intent(inout)          :: signal(:,:)
        logical*1, intent(inout)       :: mask(:,:)
        integer, intent(out)           :: status
        integer*8                      :: p, q
        integer                        :: idetector, unit, status_close
        integer, allocatable           :: imageshape(:)

        ! handle signal
        call ft_open_image(trim(obs%filename) // '_Signal.fits', unit, 3,     &
                           imageshape, status)
        if (status /= 0) return

        do idetector = 1, size(signal,2)

            p = this%pq(1,idetector)
            q = this%pq(2,idetector)
            call ft_readslice(unit, obs%first, obs%last, q+1, p+1,           &
                 imageshape, signal(dest:dest+obs%nsamples-1,idetector),status)
            if (status /= 0) go to 999

        end do

        call ft_close(unit, status)
        if (status /= 0) return

        ! handle mask
        call ft_open_image(trim(obs%filename) // '_Mask.fits', unit, 3,       &
                           imageshape, status)
        if (status /= 0) return

        do idetector = 1, size(mask,2)

            p = this%pq(1,idetector)
            q = this%pq(2,idetector)

            if (this%mask(p+1,q+1)) then
                mask(dest:dest+obs%nsamples-1,idetector) = .true.
                cycle
            end if

            call ft_readslice(unit, obs%first, obs%last, q+1, p+1,           &
                 imageshape, mask(dest:dest+obs%nsamples-1,idetector), status)
            if (status /= 0) go to 999

        end do

    999 call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine read_oldstyle


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


    !---------------------------------------------------------------------------


    subroutine destroy(this)

        class(pacsinstrument), intent(inout) :: this

        deallocate(this%ij)
        deallocate(this%pq)
        deallocate(this%flatfield)
        deallocate(this%corners_uv)
        deallocate(this%corners_yz)

    end subroutine destroy


end module module_pacsinstrument
