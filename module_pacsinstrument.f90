module module_pacsinstrument

    use ISO_FORTRAN_ENV,        only : ERROR_UNIT, OUTPUT_UNIT
    use module_fitstools,       only : ft_create_header, ft_readextension
    use module_math,            only : pInf, mInf, NaN, nint_down, nint_up
    use module_pacsobservation, only : pacsobservation
    use module_pacspointing,    only : pacspointing
    use module_pointingmatrix,  only : pointingelement
    use module_projection,      only : convex_hull,                            &
                                       intersection_polygon_unity_square
    use string,                 only : strinteger
    use module_wcs,             only : init_astrometry, ad2xy_gnomonic
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
        procedure, public :: init_scalar
        procedure         :: read_calibration_files
        procedure         :: filter_detectors
        procedure, public :: compute_mapheader
        procedure, public :: find_minmax
        procedure, public :: compute_projection_sharp_edges
        procedure, public, nopass  :: uv2yz
        procedure, public, nopass  :: yz2ad
        procedure, public, nopass  :: xy2roi

        procedure :: filter_detectors_array

    end type pacsinstrument


contains

    
    subroutine init_scalar(this, obs, fine_sampling_factor, keep_bad_detectors,&
                           status, bad_detector_mask)
        class(pacsinstrument), intent(inout) :: this
        type(pacsobservation), intent(in)    :: obs
        integer, intent(in)                  :: fine_sampling_factor
        logical, intent(in)                  :: keep_bad_detectors
        integer, intent(out)                 :: status
        logical*1, intent(in), optional      :: bad_detector_mask(:,:)

        if (all(fine_sampling_factor /= [1,2,4,8,16,32])) then
            status = 1
            write (ERROR_UNIT,'(a)') "ERROR: invalid sampling factor '" //     &
                strinteger(fine_sampling_factor) // "'. Valid values are 1,2,4,&
                &8,16,32."
            return
        end if

        this%channel = obs%channel
        this%fine_sampling_factor = fine_sampling_factor
        this%transparent_mode = obs%transparent_mode
        this%keep_bad_detectors = keep_bad_detectors

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

    end subroutine init_scalar


    !---------------------------------------------------------------------------


    subroutine init(this, obs, fine_sampling_factor, keep_bad_detectors,       &
                    status, bad_detector_mask)
        class(pacsinstrument), intent(inout) :: this
        type(pacsobservation), intent(in)    :: obs(:)
        integer, intent(in)                  :: fine_sampling_factor
        logical, intent(in)                  :: keep_bad_detectors
        integer, intent(out)                 :: status
        logical*1, intent(in), optional      :: bad_detector_mask(:,:)

        ! check that all observations have the same filter
        if (any(obs%channel /= obs(1)%channel)) then
            status = 1
            write (ERROR_UNIT,'(a)') 'ERROR: Observations have different channe&
                                     &s.'
            return
        end if

        ! check that all observations have the same transparent mode
        if (any(obs%transparent_mode .neqv. obs(1)%transparent_mode)) then
            status = 1
            write (ERROR_UNIT,'(a)') 'ERROR: Observations have different transp&
                                     &rent modes.'
            return
        end if

        call this%init_scalar(obs(1), fine_sampling_factor, keep_bad_detectors,&
                              status, bad_detector_mask)

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
              write(*,*) 'FILTER_DETECTORS: Check calibration files for the green band...'
              call this%filter_detectors_array(this%mask_blue,                 &
                   this%corners_uv_blue, this%distortion_yz_blue,              &
                   this%flatfield_green)
           case ('b')
              call this%filter_detectors_array(this%mask_blue,                &
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


    subroutine roi2pmatrix(roi, coords, nx, ny, itime, nroi, pmatrix)
        integer, intent(in)                  :: roi(:,:,:)
        real*8, intent(in)                   :: coords(:,:)
        type(pointingelement), intent(inout) :: pmatrix(:,:,:)
        integer, intent(in)                  :: nx, ny, itime
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
                        pmatrix(iroi,itime,idetector)%pixel  = ipixel
                        pmatrix(iroi,itime,idetector)%weight = weight
                    end if
                    iroi = iroi + 1
                end do
            end do
            ! fill the rest of the pointing matrix
            pmatrix(iroi:,itime,idetector)%pixel  = ipixel
            pmatrix(iroi:,itime,idetector)%weight = 0
            nroi = max(nroi, iroi-1)
        end do

    end subroutine roi2pmatrix


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

        !XXX
        if (finer_sampling) then
            stop 'COMPUTE_MAPHEADER: Finer sampling is not implemented yet.'
        end if

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
    subroutine find_minmax(this, pointing, finer_sampling, xmin, xmax,         &
                           ymin, ymax)
        class(pacsinstrument), intent(in) :: this
        class(pacspointing), intent(in)   :: pointing
        logical, intent(in)               :: finer_sampling
        real*8, intent(out)               :: xmin, xmax, ymin, ymax
        real*8                            :: ra, dec, pa, chop
        real*8, allocatable               :: hull_uv(:,:), hull(:,:)
        integer, allocatable              :: ihull(:)
        integer                           :: islice, itime, index

        xmin = pInf
        xmax = mInf
        ymin = pInf
        ymax = mInf

        call convex_hull(this%corners_uv, ihull)

        allocate(hull_uv(ndims,size(ihull)))
        allocate(hull(ndims, size(ihull)))
        hull_uv = this%corners_uv(:,ihull)

        do islice = 1, pointing%nslices

            index = 0

            !$omp parallel do default(shared) reduction(min:xmin,ymin) reduction(max:xmax,ymax) &
            !$omp private(itime,ra,dec,pa,chop,hull) firstprivate(index)
            do itime = pointing%first(islice), pointing%last(islice)
                call pointing%get_position(islice, pointing%time(itime), ra, dec, pa, chop, index)
                hull = uv2yz(hull_uv, this%distortion_yz_blue, chop)
                hull = yz2ad(hull, ra, dec, pa)
                hull = ad2xy_gnomonic(hull)
                xmin = min(xmin, minval(hull(1,:)))
                xmax = max(xmax, maxval(hull(1,:)))
                ymin = min(ymin, minval(hull(2,:)))
                ymax = max(ymax, maxval(hull(2,:)))

             end do
             !$omp end parallel do

        end do

    end subroutine find_minmax


    !---------------------------------------------------------------------------


    subroutine compute_projection_sharp_edges(this, pointing, finer_sampling,  &
                                              header, nx, ny, pmatrix, status)
        class(pacsinstrument), intent(in)  :: this
        type(pacspointing), intent(in)     :: pointing
        logical, intent(in)                :: finer_sampling
        character(len=*), intent(in)       :: header
        integer, intent(in)                :: nx, ny
        type(pointingelement), intent(out) :: pmatrix(:,:,:)
        integer, intent(out)               :: status

        real*8  :: coords(ndims,this%ndetectors*nvertices), coords_yz(ndims,this%ndetectors*nvertices)
        real*8  :: ra, dec, pa, chop, chop_old
        integer :: roi(ndims,2,this%ndetectors), islice, itime, npixels_per_sample, nroi, index

        !XXX
        if (finer_sampling) then
            stop 'COMPUTE_MAPHEADER: Finer sampling is not implemented yet.'
        end if

        call init_astrometry(header, status=status)
        if (status /= 0) return

        npixels_per_sample = -1

        do islice = 1, pointing%nslices
            chop_old = pInf
            index = 0
            !$omp parallel do default(shared) firstprivate(index, chop_old)   &
            !$omp private(itime, ra, dec, pa, chop, coords, coords_yz, roi) &
            !$omp reduction(max : npixels_per_sample)
            do itime = pointing%first(islice), pointing%last(islice)
                call pointing%get_position(islice, pointing%time(itime), ra, dec, pa, chop, index)
                if (abs(chop-chop_old) > 1.d-2) then
                    coords_yz = this%uv2yz(this%corners_uv, this%distortion_yz, chop)
                    chop_old = chop
                 end if
                 coords = this%yz2ad(coords_yz, ra, dec, pa)
                 coords = ad2xy_gnomonic(coords)
                 roi    = xy2roi(coords) ! [1=x|2=y,1=min|2=max,idetector]
                 call roi2pmatrix(roi, coords, nx, ny, itime, nroi, pmatrix)
                 npixels_per_sample = max(npixels_per_sample, nroi)
             end do
             !$omp end parallel do
        end do

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
