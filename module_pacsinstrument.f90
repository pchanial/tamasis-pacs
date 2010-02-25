module module_pacsinstrument

    use module_fitstools
    use module_pacspointing
    use module_pointingmatrix
    use module_projection
    use module_wcslib, only: WCSLEN
    use string, only : strlowcase
    implicit none

    integer, parameter :: shape_blue(2) = [32, 64]
    integer, parameter :: shape_red (2) = [16, 32]
    integer, parameter :: ndims = 2
    integer, parameter :: nvertices = 4
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

        procedure         :: read_calibration_files
        procedure         :: filter_detectors
        procedure         :: read_signal_file
        procedure         :: read_mask_file
        procedure         :: compute_mapheader
        procedure         :: find_minmax
        procedure         :: compute_projection_sharp_edges
        procedure, nopass :: get_array_color
        procedure, nopass :: uv2yz
        procedure, nopass :: yz2ad
        procedure, nopass :: xy2roi
        procedure         :: roi2pmatrix
        procedure, nopass :: multiplexing_direct
        procedure, nopass :: multiplexing_transpose

    end type pacsinstrument


contains


    subroutine read_calibration_files(this)

        class(pacsinstrument), intent(inout)  :: this

        character(len=*), parameter :: filename_saa = 'PCalPhotometer_SubArrayArray_FM_v5.fits'
        character(len=*), parameter :: filename_ai  = 'PCalPhotometer_ArrayInstrument_FM_v4.fits'
        character(len=*), parameter :: filename_bpm = 'PCalPhotometer_BadPixelMask_FM_v3.fits'
        character(len=*), parameter :: filename_ff  = 'PCalPhotometer_FlatField_FM_v1.fits'

        integer,          parameter :: hdu_blue(4) = [8, 12, 16, 20]
        integer,          parameter :: hdu_red (4) = [6, 10, 14, 18]

        integer                     :: ivertex, status
        logical*1, allocatable      :: tmplogical(:,:)
        real*8, allocatable         :: tmp2(:,:)
        real*8, allocatable         :: tmp3(:,:,:)
        character(len=100) :: value

        ! might be moved elsewhere
        call get_environment_variable('TAMASIS_DIR', value, length=tamasis_dir_len)
        if (tamasis_dir_len == 0) tamasis_dir_len = len('/home/pchanial/work/tamasis')

        status = 0

        ! read bad pixel mask

        call ft_readextension(get_calfile(filename_bpm) // '[blue]', tmplogical, status)
        call ft_printerror(status, filename_bpm // '[blue]')
        this%mask_blue = transpose(tmplogical)

        call ft_readextension(get_calfile(filename_bpm) // '[red]', tmplogical, status)
        call ft_printerror(status, filename_bpm // '[red]')
        this%mask_red = transpose(tmplogical)


        ! read flat fields

        call ft_readextension(get_calfile(filename_ff) // '+5',  tmp2, status)
        call ft_printerror(status, filename_ff // '+5')
        this%flatfield_blue = transpose(tmp2)
        call ft_readextension(get_calfile(filename_ff) // '+8', tmp2, status)
        call ft_printerror(status, filename_ff // '+8')
        this%flatfield_green = transpose(tmp2)
        call ft_readextension(get_calfile(filename_ff) // '+2',   tmp2, status)
        call ft_printerror(status, filename_ff // '+2')
        this%flatfield_red = transpose(tmp2)


        ! read detector corners in the (u,v) plane

        do ivertex=1, nvertices

            call ft_readextension(get_calfile(filename_saa), tmp2, status, hdu_blue(ivertex)  )
            call ft_printerror(status, filename_saa)
            this%corners_uv_blue(1,ivertex,:,:) = transpose(tmp2)
            call ft_readextension(get_calfile(filename_saa), tmp2, status, hdu_blue(ivertex)+1)
            call ft_printerror(status, filename_saa)
            this%corners_uv_blue(2,ivertex,:,:) = transpose(tmp2)
            call ft_readextension(get_calfile(filename_saa), tmp2, status, hdu_red (ivertex)  )
            call ft_printerror(status, filename_saa)
            this%corners_uv_red (1,ivertex,:,:) = transpose(tmp2)
            call ft_readextension(get_calfile(filename_saa), tmp2, status, hdu_red (ivertex)+1)
            call ft_printerror(status, filename_saa)
            this%corners_uv_red (2,ivertex,:,:) = transpose(tmp2)

        end do


        ! read the distortion coefficients in the (y,z) plane

        call ft_readextension(get_calfile(filename_ai) // '[ycoeffblue]', tmp3, status)
        call ft_printerror(status, filename_ai // '[ycoeffblue]')
        this%distortion_yz_blue(1,:,:,:) = tmp3
        call ft_readextension(get_calfile(filename_ai) // '[zcoeffblue]', tmp3, status)
        call ft_printerror(status, filename_ai // '[zcoeffblue]')
        this%distortion_yz_blue(2,:,:,:) = tmp3
        call ft_readextension(get_calfile(filename_ai) // '[ycoeffred]', tmp3, status)
        call ft_printerror(status, filename_ai // '[ycoeffred]')
        this%distortion_yz_red (1,:,:,:) = tmp3
        call ft_readextension(get_calfile(filename_ai) // '[zcoeffred]', tmp3, status)
        call ft_printerror(status, filename_ai // '[zcoeffred]')
        this%distortion_yz_red (2,:,:,:) = tmp3

    end subroutine read_calibration_files


    !-------------------------------------------------------------------------------


    function get_calfile(filename)
        character(len=*), intent(in)                   :: filename
        character(len=tamasis_dir_len+6+len(filename)) :: get_calfile
        integer                                        :: length

        ! Environment variables are not visible in photran...
        call get_environment_variable('TAMASIS_DIR', get_calfile, length)
        if (length == 0) get_calfile(1:tamasis_dir_len) = '/home/pchanial/work/tamasis'
        get_calfile(tamasis_dir_len+1:) = '/data/' // filename

    end function get_calfile


    !-------------------------------------------------------------------------------


    subroutine read_signal_file(this, filename, first, last, signal)

        class(pacsinstrument), intent(in)  :: this
        character(len=*), intent(in)       :: filename
        integer*8, intent(in)              :: first, last
        real*8, intent(out)                :: signal(last-first+1, this%ndetectors)
        integer*8                          :: p, q
        integer                            :: idetector, status, unit
        integer, allocatable               :: imageshape(:)

        status = 0
        !should be tested with cfitsio compiled with ./configure --enable-reentrant
        !!$omp parallel default(shared) private(idetector, p, q, unit, imageshape) firstprivate(status)
        call ft_openimage(filename, unit, 3, imageshape, status)
        if (status /= 0) then
            call ft_printerror(status, filename)
        end if

        !$omp do
        do idetector = 1, this%ndetectors

            p = this%pq(1,idetector)
            q = this%pq(2,idetector)
            call ft_readslice(unit, first, last, q+1, p+1, imageshape, signal(:,idetector), status)
            if (status /= 0) stop "ERROR in read_signal_file: ft_readslice."

        end do
        !$omp end do

        call ft_close(unit, status)
        !!$omp end parallel

    end subroutine read_signal_file


    !-------------------------------------------------------------------------------


    subroutine read_mask_file(this, filename, first, last, mask)

        class(pacsinstrument), intent(in)  :: this
        character(len=*), intent(in)       :: filename
        integer*8, intent(in)              :: first, last
        logical*1, intent(out)             :: mask(last-first+1, this%ndetectors)
        integer*8                          :: p, q
        integer                            :: idetector, status, unit
        integer, allocatable               :: imageshape(:)

        status = 0
        !should be tested with cfitsio compiled with ./configure --enable-reentrant
        !!$omp parallel default(shared) private(idetector, p, q, unit, imageshape) firstprivate(status)
        call ft_openimage(filename, unit, 3, imageshape, status)
        if (status /= 0) then
            call ft_printerror(status, filename)
        end if

        !$omp do
        do idetector = 1, this%ndetectors

            p = this%pq(1,idetector)
            q = this%pq(2,idetector)
            call ft_readslice(unit, first, last, q+1, p+1, imageshape, mask(:,idetector), status)
            if (status /= 0) stop "ERROR in read_mask_file: ft_readslice."

        end do
        !$omp end do

        call ft_close(unit, status)
        !!$omp end parallel

    end subroutine read_mask_file


    !-------------------------------------------------------------------------------


    subroutine filter_detectors(this, side, transparent_mode, keep_bad_detectors)

        class(pacsinstrument), intent(inout)   :: this
        character(len=*), intent(in), optional :: side
        logical, intent(in), optional          :: transparent_mode
        logical, intent(in), optional          :: keep_bad_detectors

        if (present(side)) then
            if (strlowcase(side) == 'red') then
                call filter_detectors_array(this, this%mask_red, this%corners_uv_red, this%distortion_yz_red, &
                                            this%flatfield_red, transparent_mode=transparent_mode,            &
                                            keep_bad_detectors=keep_bad_detectors)
                return
            else if (strlowcase(side) == 'green') then
                write(*,*) 'Check calibration files for the green band...'
                call filter_detectors_array(this, this%mask_blue, this%corners_uv_blue, this%distortion_yz_blue, &
                                            this%flatfield_green, transparent_mode=transparent_mode,             &
                                            keep_bad_detectors=keep_bad_detectors)
                return
            else if (strlowcase(side) /= 'blue') then
                write (*,*) "FILTER_DETECTORS: invalid array side ('blue', 'green' or 'red'): " // strlowcase(side)
                stop
            endif
        endif

        call filter_detectors_array(this, this%mask_blue, this%corners_uv_blue, this%distortion_yz_blue, &
                                    this%flatfield_blue, transparent_mode=transparent_mode,              &
                                    keep_bad_detectors=keep_bad_detectors)

    end subroutine filter_detectors


    !-------------------------------------------------------------------------------


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


    !-------------------------------------------------------------------------------


    function get_array_color(filename) result(color)
        character(len=*), intent(in)                 :: filename
        character(len=get_array_color_len(filename)) :: color

        if (len(color) == 4) then
            color = 'blue'
        else if (len(color) == 5) then
            color = 'green'
        else if (len(color) == 3) then
            color = 'red'
        else
            stop "File name must contain the array identifier (blue, green, red)."
        end if

    end function get_array_color


    !-------------------------------------------------------------------------------


    pure function get_array_color_len(filename) result(length)
        character(len=*), intent(in) :: filename
        integer                      :: length

        if (strlowcase(filename(len(filename)-2:len(filename))) == 'lue') then
            length = 4
        else if (strlowcase(filename(len(filename)-2:len(filename))) == 'een') then
            length = 5
        else if (strlowcase(filename(len(filename)-2:len(filename))) == 'red') then
            length = 3
        else
            length = 0
        end if
    end function get_array_color_len


    !-------------------------------------------------------------------------------


    recursive function uv2yz(uv, distortion_yz, chop) result(yz)

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
 

    !-------------------------------------------------------------------------------


    recursive function yz2ad(yz, ra0, dec0, pa0) result (ad)

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


    !-------------------------------------------------------------------------------


   recursive function xy2roi(xy) result(roi)
        real*8, intent(in) :: xy(:,:)
        real*8             :: roi(ndims,2,size(xy,2)/nvertices)
        integer            :: idetector

        do idetector = 1, size(xy,2) / nvertices
            roi(:,1,idetector) = nint(minval(xy(:,nvertices * (idetector-1)+1:nvertices*idetector),2))
            roi(:,2,idetector) = nint(maxval(xy(:,nvertices * (idetector-1)+1:nvertices*idetector),2))
        end do

    end function xy2roi


    !------------------------------------------------------------------------------


    recursive subroutine roi2pmatrix(this, roi, coords, nx, ny, isample, nroi, pmatrix)
        class(pacsinstrument), intent(in)    :: this
        integer, intent(in)                  :: roi(ndims,2,this%ndetectors)
        real*8, intent(in)                   :: coords(ndims,this%ndetectors*nvertices)
        type(pointingelement), intent(inout) :: pmatrix(:,:,:)
        integer, intent(in)                  :: nx, ny, isample
        integer, intent(out)                 :: nroi
        real*8                               :: polygon(2,nvertices)

        integer                              :: npixels_per_sample, idetector, ix, iy, iroi, ipixel, ibad
        real*8                               :: weight

        ipixel = 0
        ibad = 0
        nroi = 0
        npixels_per_sample = size(pmatrix, 1)
        do idetector = 1, this%ndetectors
            iroi = 1
            do_roi: do iy = max(roi(2,1,idetector),1), min(roi(2,2,idetector),ny)
                do ix = max(roi(1,1,idetector),1), min(roi(1,2,idetector),nx)
                    ipixel = ix - 1 + (iy - 1) * nx
                    polygon(1,:) = coords(1,(idetector-1)*nvertices+1:idetector*nvertices) - (ix-0.5d0)
                    polygon(2,:) = coords(2,(idetector-1)*nvertices+1:idetector*nvertices) - (iy-0.5d0)
                    weight = abs(intersection_polygon_unity_square(polygon, nvertices))
                    !XXX SHOULD DO: if weight less than threshold skip it
                    if (weight <= 0) cycle
                    if (iroi <= npixels_per_sample) then
                        pmatrix(iroi,isample,idetector)%pixel  = ipixel
                        pmatrix(iroi,isample,idetector)%weight = weight
                    end if
                    iroi = iroi + 1
                end do
            end do do_roi
            ! fill the rest of the pointing matrix
            pmatrix(iroi:,isample,idetector)%pixel  = ipixel
            pmatrix(iroi:,isample,idetector)%weight = 0
            nroi = max(nroi, iroi-1)
        end do

    end subroutine roi2pmatrix


    !------------------------------------------------------------------------------


    subroutine compute_mapheader(this, pointing, times, resolution, header)

        use module_fitstools, only : ft_create_header, ft_header2wcs
        use module_wcslib, only : wcsfree
        use module_wcs, only : init_wcslib

        class(pacsinstrument), intent(in) :: this
        class(pacspointing), intent(in)   :: pointing
        real*8, intent(in)                :: times(:)  ! times at which the PACS array should be considered as inside the map
        real*8, intent(in)                :: resolution
        character(len=2880), intent(out)  :: header

        integer                           :: nx, ny, status
        integer                           :: ixmin, ixmax, iymin, iymax
        real*8                            :: ra0, dec0, xmin, xmax, ymin, ymax

        call pointing%compute_center(times, ra0, dec0)

        call ft_create_header(0, 0, -resolution/3600.d0, resolution/3600.d0, 0.d0, ra0, dec0, 1.d0, 1.d0, header)
        call init_wcslib(header)

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


    !------------------------------------------------------------------------------


    ! find minimum and maximum pixel coordinates in maps
    subroutine find_minmax(this, pointing, times, xmin, xmax, ymin, ymax)

        use module_projection, only: convex_hull
        use module_wcs, only : ad2xy_wcslib, init_wcslib
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
           hull = ad2xy_wcslib(hull)
           xmin = min(xmin, minval(hull(1,:)))
           xmax = max(xmax, maxval(hull(1,:)))
           ymin = min(ymin, minval(hull(2,:)))
           ymax = max(ymax, maxval(hull(2,:)))
        end do
        !!$omp end parallel do

    end subroutine find_minmax


    !------------------------------------------------------------------------------


    subroutine compute_projection_sharp_edges(this, pointing, time, header, nx, ny, pmatrix)

        use module_wcs, only : init_wcslib, ad2xy_wcslib, free_wcslib
        class(pacsinstrument), intent(in)    :: this
        class(pacspointing), intent(in)      :: pointing
        character(len=*), intent(in)         :: header
        real*8, intent(in)                   :: time(:)
        integer, intent(in)                  :: nx, ny
        type(pointingelement), intent(out)   :: pmatrix(:,:,:)

        real*8  :: coords(ndims,this%ndetectors*nvertices), ra, dec, pa, chop
        integer :: roi(ndims,2,this%ndetectors), isample, nsamples, npixels_per_sample, nroi, index

        nsamples = size(time)
        npixels_per_sample = -1

        call init_wcslib(header)

        index = 2
        !$omp parallel do default(shared) firstprivate(index) private(isample, ra, dec, pa, chop, coords, roi) &
        !$omp reduction(max : npixels_per_sample)
        do isample = 1, nsamples
            call pointing%get_position(time(isample), ra, dec, pa, chop, index)
            coords = this%uv2yz(this%corners_uv, this%distortion_yz_blue, chop)
            coords = this%yz2ad(coords, ra, dec, pa)
            coords = ad2xy_wcslib(coords)
            roi    = this%xy2roi(coords) ! [1=x|2=y,1=min|2=max,idetector]
            call this%roi2pmatrix(roi, coords, nx, ny, isample, nroi, pmatrix)
            npixels_per_sample = max(npixels_per_sample, nroi)
        end do
        !$omp end parallel do
        if (npixels_per_sample /= size(pmatrix,1)) then
            write(*,'(a,i0,a)') 'Warning: For the Pointing Matrix, npixels_per_sample may be updated to ', npixels_per_sample, '.'
        end if

        call free_wcslib()

    end subroutine compute_projection_sharp_edges


    !------------------------------------------------------------------------------


    subroutine multiplexing_direct(signal, sampled_signal, sampling, ij)
        real*8, intent(in)  :: signal(:,:)
        real*8, intent(out) :: sampled_signal(size(signal,1)/sampling,size(signal,2))
        integer, intent(in) :: sampling
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
                frac = fraction(j * sampling / 16.d0)
                sampled_signal(:,idetector) = (1-frac) * signal(isample  ::sampling,idetector) + &
                                                 frac  * signal(isample+1::sampling,idetector)
            end do
            !$omp end parallel do
        end if

    end subroutine multiplexing_direct


    !------------------------------------------------------------------------------


    subroutine multiplexing_transpose(sampled_signal, signal, sampling, ij)
        real*8, intent(in)  :: sampled_signal(:,:)
        real*8, intent(out) :: signal(size(sampled_signal,1)*sampling,size(sampled_signal,2))
        integer, intent(in) :: sampling
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
                frac = fraction(j * sampling / 16.d0)
                signal(isample  ::sampling,idetector) = (1-frac) * sampled_signal(isample::sampling,idetector)
                signal(isample+1::sampling,idetector) =    frac  * sampled_signal(isample::sampling,idetector)
            end do
            !$omp end do
            !$omp end parallel
        end if

    end subroutine multiplexing_transpose


end module module_pacsinstrument
