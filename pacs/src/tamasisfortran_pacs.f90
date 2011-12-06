! Copyright 2010-2011 Pierre Chanial
! All rights reserved
!
! Tamasis interface for f2py
!
! Author: P. Chanial

subroutine pacs_info_instrument_init(filename, nfilenames, band, transparent_mode, status)

    use iso_fortran_env,        only : ERROR_UNIT
    use module_fitstools,       only : FLEN_VALUE, ft_close, ft_open, ft_read_keyword, ft_read_keyword_hcss
    use module_pacsobservation, only : PacsObservationSlice
    implicit none

    character(len=*), intent(in)  :: filename
    integer, intent(in)           :: nfilenames
    character(len=7), intent(out) :: band
    logical, intent(out)          :: transparent_mode
    integer, intent(out)          :: status

    class(PacsObservationSlice), allocatable :: slice
    integer                                  :: iobs, unit, first, last
    character(len=FLEN_VALUE)                :: camname, blue, compmode
    character(len=7)                         :: band_
    logical                                  :: found

    ! split input filename
    if (mod(len(filename), nfilenames) /= 0) then
        stop 'PACS_INFO_INIT: Invalid filename length.'
    end if

    band = 'unknown'
    transparent_mode = .true.
    allocate (slice)

    do iobs = 1, nfilenames

        call slice%set_filename(filename((iobs-1)*len(filename)/nfilenames+1:iobs*len(filename)/nfilenames), first, last, status)
        if (status /= 0) return

        call ft_open(trim(slice%filename), unit, status)
        if (status /= 0) return

        call ft_read_keyword_hcss(unit, 'camName', camname, status=status)
        if (status /= 0) return

        call ft_read_keyword_hcss(unit, 'blue', blue, found, status)
        if (status /= 0) return
        if (.not. found) then
            call ft_read_keyword(unit, 'FILTER', blue, found, status=status)
            if (status /= 0) return
        end if

        call ft_read_keyword_hcss(unit, 'compMode', compmode, status=status)
        if (status /= 0) return

        call ft_close(unit, status)
        if (status /= 0) return

        if (camname == 'Blue Photometer') then
            if (blue == 'blue1' .or. blue == 'blue70um') then
                band_ = 'blue'
            else if (blue == 'blue2' .or. blue == 'green100um') then
                band_ = 'green'
            else
                write (ERROR_UNIT, '(a)') "Error: Invalid blue camera identifier '" // trim(blue) // "' in file '" //              &
                      trim(slice%filename) // "'."
                status = 1
                return
            end if
        else if (camname == 'Red Photometer') then
            band_ = 'red'
        else
            write (ERROR_UNIT, '(a)') "Error: Invalid camera name '" // trim(camname) // "' in file '" // trim(slice%filename)     &
                  // "'."
            status = 1
            return
        end if

        if (iobs == 1) then
            band = band_
        else if (band_ /= band) then
            write (ERROR_UNIT, '(a)') "Error: Observations are not performed with the same band."
            status = 1
            return
        end if

        if (compmode /= 'Photometry Lossless Compression Mode') then
            transparent_mode = .false.
        end if

    end do

end subroutine pacs_info_instrument_init


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_info_instrument(band, nrows, ncolumns, input_active_fraction, detector_bad, detector_center,        &
                                detector_corner, detector_area, distortion_yz, flatfield_optical, flatfield_detector, responsivity,&
                                output_active_fraction, status)

    use module_pacsinstrument, only : read_calibration_files
    use module_tamasis,        only : p
    implicit none

    character(len=*), intent(in) :: band
    integer, intent(in)          :: nrows, ncolumns
    real(p), intent(in)          :: input_active_fraction
    logical*1, intent(out)       :: detector_bad(nrows,ncolumns)
    real(p), intent(out)         :: detector_center(2,nrows,ncolumns)
    real(p), intent(out)         :: detector_corner(2,4,nrows,ncolumns)
    real(p), intent(out)         :: detector_area(nrows,ncolumns)
    real(p), intent(out)         :: distortion_yz(2,3,3,3)
    real(p), intent(out)         :: flatfield_optical(nrows,ncolumns)
    real(p), intent(out)         :: flatfield_detector(nrows,ncolumns)
    real(p), intent(out)         :: responsivity
    real(p), intent(out)         :: output_active_fraction
    integer, intent(out)         :: status

    logical*1, allocatable :: detector_bad_all(:,:)
    real(p), allocatable   :: detector_center_all(:,:,:)
    real(p), allocatable   :: detector_corner_all(:,:,:,:)
    real(p), allocatable   :: detector_area_all(:,:)
    real(p), allocatable   :: flatfield_optical_all(:,:)
    real(p), allocatable   :: flatfield_detector_all(:,:)

    ! read calibration files
    output_active_fraction = input_active_fraction
    call read_calibration_files(band, detector_bad_all, detector_center_all, detector_corner_all, detector_area_all,               &
         flatfield_optical_all, flatfield_detector_all, distortion_yz, responsivity, output_active_fraction, status)
    if (status /= 0) return

    ! copy the allocatable arrays
    detector_bad       = detector_bad_all
    detector_center    = detector_center_all
    detector_corner    = detector_corner_all
    detector_area      = detector_area_all
    flatfield_optical  = flatfield_optical_all
    flatfield_detector = flatfield_detector_all

end subroutine pacs_info_instrument


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_info_observation_init(filename, nfilenames, nsamples, status)

    use module_fitstools,       only : ft_close, ft_open, ft_read_keyword
    use module_pacsobservation, only : PacsObservationSlice
    use module_tamasis,         only : p
    implicit none

    character(len=*), intent(in) :: filename
    integer, intent(in)          :: nfilenames
    integer, intent(out)         :: nsamples(nfilenames)
    integer, intent(out)         :: status

    class(PacsObservationSlice), allocatable :: slice
    integer                                  :: iobs, unit, first, last
    allocate (slice)

    do iobs = 1, nfilenames

        call slice%set_filename(filename((iobs-1)*len(filename)/nfilenames+1:iobs*len(filename)/nfilenames), first, last, status)
        if (status /= 0) return

        call ft_open(trim(slice%filename) // '[Signal]', unit, status)
        if (status /= 0) return

        call ft_read_keyword(unit, 'NAXIS1', nsamples(iobs), status=status)
        if (status /= 0) return

        call ft_close(unit, status)
        if (status /= 0) return

    end do

end subroutine pacs_info_observation_init


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_info_observation(filename, nfilenames, policy, nsamples_tot, masktime_calblock, obsid, cusmode, compression, unit, &
                                 ra, dec, cam_angle, scan_angle, scan_length, scan_step, scan_nlegs, frame_time, frame_ra,         &
                                 frame_dec, frame_pa, frame_chop, frame_info, frame_masked, frame_removed, nmasks, mask_name,      &
                                 mask_activated, status)

    use module_observation,     only : POINTING_INSCAN, POINTING_TURNAROUND, POINTING_OTHER, MaskPolicy
    use module_pacsobservation, only : PacsObservation
    use module_tamasis,         only : POLICY_KEEP, POLICY_MASK, POLICY_REMOVE, p
    implicit none

    !f2py character(len=70*nfilenames), intent(out), depend(nfilenames) :: cusmode, unit
    !f2py character(len=70*32*nfilenames), intent(out), depend(nfilenames) :: mask_name

    character(len=*), intent(in)                     :: filename
    integer, intent(in)                              :: nfilenames
    integer, intent(in)                              :: policy(4)
    integer, intent(in)                              :: nsamples_tot
    real(p), intent(in)                              :: masktime_calblock
    character(len=70*nfilenames), intent(out)        :: cusmode, unit !XXX hack around f2py bug with array of characters
    real(p), intent(out), dimension(nfilenames)      :: ra, dec, cam_angle, scan_angle, scan_length, scan_step
    integer, intent(out), dimension(nfilenames)      :: obsid, compression, scan_nlegs, nmasks
    real(p), intent(out), dimension(nsamples_tot)    :: frame_time, frame_ra, frame_dec, frame_pa, frame_chop, frame_info
    logical*1, intent(out), dimension(nsamples_tot)  :: frame_masked, frame_removed
    character(len=70*32*nfilenames), intent(out)     :: mask_name
    logical*1, intent(out), dimension(32,nfilenames) :: mask_activated
    integer, intent(out)                             :: status

    integer, parameter                  :: NMASKS_MAX = 32
    integer, parameter                  :: FLEN_VALUE = 70
    class(PacsObservation), allocatable :: obs
    type(MaskPolicy)                    :: mask_policy
    integer                             :: iobs, isample, imask, nsamples, dest
    character(len=len(filename)/nfilenames), allocatable :: filename_(:)

    ! split input filename
    if (mod(len(filename), nfilenames) /= 0) then
        stop 'PACS_INFO: Invalid filename length.'
    end if
    allocate(filename_(nfilenames))
    do iobs = 1, nfilenames
        filename_(iobs) = filename((iobs-1)*len(filename_)+1:iobs*len(filename_))
    end do

    ! initialise policy
    mask_policy = MaskPolicy(inscan=policy(1), turnaround=policy(2), other=policy(3), invalid=policy(4))

    ! read observations
    allocate (obs)
    call obs%init(filename_, mask_policy, status, verbose=.true., masktime_calblock=masktime_calblock)
    if (status /= 0) return

    cusmode = ''
    unit = ''
    mask_name = ''
    mask_activated = .false.
    do iobs = 1, nfilenames
        cusmode((iobs-1)*FLEN_VALUE+1:iobs*FLEN_VALUE) = obs%slice(iobs)%observing_mode
        unit   ((iobs-1)*FLEN_VALUE+1:iobs*FLEN_VALUE) = obs%slice(iobs)%unit
        nmasks(iobs) = obs%slice(iobs)%nmasks
        do imask = 1, nmasks(iobs)
            dest = ((iobs-1)*NMASKS_MAX+(imask-1))*FLEN_VALUE
            mask_name(dest+1:dest+FLEN_VALUE) = obs%slice(iobs)%mask_name(imask)
            mask_activated(imask,iobs) = obs%slice(iobs)%mask_activated(imask)
        end do
    end do
    obsid       = obs%slice%obsid
    compression = obs%slice%compression_factor
    ra          = obs%slice%ra
    dec         = obs%slice%dec
    cam_angle   = obs%slice%cam_angle
    scan_angle  = obs%slice%scan_angle
    scan_length = obs%slice%scan_length
    scan_step   = obs%slice%scan_step
    scan_nlegs  = obs%slice%scan_nlegs
    
    dest = 1
    do iobs = 1, obs%nslices
        nsamples = obs%slice(iobs)%nsamples
        frame_time   (dest:dest+nsamples-1) = obs%slice(iobs)%p%time
        frame_ra     (dest:dest+nsamples-1) = obs%slice(iobs)%p%ra
        frame_dec    (dest:dest+nsamples-1) = obs%slice(iobs)%p%dec
        frame_pa     (dest:dest+nsamples-1) = obs%slice(iobs)%p%pa
        frame_chop   (dest:dest+nsamples-1) = obs%slice(iobs)%p%chop
        frame_masked (dest:dest+nsamples-1) = obs%slice(iobs)%p%masked
        frame_removed(dest:dest+nsamples-1) = obs%slice(iobs)%p%removed
        do isample = 1, nsamples
            if (obs%slice(iobs)%p(isample)%inscan) then
                frame_info(dest) = POINTING_INSCAN
            else if (obs%slice(iobs)%p(isample)%turnaround) then
                frame_info(dest) = POINTING_TURNAROUND
            else
                frame_info(dest) = POINTING_OTHER
            end if
            dest = dest + 1
        end do
    end do

end subroutine pacs_info_observation


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_read_filter_uncorrelated(filename, ncorrelations, ndetectors, mask, nrows, ncolumns, data, status)

    use iso_fortran_env,       only : ERROR_UNIT
    use module_filtering,      only : FilterUncorrelated
    use module_pacsinstrument, only : read_filter_uncorrelated
    use module_tamasis,        only : p
    implicit none

    character(len=*), intent(in) :: filename
    integer, intent(in)          :: ncorrelations
    integer, intent(in)          :: ndetectors
    logical*1, intent(in)        :: mask(nrows,ncolumns)
    integer, intent(in)          :: nrows
    integer, intent(in)          :: ncolumns
    real(p), intent(out)         :: data(ncorrelations+1,ndetectors)
    integer, intent(out)         :: status

    type(FilterUncorrelated)     :: filter

    call read_filter_uncorrelated(filename, mask, filter, status)
    if (status /= 0) return

    if (filter%ndetectors /= ndetectors) then
        status = 1
        write (ERROR_UNIT,'(a,2(i0,a))') "The specified number of detectors '", ndetectors, "' does not match that in the calibrati&
              &on file '", filter%ndetectors, "'."
    end if
    if (filter%ncorrelations /= ncorrelations) then
        status = 1
        write (ERROR_UNIT,'(a,2(i0,a))') "The specified number of correlations '", ncorrelations, "' does not match that in the cal&
              &ibration file '", filter%ncorrelations, "'."
    end if
    if (status /= 0) return

    data = filter%data

end subroutine pacs_read_filter_uncorrelated


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_map_header(band, compression_factor, delay, fine_sampling_factor, oversampling, time, ra, dec, pa, chop,           &
                           masked, removed, npointings, detector_masked, detector_removed, nrows, ncolumns, detector_corner,       &
                           distortion_yz, resolution, header, status)

    use module_pacsinstrument, only : SAMPLING_PERIOD, PacsInstrument
    use module_observation,    only : Observation
    use module_tamasis,        only : p
    implicit none

    character(len=*), intent(in)                   :: band
    integer, intent(in)                            :: compression_factor
    real(p), intent(in)                            :: delay
    integer, intent(in)                            :: fine_sampling_factor
    logical, intent(in)                            :: oversampling
    real(p), intent(in), dimension(npointings)     :: time, ra, dec, pa, chop
    logical*1, intent(in), dimension(npointings)   :: masked, removed
    integer, intent(in)                            :: npointings
    logical*1, intent(in)                          :: detector_masked(nrows,ncolumns)
    logical*1, intent(in)                          :: detector_removed(nrows,ncolumns)
    integer, intent(in)                            :: nrows
    integer, intent(in)                            :: ncolumns
    real(p), intent(in)                            :: detector_corner(2,4,nrows,ncolumns)
    real(p), intent(in)                            :: distortion_yz(2,3,3,3)
    real(p), intent(in)                            :: resolution
    character(len=2880), intent(out)               :: header
    integer, intent(out)                           :: status

    class(Observation), allocatable    :: obs
    class(PacsInstrument), allocatable :: pacs

    ! initialise observations
    allocate(obs)
    ! offset of the fine time of the compressed frame wrt the fine time of the first frame
    ! the fine time of a compressed frame is the average of the fine times of the contributing frames
    ! the offset is in fraction of sampling interval
    call obs%init(time, ra, dec, pa, chop, masked, removed, compression_factor,                                                    &
         (compression_factor - 1) / (2._p * compression_factor) +  delay / (1000 * SAMPLING_PERIOD * compression_factor), status)
    if (status /= 0) return

    ! initialise pacs instrument
    allocate (pacs)
    call pacs%init_with_variables(band, detector_removed, fine_sampling_factor, status, detector_bad_all=detector_masked,          &
                                  detector_corner_all=detector_corner, distortion_yz=distortion_yz)
    if (status /= 0) return

    ! compute the map header
    call pacs%compute_map_header(obs, oversampling, resolution, header, status)

end subroutine pacs_map_header


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_read_tod(signal, mask, first, last, band, filename, nvalids_tot, nsamples_all, masked, removed, detector_masked,   &
                         detector_removed, flatfield_detector, do_flatfielding, do_subtraction_mean, ndetectors, selected_mask,    &
                         nrows, ncolumns, status)

    use iso_fortran_env,        only : ERROR_UNIT
    use module_observation,     only : Observation
    use module_pacsinstrument,  only : SAMPLING_PERIOD, PacsInstrument
    use module_preprocessor,    only : divide_vectordim2, subtract_meandim1
    use module_tamasis,         only : p
    use module_string,          only : strsplit
    implicit none

    real(p), intent(inout)       :: signal(nvalids_tot, ndetectors)
    logical*1, intent(inout)     :: mask(nvalids_tot, ndetectors)
    integer*8, intent(in)        :: first, last
    character(len=*), intent(in) :: band
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: nvalids_tot, nsamples_all
    logical*1, intent(in)        :: masked(nsamples_all), removed(nsamples_all)
    logical*1, intent(in)        :: detector_masked(nrows,ncolumns), detector_removed(nrows,ncolumns)
    real(p), intent(in)          :: flatfield_detector(nrows,ncolumns)
    logical, intent(in)          :: do_flatfielding, do_subtraction_mean
    integer, intent(in)          :: ndetectors
    integer, intent(in)          :: nrows
    integer, intent(in)          :: ncolumns
    character(len=*), intent(in) :: selected_mask
    integer, intent(out)         :: status

    class(PacsInstrument), allocatable             :: pacs
    character(len=len(selected_mask)), allocatable :: selected_mask_(:)

    ! initialise pacs instrument
    allocate (pacs)
    call pacs%init_with_variables(band, detector_removed, 1, status, detector_bad_all=detector_masked,                             &
                                  flatfield_detector_all=flatfield_detector)
    if (status /= 0) return

    ! check number of detectors
    if (pacs%ndetectors /= ndetectors) then
        status = 1
        write (ERROR_UNIT,'(a,2(i0,a))') "The specified number of detectors '", ndetectors, "' is incompatible with that from the i&
              &nput detector mask '", pacs%ndetectors, "'."
        return
    end if

    ! split masks
    call strsplit(selected_mask, ',', selected_mask_)

    ! read timeline
    call pacs%read(filename, masked, removed, selected_mask_, signal(first:last,:), mask(first:last,:), status)
    if (status /= 0) return

    ! detector flat fielding
    if (do_flatfielding) then
        call divide_vectordim2(signal(first:last,:), pacs%flatfield_detector)
    end if

    ! subtract the mean of each detector timeline
    if (do_subtraction_mean) then
        call subtract_meandim1(signal(first:last,:), mask(first:last,:))
    end if
    
end subroutine pacs_read_tod


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_pointing_matrix(band, nvalids, nsamples, compression_factor, delay, fine_sampling_factor, oversampling,            &
                                time, ra, dec, pa, chop, masked, removed, method, detector_mask, nrows, ncolumns, ndetectors,      &
                                detector_center, detector_corner, detector_area, distortion_yz, npixels_per_sample, header,        &
                                pmatrix, new_npixels_per_sample, status)

    use iso_fortran_env,        only : ERROR_UNIT
    use module_fitstools,       only : ft_read_keyword
    use module_pacsinstrument,  only : NEAREST_NEIGHBOUR, SHARP_EDGES, SAMPLING_PERIOD, PacsInstrument
    use module_observation,     only : Observation
    use module_pointingmatrix,  only : PointingElement
    use module_tamasis,         only : p
    implicit none

    !f2py integer*8, depend(npixels_per_sample,nvalids,ndetectors) :: pmatrix(npixels_per_sample*nvalids*ndetectors)

    character(len=*), intent(in)                   :: band
    integer*8, intent(in)                          :: nvalids, nsamples
    integer, intent(in)                            :: compression_factor
    real(p), intent(in)                            :: delay
    integer, intent(in)                            :: fine_sampling_factor
    logical, intent(in)                            :: oversampling
    real(p), intent(in), dimension(nsamples)       :: time, ra, dec, pa, chop
    logical*1, intent(in), dimension(nsamples)     :: masked, removed
    character(len=*), intent(in)                   :: method
    logical*1, intent(in)                          :: detector_mask(nrows,ncolumns)
    integer, intent(in)                            :: nrows
    integer, intent(in)                            :: ncolumns
    integer, intent(in)                            :: ndetectors
    real(p), intent(in)                            :: detector_center(2,nrows,ncolumns)
    real(p), intent(in)                            :: detector_corner(2,4,nrows,ncolumns)
    real(p), intent(in)                            :: detector_area(nrows,ncolumns)
    real(p), intent(in)                            :: distortion_yz(2,3,3,3)
    integer, intent(in)                            :: npixels_per_sample
    character(len=*), intent(in)                   :: header
    type(PointingElement), intent(inout)           :: pmatrix(npixels_per_sample,nvalids,ndetectors)
    integer, intent(out)                           :: new_npixels_per_sample
    integer, intent(out)                           :: status

    class(Observation), allocatable     :: obs
    class(PacsInstrument), allocatable  :: pacs
    integer                             :: nx, ny, nvalids_expected, method_

    ! initialise observations
    allocate(obs)
    call obs%init(time, ra, dec, pa, chop, masked, removed, compression_factor,                                                    &
         (compression_factor - 1) / (2._p * compression_factor) +  delay / (1000 * SAMPLING_PERIOD * compression_factor), status)
    if (status /= 0) return

    ! initialise pacs instrument
    allocate (pacs)
    call pacs%init_with_variables(band, detector_mask, fine_sampling_factor, status, detector_center_all=detector_center,          &
                                  detector_corner_all=detector_corner, detector_area_all=detector_area, distortion_yz=distortion_yz)
    if (status /= 0) return

    ! check number of detectors
    if (pacs%ndetectors /= ndetectors) then
        status = 1
        write (ERROR_UNIT,'(a,2(i0,a))') "The specified number of detectors '", ndetectors, "' is incompatible with that from the i&
              &nput detector mask '", pacs%ndetectors, "'."
        return
    end if

    ! check number of fine samples
    if (oversampling) then
        nvalids_expected = obs%nvalids * compression_factor * fine_sampling_factor
    else
        nvalids_expected = obs%nvalids
    end if

    if (nvalids /= nvalids_expected) then
        status = 1
        write (ERROR_UNIT,'(a,2(i0,a))') "The specified total number of samples '", nvalids, "' is incompatible with that from the&
              & observations '", nvalids_expected, "'."
        return
    end if

    ! get the size of the map
    call ft_read_keyword(header, 'naxis1', nx, status=status)
    if (status /= 0) return
    call ft_read_keyword(header, 'naxis2', ny, status=status)
    if (status /= 0) return

    ! get method id
    select case (method)

        case ('nearest')
            method_ = NEAREST_NEIGHBOUR

        case ('sharp')
            method_ = SHARP_EDGES

        case default
            write (ERROR_UNIT, '(a)') 'Unknown projection method ' // method // '.'
            status = 1
            return
    
    end select

    ! compute the projector
    call pacs%compute_projection(method_, obs, oversampling, header, nx, ny, pmatrix, new_npixels_per_sample, status)

end subroutine pacs_pointing_matrix


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_uv2yz(uv, ncoords, distortion_yz, chop, yz)
    ! Convert coordinates in the (u,v) plane into the (y,z) plane, assuming a chop angle.
    ! Input units are in millimeters, output units in degrees.

    use module_pacsinstrument, only : uv2yz
    use module_tamasis,        only : p
    implicit none

    real(p), intent(in)    :: uv(2,ncoords)
    integer, intent(in)    :: ncoords
    real(p), intent(in)    :: distortion_yz(2,3,3,3)
    real(p), intent(in)    :: chop
    real(p), intent(inout) :: yz(2,ncoords)

    yz = uv2yz(uv, distortion_yz, chop)

end subroutine pacs_uv2yz


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_uv2ad(uv, ncoords, ra, dec, pa, chop, nsamples, distortion_yz, ad)
    ! Convert coordinates in the (u,v) plane into celestial coordinates,
    ! assuming a pointing direction and a position angle.
    ! The routine is not accurate at the poles.
    ! Input units are in millimeters, output units in degrees.

    use module_pacsinstrument, only : uv2yz, yz2ad
    use module_tamasis,        only : p
    implicit none

    integer, intent(in)    :: ncoords, nsamples
    real(p), intent(in)    :: uv(2,ncoords), ra(nsamples), dec(nsamples), pa(nsamples), chop(nsamples)
    real(p), intent(in)    :: distortion_yz(2,3,3,3)
    real(p), intent(inout) :: ad(2,ncoords,nsamples)

    integer :: isample

    !$omp parallel do
    do isample = 1, nsamples
        ad(:,:,isample) = uv2yz(uv, distortion_yz, chop(isample))
        ad(:,:,isample) = yz2ad(ad(:,:,isample), ra(isample), dec(isample), pa(isample))
    end do
    !$omp end parallel do

end subroutine pacs_uv2ad


!-----------------------------------------------------------------------------------------------------------------------------------


module pacs

    use module_fitstools,      only : ft_read_keyword
    use module_math,           only : mInf, pInf
    use module_pacsinstrument, only : uv2yz, yz2ad
    use module_pointingmatrix, only : PointingElement, xy2pmatrix, xy2roi, roi2pmatrix
    use module_projection,     only : convex_hull
    use module_tamasis,        only : p
    use module_wcs,            only : ad2xy_gnomonic, ad2xys_gnomonic, init_astrometry
    implicit none
    private

    public :: uv2xy_minmax
    public :: uv2pmatrix_nearest_neighbour
    public :: uv2pmatrix_sharp_edges

contains

    subroutine uv2xy_minmax(coords, ncoords, ra, dec, pa, chop, npointings, distortion_yz, header, xmin, ymin, xmax, ymax, status)
        ! Return the minimum and maximum sky pixel coordinate values for a set of coordinates in the (u,v) plane.

        real(p), intent(in)                        :: coords(2,ncoords)       ! instrument frame coordinates
        real(p), intent(in), dimension(npointings) :: ra, dec, pa, chop       ! input pointings in celestial coordinates
        integer, intent(in)                        :: ncoords, npointings     ! #coordinates, #detectors
        real(p), intent(in)                        :: distortion_yz(2,3,3,3)  ! distortion coefficients
        character(len=2880), intent(in)            :: header                  ! input FITS header
        real(p), intent(out)                       :: xmin, ymin, xmax, ymax  ! min and max values of the map coordinates
        integer, intent(out)                       :: status                  ! status flag

        real(p), allocatable :: hull(:,:)
        integer, allocatable :: ihull(:)
        integer              :: ipointing
        external :: pacs_uv2ad

        call init_astrometry(header, status=status)
        if (status /= 0) return

        call convex_hull(coords, ihull)
        allocate (hull(2, size(ihull)))

        xmin = pInf
        xmax = mInf
        ymin = pInf
        ymax = mInf

        !$omp parallel do reduction(min:xmin,ymin) reduction(max:xmax,ymax) private(hull)
        do ipointing = 1, npointings

            hull = uv2yz(coords(:,ihull), distortion_yz, chop(ipointing))
            hull = yz2ad(hull, ra(ipointing), dec(ipointing), pa(ipointing))
            hull = ad2xy_gnomonic(hull)
            xmin = min(xmin, minval(hull(1,:)))
            ymin = min(ymin, minval(hull(2,:)))
            xmax = max(xmax, maxval(hull(1,:)))
            ymax = max(ymax, maxval(hull(2,:)))

        end do
        !$omp end parallel do

    end subroutine uv2xy_minmax


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine uv2pmatrix_nearest_neighbour(coords, ncoords, area, ra, dec, pa, chop, masked, npointings, distortion_yz, header,   &
                                            pmatrix, out, status)
        !f2py integer*8, depend(npointings,ncoords) :: pmatrix(npointings*ncoords)
        real(p), intent(in)                          :: coords(2,ncoords)      ! instrument frame coordinates
        real(p), intent(in)                          :: area(ncoords)          ! detector area / reference_area
        real(p), intent(in), dimension(npointings)   :: ra, dec, pa, chop      ! input pointings in celestial coordinates
        logical*1, intent(in), dimension(npointings) :: masked                 ! pointing flags: true if masked, removed
        integer, intent(in)                          :: ncoords, npointings    ! #coordinates, #pointings
        real(p), intent(in)                          :: distortion_yz(2,3,3,3) ! distortion coefficients
        character(len=*), intent(in)                 :: header
        type(PointingElement), intent(inout)         :: pmatrix(1,npointings,ncoords)
        logical, intent(out)                         :: out
        integer, intent(out)                         :: status

        real(p) :: coords2(2,ncoords), x(ncoords), y(ncoords), s(ncoords)
        integer :: ipointing, nx, ny

        out = .false.
        call init_astrometry(header, status=status)
        if (status /= 0) return
        ! get the size of the map
        call ft_read_keyword(header, 'naxis1', nx, status=status)
        if (status /= 0) return
        call ft_read_keyword(header, 'naxis2', ny, status=status)
        if (status /= 0) return

        !$omp parallel do private(coords2, x, y, s) reduction(.or. : out)
        ! loop over the samples which have not been removed
        do ipointing = 1, npointings

            if (masked(ipointing)) then
                pmatrix(1,ipointing,:)%pixel = -1
                pmatrix(1,ipointing,:)%weight = 0
                cycle
            end if

            coords2 = uv2yz(coords, distortion_yz, chop(ipointing))
            coords2 = yz2ad(coords2, ra(ipointing), dec(ipointing), pa(ipointing))
            call ad2xys_gnomonic(coords2, x, y, s)

            ! the input map has flux densities, not surface brightness
            ! f_idetector = f_imap * weight
            ! with weight = detector_area / pixel_area
            ! and pixel_area = reference_area / s
            call xy2pmatrix(x, y, nx, ny, out, pmatrix(1,ipointing,:))
            pmatrix(1,ipointing,:)%weight = s * area

        end do
        !$omp end parallel do

    end subroutine uv2pmatrix_nearest_neighbour


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine uv2pmatrix_sharp_edges(coords, ncoords, ra, dec, pa, chop, masked, npointings, distortion_yz, header, pmatrix,      &
                                      npixels_per_sample, new_npixels_per_sample, out, status)
        !f2py integer*8, depend(npixels_per_sample,npointings,ncoords) :: pmatrix(npixels_per_sample*npointings*ncoords/4)
        real(p), intent(in)                          :: coords(2,ncoords)      ! instrument frame coordinates
        real(p), intent(in), dimension(npointings)   :: ra, dec, pa, chop      ! input pointings in celestial coordinates
        logical*1, intent(in), dimension(npointings) :: masked                 ! pointing flags: true if masked, removed
        integer, intent(in)                          :: ncoords, npointings    ! #coordinates, #pointings
        real(p), intent(in)                          :: distortion_yz(2,3,3,3) ! distortion coefficients
        character(len=*), intent(in)                 :: header                 ! sky map FITS header
        type(PointingElement), intent(inout)         :: pmatrix(npixels_per_sample,npointings,ncoords/4) ! the pointing matrix
        integer, intent(in)  :: npixels_per_sample     ! input maximum number of sky pixels intersected by a detector
        integer, intent(out) :: new_npixels_per_sample ! actual maximum number of sky pixels intersected by a detector
        logical, intent(out) :: out                    ! true if some coordinates fall outside of the map
        integer, intent(out) :: status

        real(p) :: coords2(2,ncoords)
        integer :: roi(2,2,ncoords/4)
        integer :: ipointing, nx, ny

        new_npixels_per_sample = 0
        out = .false.
        call init_astrometry(header, status=status)
        if (status /= 0) return
        ! get the size of the map
        call ft_read_keyword(header, 'naxis1', nx, status=status)
        if (status /= 0) return
        call ft_read_keyword(header, 'naxis2', ny, status=status)
        if (status /= 0) return

        !$omp parallel do private(ipointing, coords2, roi) &
        !$omp reduction(max : npixels_per_sample) reduction(.or. : out)
        do ipointing = 1, npointings

            if (masked(ipointing)) then
                pmatrix(:,ipointing,:)%pixel = -1
                pmatrix(:,ipointing,:)%weight = 0
                cycle
            end if

            coords2 = uv2yz(coords, distortion_yz, chop(ipointing))
            coords2 = yz2ad(coords2, ra(ipointing), dec(ipointing), pa(ipointing))
            coords2 = ad2xy_gnomonic(coords2)
            roi = xy2roi(coords2, 4)
            call roi2pmatrix(roi, 4, coords2, nx, ny, new_npixels_per_sample, out, pmatrix(:,ipointing,:))

        end do
        !$omp end parallel do

    end subroutine uv2pmatrix_sharp_edges

end module pacs

!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_multiplexing_direct(signal, multiplexed, fine_sampling_factor, ij, nsamples, ndetectors)

    use module_pacsinstrument, only : multiplexing_direct
    use module_tamasis,        only : p
    implicit none

    real(p), intent(in)    :: signal(nsamples, ndetectors)
    integer*8, intent(in)  :: nsamples
    integer, intent(in)    :: ndetectors, fine_sampling_factor, ij(2, ndetectors)
    real(p), intent(inout) :: multiplexed(nsamples/fine_sampling_factor, ndetectors)

    call multiplexing_direct(signal, multiplexed, fine_sampling_factor, ij)

end subroutine pacs_multiplexing_direct


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_multiplexing_transpose(multiplexed, signal, fine_sampling_factor, ij, nsamples, ndetectors)

    use module_pacsinstrument, only : multiplexing_transpose
    use module_tamasis,        only : p
    implicit none
    
    integer*8, intent(in)  :: nsamples
    integer, intent(in)    :: ndetectors, fine_sampling_factor, ij(2, ndetectors)
    real(p), intent(inout) :: signal(nsamples, ndetectors)
    real(p), intent(in)    :: multiplexed(nsamples/fine_sampling_factor, ndetectors)

    call multiplexing_transpose(multiplexed, signal, fine_sampling_factor, ij)

end subroutine pacs_multiplexing_transpose


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_bitmask(mask, nsamples, ndetectors, bitmask)

    implicit none

    logical*1, intent(in)  :: mask(nsamples, ndetectors)
    integer, intent(in)    :: nsamples, ndetectors
    integer*4, intent(out) :: bitmask((nsamples-1)/32+1,ndetectors)

    integer :: idetector, isample

    bitmask = 0

    !$omp parallel do
    do idetector = 1, ndetectors
        do isample = 1, nsamples
            if (mask(isample,idetector)) then
                bitmask((isample-1)/32+1,idetector) = ibset(bitmask((isample-1)/32+1,idetector), modulo(isample-1,32))
            end if
        end do
    end do
    !$omp end parallel do

end subroutine pacs_bitmask
