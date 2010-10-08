! Tamasis interface for f2py
!
! Author: P. Chanial

subroutine pacs_info_band(filename, nfilenames, channel, transparent_mode, status)

    use iso_fortran_env,        only : ERROR_UNIT
    use module_fitstools,       only : ft_close, ft_open, ft_read_keyword
    use module_pacsobservation, only : PacsObservation, MaskPolicy
    use module_string,          only : strlowcase
    implicit none

    !f2py threadsafe
    !f2py intent(in)   filename
    !f2py intent(in)   nfilenames
    !f2py intent(out)  channel
    !f2py intent(out)  status

    character(len=*), intent(in)  :: filename
    integer, intent(in)           :: nfilenames
    character(len=7), intent(out) :: channel
    logical, intent(out)          :: transparent_mode
    integer, intent(out)          :: status

    character(len=len(filename)/nfilenames) :: filename_
    integer                                 :: iobs, length, pos, unit
    character(len=10)                       :: obstype
    character(len=7)                        :: channel_

    ! split input filename
    if (mod(len(filename), nfilenames) /= 0) then
        stop 'PACS_INFO_CHANNEL: Invalid filename length.'
    end if

    channel = 'unknown'
    transparent_mode = .true.

    do iobs = 1, nfilenames

        filename_ = filename((iobs-1)*len(filename_)+1:iobs*len(filename_))

        ! remove trailing section
        pos = index(filename_, '[', back=.true.)
        if (pos > 0) filename_(pos:) = ' '

        ! old style file format
        length = len_trim(filename_)
        if (filename_(length-4:length) /= '.fits') then
            status = 0
            if (strlowcase(filename_(length-3:length)) == 'blue') then
               channel_ = 'blue'
            else if (strlowcase(filename_(length-4:length)) == 'green') then
               channel_ = 'green'
            else if (strlowcase(filename_(length-2:length)) == 'red') then
               channel_ = 'red'
            else
               status = 1
               write (ERROR_UNIT,'(a)') 'File name does not contain the array channel identifier (blue, green, red).'
            end if
            return
        endif

        call ft_open(trim(filename_), unit, status)
        if (status /= 0) return

        call ft_read_keyword(unit, 'TYPE', obstype, status=status)
        if (status /= 0) return

        call ft_close(unit, status)
        if (status /= 0) return

        select case (trim(obstype))
            case ('HPPRAWBS')
                channel_ = 'blue'
            case ('HPPRAWBL')
                channel_ = 'green'
            case ('HPPRAWRS')
                channel_ = 'red'
            case ('HPPAVGBS')
                channel_ = 'blue'
                transparent_mode = .false.
            case ('HPPAVGBL')
                channel_ = 'green'
                transparent_mode = .false.
            case ('HPPAVGRS')
                channel_ = 'red'
                transparent_mode = .false.
            case default
                write (ERROR_UNIT, '(a)') "Unknown observation type '" // trim(obstype) // "' in file '" // filename_ // "'."
                status = 1
                return
        end select
        
        if (iobs == 1) then
            channel = channel_
            cycle
        end if

        if (channel_ /= channel) then
            write (ERROR_UNIT, '(a)') "Error: Observations are not performed with the same channel."
            status = 1
            return
        end if

    end do

end subroutine pacs_info_band


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_info_bad_detector_mask(tamasis_dir, channel, transparent_mode, reject_bad_line, nrows, ncolumns, detector_mask,    &
                                       status)

    use iso_fortran_env,       only : ERROR_UNIT
    use module_fitstools,      only : ft_read_image
    use module_pacsinstrument, only : FILENAME_BPM, get_calfile
    use module_tamasis,        only : init_tamasis

    !f2py threadsafe
    !f2py intent(in)   tamasis_dir
    !f2py intent(in)   channel
    !f2py intent(in)   transparent_mode
    !f2py intent(in)   reject_bad_line
    !f2py intent(in)   nrows, ncolumns
    !f2py intent(out)  detector_mask
    !f2py intent(out)  status

    character(len=*), intent(in) :: tamasis_dir
    character(len=*), intent(in) :: channel
    logical, intent(in)          :: transparent_mode
    logical, intent(in)          :: reject_bad_line
    integer, intent(in)          :: nrows, ncolumns
    logical*1, intent(out)       :: detector_mask(nrows,ncolumns)
    integer, intent(out)         :: status

    logical*1, allocatable       :: tmp(:,:)

    ! initialise tamasis
    call init_tamasis(tamasis_dir)

    ! read bad pixel mask
    call ft_read_image(get_calfile(FILENAME_BPM) // '[' // channel // ']', tmp, status)
    if (status /= 0) return

    if (size(tmp,1) /= size(detector_mask,2) .or. size(tmp,2) /= size(detector_mask,1)) then
        status = 1
        write (ERROR_UNIT, '(a,4(i0,a))') 'Invalid shape of the detector mask (', size(detector_mask,1), ',',                      &
              size(detector_mask,2), ') instead of (', size(tmp,2), ',', size(tmp,1), ')'
        return
    end if
    detector_mask = transpose(tmp)

    ! mask detectors rejected in transparent mode
    if (transparent_mode) then
        if (channel /= 'red') then
            detector_mask(1:16,1:16) = .true.
            detector_mask(1:16,33:)  = .true.
            detector_mask(17:,:)     = .true.
        else
            detector_mask(1:8,1:8) = .true.
            detector_mask(1:8,17:) = .true.
            detector_mask(9:,:)    = .true.
        end if
    end if
    
    ! mask erratic line
    if (reject_bad_line .and. channel /= 'red') then
        detector_mask(12,17:32) = .true.
    end if

end subroutine pacs_info_bad_detector_mask


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_info(tamasis_dir, filename, nfilenames, transparent_mode, fine_sampling_factor, frame_policy, detector_mask,       &
                     nrows, ncolumns, compression_factor, nsamples, unit, responsivity, detector_area, flatfield_detector,         &
                     flatfield_optical, status)

    use iso_fortran_env,        only : OUTPUT_UNIT
    use module_pacsinstrument,  only : PacsInstrument
    use module_pacsobservation, only : PacsObservation, MaskPolicy
    use module_tamasis,         only : init_tamasis, p
    implicit none

    !f2py threadsafe
    !f2py intent(in)   tamasis_dir
    !f2py intent(in)   filename
    !f2py intent(in)   nfilenames
    !f2py intent(in)   transparent_mode
    !f2py intent(in)   fine_sampling_factor
    !f2py intent(in)   frame_policy
    !f2py intent(in)   detector_mask
    !f2py intent(hide) nrows = shape(bad_detector_mask,0)
    !f2py intent(hide) ncolumns = shape(bad_detector_mask,1)
    !f2py intent(out)  compression_factor
    !f2py intent(out)  nsamples
    !f2py intent(out)  unit
    !f2py intent(out)  responsivity
    !f2py intent(out)  detector_area
    !f2py intent(out)  flatfield_detector
    !f2py intent(out)  flatfield_optical
    !f2py intent(out)  status

    character(len=*), intent(in)   :: tamasis_dir
    character(len=*), intent(in)   :: filename
    integer, intent(in)            :: nfilenames
    logical, intent(in)            :: transparent_mode
    integer, intent(in)            :: fine_sampling_factor
    integer, intent(in)            :: frame_policy(4)
    logical*1, intent(in)          :: detector_mask(nrows,ncolumns)
    integer, intent(in)            :: nrows, ncolumns
    integer, intent(out)           :: compression_factor(nfilenames)
    integer*8, intent(out)         :: nsamples(nfilenames)
    character(len=70), intent(out) :: unit
    real(p), intent(out)           :: responsivity
    real(p), intent(out)           :: detector_area(nrows,ncolumns)
    real(p), intent(out)           :: flatfield_detector(nrows,ncolumns)
    real(p), intent(out)           :: flatfield_optical(nrows,ncolumns)
    integer, intent(out)           :: status

    character(len=len(filename)/nfilenames), allocatable :: filename_(:)
    class(PacsObservation), allocatable :: obs
    class(PacsInstrument), allocatable  :: pacs
    type(MaskPolicy)                    :: policy
    integer                             :: iobs

    ! initialise tamasis
    call init_tamasis(tamasis_dir)

    ! split input filename
    if (mod(len(filename), nfilenames) /= 0) then
        stop 'PACS_INFO: Invalid filename length.'
    end if
    allocate(filename_(nfilenames))
    do iobs = 1, nfilenames
        filename_(iobs) = filename((iobs-1)*len(filename_)+1:iobs*len(filename_))
    end do

    ! initialise policy
    policy = MaskPolicy(inscan=frame_policy(1), turnaround=frame_policy(2), other=frame_policy(3), invalid=frame_policy(4))

    allocate(obs)
    call obs%init(filename_, policy, status, verbose=.true.)
    if (status /= 0) return

    allocate(pacs)
    call pacs%init(obs%channel, transparent_mode, fine_sampling_factor, detector_mask, status)
    if (status /= 0) return

    call obs%print()
    write (OUTPUT_UNIT,*)

    compression_factor = obs%slice%compression_factor
    nsamples = obs%slice%nvalids
    unit = obs%unit
    responsivity = pacs%responsivity
    detector_area = pacs%detector_area_all
    flatfield_detector = pacs%flatfield_detector
    flatfield_optical = pacs%flatfield_optical

end subroutine pacs_info


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_read_filter_calibration_ncorrelations(tamasis_dir, channel, ncorrelations, status)

    use module_pacsinstrument, only : read_filter_calibration_ncorrelations
    use module_tamasis,        only : init_tamasis
    implicit none

    !f2py intent(in)   tamasis_dir
    !f2py intent(in)   channel
    !f2py intent(out)  ncorrelations
    !f2py intent(out)  status
    
    character(len=*), intent(in) :: tamasis_dir
    character(len=*), intent(in) :: channel
    integer, intent(out)         :: ncorrelations
    integer, intent(out)         :: status
    
    call init_tamasis(tamasis_dir)
    ncorrelations = read_filter_calibration_ncorrelations(channel, status)

end subroutine pacs_read_filter_calibration_ncorrelations


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_read_filter_calibration(tamasis_dir, channel, ncorrelations, ndetectors, mask, nrows, ncolumns, data, status)

    use iso_fortran_env,       only : ERROR_UNIT
    use module_filtering,      only : FilterUncorrelated
    use module_pacsinstrument, only : read_filter_calibration
    use module_tamasis,        only : init_tamasis, p
    implicit none

    !f2py threadsafe
    !f2py intent(in)   tamasis_dir
    !f2py intent(in)   channel
    !f2py intent(in)   ncorrelations
    !f2py intent(in)   ndetectors
    !f2py intent(in)   mask
    !f2py intent(hide) nrows = shape(mask,0)
    !f2py intent(hide) ncolumns = shape(mask,1)
    !f2py intent(out)  data
    !f2py intent(out)  status

    character(len=*), intent(in) :: tamasis_dir
    character(len=*), intent(in) :: channel
    integer, intent(in)          :: ncorrelations
    integer, intent(in)          :: ndetectors
    logical*1, intent(in)        :: mask(nrows,ncolumns)
    integer, intent(in)          :: nrows
    integer, intent(in)          :: ncolumns
    real(p), intent(out)         :: data(ncorrelations+1,ndetectors)
    integer, intent(out)         :: status

    type(FilterUncorrelated)     :: filter

    call init_tamasis(tamasis_dir)
    call read_filter_calibration(channel, mask, filter, status)
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

end subroutine pacs_read_filter_calibration


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_map_header(tamasis_dir, filename, nfilenames, oversampling, fine_sampling_factor, frame_policy, detector_mask,     &
                           nrows, ncolumns, resolution, header, status)

    use module_pacsinstrument,  only : PacsInstrument
    use module_pacsobservation, only : PacsObservation, MaskPolicy
    use module_tamasis,         only : p, init_tamasis
    implicit none

    !f2py threadsafe
    !f2py intent(in)   tamasis_dir
    !f2py intent(in)   filename
    !f2py intent(in)   nfilenames
    !f2py intent(in)   oversampling
    !f2py intent(in)   fine_sampling_factor
    !f2py intent(in)   frame_policy
    !f2py intent(in)   detector_mask
    !f2py intent(hide) nrows = shape(bad_detector_mask,0)
    !f2py intent(hide) ncolumns = shape(bad_detector_mask,1)
    !f2py intent(in)   resolution
    !f2py intent(out)  header
    !f2py intent(out)  status

    character(len=*), intent(in) :: tamasis_dir
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: nfilenames
    logical, intent(in)          :: oversampling
    integer, intent(in)          :: fine_sampling_factor
    integer, intent(in)          :: frame_policy(4)
    logical*1, intent(in)        :: detector_mask(nrows,ncolumns)
    integer, intent(in)          :: nrows, ncolumns
    real(p), intent(in)          :: resolution
    character(len=2880), intent(out) :: header
    integer, intent(out)             :: status

    character(len=len(filename)/nfilenames), allocatable :: filename_(:)
    class(PacsObservation), allocatable :: obs
    class(PacsInstrument), allocatable  :: pacs
    type(MaskPolicy)                    :: policy
    integer                             :: iobs
    
    ! initialise tamasis
    call init_tamasis(tamasis_dir)

    ! split input filename
    if (mod(len(filename), nfilenames) /= 0) then
        stop 'PACS_MAP_HEADER: Invalid filename length.'
    end if
    allocate(filename_(nfilenames))
    do iobs = 1, nfilenames
        filename_(iobs) = filename((iobs-1)*len(filename_)+1:iobs*len(filename_))
    end do

    ! initialise policy
    policy = MaskPolicy(inscan=frame_policy(1), turnaround=frame_policy(2), other=frame_policy(3), invalid=frame_policy(4))

    allocate(obs)
    call obs%init(filename_, policy, status)
    if (status /= 0) go to 999

    allocate(pacs)
    call pacs%init(obs%channel, obs%slice(1)%observing_mode == 'Transparent', fine_sampling_factor, detector_mask, status)
    if (status /= 0) go to 999
    
    call pacs%compute_map_header(obs, oversampling, resolution, header, status)
    
999 continue

end subroutine pacs_map_header


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_timeline(tamasis_dir, filename, nfilenames, nsamples, ndetectors, frame_policy, detector_mask, nrow, ncol,         &
                         do_flatfielding, do_subtraction_mean, signal, mask, status)

    use iso_fortran_env,        only : ERROR_UNIT
    use module_pacsinstrument,  only : PacsInstrument
    use module_pacsobservation, only : PacsObservation, MaskPolicy
    use module_preprocessor,    only : divide_vectordim2, subtract_meandim1
    use module_tamasis,         only : init_tamasis, p
    implicit none

    !f2py threadsafe
    !f2py intent(in)   :: tamasis_dir
    !f2py intent(in)   :: filename
    !f2py intent(in)   :: nfilenames
    !f2py intent(in)   :: nsamples
    !f2py intent(in)   :: ndetectors
    !f2py intent(in)   :: frame_policy
    !f2py intent(in)   :: detector_mask
    !f2py intent(hide) :: nrow = shape(bad_detector_mask,0)
    !f2py intent(hide) :: ncol = shape(bad_detector_mask,1)
    !f2py intent(in)   :: do_flatfielding
    !f2py intent(in)   :: do_subtraction_mean
    !f2py intent(out)  :: signal
    !f2py intent(out)  :: mask
    !f2py intent(out)  :: status

    character(len=*), intent(in) :: tamasis_dir
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: nfilenames
    integer, intent(in)          :: nsamples
    integer, intent(in)          :: ndetectors
    integer, intent(in)          :: frame_policy(4)
    logical*1, intent(in)        :: detector_mask(nrow,ncol)
    integer, intent(in)          :: nrow, ncol
    logical, intent(in)          :: do_flatfielding, do_subtraction_mean
    real(p), intent(out)         :: signal(nsamples, ndetectors)
    logical*1, intent(out)       :: mask(nsamples, ndetectors)
    integer, intent(out)         :: status

    character(len=len(filename)/nfilenames), allocatable :: filename_(:)
    class(PacsObservation), allocatable :: obs
    class(PacsInstrument), allocatable  :: pacs
    type(MaskPolicy)                    :: policy
    integer                             :: iobs, destination

    ! initialise tamasis
    call init_tamasis(tamasis_dir)

    ! split input filename
    if (mod(len(filename), nfilenames) /= 0) then
        stop 'PACS_TIMELINE: Invalid filename length.'
    end if
    allocate(filename_(nfilenames))
    do iobs = 1, nfilenames
        filename_(iobs) = filename((iobs-1)*len(filename_)+1:iobs*len(filename_))
    end do

    ! initialise policy
    policy = MaskPolicy(inscan=frame_policy(1), turnaround=frame_policy(2), other=frame_policy(3), invalid=frame_policy(4))

    ! initialise observations
    allocate(obs)
    call obs%init(filename_, policy, status)
    if (status /= 0) go to 999

    ! initialise pacs instrument
    allocate(pacs)
    call pacs%init(obs%channel, obs%slice(1)%observing_mode == 'Transparent', 1, detector_mask, status)
    if (status /= 0) go to 999

    ! read timeline
    call pacs%read(obs, signal, mask, status, verbose=.true.)
    if (status /= 0) go to 999

    ! flat fielding
    if (do_flatfielding) then
        call divide_vectordim2(signal, pack(pacs%flatfield_detector, .not. pacs%mask))
    end if
    
    ! subtract the mean of each detector timeline
    if (do_subtraction_mean) then
        destination = 1
        do iobs=1, nfilenames
            call subtract_meandim1(signal(destination:destination+obs%slice(iobs)%nvalids-1,:),                                    &
                                   mask(destination:destination+obs%slice(iobs)%nvalids-1,:))
            destination = destination + obs%slice(iobs)%nvalids
         end do
    end if

999 continue

end subroutine pacs_timeline


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_pointing_matrix_filename(tamasis_dir, filename, nfilenames, method, oversampling, fine_sampling_factor,            &
                                         npixels_per_sample, nsamples, ndetectors, frame_policy, detector_mask, nrow, ncol, header,&
                                         pmatrix, status)

    use iso_fortran_env,        only : ERROR_UNIT
    use module_fitstools,       only : ft_read_keyword
    use module_pacsinstrument,  only : PacsInstrument, NEAREST_NEIGHBOUR, SHARP_EDGES
    use module_pacsobservation, only : PacsObservation, MaskPolicy
    use module_pointingmatrix,  only : PointingElement
    use module_tamasis,         only : init_tamasis
    implicit none

    !f2py threadsafe
    !f2py intent(in)   :: tamasis_dir
    !f2py intent(in)   :: filename
    !f2py intent(in)   :: nfilenames
    !f2py intent(in)   :: method
    !f2py intent(in)   :: oversampling
    !f2py intent(in)   :: fine_sampling_factor
    !f2py intent(in)   :: npixels_per_sample
    !f2py intent(in)   :: nsamples
    !f2py intent(in)   :: ndetectors
    !f2py intent(in)   :: frame_policy
    !f2py intent(in)   :: detector_mask
    !f2py intent(hide) :: nrow = shape(bad_detector_mask,0)
    !f2py intent(hide) :: ncol = shape(bad_detector_mask,1)
    !f2py intent(in)   :: header
    !f2py integer*8, intent(inout) :: pmatrix(npixels_per_sample*nsamples*ndetectors)
    !f2py intent(out)  :: status

    character(len=*), intent(in) :: tamasis_dir
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: nfilenames
    character(len=*), intent(in) :: method
    logical, intent(in)          :: oversampling
    integer, intent(in)          :: fine_sampling_factor
    integer, intent(in)          :: npixels_per_sample
    integer*8, intent(in)        :: nsamples
    integer, intent(in)          :: ndetectors
    integer, intent(in)          :: frame_policy(4)
    logical*1, intent(in)        :: detector_mask(nrow,ncol)
    integer, intent(in)          :: nrow, ncol
    character(len=*), intent(in) :: header
    type(PointingElement), intent(inout) :: pmatrix(npixels_per_sample,nsamples,ndetectors)
    integer, intent(out)         :: status

    character(len=len(filename)/nfilenames), allocatable :: filename_(:)
    class(PacsObservation), allocatable :: obs
    class(PacsInstrument), allocatable  :: pacs
    type(MaskPolicy)                    :: policy
    integer                             :: iobs, nx, ny, nsamples_expected, method_

    ! initialise tamasis
    call init_tamasis(tamasis_dir)

    ! split input filename
    if (mod(len(filename), nfilenames) /= 0) then
        stop 'PACS_POINTING_MATRIX_FILENAME: Invalid filename length.'
    end if
    allocate(filename_(nfilenames))
    do iobs = 1, nfilenames
        filename_(iobs) = filename((iobs-1)*len(filename_)+1:iobs*len(filename_))
    end do

    ! initialise policy
    policy = MaskPolicy(inscan=frame_policy(1), turnaround=frame_policy(2), other=frame_policy(3), invalid=frame_policy(4))

    ! initialise observations
    allocate(obs)
    call obs%init(filename_, policy, status)
    if (status /= 0) return

    ! initialise pacs instrument
    allocate(pacs)
    call pacs%init(obs%channel, obs%slice(1)%observing_mode == 'Transparent', fine_sampling_factor, detector_mask, status)
    if (status /= 0) return

    ! check number of detectors
    if (pacs%ndetectors /= ndetectors) then
        status = 1
        write (ERROR_UNIT,'(a,2(i0,a))') "The specified number of detectors '", ndetectors, "' is incompatible with that from the i&
              &nput bad detector mask '", pacs%ndetectors, "'."
        return
    end if

    ! check number of fine samples
    if (oversampling) then
        nsamples_expected = sum(obs%slice%nvalids * obs%slice%compression_factor)*fine_sampling_factor
    else
        nsamples_expected = sum(obs%slice%nvalids)
    end if
    if (nsamples /= nsamples_expected) then
        status = 1
        write (ERROR_UNIT,'(a,2(i0,a))') "The specified total number of samples '", nsamples, "' is incompatible with that from the&
              & observations '", nsamples_expected, "'."
        return
    end if

    ! get the size of the map
    call ft_read_keyword(header, 'naxis1', nx, status=status)
    if (status /= 0) return
    call ft_read_keyword(header, 'naxis2', ny, status=status)
    if (status /= 0) return

    ! get method id
    select case (method)

        case ('nearest neighbour')
            method_ = NEAREST_NEIGHBOUR

        case ('sharp edges')
            method_ = SHARP_EDGES

        case default
            write (ERROR_UNIT, '(a)') 'Unknown projection method ' // method
            status = 1
            return
    
    end select

    ! compute the projector
    call pacs%compute_projection(method_, obs, oversampling, header, nx, ny, pmatrix, status)

end subroutine pacs_pointing_matrix_filename


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine pacs_multiplexing_direct(signal, multiplexed, fine_sampling_factor, ij, nsamples, ndetectors)

    use module_pacsinstrument, only : multiplexing_direct
    use module_tamasis,        only : p
    implicit none

    !f2py threadsafe
    !f2py intent(in)       :: signal
    !f2py intent(inout)    :: multiplexed
    !f2py intent(in)       :: fine_sampling_factor
    !f2py intent(in)       :: ij
    !f2py intent(hide)     :: nsamples = shape(signal,0)
    !f2py intent(hide)     :: ndetectors = shape(signal,1)

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

    !f2py threadsafe
    !f2py intent(in)       :: multiplexed
    !f2py intent(inout)    :: signal
    !f2py intent(in)       :: fine_sampling_factor
    !f2py intent(in)       :: ij
    !f2py intent(hide)     :: nsamples = shape(signal,0)
    !f2py intent(hide)     :: ndetectors = shape(signal,1)
    
    integer*8, intent(in)  :: nsamples
    integer, intent(in)    :: ndetectors, fine_sampling_factor, ij(2, ndetectors)
    real(p), intent(inout) :: signal(nsamples, ndetectors)
    real(p), intent(in)    :: multiplexed(nsamples/fine_sampling_factor, ndetectors)

    call multiplexing_transpose(multiplexed, signal, fine_sampling_factor, ij)

end subroutine pacs_multiplexing_transpose
