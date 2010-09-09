! Tamasis interface for f2py
!
! Author: P. Chanial

subroutine madmap1_nslices(invnttfile, ndetectors, nslices, status)

    use module_madcap, only : get_nfilters
    implicit none

    !f2py threadsafe
    !f2py intent(in)  :: invnttfile
    !f2py intent(in)  :: ndetectors
    !f2py intent(out) :: nslices
    !f2py intent(out) :: status

    character(len=*), intent(in) :: invnttfile
    integer, intent(in)          :: ndetectors
    integer, intent(out)         :: nslices
    integer, intent(out)         :: status

    integer                      :: nfilterfiles

    nfilterfiles = get_nfilters(invnttfile, ndetectors, status)
    if (status /= 0) return
    nslices = nfilterfiles / ndetectors

end subroutine madmap1_nslices


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine madmap1_info(todfile, invnttfile, convert, ndetectors, nslices, npixels_per_sample, nsamples, ncorrelations, status)

    use iso_fortran_env, only : ERROR_UNIT
    use module_madcap,   only : read_filter_headers, read_tod_header
    use module_string,   only : strinteger
    implicit none

    !f2py threadsafe
    !f2py intent(in)  :: todfile
    !f2py intent(in)  :: invnttfile
    !f2py intent(in)  :: convert
    !f2py intent(in)  :: ndetectors
    !f2py intent(in)  :: nslices
    !f2py intent(out) :: npixels_per_sample
    !f2py intent(out) :: nsamples
    !f2py intent(out) :: ncorrelations
    !f2py intent(out) :: status

    character(len=*), intent(in) :: todfile
    character(len=*), intent(in) :: invnttfile
    character(len=*), intent(in) :: convert
    integer, intent(in)          :: ndetectors
    integer, intent(in)          :: nslices
    integer, intent(out)         :: npixels_per_sample
    integer, intent(out)         :: nsamples(nslices)
    integer, intent(out)         :: ncorrelations
    integer, intent(out)         :: status

    integer*8, allocatable       :: first(:), last(:)
    integer, allocatable         :: ncorrelations_(:)
    integer                      :: idetector, ifilter, islice
    integer*8                    :: nsamples_tot

    call read_tod_header(todfile, convert, ndetectors, nsamples_tot, npixels_per_sample, status)
    if (status /= 0) return

    call read_filter_headers(invnttfile, convert, ndetectors, first, last, ncorrelations_, status)
    if (status /= 0) return

    ncorrelations = ncorrelations_(1)
    ! check number of correlations
    if (any(ncorrelations_ /= ncorrelations)) then
        status = 1
        write (ERROR_UNIT,'(a)') 'Error: Filter files do not have the same correlation length.'
        return
    end if

    ! check number of samples per slice
    do islice=1, nslices
        
        nsamples(islice) = last((islice-1)*ndetectors+1) - first((islice-1)*ndetectors+1) + 1
        
        do idetector = 1, ndetectors
            
            ifilter = (islice-1)*ndetectors + idetector
            if (last(ifilter) - first(ifilter) + 1 /= nsamples(islice)) then
                write (ERROR_UNIT,'(a)') "Error: Filter '" // invnttfile // '.' // strinteger(ifilter-1)// "' does not apply to th&
                     &e same number of samples than filter '" // invnttfile // '.' // strinteger((islice-1)*ndetectors) // "'."
                status = 1
                return
            end if
            
        end do

    end do

    ! check number of samples
    if (nsamples_tot /= sum(nsamples)) then
        status = 1
        write (ERROR_UNIT,'(a,2(i0,a))') "Error: Invalid number of samples in tod ('", nsamples_tot, "' instead of '",    &
             sum(nsamples), "')."
        return
    end if

end subroutine madmap1_info


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine madmap1_read_tod(todfile, invnttfile, convert, npixels_per_sample, nsamples_tot, ndetectors, tod, pmatrix, status)

    use module_filtering,      only : FilterUncorrelated
    use module_madcap,         only : read_filter, read_tod
    use module_pointingmatrix, only : PointingElement
    implicit none

    !f2py threadsafe
    !f2py intent(in)               :: todfile
    !f2py intent(in)               :: invnttfile
    !f2py intent(in)               :: convert
    !f2py intent(in)               :: npixels_per_sample
    !f2py intent(hide)             :: nsamples_tot = shape(tod,0)
    !f2py intent(hide)             :: ndetectors = shape(tod,1)
    !f2py intent(inout)            :: tod
    !f2py integer*8, intent(inout) :: pmatrix(npixels_per_sample*nsamples_tot*ndetectors)
    !f2py intent(out)              :: status

    character(len=*), intent(in)         :: todfile
    character(len=*), intent(in)         :: invnttfile
    character(len=*), intent(in)         :: convert
    integer, intent(in)                  :: npixels_per_sample
    integer, intent(in)                  :: nsamples_tot
    integer, intent(in)                  :: ndetectors
    real*8, intent(inout)                :: tod(nsamples_tot,ndetectors)
    type(PointingElement), intent(inout) :: pmatrix(npixels_per_sample,nsamples_tot,ndetectors)
    integer, intent(out)                 :: status

    type(FilterUncorrelated), allocatable :: filter(:)
    integer, allocatable                  :: nsamples(:)
    
    call read_filter(invnttfile, convert, ndetectors, filter, nsamples, status)
    if (status /= 0) return

    call read_tod(todfile, convert, nsamples, tod, pmatrix, status)

end subroutine madmap1_read_tod


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine madmap1_read_filter(filename, convert, ncorrelations, ndetectors, nslices, data, nsamples, status)

    use iso_fortran_env,  only : ERROR_UNIT
    use module_filtering, only : FilterUncorrelated, create_filter_uncorrelated
    use module_madcap,    only : read_filter
    implicit none

    !f2py threadsafe
    !f2py intent(in)   :: filename
    !f2py intent(in)   :: convert
    !f2py intent(in)   :: ncorrelations
    !f2py intent(in)   :: ndetectors
    !f2py intent(in)   :: nslices
    !f2py intent(out)  :: data
    !f2py intent(out)  :: nsamples
    !f2py intent(out)  :: status

    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: convert
    integer, intent(in)          :: ncorrelations
    integer, intent(in)          :: ndetectors
    integer, intent(in)          :: nslices
    real*8, intent(out)          :: data(ncorrelations+1, ndetectors, nslices)
    integer, intent(out)         :: nsamples(nslices)
    integer, intent(out)         :: status

    type(FilterUncorrelated), allocatable :: filter(:)
    integer, allocatable                  :: nsamples_(:)
    integer                               :: islice

    call read_filter(filename, convert, ndetectors, filter, nsamples_, status)
    if (status /= 0) return

    if (size(filter) /= nslices) then
        status = 1
        write (ERROR_UNIT,'(a)') 'READ_FILTER_MADMAP1: Incompatible number of slices.'
    else
        nsamples = nsamples_
    end if

    if (any(filter%ncorrelations /= ncorrelations)) then
        status = 1
        write (ERROR_UNIT,'(a)') 'READ_FILTER_MADMAP1: filters have an invalid correlation length.'
    end if

    if (any(filter%ndetectors /= ndetectors)) then
        status = 1
        write (ERROR_UNIT,'(a)') 'READ_FILTER_MADMAP1: filters have an invalid number of detectors.'
    end if

    if (status /= 0) return

    do islice = 1, nslices
        data(:,:,islice) = filter(islice)%data
    end do

end subroutine madmap1_read_filter
