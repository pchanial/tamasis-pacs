! Copyright 2010-2011 Pierre Chanial
! All rights reserved
!
module module_filtering

    use iso_fortran_env, only : ERROR_UNIT
    use module_tamasis,  only : p
    implicit none
    private

    public :: FilterUncorrelated
    public :: create_filter_uncorrelated
    public :: fft_tod
    public :: convolution_trexp_direct
    public :: convolution_trexp_transpose

    include 'fftw3.f'

    type FilterUncorrelated
        integer                :: ncorrelations
        integer                :: bandwidth
        integer                :: ndetectors
        real(p), allocatable   :: data(:,:) ! (filter values, detector)
    end type FilterUncorrelated


contains


    subroutine create_filter_uncorrelated(filter, nsamples, ndetectors, tod, status)

        type(FilterUncorrelated), intent(in) :: filter(:)
        integer, intent(in)                  :: nsamples(:)
        integer, intent(in)                  :: ndetectors
        real(p), intent(out)                 :: tod(sum(nsamples),ndetectors)
        integer, intent(out)                 :: status

        real(p), allocatable :: in(:), out(:)
        integer              :: dest, idetector, ifilter, islice, nslices, isample
        integer*8            :: plan

        status = 1

        nslices = size(nsamples)
        if (nslices /= size(filter) .and. nslices /= 1) then
            write (ERROR_UNIT,'(a,2(i0,a))') "CREATE_FILTER_UNCORRELATED: The timeline's slice number '", nslices, "' is incompatib&
                  &le with the number of filter slices '", size(filter), "'."
            return
        end if

        if (any(filter%ndetectors /= filter(1)%ndetectors)) then
            write (ERROR_UNIT,'(a)') 'CREATE_FILTER_UNCORRELATED: Slices do not have the same number of detectors.'
            return
        end if

        if (filter(1)%ndetectors /= ndetectors) then
            write (ERROR_UNIT,'(a,2(i0,a))') "CREATE_FILTER_UNCORRELATED: The specified number of detectors '", ndetectors, "'does &
                  &not match that in the filter '", filter(1)%ndetectors, "'."
            return
        end if

        dest = 1
        do islice = 1, nslices

            if (size(filter) /= 1) then
                ifilter = islice
            else
                ifilter = 1
            end if

            if (nsamples(islice) < filter(ifilter)%bandwidth) then
                write (ERROR_UNIT,'(a,3(i0,a))') "CREATE_FILTER_UNCORRELATED: Slice '", islice, "' cannot be filtered because it ha&
                      &s a size '", nsamples(islice), "' smaller than the filter bandwidth '", filter(ifilter)%bandwidth, "'."
                return
            end if

            allocate (in(nsamples(islice)), out(nsamples(islice)))

            call dfftw_plan_r2r_1d(plan, nsamples(islice), in, out, FFTW_R2HC, FFTW_ESTIMATE)

            do idetector = 1, ndetectors

                in(1:filter(ifilter)%ncorrelations+1) = filter(ifilter)%data(:,idetector)
                in(filter(ifilter)%ncorrelations+2:nsamples(islice)-filter(ifilter)%ncorrelations) = 0
                in(nsamples(islice)-filter(ifilter)%ncorrelations+1:nsamples(islice)) = &
                     filter(ifilter)%data(filter(ifilter)%ncorrelations+1:2:-1,idetector)

                call dfftw_execute_r2r(plan, in, out)

                out(1) = out(1) / nsamples(islice)
                do isample = 2, nsamples(islice) / 2 + 1
                    out(isample) = out(isample) / nsamples(islice)
                    out(nsamples(islice)-isample+2) = out(isample)
                end do
                tod(dest:dest+nsamples(islice)-1,idetector) = out

            end do

            dest = dest + nsamples(islice)
            deallocate (in, out)
            
            call dfftw_destroy_plan(plan)

        end do

        status = 0

    end subroutine create_filter_uncorrelated


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine fft_tod(plan, data)
        integer*8, intent(in)  :: plan
        real(p), intent(inout) :: data(:,:)

        integer :: idetector, ndetectors
        real(p) :: out(size(data,1))

        ndetectors = size(data, 2)

        !$omp parallel do default(shared) private(idetector, out)
        do idetector = 1, ndetectors
        
            call dfftw_execute_r2r(plan, data(:,idetector), out)
            data(:,idetector) = out
        
        end do
        !$omp end parallel do
        
    end subroutine fft_tod


    !-------------------------------------------------------------------------------------------------------------------------------
    
    
    subroutine convolution_trexp_direct(data, tau)

        real(p), intent(inout) :: data(:,:)
        real(p), intent(in)    :: tau(size(data,2))

        integer :: i, j
        real(p) :: v, w

        !$omp parallel do private(i,j,v,w)
        do j = 1, size(data,2)
            if (tau(j) <= 0 .or. tau(j) /= tau(j)) then
                data(:,j) = 0
                cycle
            end if
            w = exp(-1/tau(j))
            v = 1._p - w
            do i = 2, size(data,1)
                data(i,j) = w * data(i-1,j) + v * data(i,j)
            end do
        end do
        !$omp end parallel do
        
    end subroutine convolution_trexp_direct


    !-------------------------------------------------------------------------------------------------------------------------------
    
    
    subroutine convolution_trexp_transpose(data, tau)

        real(p), intent(inout) :: data(:,:)
        real(p), intent(in)    :: tau(size(data,2))

        integer :: m, n, i, j
        real(p) :: v, w

        m = size(data,1)
        n = size(data,2)
        !$omp parallel do private(i,j,v,w)
        do j = 1, n
            if (tau(j) <= 0 .or. tau(j) /= tau(j)) then
                data(:,j) = 0
                cycle
            end if
            w = exp(-1/tau(j))
            v = 1._p - w
            data(m,j) = v * data(m,j)
            do i = m-1, 1, -1
                data(i,j) = w * data(i+1,j) + v * data(i,j)
            end do
            data(1,j) = data(1,j) / v
        end do
        !$omp end parallel do
        
    end subroutine convolution_trexp_transpose


end module module_filtering


