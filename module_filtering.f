module module_filtering

    use iso_fortran_env,  only : ERROR_UNIT
    use module_precision, only : dp, p
    implicit none
    private

    public :: filterset
    public :: create_filter_uncorrelated

    include 'fftw3.f'

    type filterset
        integer                :: ncorrelations
        integer                :: ndetectors
        integer                :: nslices
        integer                :: bandwidth
        integer*8, allocatable :: first(:), last(:)
        real(p), allocatable   :: data(:,:,:) ! (filter values, detector, slice)
    end type filterset

contains

    subroutine create_filter_uncorrelated(filter, nsamples, nsamples_tot, tod, status)

        type(filterset), intent(in) :: filter
        integer, intent(in)         :: nsamples(:), nsamples_tot
!        real(dp), intent(out)       :: tod(sum(nsamples),filter%ndetectors)
        real(dp), intent(out)       :: tod(nsamples_tot,filter%ndetectors)
        integer, intent(out)        :: status

        integer :: ndetectors, nslices, islice, islicefilter, idetector, dest
        real(p), allocatable :: in(:), out(:)
        integer*8 :: plan

        status = 1

        ndetectors = filter%ndetectors
        nslices = size(nsamples)
        if (nslices /= filter%nslices .and. nslices /= 1) then
            write (ERROR_UNIT,'(a,2(i0,a))') "CREATE_FILTER_UNCORRELATED: The timeline's slice number '", nslices, "' is incompatib&
                  &le with the number of filter slices '", filter%nslices, "'."
            return
        end if

        dest = 1
        do islice = 1, nslices

            if (filter%nslices /= 1) then
                islicefilter = islice
            else
                islicefilter = 1
            end if

            if (nsamples(islice) < filter%bandwidth) then
                write (ERROR_UNIT,'(a,3(i0,a))') "CREATE_FILTER_UNCORRELATED: Slice '", islice, "' cannot be filtered because it ha&
                      &s a size '", nsamples(islice), "' smaller than the filter bandwidth '", filter%bandwidth, "'."
                return
            end if

            allocate (in(nsamples(islice)), out(nsamples(islice)))

            call dfftw_plan_r2r_1d(plan, nsamples(islice), in, out, FFTW_R2HC, FFTW_ESTIMATE)

            do idetector = 1, ndetectors

                in(1:filter%ncorrelations+1) = filter%data(:,idetector,islicefilter)
                in(filter%ncorrelations+2:nsamples(islice)-filter%ncorrelations) = 0
                in(nsamples(islice)-filter%ncorrelations+1:nsamples(islice)) = &
                     filter%data(filter%ncorrelations+1:2:-1,idetector,islicefilter)

                call dfftw_execute_r2r(plan, in, out)

                tod(dest:dest+nsamples(islice)-1,idetector) = out

            end do

            dest = dest + nsamples(islice)
            deallocate (in, out)
            
            call dfftw_destroy_plan(plan)

        end do

        status = 0

    end subroutine create_filter_uncorrelated


end module module_filtering
