!-------------------------------------------------------------------------------
!
! This module defines the pacsobservation derived type.
! It contains all the information of a PACS observation relevant to its data
! processing. In particular:
!    - file name of the observation
!    - channel ('b', 'g' or 'r')
!    - mode ('prime', 'parallel', 'transparent')
!    - section to be processed
!
! Author: Pierre Chanial
!
!-------------------------------------------------------------------------------

module module_pacsobservation

    use iso_fortran_env,    only : ERROR_UNIT, OUTPUT_UNIT
    use module_fitstools,   only : ft_open, ft_open_bintable, ft_read_column, ft_read_extension, ft_close
    use module_math,        only : median_nocopy, neq_real
    use module_observation, only : maskarray, observation, pointing, pacsobservationslice
    use module_precision,   only : dp, p
    use module_string,      only : strinteger, strlowcase, strreal, strsection, strternary
    implicit none
    private

    public :: maskarray
!!$    public :: observationslice
    public :: pacsobservationslice
    public :: pacsobservation

    type, extends(observation) :: pacsobservation

    contains

        private
        procedure, public :: init
        procedure, public :: init_sim
        procedure, public :: print

    end type pacsobservation


contains


    subroutine init(this, filename, maskarray_policy, status, verbose)

        class(pacsobservation), intent(inout) :: this
        character(len=*), intent(in)          :: filename(:)
        type(maskarray), intent(in)           :: maskarray_policy
        integer, intent(out)                  :: status
        logical, intent(in), optional         :: verbose
        integer*8                             :: first, last
        integer                               :: iobs
        logical                               :: verbose_

        ! parameter checking
        if (size(filename) == 0) then 
            status = 1
            write (ERROR_UNIT,'(a)') 'INIT_PACSOBSERVATION: Array has zero size.'
            return
        end if

        if (present(verbose)) then
            verbose_ = verbose
        else
            verbose_ = .false.
        end if

        ! set number of observations
        this%nslices = size(filename)

        ! set the mask policy '.true.' means 'taken into account'
        this%maskarray_policy = maskarray_policy        
        
        if (allocated(this%slice)) then
            do iobs = 1, size(this%slice)
                if (allocated(this%slice(iobs)%p)) deallocate (this%slice(iobs)%p)
            end do
            deallocate(this%slice)
        end if
        allocate(this%slice(this%nslices))

        do iobs = 1, this%nslices

            call this%slice(iobs)%set_filename(filename(iobs), first, last, status)
            if (status /= 0) return

            call this%slice(iobs)%set_valid_slice(first, last, this%maskarray_policy, status)
            if (status /= 0) return

            call this%slice(iobs)%set_channel(status)
            if (status /= 0) return
                
            call this%slice(iobs)%set_compression_mode(status)
            if (status /= 0) return
        
            call this%slice(iobs)%set_observing_mode(status)
            if (status /= 0) return

            call this%slice(iobs)%set_pointing(status, verbose_)
            if (status /= 0) return
                
            call this%slice(iobs)%validate_pointing(this%maskarray_policy, status, verbose_)
            if (status /= 0) return

        end do

        ! make sure the channel is the same for all observations
        this%channel = this%slice(1)%channel
        if (any(this%slice%channel /= this%channel)) then
            status = 1
            write (ERROR_UNIT,'(a)') 'ERROR: Observations do not have the same channel: ', this%slice%channel, '.'
            return
        end if

        ! make sure the observing mode is the same for all observations
        ! that could be relaxed, but we would need to properly handle the bad detector mask for transparent observations
        this%observing_mode = this%slice(1)%observing_mode
        if (any(this%slice%observing_mode /= this%observing_mode)) then
            status = 1
            write (ERROR_UNIT,'(a)') 'ERROR: Observations do not have the same observing mode.'
            return
        end if

        this%nsamples_tot = sum(this%slice%nsamples)

    end subroutine init


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine init_sim(this, time, ra, dec, pa, chop, status, verbose)
        class(pacsobservation), intent(inout) :: this
        real*8, intent(in)                    :: time(:),ra(:),dec(:),pa(:),chop(:)
        integer, intent(out)                  :: status
        logical, intent(in), optional         :: verbose
        integer*8                             :: nsamples
        logical                               :: verbose_

        if (present(verbose)) then
            verbose_ = verbose
        else
            verbose_ = .false.
        end if

        status = 1

        ! check conformity of time, ra, dec, pa
        nsamples = size(time)

        if (size(ra) /= nsamples) then
            write (ERROR_UNIT,'(a)') "Input R.A. has an invalid number of samples."
            return
        endif
        if (size(dec) /= nsamples) then
            write (ERROR_UNIT,'(a)') "Input declination has an invalid number of samples."
            return
        endif
        if (size(pa) /= nsamples) then
            write (ERROR_UNIT,'(a)') "Input P.A. has an invalid number of samples."
            return
        endif
        if (size(chop) /= nsamples) then
            write (ERROR_UNIT,'(a)') "Input chop angle has an invalid number of samples."
            return
        endif

        this%nslices      = 1
        this%nsamples_tot = nsamples
        this%channel      = ' '
        this%observing_mode = 'Unknown'
        this%maskarray_policy = maskarray(wrong_time=.true., off_target=.true.)

        allocate (this%slice(1))
        this%slice(1)%first    = 1
        this%slice(1)%last     = nsamples
        this%slice(1)%nsamples = nsamples
        this%slice(1)%filename = 'simulation'
        this%slice(1)%channel  = ' '
        this%slice(1)%observing_mode = 'Unknown'
        this%slice(1)%compression_factor = nint((time(2)-time(1)) / 0.024996_dp)
        this%slice(1)%delta = time(2)-time(1)
        allocate (this%slice(1)%p(nsamples))

        this%slice(1)%p%time       = time
        this%slice(1)%p%ra         = ra
        this%slice(1)%p%dec        = dec
        this%slice(1)%p%pa         = pa
        this%slice(1)%p%chop       = chop
        this%slice(1)%p%invalid    = .false.
        this%slice(1)%p%wrong_time = .false.
        this%slice(1)%p%off_target = .false.

        call this%slice(1)%validate_pointing(this%maskarray_policy, status, verbose_)

    end subroutine init_sim


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine print(this)

        class(pacsobservation), intent(in) :: this

        integer :: iobs, nrejected, nsamples

        write (*,*)

        do iobs = 1, this%nslices

            nsamples = this%slice(iobs)%nsamples

            ! observation number & file name
            write (OUTPUT_UNIT,'(a)') 'Info: Observation' // strternary(this%nslices>1, ' ' // strinteger(iobs), '') // ': ' //    &
                trim(this%slice(iobs)%filename)
            write (OUTPUT_UNIT,'(a)') '      Section: [' // strsection(this%slice(iobs)%first,this%slice(iobs)%last) // ']'
            
            ! channel
            write (OUTPUT_UNIT,'(a,$)') "      Channel: "
            select case (this%slice(iobs)%channel)
                case ('b')
                    write (OUTPUT_UNIT,'(a)') 'Blue'
                case ('g')
                    write (OUTPUT_UNIT,'(a)') 'Green'
                case ('r')
                    write (OUTPUT_UNIT,'(a)') 'Red'
                case default
                    write (OUTPUT_UNIT,'(a)') 'Unknown'
            end select

            ! observing mode
            write (OUTPUT_UNIT,'(a)') '      Observing mode: ' // trim(this%observing_mode)

            ! compression factor
            write (OUTPUT_UNIT,'(a,i0)') '      Compression factor: ', this%slice(iobs)%compression_factor

            ! print maskarray information
            if (this%maskarray_policy%off_target) then
                nrejected = count(this%slice(iobs)%p%off_target)
                if (nrejected > 0) then
                    write (OUTPUT_UNIT,'(a,2(i0,a))') "      Mask 'off target': Rejecting ", nrejected, ' / ', nsamples, '.'
                end if
            end if

            if (this%maskarray_policy%wrong_time) then
                nrejected = count(this%slice(iobs)%p%wrong_time)
                if (nrejected > 0) then
                    write (OUTPUT_UNIT,'(a,2(i0,a))') "      Mask 'wrong time': Rejecting ", nrejected, ' / ', nsamples, '.'
                end if
            end if

            write (*,*)

        end do

    end subroutine print

 end module module_pacsobservation
