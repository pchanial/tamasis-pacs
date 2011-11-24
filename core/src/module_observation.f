! Copyright 2010-2011 Pierre Chanial
! All rights reserved

module module_observation

    use iso_fortran_env,  only : ERROR_UNIT, OUTPUT_UNIT
    use module_math,      only : NaN, barycenter_lonlat
    use module_string,    only : strinteger, strlowcase, strreal, strsection, strternary
    use module_tamasis,   only : p, POLICY_KEEP, POLICY_MASK, POLICY_REMOVE
    implicit none
    private

    public :: MaskPolicy
    public :: Observation
    public :: ObservationPointing
    public :: POINTING_INSCAN, POINTING_TURNAROUND, POINTING_OTHER

    integer, parameter :: POINTING_INSCAN = 1
    integer, parameter :: POINTING_TURNAROUND = 2
    integer, parameter :: POINTING_OTHER = 3

    type MaskPolicy
        integer :: inscan = POLICY_KEEP
        integer :: turnaround = POLICY_KEEP
        integer :: other = POLICY_REMOVE
        integer :: invalid = POLICY_MASK
    end type MaskPolicy

    type ObservationPointing
        real(p)   :: time
        real(p)   :: ra
        real(p)   :: dec
        real(p)   :: pa
        real(p)   :: chop

        logical*1 :: inscan
        logical*1 :: turnaround
        logical*1 :: other      ! not inscan and not turnaround
        logical*1 :: invalid    ! wrong finetime, undef. band

        logical*1 :: removed
        logical*1 :: masked

        logical*1 :: notused1
        logical*1 :: notused2
    end type ObservationPointing

    type Observation
        character(len=512) :: id
        integer            :: nsamples
        integer            :: nvalids
        integer            :: factor
        real(p)            :: offset
        type(ObservationPointing), allocatable :: pointing(:)
    contains
        procedure :: init
        procedure :: get_position_time
        procedure :: get_position_index
        procedure :: compute_center
    end type Observation

    public :: strpolicy


contains


    subroutine init(this, time, ra, dec, pa, chop, masked, removed, factor, offset, status)
        ! delay is the instrument lag wrt the telescope

        class(Observation), intent(inout) :: this
        real(p), intent(in)               :: time(:), ra(:), dec(:), pa(:), chop(:)
        logical*1, intent(in)             :: masked(:), removed(:)
        integer, intent(in)               :: factor
        real(p), intent(in)               :: offset
        integer, intent(out)              :: status

        integer :: nsamples

        status = 1

        ! check conformity of time, ra, dec, pa, etc.
        nsamples = size(time)
        if (size(ra) /= nsamples) then
            write (ERROR_UNIT,'(a)') "INIT: Input argument 'ra' has an invalid number of samples."
            return
        end if
        if (size(dec) /= nsamples) then
            write (ERROR_UNIT,'(a)') "INIT: Input argument 'dec' has an invalid number of samples."
            return
        end if
        if (size(pa) /= nsamples) then
            write (ERROR_UNIT,'(a)') "INIT: Input argument 'pa' has an invalid number of samples."
            return
        end if
        if (size(chop) /= nsamples) then
            write (ERROR_UNIT,'(a)') "INIT: Input argument 'chop' has an invalid number of samples."
            return
        end if
        if (size(masked) /= nsamples) then
            write (ERROR_UNIT,'(a)') "INIT: Input argument 'masked' has an invalid number of samples."
            return
        end if
        if (size(removed) /= nsamples) then
            write (ERROR_UNIT,'(a)') "INIT: Input argument 'removed' has an invalid number of samples."
            return
        end if

        allocate (this%pointing(nsamples))

        this%nsamples = nsamples
        this%factor = factor
        this%offset = offset
        this%nvalids = count(.not. removed)
        this%pointing%time    = time
        this%pointing%ra      = ra
        this%pointing%dec     = dec
        this%pointing%pa      = pa
        this%pointing%chop    = chop
        this%pointing%masked  = masked
        this%pointing%removed = removed

        status = 0

    end subroutine init


    !-------------------------------------------------------------------------------------------------------------------------------


    ! linear interpolation based on time values, useful if time samples are not evenly spaced.
    ! for the first call, index should be set to zero
    subroutine get_position_time(this, time, ra, dec, pa, chop, index)

        class(Observation), intent(in) :: this
        real(p), intent(in)            :: time
        real(p), intent(out)           :: ra, dec, pa, chop
        integer, intent(inout)         :: index

        integer                        :: i
        real(p)                        :: frac

        if (time > 2._p * this%pointing(this%nsamples)%time - this%pointing(this%nsamples-1)%time .or.                             &
            time < 2._p * this%pointing(1)%time - this%pointing(2)%time) then
            ra   = NaN
            dec  = NaN
            pa   = NaN
            chop = NaN
            return
        end if

        if (index < 1 .or. index-1 > this%nsamples) then
            index = 2
        else if (time <= this%pointing(index-1)%time) then
            index = 2
        end if

        do i = index, this%nsamples
            if (time <= this%pointing(i)%time) go to 100
        end do
        i = i - 1

    100 frac = (time - this%pointing(i-1)%time) / (this%pointing(i)%time - this%pointing(i-1)%time)

        ra   = this%pointing(i-1)%ra   * (1 - frac) + this%pointing(i)%ra   * frac
        dec  = this%pointing(i-1)%dec  * (1 - frac) + this%pointing(i)%dec  * frac
        pa   = this%pointing(i-1)%pa   * (1 - frac) + this%pointing(i)%pa   * frac
        chop = this%pointing(i-1)%chop * (1 - frac) + this%pointing(i)%chop * frac

        index = i

    end subroutine get_position_time


   !-------------------------------------------------------------------------------------------------------------------------------


    ! The sampling factor is the compression factor times the fine sampling
    ! factor. itime is a fine sampling index
    ! if there is a gap in the timeline (which should not happen in one slice, and should be
    ! corrected beforehand), it will not be taken into account so that finer
    ! sampling will be interpolated using the start and end of the gap.
    subroutine get_position_index(this, itime, sampling_factor, offset, ra, dec, pa, chop)

        class(Observation), intent(in) :: this
        integer, intent(in)            :: itime, sampling_factor
        real(p), intent(in)            :: offset
        real(p), intent(out)           :: ra, dec, pa, chop

        integer                        :: i, itime_max
        real(p)                        :: frac

        itime_max = this%nsamples * sampling_factor
        if (itime < 1 .or. itime > itime_max) then
            write (ERROR_UNIT,'(a,i0,a,i0,a)') "GET_POSITION_INDEX: Invalid time index '", itime, "'. Valid range is [1:",         &
                  itime_max,"]."
            stop
        end if

        i    = min((itime - 1) / sampling_factor, this%nsamples-2)
        frac = real(itime - 1 - sampling_factor * i, p) / sampling_factor - offset
        i    = i + 1

        ra   = this%pointing(i)%ra   * (1 - frac) + this%pointing(i+1)%ra   * frac
        dec  = this%pointing(i)%dec  * (1 - frac) + this%pointing(i+1)%dec  * frac
        pa   = this%pointing(i)%pa   * (1 - frac) + this%pointing(i+1)%pa   * frac
        chop = this%pointing(i)%chop * (1 - frac) + this%pointing(i+1)%chop * frac

    end subroutine get_position_index


    !-------------------------------------------------------------------------------------------------------------------------------


    ! compute the mean value of R.A. and dec, taking into account RA's singularity at 0
    subroutine compute_center(this, ra0, dec0)

        class(Observation), intent(in) :: this
        real(p), intent(out)           :: ra0, dec0

        call barycenter_lonlat(pack(this%pointing%ra,  .not. this%pointing%removed .and. .not. this%pointing%masked),              &
                               pack(this%pointing%dec, .not. this%pointing%removed .and. .not. this%pointing%masked), ra0, dec0)
        
    end subroutine compute_center


    !-------------------------------------------------------------------------------------------------------------------------------


    pure function strpolicy_len(policy)
        
        integer, intent(in) :: policy
        integer             :: strpolicy_len

        select case (policy)
            case (POLICY_KEEP)
                strpolicy_len = 4
            case (POLICY_MASK)
                strpolicy_len = 6
            case (POLICY_REMOVE)
                strpolicy_len = 7
            case default
                strpolicy_len = 0
        end select

    end function strpolicy_len


    !-------------------------------------------------------------------------------------------------------------------------------


    function strpolicy(policy)

        integer, intent(in)                  :: policy
        character(len=strpolicy_len(policy)) :: strpolicy

        select case (policy)
            case (POLICY_KEEP)
                strpolicy = 'kept'
            case (POLICY_MASK)
                strpolicy = 'masked'
            case (POLICY_REMOVE)
                strpolicy = 'removed'
        end select

    end function strpolicy


end module module_observation
