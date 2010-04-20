module module_pacspointing

    use iso_fortran_env,        only : OUTPUT_UNIT, ERROR_UNIT
    use module_fitstools,       only : ft_read_column, ft_read_extension, ft_open_bintable, ft_close
    use module_math,            only : NaN, median_nocopy, neq_real
    use module_pacsobservation, only : pacsobservation, pacsobsinfo
    use module_precision,       only : p, dp
    use module_string,          only : strinteger, strreal, strternary
    implicit none
    private

    public :: pacspointing

    type pacspointing

        integer                :: nslices, nsamples_tot
        integer*8, allocatable :: nsamples(:), first(:), last(:)
        integer, allocatable   :: compression_factor(:)
        real*8, allocatable    :: time(:), ra(:), dec(:), pa(:), chop(:)
        logical*1, allocatable :: invalid_time(:)
        real(dp), allocatable  :: delta(:)

    contains

        private
        procedure, public :: init
        procedure         :: init_oldstyle
        procedure, public :: init_sim
        procedure         :: init2
        procedure, public :: get_position_time
        procedure, public :: get_position_index
        procedure, public :: compute_center
        procedure, public :: print
        procedure, public :: destructor ! XXX should be final, but not handled by gfortran 4.5

    end type pacspointing


contains


    subroutine init(this, obs, status, verbose)

        class(pacspointing), intent(inout)    :: this
        class(pacsobservation), intent(inout) :: obs
        integer, intent(out)                  :: status
        logical, intent(in), optional         :: verbose

        real(kind=dp), allocatable :: buffer(:)
        integer*8, allocatable     :: timeus(:)
        integer                    :: nsamples
        integer :: first, last
        integer :: length, unit, status_close, iobs

        this%nslices      = size(obs%info)
        this%nsamples_tot = sum(obs%info%nsamples)
        allocate(this%nsamples(this%nslices))
        allocate(this%first   (this%nslices))
        allocate(this%last    (this%nslices))
        allocate(this%delta   (this%nslices))
        allocate(this%compression_factor(this%nslices))
        allocate(this%time    (this%nsamples_tot))
        allocate(this%ra      (this%nsamples_tot))
        allocate(this%dec     (this%nsamples_tot))
        allocate(this%pa      (this%nsamples_tot))
        allocate(this%chop    (this%nsamples_tot))
        allocate(this%invalid_time(this%nsamples_tot))
        this%nsamples = obs%info%nsamples

        do iobs = 1, this%nslices

            ! compute range of slices. They are contiguous in pacspointing
            this%first(iobs) = sum(obs%info(1:iobs-1)%nsamples) + 1
            this%last (iobs) = this%first(iobs) + this%nsamples(iobs) - 1

            length = len_trim(obs%info(iobs)%filename)
            if (obs%info(iobs)%filename(length-4:length) /= '.fits') then
                 call this%init_oldstyle(obs%info(iobs), iobs, status)
                 if (status /= 0) return
                 cycle
            end if

            allocate(timeus(this%nsamples(iobs)))
            allocate(buffer(this%nsamples(iobs)))

            call ft_open_bintable(trim(obs%info(iobs)%filename) // '[Status]', unit, nsamples, status)
            if (status /= 0) return
            if (nsamples < this%nsamples(iobs)) then
                status = 1
                write (ERROR_UNIT,'(a)') 'ERROR: There is not pointing information for every frame.'
                go to 999
            end if

            first = obs%info(iobs)%first
            last  = obs%info(iobs)%last

            call ft_read_column(unit, 'FINETIME', first, last, timeus, status)
            if (status /= 0) go to 999
            buffer = timeus * 1.d-6
            this%time(this%first(iobs):this%last(iobs)) = buffer

            call ft_read_column(unit, 'RaArray', first, last, buffer, status)
            if (status /= 0) go to 999
            this%ra(this%first(iobs):this%last(iobs)) = buffer

            call ft_read_column(unit, 'DecArray', first, last, buffer, status)
            if (status /= 0) go to 999
            this%dec(this%first(iobs):this%last(iobs)) = buffer

            call ft_read_column(unit, 'PaArray', first, last, buffer, status)
            if (status /= 0) go to 999
            this%pa(this%first(iobs):this%last(iobs)) = buffer

            call ft_read_column(unit, 'CHOPFPUANGLE', first, last,buffer,status)
            if (status /= 0) go to 999
            this%chop(this%first(iobs):this%last(iobs)) = buffer

            deallocate(timeus)
            deallocate(buffer)

            call ft_close(unit, status)
            if (status /= 0) return

        end do

        this%invalid_time = .false.

        call this%init2(status)
        if (status /= 0) return

        ! some post processing
        do iobs = 1, this%nslices

            ! check that the compression factor taken from the observation headers and from the pointing is the same
            if (this%compression_factor(iobs) /= obs%info(iobs)%compression_factor) then
                status = 1
                write (ERROR_UNIT, '(a)') 'Error: ' // strternary(this%nslices>1,'In observation '//strinteger(iobs)//', i','I') //&
                      "ncompatible compression factor from header '" // strinteger(obs%info(iobs)%compression_factor) // "' and poi&
                      &nting '" // strinteger(this%compression_factor(iobs)) // "'."
            end if
            
            ! propagate the mask up the pacs observation
            call obs%info(iobs)%set_maskarray('invalid time', obs%maskarray_policy, this%invalid_time(this%first(iobs):            &
                 this%last(iobs)), verbose)

        end do
        return

    999 call ft_close(unit, status_close)
        if (allocated(buffer)) deallocate(buffer)
        if (allocated(timeus)) deallocate(timeus)

    end subroutine init


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine init_oldstyle(this, obs, iobs, status)

        class(pacspointing), intent(inout) :: this
        class(pacsobsinfo), intent(in)     :: obs
        integer, intent(in)                :: iobs
        integer, intent(out)               :: status
        real*8, allocatable                :: buffer(:)
        integer*8, allocatable             :: timeus(:)

        call ft_read_extension(trim(obs%filename)//'_Time.fits', timeus, status)
        if (status /= 0) return
        allocate(buffer(size(timeus)))
        buffer = timeus * 1.0d-6
        this%time(this%first(iobs):this%last(iobs)) = buffer(obs%first:obs%last)

        call ft_read_extension(trim(obs%filename) // '_RaArray.fits', buffer, status)
        if (status /= 0) return
        this%ra(this%first(iobs):this%last(iobs)) = buffer(obs%first:obs%last)

        call ft_read_extension(trim(obs%filename) // '_DecArray.fits', buffer, status)
        if (status /= 0) return
        this%dec(this%first(iobs):this%last(iobs)) = buffer(obs%first:obs%last)

        call ft_read_extension(trim(obs%filename) // '_PaArray.fits', buffer, status)
        if (status /= 0) return
        this%pa(this%first(iobs):this%last(iobs)) = buffer(obs%first:obs%last)

        call ft_read_extension(trim(obs%filename) // '_ChopFpuAngle.fits', buffer, status)
        if (status /= 0) return
        this%chop(this%first(iobs):this%last(iobs)) = buffer(obs%first:obs%last)

    end subroutine init_oldstyle


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine init_sim(this, time, ra, dec, pa, chop, status)
        class(pacspointing), intent(inout) :: this
        real*8, intent(in)                 :: time(:),ra(:),dec(:),pa(:),chop(:)
        integer, intent(out)               :: status
        integer*8                          :: nsamples

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

        allocate(this%nsamples(1))
        allocate(this%first   (1))
        allocate(this%last    (1))
        allocate(this%delta   (1))
        allocate(this%compression_factor(1))
        allocate(this%time(nsamples))
        allocate(this%ra  (nsamples))
        allocate(this%dec (nsamples))
        allocate(this%pa  (nsamples))
        allocate(this%chop(nsamples))
        allocate(this%invalid_time(nsamples))

        this%nslices      = 1
        this%nsamples_tot = nsamples
        this%nsamples(1)  = nsamples
        this%first(1)     = 1
        this%last (1)     = nsamples
        this%compression_factor = 1
        this%time         = time
        this%ra           = ra
        this%dec          = dec
        this%pa           = pa
        this%chop         = chop
        this%invalid_time = .false.

        call this%init2(status)

    end subroutine init_sim


    !-------------------------------------------------------------------------------------------------------------------------------


    ! in case of a negative pointing time jump (should affect only one frame, otherwise an error is triggered)
    ! the position is taken as the previous one and the frame is masked
    subroutine init2(this, status)
        class(pacspointing), intent(inout) :: this
        integer, intent(out) :: status
        integer   :: islice, isample, njumps
        integer*8 :: nsamples, first, last
        real(dp), allocatable :: delta(:)
        real(dp)              :: delta_max

        status = 0

        do islice = 1, this%nslices

            nsamples = this%nsamples(islice)
            first    = this%first(islice)
            last     = this%last(islice)
           
            if (nsamples <= 1) then
                status = 1
                write (ERROR_UNIT,'(a)') 'Input time has less than 2 samples.'
                return
            end if

            ! check that the input time is monotonous (and increasing)
            allocate(delta(this%nsamples(islice)-1))
            delta = this%time(first+1:last) - this%time(first:last-1)
            this%delta(islice) = median_nocopy(delta)
            njumps = 0
            do isample = 2, nsamples
                if (delta(isample-1) <= 0) then
                    if (isample < nsamples) then
                        if (delta(isample) <= 0) then
                            status = 1
                            deallocate(delta)
                            write (ERROR_UNIT,'(a)') "ERROR: The pointing time is not strictly increasing."
                        end if
                    end if
                    njumps = njumps + 1
                    this%invalid_time(first+isample-1) = .true.
                    this%time(first+isample-1) = this%time(first+isample-2) + this%delta(islice)
                    this%ra  (first+isample-1) = this%ra  (first+isample-2)
                    this%dec (first+isample-1) = this%dec (first+isample-2)
                    this%pa  (first+isample-1) = this%pa  (first+isample-2)
                    this%chop(first+isample-1) = this%chop(first+isample-2)
                end if
            end do

            if (njumps > 0) then
                write (OUTPUT_UNIT,'(a,i0,a)') 'Warning: The pointing fine time has ', njumps, ' negative jump(s). The affected fra&
                      &mes have been masked.'
                delta = this%time(first+1:last) - this%time(first:last-1)
            end if

            ! check if there are gaps
            delta_max = maxval(abs(delta))
            if (any(neq_real(delta, this%delta(islice), 3))) then
                write (*,'(a)') 'Warning: ' // strternary(this%nslices>1, ' In observation '//strinteger(islice)//', t','T') //    &
                      'he pointing time is not evenly spaced.'
                if (delta_max > 1.5_p * this%delta(islice)) then
                    write (*,'(a)') '         Largest gap is ' // strreal(delta_max*1000._p,1) // 'ms.'
                end if
            end if

            ! check the compression factor from the data themselves
            this%compression_factor(islice) = nint(this%delta(islice) / 0.024996_dp)
            if (neq_real(this%compression_factor(islice) * 0.024996_dp, this%delta(islice), 2)) then
                status = 1
                write (*,'(a)') 'Error: The sampling time is not an integer number of PACS sampling time (40Hz).'
                return
            end if

            deallocate(delta)
          
        end do

    end subroutine init2


    !-------------------------------------------------------------------------------------------------------------------------------


    ! linear interpolation if time samples are not evenly spaced.
    ! for the first call, index should be set to zero
    subroutine get_position_time(this, islice, time, ra, dec, pa, chop, index)
        class(pacspointing), intent(in) :: this
        integer, intent(in)             :: islice
        real*8, intent(in)              :: time
        real*8, intent(out)             :: ra, dec, pa, chop
        integer, intent(inout)          :: index
        integer                         :: i, first, last
        real*8                          :: frac

        if (islice < 1 .or. islice > this%nslices) then
            write (ERROR_UNIT,'(a,i0,a)') "GET_POSITION: Invalid slice number '&
                &", islice, "'."
            stop
        end if

        first = this%first(islice)
        last  = this%last (islice)

        if (time > 2.d0 * this%time(last) - this%time(last-1) .or. time < 2.d0 * this%time(first) - this%time(first+1)) then
            ra   = NaN
            dec  = NaN
            pa   = NaN
            chop = NaN
            return
        end if

        if (index == 0 .or. index-1 > size(this%time)) then
            index = first + 1
        else if (time <= this%time(index-1)) then
            index = first + 1
        end if

        do i = index, last
            if (time <= this%time(i)) go to 100
        end do
        i = i - 1

    100 frac = (time - this%time(i-1)) / (this%time(i) - this%time(i-1))

        ra   = this%ra  (i-1) * (1 - frac) + this%ra  (i) * frac
        dec  = this%dec (i-1) * (1 - frac) + this%dec (i) * frac
        pa   = this%pa  (i-1) * (1 - frac) + this%pa  (i) * frac
        chop = this%chop(i-1) * (1 - frac) + this%chop(i) * frac

        index = i

    end subroutine get_position_time


    !-------------------------------------------------------------------------------------------------------------------------------


    ! The sampling factor is the compression factor times the fine sampling
    ! factor. itime is a fine sampling index
    ! if there is a gap in the timeline (which should not happen, and should be
    ! corrected beforehand), it will not be taken into account so that finer
    ! sampling will be interpolated using the start and end of the gap.
    subroutine get_position_index(this, islice, itime, sampling_factor, ra, dec, pa, chop)
        class(pacspointing), intent(in) :: this
        integer, intent(in)             :: islice, itime, sampling_factor
        real*8, intent(out)             :: ra, dec, pa, chop
        integer                         :: i, itime_max
        real*8                          :: frac

        if (islice < 1 .or. islice > this%nslices) then
            write (ERROR_UNIT,'(a,i0,a)') "GET_POSITION: Invalid slice number '", islice, "'."
            stop
        end if

        itime_max = this%nsamples(islice) * sampling_factor
        if (itime < 1 .or. itime > itime_max) then
            write (ERROR_UNIT,'(a,i0,a,i0,a)') "GET_POSITION: Invalid time index '", itime, "'. Valid range is [1:", itime_max,"]."
            stop
        end if

        i    = min((itime - 1) / sampling_factor, this%nsamples(islice)-2)
        frac = real(itime - 1 - sampling_factor * i, p) / sampling_factor
        i    = i + this%first(islice)

        ra   = this%ra  (i) * (1 - frac) + this%ra  (i+1) * frac
        dec  = this%dec (i) * (1 - frac) + this%dec (i+1) * frac
        pa   = this%pa  (i) * (1 - frac) + this%pa  (i+1) * frac
        chop = this%chop(i) * (1 - frac) + this%chop(i+1) * frac

    end subroutine get_position_index


    !-------------------------------------------------------------------------------------------------------------------------------


    ! compute the mean value of R.A. and dec, taking into account RA's singularity at 0
    subroutine compute_center(this, ra0, dec0)
        class(pacspointing), intent(in) :: this
        real*8, intent(out)             :: ra0, dec0
        integer                         :: isample, n180
        real*8                          :: ra, dec
        logical                         :: zero_minus, zero_plus

        ra0  = 0.0d0
        dec0 = 0.0d0
        zero_minus = .false.
        zero_plus  = .false.
        n180 = 0

        !$omp parallel do default(shared) reduction(+:n180,ra0,dec0)           &
        !$omp reduction(.or.:zero_minus,zero_plus) private(isample, ra, dec)
        do isample=1, this%nsamples_tot
           ra  = this%ra(isample)
           dec = this%dec(isample)
           zero_minus = zero_minus .or. ra > 270.d0
           zero_plus  = zero_plus  .or. ra <= 90.d0
           if (ra >= 180.d0) n180 = n180 + 1
           ra0  = ra0  + ra
           dec0 = dec0 + dec
        end do
        !$omp end parallel do

        if (zero_minus .and. zero_plus) ra0 = ra0 - 360.d0 * n180

        ra0  = ra0  / this%nsamples_tot
        dec0 = dec0 / this%nsamples_tot

    end subroutine compute_center


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine print(this)
        class(pacspointing), intent(in) :: this

        write (*,*) 'Time: ', minval(this%time), '...', maxval(this%time)
        write (*,*) 'RA  : ', minval(this%ra)  , '...', maxval(this%ra)
        write (*,*) 'Dec : ', minval(this%dec) , '...', maxval(this%dec)
        write (*,*) 'PA  : ', minval(this%pa)  , '...', maxval(this%pa)
        write (*,*) 'Chop: ', minval(this%chop), '...', maxval(this%chop)

    end subroutine print


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine destructor(this)

        class(pacspointing), intent(inout) :: this

        if (allocated(this%time))     deallocate (this%time)
        if (allocated(this%ra))       deallocate (this%ra)
        if (allocated(this%dec))      deallocate (this%dec)
        if (allocated(this%pa))       deallocate (this%pa)
        if (allocated(this%chop))     deallocate (this%chop)
        if (allocated(this%invalid_time)) deallocate (this%invalid_time)
        if (allocated(this%nsamples)) deallocate (this%nsamples)
        if (allocated(this%first))    deallocate (this%first)
        if (allocated(this%last ))    deallocate (this%last)
        if (allocated(this%delta))    deallocate (this%delta)
        if (allocated(this%compression_factor)) deallocate (this%compression_factor)

    end subroutine


end module module_pacspointing
