module module_pacspointing
    use ISO_FORTRAN_ENV,        only : OUTPUT_UNIT, ERROR_UNIT
    use module_fitstools,       only : ft_read_column, ft_readextension, ft_open_bintable, ft_close
    use module_math,            only : mean, NaN
    use module_pacsobservation, only : pacsobservation, pacsobsinfo
    use precision,              only : p, dp
    use string,                 only : strinteger, strternary
    implicit none
    private

    public :: pacspointing

    type pacspointing

        integer                :: nslices, nsamples_tot
        integer*8, allocatable :: nsamples(:), first(:), last(:)
        real*8, allocatable    :: time(:), ra(:), dec(:), pa(:), chop(:)
        real(dp), allocatable  :: delta(:)

    contains

        private
        procedure, public :: init
        procedure         :: init_oldstyle
        procedure, public :: init_sim
        procedure         :: init2
        procedure, public :: get_position
        procedure, public :: compute_center
        procedure, public :: print
        procedure, public :: destructor ! XXX should be final, but not handled by gfortran 4.5

    end type pacspointing


contains


    subroutine init(this, obs, status)

        class(pacspointing), intent(inout) :: this
        class(pacsobservation), intent(in) :: obs
        integer, intent(out)               :: status

        real(kind=dp), allocatable         :: buffer(:)
        integer*8, allocatable             :: timeus(:)
        integer                            :: nsamples
        integer :: first, last
        integer :: length, unit, status_close, iobs

        this%nslices      = size(obs%info)
        this%nsamples_tot = sum(obs%info%nsamples)
        allocate(this%nsamples(this%nslices))
        allocate(this%first   (this%nslices))
        allocate(this%last    (this%nslices))
        allocate(this%delta   (this%nslices))
        allocate(this%time    (this%nsamples_tot))
        allocate(this%ra      (this%nsamples_tot))
        allocate(this%dec     (this%nsamples_tot))
        allocate(this%pa      (this%nsamples_tot))
        allocate(this%chop    (this%nsamples_tot))
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

            call ft_open_bintable(trim(obs%info(iobs)%filename) // '[Status]', unit,&
                                  nsamples, status)
            if (status /= 0) return
            if (nsamples < this%nsamples(iobs)) then
                status = 1
                write (ERROR_UNIT,'(a)') 'ERROR: There is not pointing informat&
                                         &ion for every frame.'
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

        call this%init2(status)
        return

    999 call ft_close(unit, status_close)
        if (allocated(buffer)) deallocate(buffer)
        if (allocated(timeus)) deallocate(timeus)

    end subroutine init


    !---------------------------------------------------------------------------


    subroutine init_oldstyle(this, obs, iobs, status)
        class(pacspointing), intent(inout) :: this
        type(pacsobsinfo), intent(in)      :: obs
        integer, intent(in)                :: iobs
        integer, intent(out)               :: status
        real*8, allocatable                :: buffer(:)
        integer*8, allocatable             :: timeus(:)

        call ft_readextension(trim(obs%filename)//'_Time.fits', timeus, status)
        if (status /= 0) return
        allocate(buffer(size(timeus)))
        buffer = timeus * 1.0d-6
        this%time(this%first(iobs):this%last(iobs)) = buffer(obs%first:obs%last)

        call ft_readextension(trim(obs%filename) // '_RaArray.fits', buffer, status)
        if (status /= 0) return
        this%ra(this%first(iobs):this%last(iobs)) = buffer(obs%first:obs%last)

        call ft_readextension(trim(obs%filename) // '_DecArray.fits', buffer, status)
        if (status /= 0) return
        this%dec(this%first(iobs):this%last(iobs)) = buffer(obs%first:obs%last)

        call ft_readextension(trim(obs%filename) // '_PaArray.fits', buffer, status)
        if (status /= 0) return
        this%pa(this%first(iobs):this%last(iobs)) = buffer(obs%first:obs%last)

        call ft_readextension(trim(obs%filename) // '_ChopFpuAngle.fits', buffer, status)
        if (status /= 0) return
        this%chop(this%first(iobs):this%last(iobs)) = buffer(obs%first:obs%last)

    end subroutine init_oldstyle


    !---------------------------------------------------------------------------


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
        allocate(this%time(nsamples))
        allocate(this%ra  (nsamples))
        allocate(this%dec (nsamples))
        allocate(this%pa  (nsamples))
        allocate(this%chop(nsamples))

        this%nslices      = 1
        this%nsamples_tot = nsamples
        this%nsamples(1)  = nsamples
        this%first(1)     = 1
        this%last (1)     = nsamples
        this%time         = time
        this%ra           = ra
        this%dec          = dec
        this%pa           = pa
        this%chop         = chop

        call this%init2(status)

    end subroutine init_sim


    !---------------------------------------------------------------------------


    subroutine init2(this, status)
        class(pacspointing), intent(inout) :: this
        integer, intent(out) :: status
        integer   :: islice
        integer*8 :: isample, nsamples, first, last
        real(dp), allocatable :: delta(:)
        real(dp), parameter :: tol = 0.01d0 ! fractional threshold of sampling above which sampling is considered uneven

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
            if (any(delta <= 0)) then
                status = 1
                deallocate(delta)
                write (ERROR_UNIT,'(a)') "ERROR: The pointing time is not stric&
                    &tly increasing"
                return
            end if

            ! check there is no time drifts
            this%delta(islice) = mean(delta)
            if (any(abs([(this%delta(islice)*(isample-1),isample=1,nsamples)] -&
                (this%time(first:last)-this%time(first))) > tol *              &
                this%delta(islice))) then
                write (*,'(a)') 'Warning: the pointing time is not evenly space&
                    &d or is drifting' // strternary(this%nslices>1, ' in obser&
                    &vation '//strinteger(islice),'') // '.'
            end if

            deallocate(delta)
          
        end do

    end subroutine init2


    !---------------------------------------------------------------------------


    ! time samples are not evenly spaced, slow linear interpolation
    ! for the first call, index should be set to zero
    subroutine get_position(this, islice, time, ra, dec, pa, chop, index)
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

        if (time > this%time(last)  + this%delta(islice) .or.                  &
            time < this%time(first) - this%delta(islice)) then
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

    end subroutine get_position


    !---------------------------------------------------------------------------


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


    !---------------------------------------------------------------------------


    subroutine print(this)
        class(pacspointing), intent(in) :: this

        write (*,*) 'Time: ', minval(this%time), '...', maxval(this%time)
        write (*,*) 'RA  : ', minval(this%ra)  , '...', maxval(this%ra)
        write (*,*) 'Dec : ', minval(this%dec) , '...', maxval(this%dec)
        write (*,*) 'PA  : ', minval(this%pa)  , '...', maxval(this%pa)
        write (*,*) 'Chop: ', minval(this%chop), '...', maxval(this%chop)

    end subroutine print


    !---------------------------------------------------------------------------

    subroutine destructor(this)
        class(pacspointing), intent(inout) :: this

        if (allocated(this%time)) deallocate(this%time)
        if (allocated(this%ra))   deallocate(this%ra)
        if (allocated(this%dec))  deallocate(this%dec)
        if (allocated(this%pa))   deallocate(this%pa)
        if (allocated(this%chop)) deallocate(this%chop)
        if (allocated(this%nsamples)) deallocate(this%nsamples)
        if (allocated(this%first)) deallocate(this%first)
        if (allocated(this%last )) deallocate(this%last)
        if (allocated(this%delta)) deallocate(this%delta)

    end subroutine


end module module_pacspointing
