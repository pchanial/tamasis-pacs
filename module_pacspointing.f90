module module_pacspointing
    implicit none
    private

    public :: pacspointing

    type pacspointing

        private
        integer                           :: nsamples
        real*8                            :: delta
        real*8, allocatable, dimension(:), public :: time, ra, dec, pa, chop
        procedure(get_position_ev), pointer, public :: get_position => null()

    contains

        private
!        generic          :: load => load_filename, load_array
        procedure, public :: load_filename
        procedure         :: load_filename_oldstyle
        procedure, public :: load_array
        procedure, public :: compute_center
        procedure, public :: print
        procedure, public :: get_position_gen
        procedure, public :: destructor ! XXX should be final, but not handled by gfortran 4.5

    end type pacspointing


contains


    subroutine load_filename(this, filename, status)
        use, intrinsic :: ISO_FORTRAN_ENV
        use            :: precision, only : p, dp
        use            :: module_fitstools, only : ft_read_column,   &
                                                   ft_open_bintable, ft_close

        class(pacspointing), intent(inout) :: this
        character(len=*), intent(in)       :: filename
        integer, intent(out)               :: status
        integer                            :: pos

        real(kind=dp), allocatable  :: time(:), ra(:), dec(:), pa(:), chop(:)
        integer*8, allocatable      :: timeus(:)
        integer unit
        integer ncolumns, nrecords

        pos = len_trim(filename)
        if (filename(pos-4:pos) /= '.fits') then
            write (OUTPUT_UNIT,'(a)') 'LOAD_FILENAME: obsolete file format.'
            call this%load_filename_oldstyle(filename, status)
            return
        endif

        call ft_open_bintable(filename // '[Status]', unit, ncolumns, nrecords,&
                              status)
        if (status /= 0) return

        allocate(time  (nrecords))
        allocate(timeus(nrecords))
        allocate(ra    (nrecords))
        allocate(dec   (nrecords))
        allocate(pa    (nrecords))
        allocate(chop  (nrecords))

        call ft_read_column(unit, 'FINETIME', timeus, status)
        if (status /= 0) return
        time = timeus * 1.d-6
        call ft_read_column(unit, 'RaArray', ra, status)
        if (status /= 0) return
        call ft_read_column(unit, 'DecArray', dec, status)
        if (status /= 0) return
        call ft_read_column(unit, 'PaArray', pa, status)
        if (status /= 0) return
        call ft_read_column(unit, 'CHOPFPUANGLE', chop, status)
        if (status /= 0) return

        call ft_close(unit, status)
        if (status /= 0) return

        call this%load_array(time, ra, dec, pa, chop, status)

    end subroutine load_filename


    !---------------------------------------------------------------------------


    subroutine load_filename_oldstyle(this, filename, status)
        use module_fitstools, only : ft_readextension

        class(pacspointing), intent(inout) :: this
        character(len=*), intent(in)       :: filename
        integer, intent(out)               :: status
        real*8, allocatable                :: time(:), ra(:), dec(:), pa(:), chop(:)
        integer*8, allocatable             :: timeus(:)

        call ft_readextension(filename // '_Time.fits', timeus, status)
        if (status /= 0) return
        allocate(time(size(timeus)))
        time = timeus * 1.0d-6

        call ft_readextension(filename // '_RaArray.fits', ra, status)
        if (status /= 0) return

        call ft_readextension(filename // '_DecArray.fits', dec, status)
        if (status /= 0) return

        call ft_readextension(filename // '_PaArray.fits', pa, status)
        if (status /= 0) return

        call ft_readextension(filename // '_ChopFpuAngle.fits', chop, status)
        if (status /= 0) return

        call this%load_array(time, ra, dec, pa, chop, status)

    end subroutine load_filename_oldstyle


    !---------------------------------------------------------------------------


    subroutine load_array(this, time, ra, dec, pa, chop, status)
        use, intrinsic :: ISO_FORTRAN_ENV
        class(pacspointing), intent(inout) :: this
        real*8, intent(in)                 :: time(:), ra(:), dec(:), pa(:), chop(:)
        integer, intent(out)               :: status
        real*8, parameter                  :: tol = 0.001d0 ! fractional threshold of sampling above which sampling is considered uneven
        integer                            :: isample

        status = 1

        ! check conformity of time, ra, dec, pa
        this%nsamples = size(time)
        if (this%nsamples <= 1) then
            write (ERROR_UNIT,'(a)') "Input time has less than two time samples."
            return
        endif
        if (size(ra) /= this%nsamples) then
            write (ERROR_UNIT,'(a)') "Input R.A. has an invalid number of samples."
            return
        endif
        if (size(dec) /= this%nsamples) then
            write (ERROR_UNIT,'(a)') "Input declination has an invalid number of samples."
            return
        endif
        if (size(pa) /= this%nsamples) then
            write (ERROR_UNIT,'(a)') "Input P.A. has an invalid number of samples."
            return
        endif
        if (size(chop) /= this%nsamples) then
            write (ERROR_UNIT,'(a)') "Input chop angle has an invalid number of samples."
            return
        endif
        this%delta = time(2) - time(1)

        ! check that the input time is monotonous (and increasing)
        if (any(time(2:)-time(1:this%nsamples-1) <= 0)) then
            write (ERROR_UNIT,'(a)') "Error: the pointing time is not strictly increasing"
            return
        endif

        status = 0

        ! check there is no time drifts to ensure that fast interpolation can be relied upon
        if (any(abs([(this%delta * (isample-1), isample = 1, this%nsamples)] - (time-time(1))) > tol * this%delta)) then
            write (*,'(a)') "Warning: the pointing time is not evenly spaced or is drifting."
            this%get_position => get_position_gen
        else
            this%get_position => get_position_ev
        end if

        allocate(this%time(this%nsamples))
        allocate(this%ra  (this%nsamples))
        allocate(this%dec (this%nsamples))
        allocate(this%pa  (this%nsamples))
        allocate(this%chop(this%nsamples))

        this%time = time
        this%ra   = ra
        this%dec  = dec
        this%pa   = pa
        this%chop = chop

    end subroutine


    !---------------------------------------------------------------------------


    ! time samples are not evenly spaced, slow linear interpolation
    subroutine get_position_gen(this, time, ra, dec, pa, chop, index)
        !use, intrinsic :: ieee_arithmetic (not implemented in gfortran 4.5)
        use, intrinsic :: ISO_FORTRAN_ENV, only : OUTPUT_UNIT
        class(pacspointing), intent(in) :: this
        real*8, intent(in)              :: time
        real*8, intent(out)             :: ra, dec, pa, chop
        integer, intent(inout)          :: index
        integer                         :: i
        real*8                          :: frac, zero

        if (time > this%time(this%nsamples) + this%delta .or. time < this%time(1) - this%delta) then
            zero = 0.d0
            ra   = 0.d0 / zero
            dec  = 0.d0 / zero
            pa   = 0.d0 / zero
            chop = 0.d0 / zero
            return
        endif

        if (time <= this%time(index-1)) index = 2

        do i = index, this%nsamples
            if (time <= this%time(i)) go to 100
        end do
        i = i - 1

100     frac = (time - this%time(i-1)) / (this%time(i) - this%time(i-1))

        ra   = this%ra  (i-1) * (1 - frac) + this%ra  (i) * frac
        dec  = this%dec (i-1) * (1 - frac) + this%dec (i) * frac
        pa   = this%pa  (i-1) * (1 - frac) + this%pa  (i) * frac
        chop = this%chop(i-1) * (1 - frac) + this%chop(i) * frac

        index = i

    end subroutine get_position_gen


    !---------------------------------------------------------------------------


    ! time samples are assumed to be equally spaced (hence a fast linear interpolation)
    subroutine get_position_ev(this, time, ra, dec, pa, chop, index)
        !use, intrinsic :: ieee_arithmetic (not implemented by gfortran 4.5)
        class(pacspointing), intent(in)  :: this
        real*8, intent(in)               :: time
        real*8, intent(out)              :: ra, dec, pa, chop
        integer, intent(inout)           :: index
        integer                          :: i
        real*8                           :: frac, dsample, zero

        if (time > this%time(this%nsamples) + this%delta .or. time < this%time(1) - this%delta) then
            zero = 0.d0
            ra   = 0.d0 / zero
            dec  = 0.d0 / zero
            pa   = 0.d0 / zero
            chop = 0.d0 / zero
            return
        endif

        dsample = (time-this%time(1)) / this%delta

        i = floor(dsample)
        frac = dsample - i
        if (i < 0) then
            frac = frac + i
            i = 1
        else if (i + 1 >= this%nsamples) then
            frac = frac + i + 2 - this%nsamples
            i = this%nsamples - 1
        else
            i = i + 1
        end if

        ra   = this%ra  (i) * (1 - frac) + this%ra  (i+1) * frac
        dec  = this%dec (i) * (1 - frac) + this%dec (i+1) * frac
        pa   = this%pa  (i) * (1 - frac) + this%pa  (i+1) * frac
        chop = this%chop(i) * (1 - frac) + this%chop(i+1) * frac

    end subroutine get_position_ev


    !------------------------------------------------------------------------------


    ! compute the mean value of R.A. and dec, taking into account RA's singularity at 0
    subroutine compute_center(this, times, ra0, dec0)
        class(pacspointing), intent(in) :: this
        real*8, intent(in)              :: times(:)  ! times at which the PACS array should be considered as inside the map
        real*8, intent(out)             :: ra0, dec0
        integer                         :: ntimes, isample, n180, index
        real*8                          :: ra, dec, pa, chop
        logical                         :: zero_minus, zero_plus

        ntimes = size(times)
        ra0  = 0.0d0
        dec0 = 0.0d0
        zero_minus = .false.
        zero_plus  = .false.
        n180 = 0
        index = 2

        !$omp parallel do default(shared) reduction(+:n180,ra0,dec0) reduction(.or.:zero_minus,zero_plus) &
        !$omp private(isample, ra, dec, pa, chop) firstprivate(index)
        do isample=1, size(times)
           call this%get_position(times(isample), ra, dec, pa, chop, index)
           zero_minus = zero_minus .or. ra > 270.d0
           zero_plus  = zero_plus  .or. ra <= 90.d0
           if (ra >= 180) n180 = n180 + 1
           ra0  = ra0  + ra
           dec0 = dec0 + dec
        end do
        !$omp end parallel do

        if (zero_minus .and. zero_plus) ra0 = ra0 - 360.d0 * n180

        ra0  = ra0  / ntimes
        dec0 = dec0 / ntimes

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
    end subroutine


end module module_pacspointing
