module module_wcs

    use, intrinsic :: ISO_FORTRAN_ENV
    implicit none

    private
    real*8, parameter :: pi = 4.0d0 * atan(1.0d0)
    real*8, parameter :: radeg = 180.d0 / pi

    type, public :: astrometry
        integer :: naxis(2)
        real*8  :: cdelt(2), crpix(2), crval(2), cd(2,2)
        character(len=8) :: ctype(2), cunit(2)
    end type astrometry

    public :: init_astrometry
    public :: print_astrometry
    public :: init_gnomonic
    public :: ad2xy_gnomonic
    public :: ad2xy_gnomonic_vect
    public :: init_rotation
    public :: xy2xy_rotation
    public :: init_wcslib
    public :: ad2xy_wcslib
    public :: free_wcslib
    public :: header2wcslib


contains


    subroutine init_astrometry(header, astr, status)
        use string, only : strlowcase, strupcase
        use module_fitstools, only : ft_readparam
        character(len=*), intent(in)  :: header
        type(astrometry), intent(out), optional :: astr
        integer, intent(out)          :: status

        type(astrometry)              :: myastr
        integer                       :: count
        integer                       :: has_cd, has_cdelt
        real*8                        :: crota2
        character(len=70)             :: buffer

        has_cd = 0
        has_cdelt = 0

        call ft_readparam(header, 'naxis1', count, myastr%naxis(1), must_exist=.true., status=status)
        if (status /= 0) return

        call ft_readparam(header, 'naxis2', count, myastr%naxis(2), must_exist=.true., status=status)
        if (status /= 0) return

        call ft_readparam(header, 'crpix1', count, myastr%crpix(1), must_exist=.true., status=status)
        if (status /= 0) return

        call ft_readparam(header, 'crpix2', count, myastr%crpix(2), must_exist=.true., status=status)
        if (status /= 0) return

        call ft_readparam(header, 'crval1', count, myastr%crval(1), must_exist=.true., status=status)
        if (status /= 0) return

        call ft_readparam(header, 'crval2', count, myastr%crval(2), must_exist=.true., status=status)
        if (status /= 0) return

        call ft_readparam(header, 'cdelt1', count, myastr%cdelt(1), must_exist=.true., status=status)
        if (status /= 0) return
        if (count /= 0) has_cdelt = has_cdelt + 1

        call ft_readparam(header, 'cdelt2', count, myastr%cdelt(2), must_exist=.true., status=status)
        if (status /= 0) return
        if (count /= 0) has_cdelt = has_cdelt + 1

        call ft_readparam(header, 'cd1_1', count, myastr%cd(1,1), must_exist=.false., status=status)
        if (status /= 0) return
        if (count /= 0) has_cd = has_cd + 1

        call ft_readparam(header, 'cd2_1', count, myastr%cd(2,1), must_exist=.false., status=status)
        if (status /= 0) return
        if (count /= 0) has_cd = has_cd + 1

        call ft_readparam(header, 'cd1_2', count, myastr%cd(1,2), must_exist=.false., status=status)
        if (status /= 0) return
        if (count /= 0) has_cd = has_cd + 1

        call ft_readparam(header, 'cd2_2', count, myastr%cd(2,2), must_exist=.false., status=status)
        if (status /= 0) return
        if (count /= 0) has_cd = has_cd + 1

        call ft_readparam(header, 'cd2_2', count, myastr%cd(2,2), must_exist=.false., status=status)
        if (status /= 0) return

        call ft_readparam(header, 'ctype1', count, buffer, must_exist=.true., status=status)
        if (status /= 0) return
        myastr%ctype(1) = strupcase(buffer(1:8))

        call ft_readparam(header, 'ctype2', count, buffer, must_exist=.true., status=status)
        if (status /= 0) return
        myastr%ctype(2) = strupcase(buffer(1:8))

        call ft_readparam(header, 'cunit1', count, buffer, must_exist=.false., status=status)
        if (status /= 0) return
        myastr%cunit(1) = strlowcase(buffer(1:8))

        call ft_readparam(header, 'cunit2', count, buffer, must_exist=.false., status=status)
        if (status /= 0) return
        myastr%cunit(2) = strlowcase(buffer(1:8))

        if (has_cd /= 0 .and. has_cd /= 4) then
            write (ERROR_UNIT,'(a)') 'Header has incomplete CD matrix.'
            status = 1
            return
        end if

        if (has_cd == 0) then
            call ft_readparam(header, 'crota2', count, crota2, must_exist=.false., status=status)
            if (status /= 0) return
            if (count == 0) then
                write (ERROR_UNIT,'(a)') 'Header has no definition for the CD matrix'
                status = 1
                return
            end if
            
            crota2 = crota2 / radeg
            myastr%cd(1,1) =  cos(crota2)
            myastr%cd(2,1) = -sin(crota2)
            myastr%cd(1,2) =  sin(crota2)
            myastr%cd(2,2) =  cos(crota2)
        end if

        if (has_cdelt /=0 .and. has_cdelt /= 2) then
            write (ERROR_UNIT,'(a)') 'Header has incomplete CDELTi.'
            status = 1
            return
        end if

        if (has_cdelt == 0) then
            myastr%cdelt = 1.d0
        end if

        if (myastr%ctype(1) == 'RA---TAN' .and. myastr%ctype(2) == 'DEC--TAN') then
            call init_gnomonic(myastr)
        else
            write (ERROR_UNIT,'(a)') "Type '" // myastr%ctype(1) // "', '" // myastr%ctype(2) // "' is not implemented."
            status = 1
            return
        end if

        call init_rotation(myastr)

        if (present(astr)) astr = myastr

    end subroutine init_astrometry


    !---------------------------------------------------------------------------


    subroutine print_astrometry(astr)
        use string, only : strinteger
        type(astrometry), intent(in) :: astr
        integer :: naxis, i

        naxis = 0
        do i=1, size(astr%naxis)
            if (astr%naxis(i) /= 0) naxis = naxis + 1
        end do

        write (*,*) 'NAXIS: ', strinteger(naxis), ' (', strinteger(astr%naxis(1)), ',', strinteger(astr%naxis(2)),')'
        write (*,*) 'CDELT: ', astr%cdelt
        write (*,*) 'CRPIX: ', astr%crpix
        write (*,*) 'CRVAL: ', astr%crval
        write (*,*) 'CD   : ', astr%cd(1,:)
        write (*,*) '       ', astr%cd(2,:)
        write (*,*) 'CUNIT: ', astr%cunit(1), ', ', astr%cunit(2)
        write (*,*) 'CTYPE: ', astr%ctype(1), ', ', astr%ctype(2)

    end subroutine print_astrometry

    !---------------------------------------------------------------------------


    subroutine init_gnomonic(astr)
        type(astrometry), intent(in) :: astr
        real*8 :: lambda0                ! crval(1) in rad
        real*8 :: phi1, cosphi1, sinphi1 ! cos and sin of crval(2)
        common /gnomonic/ lambda0, cosphi1, sinphi1

        lambda0 = astr%crval(1) / radeg
        phi1 = astr%crval(2) / radeg
        cosphi1 = cos(phi1)
        sinphi1 = sin(phi1)

    end subroutine init_gnomonic


    !---------------------------------------------------------------------------


    pure function ad2xy_gnomonic(ad) result(xy)
        real*8, intent(in) :: ad(:,:)          ! R.A. and declination in degrees
        real*8             :: xy(size(ad,1),size(ad,2))
        real*8             :: lambda, phi, invcosc, xsi, eta
        integer            :: i
        real*8             :: lambda0          ! crval[0] in rad
        real*8             :: cosphi1, sinphi1 ! cos and sin of crval[1]
        common /gnomonic/ lambda0, cosphi1, sinphi1

        do i = 1, size(ad,2)
            lambda = ad(1,i) / radeg
            phi = ad(2,i) / radeg
            invcosc = radeg/(sinphi1*sin(phi)+cosphi1*cos(phi)*cos(lambda-lambda0))
            xsi = invcosc * cos(phi)*sin(lambda-lambda0)
            eta = invcosc * (cosphi1*sin(phi)-sinphi1*cos(phi)*cos(lambda-lambda0))

            call xy2xy_rotation(xsi, eta, xy(1,i), xy(2,i))
        end do

    end function ad2xy_gnomonic


    !---------------------------------------------------------------------------


    pure elemental subroutine ad2xy_gnomonic_vect(a, d, x, y)
        real*8, intent(in)  :: a, d             ! R.A. and declination in degrees
        real*8, intent(out) :: x, y
        real*8              :: lambda, phi, invcosc, xsi, eta
        real*8              :: lambda0          ! crval[0] in rad
        real*8              :: cosphi1, sinphi1 ! cos and sin of crval[1]
        common /gnomonic/ lambda0, cosphi1, sinphi1

        lambda = a / radeg
        phi = d / radeg
        invcosc = radeg /(sinphi1*sin(phi)+cosphi1*cos(phi)*cos(lambda-lambda0))
        xsi = invcosc * cos(phi)*sin(lambda-lambda0)
        eta = invcosc * (cosphi1*sin(phi)-sinphi1*cos(phi)*cos(lambda-lambda0))

        call xy2xy_rotation(xsi, eta, x, y)

    end subroutine ad2xy_gnomonic_vect


    !---------------------------------------------------------------------------


    subroutine init_rotation(astr)

        type(astrometry), intent(in) :: astr
        real*8 :: cdinv(2,2), crpix(2)
        real*8 :: cd(2,2)
        common /rotation/ cdinv,crpix

        cd(1,:) = astr%cd(1,:) * astr%cdelt(1)
        cd(2,:) = astr%cd(2,:) * astr%cdelt(2)

        cdinv = reshape([cd(2,2), -cd(2,1), -cd(1,2), cd(1,1)], [2,2]) / &
                (cd(1,1)*cd(2,2) - cd(2,1)*cd(1,2))
        crpix = astr%crpix

    end subroutine init_rotation


    !---------------------------------------------------------------------------


    elemental subroutine xy2xy_rotation(xsi, eta, x, y)

        real*8, intent(in)  :: xsi, eta
        real*8, intent(out) :: x, y
        real*8              :: cdinv(2,2), crpix(2)
        common /rotation/ cdinv, crpix
!check in place

        x = cdinv(1,1)*xsi + cdinv(1,2)*eta + crpix(1)
        y = cdinv(2,1)*xsi + cdinv(2,2)*eta + crpix(2)

    end subroutine xy2xy_rotation


    !---------------------------------------------------------------------------


    subroutine init_wcslib(header, status)
        use module_fitstools
        use module_wcslib, only : WCSLEN
        character(len=*), intent(in) :: header
        integer, intent(out)         :: status
        integer                      :: wcs(WCSLEN), nx, ny
        common /wcslib/ wcs

        call header2wcslib(header, wcs, nx, ny, status)

    end subroutine init_wcslib


    !---------------------------------------------------------------------------


    recursive function ad2xy_wcslib(ad) result(xy)
        use module_wcslib
        real*8, intent(in)   :: ad(:,:)
        integer              :: status
        integer              :: wcs(WCSLEN)
        real*8               :: xy(size(ad,1),size(ad,2))
        integer              :: stat(size(ad,2))
        real*8               :: phi(size(ad,2)), theta(size(ad,2)), imgcrd(size(ad,2),size(ad,2))
        common /wcslib/ wcs

        status = wcss2p(wcs, size(ad,2), 2, ad, phi, theta, imgcrd, xy, stat)

    end function ad2xy_wcslib


    !---------------------------------------------------------------------------


    ! translate the astrometry in a FITS header into a wcslib wcs
    ! XXX this routine currently leaks wcsp
    ! investigate how to keep wcs without wcsp
    subroutine header2wcslib(header, wcs, nx, ny, status)

        use, intrinsic :: ISO_C_BINDING
        use module_fitstools, only : ft_readparam
        use module_wcslib
        character(len=*), intent(in)     :: header
        integer(kind=C_INT), intent(out) :: wcs(WCSLEN)
        integer, intent(out)             :: nx, ny
        integer, intent(out)             :: status
        integer(kind=C_INT)              :: wcs_(WCSLEN)

        integer(kind=C_INT)              :: wcsp
        integer                          :: nreject, nwcs, statfix(WCSFIX_NWCS), count1, count2

        status = wcspih(header // C_NULL_CHAR, len(header)/80, WCSHDR_all, 0, nreject, nwcs, wcsp)
        if (status /= 0) then
            write (ERROR_UNIT,'(a)') 'WCSPIH error: ', status
            return
        end if

        if (nwcs == 0) then
            status = 1
            write (ERROR_UNIT,'(a)') 'ft_header2wcs: Header has no astrometry.'
            return
        end if

        if (wcsp < 0) then
            status = 1
            write (ERROR_UNIT,'(a)') 'ft_header2wcs: wcsp negative. Unrecoverable error.'
            return
        end if

        status = wcsvcopy(wcsp, 0, wcs_)
        if (status /= 0) then
            write (ERROR_UNIT,'(a)') 'ft_header2wcs: wcsvcopy error: ', status
            return
        end if
        status = wcsput(wcs, WCS_FLAG, -1_C_INT, 0_C_INT, 0_C_INT)
        if (status /= 0) then
            write (ERROR_UNIT,'(a)') 'ft_header2wcs: wcsput error: ', status
            return
        end if
        status = wcscopy(wcs_, wcs)
        if (status /= 0) then
            write (ERROR_UNIT,'(a)') 'ft_header2wcs: wcscopy error: ', status
            return
        end if
        status = wcsfree(wcs_)
        if (status /= 0) then
            write (ERROR_UNIT,'(a)') 'ft_header2wcs: wcsfree error: ', status
            return
        end if
        status = wcsvfree(nwcs, wcsp)
        if (status /= 0) then
            write (ERROR_UNIT,'(a)') 'ft_header2wcs: wcsvfree error: ', status
            return
        end if

        !Fix non-standard WCS keyvalues.
        status = wcsfix(ctrl=7, naxis=c_null_ptr, wcs=wcs, stat=statfix)
        if (status /= 0) then
            write (ERROR_UNIT,'(a)') 'ft_header2wcs: wcsfix error: ', status
            return
        end if

        status = wcsset(wcs)
        if (status /= 0) then
            write (ERROR_UNIT,'(a)') 'ft_header2wcs: wcsset error: ', status
            return
        end if

        ! extract NAXIS1 and NAXIS2 keywords
        call ft_readparam(header, 'naxis1', count1, nx, status=status)
        if (status /= 0) return
        call ft_readparam(header, 'naxis2', count2, ny, status=status)
        if (status /= 0) return

        if (count1 == 0 .or. count2 == 0) then
            status = 1
            write (ERROR_UNIT,'(a)') "Missing NAXISn keyword(s) in header."
            return
        end if

    end subroutine header2wcslib


    !---------------------------------------------------------------------------


    subroutine free_wcslib()
        use module_wcslib
        integer :: wcs(WCSLEN), status
        common /wcslib/ wcs

        status = wcsfree(wcs)

    end subroutine free_wcslib


 end module module_wcs
