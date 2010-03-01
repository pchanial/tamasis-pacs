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


contains


    subroutine init_astrometry(header, astr, status)
        use string, only : strlowcase, strupcase
        use module_fitstools, only : ft_readparam
        character(len=*), intent(in)  :: header
        type(astrometry), intent(out) :: astr
        integer, intent(out)          :: status
        integer                       :: count
        integer                       :: has_cd, has_cdelt
        real*8                        :: crota2
        character(len=70)             :: buffer

        has_cd = 0
        has_cdelt = 0

        call ft_readparam(header, 'naxis1', count, astr%naxis(1), must_exist=.true., status=status)
        if (status /= 0) return

        call ft_readparam(header, 'naxis2', count, astr%naxis(2), must_exist=.true., status=status)
        if (status /= 0) return

        call ft_readparam(header, 'crpix1', count, astr%crpix(1), must_exist=.true., status=status)
        if (status /= 0) return

        call ft_readparam(header, 'crpix2', count, astr%crpix(2), must_exist=.true., status=status)
        if (status /= 0) return

        call ft_readparam(header, 'crval1', count, astr%crval(1), must_exist=.true., status=status)
        if (status /= 0) return

        call ft_readparam(header, 'crval2', count, astr%crval(2), must_exist=.true., status=status)
        if (status /= 0) return

        call ft_readparam(header, 'cdelt1', count, astr%cdelt(1), must_exist=.true., status=status)
        if (status /= 0) return
        if (count /= 0) has_cdelt = has_cdelt + 1

        call ft_readparam(header, 'cdelt2', count, astr%cdelt(2), must_exist=.true., status=status)
        if (status /= 0) return
        if (count /= 0) has_cdelt = has_cdelt + 1

        call ft_readparam(header, 'cd1_1', count, astr%cd(1,1), must_exist=.false., status=status)
        if (status /= 0) return
        if (count /= 0) has_cd = has_cd + 1

        call ft_readparam(header, 'cd2_1', count, astr%cd(2,1), must_exist=.false., status=status)
        if (status /= 0) return
        if (count /= 0) has_cd = has_cd + 1

        call ft_readparam(header, 'cd1_2', count, astr%cd(1,2), must_exist=.false., status=status)
        if (status /= 0) return
        if (count /= 0) has_cd = has_cd + 1

        call ft_readparam(header, 'cd2_2', count, astr%cd(2,2), must_exist=.false., status=status)
        if (status /= 0) return
        if (count /= 0) has_cd = has_cd + 1

        call ft_readparam(header, 'cd2_2', count, astr%cd(2,2), must_exist=.false., status=status)
        if (status /= 0) return

        call ft_readparam(header, 'ctype1', count, buffer, must_exist=.true., status=status)
        if (status /= 0) return
        astr%ctype(1) = strupcase(buffer(1:8))

        call ft_readparam(header, 'ctype2', count, buffer, must_exist=.true., status=status)
        if (status /= 0) return
        astr%ctype(2) = strupcase(buffer(1:8))

        call ft_readparam(header, 'cunit1', count, buffer, must_exist=.false., status=status)
        if (status /= 0) return
        astr%cunit(1) = strlowcase(buffer(1:8))

        call ft_readparam(header, 'cunit2', count, buffer, must_exist=.false., status=status)
        if (status /= 0) return
        astr%cunit(2) = strlowcase(buffer(1:8))

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
            astr%cd(1,1) =  cos(crota2)
            astr%cd(2,1) = -sin(crota2)
            astr%cd(1,2) =  sin(crota2)
            astr%cd(2,2) =  cos(crota2)
        end if

        if (has_cdelt /=0 .and. has_cdelt /= 2) then
            write (ERROR_UNIT,'(a)') 'Header has incomplete CDELTi.'
            status = 1
            return
        end if

        if (has_cdelt == 0) then
            astr%cdelt = 1.d0
        end if

        if (astr%ctype(1) == 'RA---TAN' .and. astr%ctype(2) == 'DEC--TAN') then
            call init_gnomonic(astr)
        else
            write (ERROR_UNIT,'(a)') "Type '" // astr%ctype(1) // "', '" // astr%ctype(2) // "' is not implemented."
            status = 1
            return
        end if

        call init_rotation(astr)

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
        common /ad2xy_gnomonic/ lambda0, cosphi1, sinphi1

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
        common /ad2xy_gnomonic/ lambda0, cosphi1, sinphi1

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
        common /ad2xy_gnomonic/ lambda0, cosphi1, sinphi1

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
        common /xy2xy_rotation/ cdinv,crpix

        cd(1,:) = astr%cd(1,:) * astr%cdelt(1)
        cd(2,:) = astr%cd(2,:) * astr%cdelt(2)

        cdinv = reshape([cd(2,2), -cd(2,1), -cd(1,2), cd(1,1)], [2,2]) / &
                (cd(1,1)*cd(2,2) - cd(2,1)*cd(1,2))
        crpix = astr%crpix

    end subroutine init_rotation


    !---------------------------------------------------------------------------


    pure elemental subroutine xy2xy_rotation(xsi, eta, x, y)

        real*8, intent(in)  :: xsi, eta
        real*8, intent(out) :: x, y
        real*8              :: cdinv(2,2), crpix(2)
        common /xy2xy_rotation/ cdinv, crpix
!check in place

        x = cdinv(1,1)*xsi + cdinv(1,2)*eta + crpix(1)
        y = cdinv(2,1)*xsi + cdinv(2,2)*eta + crpix(2)

    end subroutine xy2xy_rotation


    !---------------------------------------------------------------------------


    subroutine init_wcslib(header, status)
        use module_wcslib
        use module_fitstools
        character(len=*), intent(in) :: header
        integer, intent(out)         :: status
        integer                      :: wcs(WCSLEN), nx, ny
        common /wcslib/ wcs

        call ft_header2wcs(header, wcs, nx, ny, status)

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


    subroutine free_wcslib()
        use module_wcslib
        integer :: wcs(WCSLEN), status
        common /wcslib/ wcs

        status = wcsfree(wcs)

    end subroutine free_wcslib


 end module module_wcs
