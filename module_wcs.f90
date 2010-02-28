module module_wcs

    implicit none

    private
    real*8, parameter :: pi = 4.0d0 * atan(1.0d0)
    real*8, parameter :: radeg = 180.d0 / pi

    type, public :: astrometry
        integer :: naxis(2)
        real*8  :: crval(2), crpix(2), cdelt(2), cd(2,2)
        character(len=8) :: ctype(2), cunit(2)
    end type astrometry

    public :: init_astrometry
    public :: init_gnomonic
    public :: ad2xy_gnomonic
    public :: init_rotation
    public :: xy2xy_rotation
    public :: init_wcslib
    public :: ad2xy_wcslib
    public :: free_wcslib


contains


    subroutine init_astrometry(header, astr, status)
        use module_fitstools, only : ft_readparam
        character(len=*), intent(in)  :: header
        type(astrometry), intent(out) :: astr
        integer, intent(out)          :: status
        integer                       :: count

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
        call ft_readparam(header, 'cdelt2', count, astr%cdelt(2), must_exist=.true., status=status)
        if (status /= 0) return
        call ft_readparam(header, 'cd1_1', count, astr%cd(1,1), must_exist=.false., status=status)
        if (status /= 0) return
        call ft_readparam(header, 'cd2_1', count, astr%cd(2,1), must_exist=.false., status=status)
        if (status /= 0) return
        call ft_readparam(header, 'cd1_2', count, astr%cd(1,2), must_exist=.false., status=status)
        if (status /= 0) return
        call ft_readparam(header, 'cd2_2', count, astr%cd(2,2), must_exist=.false., status=status)
        if (status /= 0) return
        call ft_readparam(header, 'ctype1', count, astr%ctype(1), must_exist=.true., status=status)
        if (status /= 0) return
        call ft_readparam(header, 'ctype2', count, astr%ctype(2), must_exist=.true., status=status)
        if (status /= 0) return
        call ft_readparam(header, 'cunit1', count, astr%cunit(1), must_exist=.false., status=status)
        if (status /= 0) return
        call ft_readparam(header, 'cunit2', count, astr%cunit(2), must_exist=.false., status=status)
        if (status /= 0) return

    end subroutine init_astrometry


    !-------------------------------------------------------------------------------


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


    !-------------------------------------------------------------------------------


    pure subroutine ad2xy_gnomonic(a, d, x, y)
        real*8, intent(in)  :: a, d             ! R.A. and declination in degrees
        real*8, intent(out) :: x, y
        real*8              :: lambda, phi, invcosc
        real*8              :: lambda0          ! crval[0] in rad
        real*8              :: cosphi1, sinphi1 ! cos and sin of crval[1]
        common /ad2xy_gnomonic/ lambda0, cosphi1, sinphi1

        lambda = a / radeg
        phi = d / radeg
        invcosc = 1.d0/(sinphi1*sin(phi)+cosphi1*cos(phi)*cos(lambda-lambda0))
        x = radeg * invcosc * cos(phi)*sin(lambda-lambda0)
        y = radeg * invcosc * (cosphi1*sin(phi)-sinphi1*cos(phi)*cos(lambda-lambda0))

    end subroutine ad2xy_gnomonic


    !-------------------------------------------------------------------------------


    subroutine init_rotation(astr)

        type(astrometry), intent(in) :: astr
        real*8 :: cdinv(2,2), crpix(2)
        real*8 :: cd(2,2), det
        common /xy2xy_rotation_params/ cdinv,crpix

        if (astr%cdelt(1) /= 1.0d0) then
            cd(1,:) = astr%cd(1,:) * astr%cdelt(1)
            cd(2,:) = astr%cd(2,:) * astr%cdelt(2)
        endif

        crpix = astr%crpix - 1
        det = 1.d0/(astr%cd(1,1)*astr%cd(2,2)-astr%cd(2,1)*astr%cd(1,2))
        cdinv = reshape([astr%cd(2,2), -astr%cd(2,1), -astr%cd(1,2), astr%cd(1,1)], [2,2])
        cdinv = cdinv / det

    end subroutine init_rotation


    !-------------------------------------------------------------------------------


    pure subroutine xy2xy_rotation(xsi, eta, x, y)

        real*8, intent(in)  :: xsi, eta
        real*8, intent(out) :: x, y
        real*8              :: cdinv(2,2), crpix(2)
        common /xy2xy_rotation/ cdinv, crpix
!check in place
        x = cdinv(1,1)*xsi + cdinv(1,2)*eta + crpix(1)
        y = cdinv(2,1)*xsi + cdinv(2,2)*eta + crpix(2)

    end subroutine xy2xy_rotation


    !-------------------------------------------------------------------------------


    subroutine init_wcslib(header, status)
        use module_wcslib
        use module_fitstools
        character(len=*), intent(in) :: header
        integer, intent(out)         :: status
        integer                      :: wcs(WCSLEN), nx, ny
        common /wcslib/ wcs

        call ft_header2wcs(header, wcs, nx, ny, status)

    end subroutine init_wcslib


    !-------------------------------------------------------------------------------


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


    !-------------------------------------------------------------------------------


    subroutine free_wcslib()
        use module_wcslib
        integer :: wcs(WCSLEN), status
        common /wcslib/ wcs

        status = wcsfree(wcs)

    end subroutine free_wcslib


 end module module_wcs
