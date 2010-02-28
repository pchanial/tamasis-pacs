program test_fitstools

    use module_fitstools
    use module_wcslib, only : wcsfree, WCSLEN
    implicit none

    character(len=*), parameter :: filename_header = 'tests/csh_header_ngc6946.fits'

    integer                     :: wcs(WCSLEN)
    integer                     :: count, status, nx, ny
    character(len=2880)         :: header
    real*8                      :: image(10,15)
    logical                     :: test
    logical                     :: bvalue
    integer*4                   :: ivalue
    integer*8                   :: lvalue
    real*8                      :: dvalue
    character(len=70)           :: cvalue

    ! test extract of wcs in FITS header
    status = 0
    call ft_header2str(filename_header, header, status=status)
    if (status /= 0) stop 'FAILED.'
    call ft_header2wcs(header, wcs, nx, ny, status=status)
    if (status /= 0) stop 'FAILED.'

    test = .true.
    call ft_readparam(header, 'donotexist', count, ivalue, status=status)
    if (status /= 0) stop 'FAILED.'
    test = test .and. count == 0
    call ft_readparam(header, 'donotexist', count, lvalue, status=status)
    if (status /= 0) stop 'FAILED.'
    test = test .and. count == 0
    call ft_readparam(header, 'donotexist', count, dvalue, status=status)
    if (status /= 0) stop 'FAILED.'
    test = test .and. count == 0
    call ft_readparam(header, 'donotexist', count, cvalue, status=status)
    if (status /= 0) stop 'FAILED.'
    test = test .and. count == 0
    call ft_readparam(header, 'donotexist', count, lvalue, status=status)
    if (status /= 0) stop 'FAILED.'
    test = test .and. count == 0
    call ft_readparam(header, 'naxis', count, ivalue, status=status)
    if (status /= 0) stop 'FAILED.'
    test = test .and. ivalue == 2
    call ft_readparam(header, 'naxis1', count, lvalue, status=status)
    if (status /= 0) stop 'FAILED.'
    test = test .and. lvalue == 487
    call ft_readparam(header, 'naxis2', count, lvalue, status=status)
    if (status /= 0) stop 'FAILED.'
    test = test .and. lvalue == 487
    call ft_readparam(header, 'cdelt1', count, dvalue, status=status)
    if (status /= 0) stop 'FAILED.'
    test = test .and. dvalue == 0.00083333333333333d0
    call ft_readparam(header, 'cdelt2', count, dvalue, status=status)
    if (status /= 0) stop 'FAILED.'
    test = test .and. dvalue == -0.00083333333333333d0
    call ft_readparam(header, 'crota2', count, dvalue, status=status)
    if (status /= 0) stop 'FAILED.'
    test = test .and. dvalue == 0
    call ft_readparam(header, 'crpix1', count, dvalue, status=status)
    if (status /= 0) stop 'FAILED.'
    test = test .and. dvalue == 1
    call ft_readparam(header, 'crpix2', count, dvalue, status=status)
    if (status /= 0) stop 'FAILED.'
    test = test .and. dvalue == 1
    call ft_readparam(header, 'crval1', count, dvalue, status=status)
    if (status /= 0) stop 'FAILED.'
    test = test .and. dvalue == 308.3404249521426d0
    call ft_readparam(header, 'crval2', count, dvalue, status=status)
    if (status /= 0) stop 'FAILED.'
    test = test .and. dvalue == 60.357773247484651d0
    call ft_readparam(header, 'ctype1', count, cvalue, status=status)
    if (status /= 0) stop 'FAILED.'
    test = test .and. cvalue == 'RA---TAN'
    call ft_readparam(header, 'ctype2', count, cvalue, status=status)
    if (status /= 0) stop 'FAILED.'
    test = test .and. cvalue == 'DEC--TAN'
    call ft_readparam(header, 'cunit1', count, cvalue, status=status)
    if (status /= 0) stop 'FAILED.'
    test = test .and. cvalue == 'deg'
    call ft_readparam(header, 'cunit2', count, cvalue, status=status)
    if (status /= 0) stop 'FAILED.'
    test = test .and. cvalue == 'deg'
    call ft_readparam(header, 'testme1', count, bvalue, status=status)
    if (status /= 0) stop 'FAILED.'
    test = test .and. (bvalue .eqv. .true.)
    call ft_readparam(header, 'testme2', count, bvalue, status=status)
    if (status == 0) stop 'FAILED.'
    call ft_readparam(header, 'testme3', count, bvalue, status=status)
    if (status /= 0) stop 'FAILED.'
    test = test .and. (bvalue .eqv. .false.)
    call ft_readparam(header, 'testme4', count, bvalue, status=status)
    if (status == 0) stop 'FAILED.'
    call ft_readparam(header, 'testme5', count, ivalue, status=status)
    if (status == 0) stop 'FAILED.'
    call ft_readparam(header, 'testme5', count, lvalue, status=status)
    if (status == 0) stop 'FAILED.'
    call ft_readparam(header, 'testme6', count, dvalue, status=status)
    if (status == 0) stop 'FAILED.'
    call ft_readparam(header, 'testme7', count, cvalue, status=status)
    if (status /= 0) stop 'FAILED.'
    test = test .and. (cvalue == "3'D")
    call ft_readparam(header, 'testme8', count, cvalue, status=status)
    if (status /= 0) stop 'FAILED.'
    test = test .and. (cvalue == "'/")
    call ft_readparam(header, 'testme9', count, cvalue, status=status)
    if (status == 0) stop 'FAILED.'
    call ft_readparam(header, 'testme10',count, cvalue, status=status)
    if (status == 0) stop 'FAILED.'

    ! test writing an image
    status = 0
    image = 0
    image(1,2) = 1
    call ft_write('/tmp/test_fitstools.fits', image, wcs, status)
    if (status /= 0) stop 'FAILED.'

    status = wcsfree(wcs)

    if (.not. test) stop 'FAILED.'
    stop 'OK.'

end program test_fitstools
