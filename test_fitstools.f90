program test_fitstools

    use module_fitstools
    implicit none

    character(len=*), parameter :: filename_header = 'tests/csh_header_ngc6946.fits'

    integer                     :: count, status
    character(len=28800)        :: header
    real*8                      :: image(10,15)
    logical                     :: bvalue
    integer*4                   :: ivalue
    integer*8                   :: lvalue
    real*8                      :: dvalue
    character(len=70)           :: cvalue

    ! test extraction of header in FITS file
    call ft_header2str(filename_header, header, status=status)
    if (status /= 0) stop 'FAILED: ft_header2str'

    call ft_readparam(header, 'donotexist', count, ivalue, status=status)
    if (status /= 0 .or. count /= 0) stop 'FAILED: donotexist1'

    call ft_readparam(header, 'donotexist', count, lvalue, status=status)
    if (status /= 0 .or. count /= 0) stop 'FAILED: donotexist2'

    call ft_readparam(header, 'donotexist', count, dvalue, status=status)
    if (status /= 0 .or. count /= 0) stop 'FAILED: donotexist3'

    call ft_readparam(header, 'donotexist', count, cvalue, status=status)
    if (status /= 0 .or. count /= 0) stop 'FAILED: donotexist4'

    call ft_readparam(header, 'donotexist', count, lvalue, status=status)
    if (status /= 0 .or. count /= 0) stop 'FAILED: donotexist5'

    call ft_readparam(header, 'naxis', count, ivalue, status=status)
    if (status /= 0 .or. ivalue /= 2) stop 'FAILED: naxis'

    call ft_readparam(header, 'naxis1', count, lvalue, status=status)
    if (status /= 0 .or. lvalue /= 487) stop 'FAILED: naxis1'

    call ft_readparam(header, 'naxis2', count, lvalue, status=status)
    if (status /= 0 .or. lvalue /= 487) stop 'FAILED: naxis2'

    call ft_readparam(header, 'cdelt1', count, dvalue, status=status)
    if (status /= 0 .or. dvalue /= -0.00083333333333333d0) stop 'FAILED: cdelt1'

    call ft_readparam(header, 'cdelt2', count, dvalue, status=status)
    if (status /= 0 .or. dvalue /= 0.00083333333333333d0) stop 'FAILED: cdelt2'

    call ft_readparam(header, 'crota2', count, dvalue, status=status)
    if (status /= 0 .or. dvalue /= 0) stop 'FAILED.'

    call ft_readparam(header, 'crpix1', count, dvalue, status=status)
    if (status /= 0 .or. dvalue /= 1) stop 'FAILED.'

    call ft_readparam(header, 'crpix2', count, dvalue, status=status)
    if (status /= 0 .or. dvalue /= 1) stop 'FAILED.'

    call ft_readparam(header, 'crval1', count, dvalue, status=status)
    if (status /= 0 .or. dvalue /= 308.3404249521426d0) stop 'FAILED: crval1'

    call ft_readparam(header, 'crval2', count, dvalue, status=status)
    if (status /= 0 .or. dvalue /= 60.357773247484651d0) stop 'FAILED: crval2'

    call ft_readparam(header, 'ctype1', count, cvalue, status=status)
    if (status /= 0 .or. cvalue /= 'RA---TAN') stop 'FAILED: ctype1'

    call ft_readparam(header, 'ctype2', count, cvalue, status=status)
    if (status /= 0 .or. cvalue /= 'DEC--TAN') stop 'FAILED: ctype2'

    call ft_readparam(header, 'cunit1', count, cvalue, status=status)
    if (status /= 0 .or. cvalue /= 'deg') stop 'FAILED: cunit1'

    call ft_readparam(header, 'cunit2', count, cvalue, status=status)
    if (status /= 0 .or. cvalue /= 'deg') stop 'FAILED: cunit2'

    call ft_readparam(header, 'testme1', count, bvalue, status=status)
    if (status /= 0 .or. bvalue .neqv. .true.) stop 'FAILED: testme1'

    call ft_readparam(header, 'testme2', count, bvalue, status=status)
    if (status == 0) stop 'FAILED: testme2'

    call ft_readparam(header, 'testme3', count, bvalue, status=status)
    if (status /= 0 .or. bvalue .neqv. .false.) stop 'FAILED: testme3'

    call ft_readparam(header, 'testme4', count, bvalue, status=status)
    if (status == 0) stop 'FAILED: testme4'

    call ft_readparam(header, 'testme5', count, ivalue, status=status)
    if (status == 0) stop 'FAILED: testme5a'

    call ft_readparam(header, 'testme5', count, lvalue, status=status)
    if (status == 0) stop 'FAILED: testme5b'

    call ft_readparam(header, 'testme6', count, dvalue, status=status)
    if (status == 0) stop 'FAILED: testme6'

    call ft_readparam(header, 'testme7', count, cvalue, status=status)
    if (status /= 0 .or. cvalue /= "3'D") stop 'FAILED: testme7'

    call ft_readparam(header, 'testme8', count, cvalue, status=status)
    if (status /= 0 .or. cvalue /= "'/") stop 'FAILED: testme8'

    call ft_readparam(header, 'testme9', count, cvalue, status=status)
    if (status == 0) stop 'FAILED: testme9'

    call ft_readparam(header, 'testme10',count, cvalue, status=status)
    if (status == 0) stop 'FAILED: testme10'

    ! test writing an image
    status = 0
    image = 0
    image(1,2) = 1
    call ft_write('/tmp/test_fitstools.fits', image, header, status)
    if (status /= 0) stop 'FAILED: ft_write'

    stop 'OK.'

end program test_fitstools
