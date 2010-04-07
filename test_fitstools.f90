program test_fitstools
    use module_fitstools
    implicit none

    character(len=*), parameter :: filename_header = 'tests/csh_header_ngc6946.fits'

    integer                     :: count, status
    character(len=2880)         :: header
    character(len=160)          :: header_too_small
    real*8                      :: image(10,15)
    logical                     :: bvalue
    integer*4                   :: ivalue
    integer*8                   :: lvalue
    real*8                      :: dvalue
    character(len=70)           :: cvalue

    ! test extraction of header in FITS file
    call ft_header2str(filename_header, header_too_small, status=status)
    if (status == 0) stop 'FAILED: ft_header2str'

    call ft_header2str(filename_header, header, status=status)
    if (status /= 0) stop 'FAILED: ft_header2str'

    call ft_read_parameter(header, 'donotexist', ivalue, count, status=status)
    if (status /= 0 .or. count /= 0) stop 'FAILED: donotexist1'

    call ft_read_parameter(header, 'donotexist', lvalue, count, status=status)
    if (status /= 0 .or. count /= 0) stop 'FAILED: donotexist2'

    call ft_read_parameter(header, 'donotexist', dvalue, count, status=status)
    if (status /= 0 .or. count /= 0) stop 'FAILED: donotexist3'

    call ft_read_parameter(header, 'donotexist', cvalue, count, status=status)
    if (status /= 0 .or. count /= 0) stop 'FAILED: donotexist4'

    call ft_read_parameter(header, 'donotexist', lvalue, count, status=status)
    if (status /= 0 .or. count /= 0) stop 'FAILED: donotexist5'

    call ft_read_parameter(header, 'naxis', ivalue, count, status=status)
    if (status /= 0 .or. ivalue /= 2) stop 'FAILED: naxis'

    call ft_read_parameter(header, 'naxis1', lvalue, count, status=status)
    if (status /= 0 .or. lvalue /= 487) stop 'FAILED: naxis1'

    call ft_read_parameter(header, 'naxis2', lvalue, count, status=status)
    if (status /= 0 .or. lvalue /= 487) stop 'FAILED: naxis2'

    call ft_read_parameter(header, 'cdelt1', dvalue, count, status=status)
    if (status /= 0 .or. dvalue /= -0.00083333333333333d0) stop 'FAILED: cdelt1'

    call ft_read_parameter(header, 'cdelt2', dvalue, count, status=status)
    if (status /= 0 .or. dvalue /= 0.00083333333333333d0) stop 'FAILED: cdelt2'

    call ft_read_parameter(header, 'crota2', dvalue, count, status=status)
    if (status /= 0 .or. dvalue /= 0) stop 'FAILED.'

    call ft_read_parameter(header, 'crpix1', dvalue, count, status=status)
    if (status /= 0 .or. dvalue /= 1) stop 'FAILED.'

    call ft_read_parameter(header, 'crpix2', dvalue, count, status=status)
    if (status /= 0 .or. dvalue /= 1) stop 'FAILED.'

    call ft_read_parameter(header, 'crval1', dvalue, count, status=status)
    if (status /= 0 .or. dvalue /= 308.3404249521426d0) stop 'FAILED: crval1'

    call ft_read_parameter(header, 'crval2', dvalue, count, status=status)
    if (status /= 0 .or. dvalue /= 60.357773247484651d0) stop 'FAILED: crval2'

    call ft_read_parameter(header, 'ctype1', cvalue, count, status=status)
    if (status /= 0 .or. cvalue /= 'RA---TAN') stop 'FAILED: ctype1'

    call ft_read_parameter(header, 'ctype2', cvalue, count, status=status)
    if (status /= 0 .or. cvalue /= 'DEC--TAN') stop 'FAILED: ctype2'

    call ft_read_parameter(header, 'cunit1', cvalue, count, status=status)
    if (status /= 0 .or. cvalue /= 'deg') stop 'FAILED: cunit1'

    call ft_read_parameter(header, 'cunit2', cvalue, count, status=status)
    if (status /= 0 .or. cvalue /= 'deg') stop 'FAILED: cunit2'

    call ft_read_parameter(header, 'testme1', bvalue, count, status=status)
    if (status /= 0 .or. bvalue .neqv. .true.) stop 'FAILED: testme1'

    call ft_read_parameter(header, 'testme2', bvalue, count, status=status)
    if (status == 0) stop 'FAILED: testme2'

    call ft_read_parameter(header, 'testme3', bvalue, count, status=status)
    if (status /= 0 .or. bvalue .neqv. .false.) stop 'FAILED: testme3'

    call ft_read_parameter(header, 'testme4', bvalue, count, status=status)
    if (status == 0) stop 'FAILED: testme4'

    call ft_read_parameter(header, 'testme5', ivalue, count, status=status)
    if (status == 0) stop 'FAILED: testme5a'

    call ft_read_parameter(header, 'testme5', lvalue, count, status=status)
    if (status == 0) stop 'FAILED: testme5b'

    call ft_read_parameter(header, 'testme6', dvalue, count, status=status)
    if (status == 0) stop 'FAILED: testme6'

    call ft_read_parameter(header, 'testme7', cvalue, count, status=status)
    if (status /= 0 .or. cvalue /= "3'D") stop 'FAILED: testme7'

    call ft_read_parameter(header, 'testme8', cvalue, count, status=status)
    if (status /= 0 .or. cvalue /= "'/") stop 'FAILED: testme8'

    call ft_read_parameter(header, 'testme9', cvalue, count, status=status)
    if (status == 0) stop 'FAILED: testme9'

    call ft_read_parameter(header, 'testme10', cvalue, count, status=status)
    if (status == 0) stop 'FAILED: testme10'

    ! test writing an image
    status = 0
    image = 0
    image(1,2) = 1
    call ft_write('/tmp/test_fitstools.fits', image, header, status)
    if (status /= 0) stop 'FAILED: ft_write'

    stop 'OK.'

end program test_fitstools
