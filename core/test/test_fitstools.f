program test_fitstools

    use iso_fortran_env, only : ERROR_UNIT
    use module_fitstools
    use module_tamasis,  only : p
    implicit none

    character(len=*), parameter :: filename_header = 'core/test/data/header.fits'

    integer                     :: status, unit
    character(len=2880*2)       :: header
    character(len=160)          :: header_too_small
    real(p)                     :: image(10,15)
    real(p), allocatable        :: image_(:,:)
    logical                     :: bvalue
    integer*4                   :: ivalue
    integer*8                   :: lvalue
    real(p)                     :: dvalue
    character(len=70)           :: cvalue
    logical                     :: found

    ! test extraction of header in FITS file
    call ft_header2str(filename_header, header_too_small, status=status)
    if (status == 0) call failure('ft_header2str1')

    call ft_header2str(filename_header, header, status=status)
    if (status /= 0) call failure('ft_header2str2')

    call ft_read_keyword(header, 'donotexist', ivalue, found, status)
    if (status /= 0 .or. found) call failure('donotexist1')

    call ft_read_keyword(header, 'donotexist', lvalue, found, status)
    if (status /= 0 .or. found) call failure('donotexist2')

    call ft_read_keyword(header, 'donotexist', dvalue, found, status)
    if (status /= 0 .or. found) call failure('donotexist3')

    call ft_read_keyword(header, 'donotexist', cvalue, found, status)
    if (status /= 0 .or. found) call failure('donotexist4')

    call ft_read_keyword(header, 'donotexist', lvalue, found, status)
    if (status /= 0 .or. found) call failure('donotexist5')

    call ft_read_keyword(header, 'donotexist', ivalue, status=status)
    if (status == 0) call failure('donotexist1b')

    call ft_read_keyword(header, 'donotexist', lvalue, status=status)
    if (status == 0) call failure('donotexist2b')

    call ft_read_keyword(header, 'donotexist', dvalue, status=status)
    if (status == 0) call failure('donotexist3b')

    call ft_read_keyword(header, 'donotexist', cvalue, status=status)
    if (status == 0) call failure('donotexist4b')

    call ft_read_keyword(header, 'donotexist', lvalue, status=status)
    if (status == 0) call failure('donotexist5b')

    call ft_read_keyword(header, 'naxis', ivalue, found, status)
    if (status /= 0 .or. ivalue /= 2) call failure('naxis')

    call ft_read_keyword(header, 'naxis1', lvalue, found, status)
    if (status /= 0 .or. lvalue /= 487) call failure('naxis1')

    call ft_read_keyword(header, 'naxis2', lvalue, found, status)
    if (status /= 0 .or. lvalue /= 487) call failure('naxis2')

    call ft_read_keyword(header, 'cdelt1', dvalue, found, status)
    if (status /= 0 .or. dvalue /= -0.00083333333333333_p) call failure('cdelt1')

    call ft_read_keyword(header, 'cdelt2', dvalue, found, status)
    if (status /= 0 .or. dvalue /= 0.00083333333333333_p) call failure('cdelt2')

    call ft_read_keyword(header, 'crota2', dvalue, found, status)
    if (status /= 0 .or. dvalue /= 0) call failure('crota2')

    call ft_read_keyword(header, 'crpix1', dvalue, found, status)
    if (status /= 0 .or. dvalue /= 1) call failure('crpix1')

    call ft_read_keyword(header, 'crpix2', dvalue, found, status)
    if (status /= 0 .or. dvalue /= 1) call failure('crpix2')

    call ft_read_keyword(header, 'crval1', dvalue, found, status)
    if (status /= 0 .or. dvalue /= 308.3404249521426_p) call failure('crval1')

    call ft_read_keyword(header, 'crval2', dvalue, found, status)
    if (status /= 0 .or. dvalue /= 60.357773247484651_p) call failure('crval2')

    call ft_read_keyword(header, 'ctype1', cvalue, found, status)
    if (status /= 0 .or. cvalue /= 'RA---TAN') call failure('ctype1')

    call ft_read_keyword(header, 'ctype2', cvalue, found, status)
    if (status /= 0 .or. cvalue /= 'DEC--TAN') call failure('ctype2')

    call ft_read_keyword(header, 'cunit1', cvalue, found, status)
    if (status /= 0 .or. cvalue /= 'deg') call failure('cunit1')

    call ft_read_keyword(header, 'cunit2', cvalue, found, status)
    if (status /= 0 .or. cvalue /= 'deg') call failure('cunit2')

    call ft_read_keyword(header, 'testme1a', bvalue, found, status)
    if (status /= 0 .or. bvalue .neqv. .true.) call failure('testme1a')

    call ft_read_keyword(header, 'testme1b', bvalue, found, status)
    if (status /= 0 .or. bvalue .neqv. .true.) call failure('testme1b')

    call ft_read_keyword(header, 'testme1c', bvalue, found, status)
    if (status /= 0 .or. bvalue .neqv. .true.) call failure('testme1c')

    call ft_read_keyword(header, 'testme2a', bvalue, found, status)
    if (status /= 0 .or. bvalue .neqv. .false.) call failure('testme2a')

    call ft_read_keyword(header, 'testme2b', bvalue, found, status)
    if (status /= 0 .or. bvalue .neqv. .false.) call failure('testme2b')

    call ft_read_keyword(header, 'testme2c', bvalue, found, status)
    if (status /= 0 .or. bvalue .neqv. .false.) call failure('testme2c')

    call ft_read_keyword(header, 'testme4', bvalue, found, status)
    if (status == 0) call failure('testme4')

    call ft_read_keyword(header, 'testme5', ivalue, found, status)
    if (status == 0) call failure('testme5a')

    call ft_read_keyword(header, 'testme5', lvalue, found, status)
    if (status == 0) call failure('testme5b')

    call ft_read_keyword(header, 'testme6', dvalue, found, status)
    if (status == 0) call failure('testme6')

    call ft_read_keyword(header, 'testme7', cvalue, found, status)
    if (status /= 0 .or. cvalue /= "3''D") call failure('testme7')

    call ft_read_keyword(header, 'testme8', cvalue, found, status)
    if (status /= 0 .or. cvalue /= "'/") call failure('testme8')

    call ft_read_keyword(header, 'testme9', cvalue, found, status)
    if (status /= 0 .or. cvalue /= "'/") call failure('testme9')

    call ft_read_keyword(header, 'testme10', cvalue, found, status)
    if (status /= 0 .or. cvalue /= "'/  / myc10") call failure('testme10')

    ! test writing an image
    status = 0
    image = 0
    image(1,2) = 1
    call ft_write_image('/tmp/test_fitstools.fits', image, header, status)
    if (status /= 0) then
        call cleanup
        call failure('ft_write_image')
    end if

    call ft_read_image('/tmp/test_fitstools.fits', image_, status)
    call cleanup
    if (status /= 0) call failure('ft_read_image')
    if (any(shape(image) /= shape(image_))) call failure('ft_read_image shape')
    if (any(image /= image_)) call failure('ft_read_image comparison')

    ! test reading keywords from unit
    call ft_open(filename_header, unit, status)
    if (status /= 0) call failure('open read_keyword')

    call ft_read_keyword(unit, 'donotexist', ivalue, found, status)
    if (status /= 0 .or. found) call failure('unit,donotexist1')

    call ft_read_keyword(unit, 'donotexist', lvalue, found, status)
    if (status /= 0 .or. found) call failure('unit,donotexist2')

    call ft_read_keyword(unit, 'donotexist', dvalue, found, status)
    if (status /= 0 .or. found) call failure('unit,donotexist3')

    call ft_read_keyword(unit, 'donotexist', cvalue, found, status)
    if (status /= 0 .or. found) call failure('unit,donotexist4')

    call ft_read_keyword(unit, 'donotexist', lvalue, found, status)
    if (status /= 0 .or. found) call failure('unit,donotexist5')

    call ft_read_keyword(unit, 'donotexist', ivalue, status=status)
    if (status == 0) call failure('unit,donotexist1b')
    call ft_open(filename_header, unit, status)

    call ft_read_keyword(unit, 'donotexist', lvalue, status=status)
    if (status == 0) call failure('unit,donotexist2b')
    call ft_open(filename_header, unit, status)

    call ft_read_keyword(unit, 'donotexist', dvalue, status=status)
    if (status == 0) call failure('unit,donotexist3b')
    call ft_open(filename_header, unit, status)

    call ft_read_keyword(unit, 'donotexist', cvalue, status=status)
    if (status == 0) call failure('unit,donotexist4b')
    call ft_open(filename_header, unit, status)

    call ft_read_keyword(unit, 'donotexist', lvalue, status=status)
    if (status == 0) call failure('unit,donotexist5b')
    call ft_open(filename_header, unit, status)

    call ft_read_keyword(unit, 'naxis', ivalue, found, status)
    if (status /= 0 .or. ivalue /= 2) call failure('unit,naxis')

    call ft_read_keyword(unit, 'naxis1', lvalue, found, status)
    if (status /= 0 .or. lvalue /= 487) call failure('unit,naxis1')

    call ft_read_keyword(unit, 'naxis2', lvalue, found, status)
    if (status /= 0 .or. lvalue /= 487) call failure('unit,naxis2')

    call ft_read_keyword(unit, 'cdelt1', dvalue, found, status)
    if (status /= 0 .or. dvalue /= -0.00083333333333333_p) call failure('unit,cdelt1')

    call ft_read_keyword(unit, 'cdelt2', dvalue, found, status)
    if (status /= 0 .or. dvalue /= 0.00083333333333333_p) call failure('unit,cdelt2')

    call ft_read_keyword(unit, 'crota2', dvalue, found, status)
    if (status /= 0 .or. dvalue /= 0) call failure('unit,crota2')

    call ft_read_keyword(unit, 'crpix1', dvalue, found, status)
    if (status /= 0 .or. dvalue /= 1) call failure('unit,crpix1')

    call ft_read_keyword(unit, 'crpix2', dvalue, found, status)
    if (status /= 0 .or. dvalue /= 1) call failure('unit,crpix2')

    call ft_read_keyword(unit, 'crval1', dvalue, found, status)
    if (status /= 0 .or. dvalue /= 308.3404249521426_p) call failure('unit,crval1')

    call ft_read_keyword(unit, 'crval2', dvalue, found, status)
    if (status /= 0 .or. dvalue /= 60.357773247484651_p) call failure('unit,crval2')

    call ft_read_keyword(unit, 'ctype1', cvalue, found, status)
    if (status /= 0 .or. cvalue /= 'RA---TAN') call failure('unit,ctype1')

    call ft_read_keyword(unit, 'ctype2', cvalue, found, status)
    if (status /= 0 .or. cvalue /= 'DEC--TAN') call failure('unit,ctype2')

    call ft_read_keyword(unit, 'cunit1', cvalue, found, status)
    if (status /= 0 .or. cvalue /= 'deg') call failure('unit,cunit1')

    call ft_read_keyword(unit, 'cunit2', cvalue, found, status)
    if (status /= 0 .or. cvalue /= 'deg') call failure('unit,cunit2')

    call ft_read_keyword(unit, 'testme1a', bvalue, found, status)
    if (status /= 0 .or. bvalue .neqv. .true.) call failure('unit,testme1a')

    call ft_read_keyword(unit, 'testme1b', bvalue, found, status)
    if (status /= 0 .or. bvalue .neqv. .true.) call failure('unit,testme1b')

    call ft_read_keyword(unit, 'testme1c', bvalue, found, status)
    if (status /= 0 .or. bvalue .neqv. .true.) call failure('unit,testme1c')

    call ft_read_keyword(unit, 'testme2a', bvalue, found, status)
    if (status /= 0 .or. bvalue .neqv. .false.) call failure('unit,testme2a')

    call ft_read_keyword(unit, 'testme2b', bvalue, found, status)
    if (status /= 0 .or. bvalue .neqv. .false.) call failure('unit,testme2b')

    call ft_read_keyword(unit, 'testme2c', bvalue, found, status)
    if (status /= 0 .or. bvalue .neqv. .false.) call failure('unit,testme2c')

    call ft_read_keyword(unit, 'testme4', bvalue, found, status)
    if (status == 0) call failure('unit,testme4')
    call ft_open(filename_header, unit, status)

    call ft_read_keyword(unit, 'testme5', ivalue, found, status)
    if (status == 0) call failure('unit,testme5a')
    call ft_open(filename_header, unit, status)

    call ft_read_keyword(unit, 'testme5', lvalue, found, status)
    if (status == 0) call failure('unit,testme5b')
    call ft_open(filename_header, unit, status)

    call ft_read_keyword(unit, 'testme6', dvalue, found, status)
    if (status == 0) call failure('unit,testme6')
    call ft_open(filename_header, unit, status)

    call ft_read_keyword(unit, 'testme7', cvalue, found, status)
    if (status /= 0 .or. cvalue /= "3''D") call failure('unit,testme7')

    call ft_read_keyword(unit, 'testme8', cvalue, found, status)
    if (status /= 0 .or. cvalue /= "'/") call failure('unit,testme8')

    call ft_read_keyword(unit, 'testme9', cvalue, found, status)
    if (status /= 0 .or. cvalue /= "'/") call failure('unit,testme9')

    call ft_read_keyword(unit, 'testme10', cvalue, found, status)
    if (status /= 0 .or. cvalue /= "'/  / myc10") call failure('unit,testme10')
    
    ! ft_read_keyword_hcss
    call ft_read_keyword_hcss(unit, 'donotexist', bvalue, found, status)
    if (status /= 0 .or. found) call failure('hcss,donotexist1a')

    call ft_read_keyword_hcss(unit, 'donotexist', ivalue, found, status)
    if (status /= 0 .or. found) call failure('hcss,donotexist2a')

    call ft_read_keyword_hcss(unit, 'donotexist', lvalue, found, status)
    if (status /= 0 .or. found) call failure('hcss,donotexist3a')

    call ft_read_keyword_hcss(unit, 'donotexist', dvalue, found, status)
    if (status /= 0 .or. found) call failure('hcss,donotexist4a')

    call ft_read_keyword_hcss(unit, 'donotexist', cvalue, found, status)
    if (status /= 0 .or. found) call failure('hcss,donotexist5a')

    call ft_read_keyword_hcss(unit, 'donotexist', bvalue, status=status)
    if (status == 0) call failure('hcss,donotexist1b')
    call ft_open(filename_header, unit, status)
    
    call ft_read_keyword_hcss(unit, 'donotexist', ivalue, status=status)
    if (status == 0) call failure('hcss,donotexist2b')
    call ft_open(filename_header, unit, status)

    call ft_read_keyword_hcss(unit, 'donotexist', lvalue, status=status)
    if (status == 0) call failure('hcss,donotexist3b')
    call ft_open(filename_header, unit, status)

    call ft_read_keyword_hcss(unit, 'donotexist', dvalue, status=status)
    if (status == 0) call failure('hcss,donotexist4b')
    call ft_open(filename_header, unit, status)

    call ft_read_keyword_hcss(unit, 'donotexist', cvalue, status=status)
    if (status == 0) call failure('hcss,donotexist5b')
    call ft_open(filename_header, unit, status)

    call ft_read_keyword_hcss(unit, 'camBool', bvalue, found, status)
    if (status /= 0 .or. .not. found .or. .not. bvalue) call failure('hcss,camBool')

    call ft_read_keyword_hcss(unit, 'detRow', ivalue, found, status)
    if (status /= 0 .or. .not. found .or. ivalue /= 32) call failure('hcss,detRow')

    call ft_read_keyword_hcss(unit, 'detRow', lvalue, found, status)
    if (status /= 0 .or. .not. found .or. lvalue /= 32) call failure('hcss,detRow')

    call ft_read_keyword_hcss(unit, 'detDouble', dvalue, found, status)
    if (status /= 0 .or. .not. found .or. dvalue /= 64.5) call failure('hcss,detDouble')

    call ft_read_keyword_hcss(unit, 'camName', cvalue, found, status)
    if (status /= 0 .or. .not. found .or. cvalue /= 'Blue Photometer') call failure('hcss,camName')

    call ft_close(unit, status)
    if (status /= 0) call failure('close read_keyword')

contains

    subroutine failure(errmsg)
        character(len=*), intent(in) :: errmsg
        write (ERROR_UNIT,'(a)'), 'FAILED: ' // errmsg
        stop 1
    end subroutine failure

    subroutine cleanup
        open (unit=10, file='/tmp/test_fitstools.fits')
        close (10, status='delete')
    end subroutine cleanup

end program test_fitstools
