program test_cfitsio

    use, intrinsic :: ISO_C_BINDING
    use module_stdio
    use module_cfitsio
    implicit none

    character(len=*), parameter    :: filename = '/home/pchanial/software/wcslib-4.4.4/C/pih.fits'
    character(len=80*1000+1,kind=C_CHAR), pointer :: header
    integer(kind=C_INT)            :: nkeyrec, status
    type(C_PTR)                    :: fptr, c_header
    integer                        :: i

    call init_stdio()
    status = 0

    call fits_open_file(fptr, filename // C_NULL_CHAR, CFITSIO_READONLY, status)
    if (status /= 0) then
        call fits_report_error(stderr, status)
        stop
    end if

    call fits_hdr2str(fptr, 1, C_NULL_PTR, 0, c_header, nkeyrec, status)
    if (status /= 0) then
        call fits_report_error(stderr, status)
        stop
    end if

    call c_f_pointer(c_header, header)
    do i=1, nkeyrec
       write(*,'(i3,a,a)') i, ' : ', header((i-1)*80+1:i*80)
    end do

    call fits_close_file(fptr, status)
    if (status /= 0) then
        call fits_report_error(stderr, status)
        stop
    end if

    write (*,*) 'Test CFITSIO: done.'

end program test_cfitsio
