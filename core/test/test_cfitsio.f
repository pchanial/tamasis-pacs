program test_cfitsio

    use, intrinsic :: ISO_C_BINDING
    use, intrinsic :: ISO_FORTRAN_ENV
    use            :: module_cfitsio
    implicit none

    character(len=*), parameter    :: filename = 'core/test/data/pih.fits'
    character(len=80*1000+1,kind=C_CHAR), pointer :: header
    integer(kind=C_INT)            :: nkeyrec, status
    type(C_PTR)                    :: fptr, c_header
    integer                        :: i

    status = 0
    call fits_open_file(fptr, filename // C_NULL_CHAR, CFITSIO_READONLY, status)
    if (status /= 0) then
        call cfitsio_report_error(status)
        call failure('fits_open_file')
    end if

    call fits_hdr2str(fptr, 1, C_NULL_PTR, 0, c_header, nkeyrec, status)
    if (status /= 0) then
        call cfitsio_report_error(status)
        call failure('fits_hdr2str')
    end if

    call c_f_pointer(c_header, header)
    do i=1, 5
       write(*,'(i3,a,a)') i, ' : ', trim(header((i-1)*80+1:i*80))
    end do

    call fits_close_file(fptr, status)
    if (status /= 0) then
        call cfitsio_report_error(status)
        call failure('fits_close_file')
    end if

contains

    subroutine failure(errmsg)
        character(len=*), intent(in) :: errmsg
        write (ERROR_UNIT,'(a)'), 'FAILED: ' // errmsg
        stop 1
    end subroutine failure

end program test_cfitsio
