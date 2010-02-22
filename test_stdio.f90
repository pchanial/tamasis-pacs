program test_stdio

    use, intrinsic :: ISO_C_BINDING
    use module_stdio
    use module_cfitsio, only : fits_report_error
    implicit none

    integer(kind=C_INT) :: status

    ! opens stdin, stdout, stderr streams and make them available for C programs
    call init_stdio()

    ! test stdout & stderr
    status = 1
    call fits_report_error(stdout, status)
    call fits_report_error(stderr, status)
    status = 301
    call fits_report_error(stdout, status)
    call fits_report_error(stderr, status)

end program test_stdio
