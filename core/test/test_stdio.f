program test_stdio

    use iso_c_binding
    use iso_fortran_env, only : OUTPUT_UNIT
    use module_stdio
    use module_cfitsio, only : fits_report_error
    implicit none

    integer(kind=C_INT) :: status

    ! opens stdin, stdout, stderr streams and make them available for C programs
    call init_stdio()

    ! test stdout & stderr
    status = 1
    write (OUTPUT_UNIT,'(a)') 'Testing stdout/stderr...'
    call fits_report_error(stdout, status)
    call fits_report_error(stderr, status)
    write (OUTPUT_UNIT,'(a)') 'OK.'

end program test_stdio
