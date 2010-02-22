program test_fitstools

    use module_fitstools
    use module_wcslib, only : wcsfree, WCSLEN
    implicit none

    character(len=*), parameter :: filename_header = '/home/pchanial/work/tamasis/csh_header_ngc6946.fits'

    integer                     :: wcs(WCSLEN)
    integer                     :: status, nx, ny
    character(len=2880)         :: header
    real*8                      :: image(10,15)

    status = 0
    call ft_header2str(filename_header, header, status)
    call ft_printerror(status, filename_header)

    call ft_header2wcs(header, wcs, nx, ny)

    ! test writing an image
    status = 0
    image = 0
    image(1,2) = 1
    call ft_write('/tmp/test_fitstools.fits', image, wcs, status)
    call ft_printerror(status)

    status = wcsfree(wcs)
    stop "OK."

end program test_fitstools
