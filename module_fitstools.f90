module module_fitstools

    use, intrinsic :: ISO_C_BINDING
    use, intrinsic :: ISO_FORTRAN_ENV
    use module_cfitsio, only: CFITSIO_READONLY, CFITSIO_READWRITE
    use module_wcslib, only: WCSLEN, wcshdo
    implicit none

    private
    integer, private, parameter :: GROUP = 1
    integer, private, parameter :: NULLVAL = 0
    integer, private, parameter :: BUFFERSIZE = 1024

    public :: ft_openimage
    public :: ft_close
    public :: ft_readextension
    public :: ft_readslice
    public :: ft_create_header
    public :: ft_write
    public :: ft_header2str
    public :: ft_header2wcs
    public :: ft_printerror

    interface ft_readextension
        module procedure readext_logical_1d, readext_int64_1d, readext_double_1d, &
                         readext_logical_2d, readext_double_2d,                   &
                         readext_logical_3d, readext_double_3d
    end interface ft_readextension

    interface ft_readslice
        module procedure readslice_logical_filename, readslice_logical_unit,      &
                         readslice_int64_filename,   readslice_int64_unit,        &
                         readslice_double_filename,  readslice_double_unit,       &
                         readslice_logical_3d, readslice_double_3d
    end interface ft_readslice

    interface ft_write
        module procedure writefits_double_1d, writefits_double_2d
    end interface ft_write


contains


    recursive subroutine readext_logical_1d(filename, output, status, hdu)

        character(len=*), intent(in)        :: filename
        logical*1, allocatable, intent(out) :: output(:)
        integer, intent(inout)              :: status
        integer, optional, intent(in)       :: hdu

        integer                             :: unit, firstpix, anynull
        integer, allocatable                :: imageshape(:)

        if (status /= 0) return

        call ft_openimage(filename, unit, 1, imageshape, status, hdu=hdu)
        if (status /= 0) return

        !  Initialize variables
        firstpix=1
        allocate(output(imageshape(1)))

        call ftgpvb(unit, GROUP, firstpix, imageshape(1), NULLVAL, output, anynull, status)
        if (status /= 0) return

        call ft_close(unit, status)

    end subroutine readext_logical_1d


    !-------------------------------------------------------------------------------


    recursive subroutine readext_int64_1d(filename, output, status, hdu)

        character(len=*), intent(in)        :: filename
        integer*8, allocatable, intent(out) :: output(:)
        integer, intent(inout)              :: status
        integer, optional, intent(in)       :: hdu

        integer                             :: unit, firstpix, anynull
        integer, allocatable                :: imageshape(:)

        if (status /= 0) return

        call ft_openimage(filename, unit, 1, imageshape, status, hdu=hdu)
        if (status /= 0) return

        !  Initialize variables
        firstpix=1
        allocate(output(imageshape(1)))

        call ftgpvk(unit, GROUP, firstpix, imageshape(1), NULLVAL, output, anynull, status)
        if (status /= 0) return

        call ft_close(unit, status)

    end subroutine readext_int64_1d


    !-------------------------------------------------------------------------------


    recursive subroutine readext_double_1d(filename, output, status, hdu)

        character(len=*), intent(in)     :: filename
        real*8, allocatable, intent(out) :: output(:)
        integer, intent(inout)           :: status
        integer, optional, intent(in)    :: hdu

        integer                          :: unit, firstpix, anynull
        integer, allocatable             :: imageshape(:)

        if (status /= 0) return

        call ft_openimage(filename, unit, 1, imageshape, status, hdu=hdu)
        if (status /= 0) return

        !  Initialize variables
        firstpix=1
        allocate(output(imageshape(1)))

        call ftgpvd(unit, GROUP, firstpix, imageshape(1), NULLVAL, output, anynull, status)
        if (status /= 0) return

        call ft_close(unit, status)

    end subroutine readext_double_1d


    !-------------------------------------------------------------------------------


    recursive subroutine readext_logical_2d(filename, output, status, hdu)

        character(len=*), intent(in)        :: filename
        integer, intent(inout)              :: status
        logical*1, allocatable, intent(out) :: output(:,:)
        integer, optional, intent(in)       :: hdu

        integer                             :: unit, BUFFERSIZE, firstpix, anynull, j
        integer, allocatable                :: imageshape(:)

        if (status /= 0) return

        call ft_openimage(filename, unit, 2, imageshape, status, hdu=hdu)
        if (status /= 0) return

        ! Initialize variables
        firstpix=1
        BUFFERSIZE = imageshape(1)
        allocate(output(imageshape(1), imageshape(2)))

        do j = 1, imageshape(2)

           call ftgpvb(unit, GROUP, firstpix, BUFFERSIZE, NULLVAL, output(:, j), anynull, status)
           if (status /= 0) exit

           firstpix = firstpix + BUFFERSIZE

        end do

        call ft_close(unit, status)

    end subroutine readext_logical_2d


    !-------------------------------------------------------------------------------


    recursive subroutine readext_double_2d(filename, output, status, hdu)
    
        character(len=*), intent(in)     :: filename
        real*8, allocatable, intent(out) :: output(:,:)
        integer, intent(inout)           :: status
        integer, optional, intent(in)    :: hdu

        integer                          :: unit, BUFFERSIZE, firstpix, anynull, j
        integer, allocatable             :: imageshape(:)

        if (status /= 0) return
    
        call ft_openimage(filename, unit, 2, imageshape, status, hdu=hdu)
        if (status /= 0) return

        !  Initialize variables
        firstpix=1
        BUFFERSIZE = imageshape(1)
        allocate(output(imageshape(1),imageshape(2)))

        do j = 1, imageshape(2)
    
            call ftgpvd(unit, GROUP, firstpix, BUFFERSIZE, NULLVAL, output(:,j), anynull, status)
            if (status /= 0) exit

            firstpix = firstpix + BUFFERSIZE

        end do

        call ft_close(unit, status)

    end subroutine readext_double_2d


    !-------------------------------------------------------------------------------


    recursive subroutine readext_logical_3d(filename, output, status)

        character(len=*), intent(in)        :: filename
        logical*1, allocatable, intent(out) :: output(:,:,:)
        integer, intent(inout)              :: status

        integer                             :: unit, BUFFERSIZE, firstpix, anynull, j, k
        integer, allocatable                :: imageshape(:)

        if (status /= 0) return

        call ft_openimage(filename, unit, 3, imageshape, status)
        if (status /= 0) return

        !  Initialize variables
        firstpix = 1
        BUFFERSIZE = imageshape(1)
        allocate(output(imageshape(1), imageshape(2), imageshape(3)))

        outer: do k = 1, imageshape(3)
            do j = 1, imageshape(2)

                call ftgpvb(unit, GROUP, firstpix, BUFFERSIZE, NULLVAL, output(:,j,k), anynull, status)
                if (status /= 0) exit outer

                firstpix = firstpix + BUFFERSIZE

            end do
        end do outer

        call ft_close(unit, status)

    end subroutine readext_logical_3d


    !-------------------------------------------------------------------------------


    recursive subroutine readext_double_3d(filename, output, status)

        character(len=*), intent(in)     :: filename
        real*8, allocatable, intent(out) :: output(:,:,:)
        integer, intent(inout)           :: status

        integer                          :: unit, BUFFERSIZE, firstpix, anynull, j, k
        integer, allocatable             :: imageshape(:)

        if (status /= 0) return

        call ft_openimage(filename, unit, 3, imageshape, status)
        if (status /= 0) return

        !  Initialize variables
        firstpix = 1
        BUFFERSIZE = imageshape(1)
        allocate(output(imageshape(1), imageshape(2), imageshape(3)))

        outer: do k = 1, imageshape(3)

            do j = 1, imageshape(2)

                call ftgpvd(unit, GROUP, firstpix, BUFFERSIZE, NULLVAL, output(:,j,k), anynull, status)
                if (status /= 0) exit outer

                firstpix = firstpix + BUFFERSIZE

            end do

        end do outer

        call ft_close(unit, status)

    end subroutine readext_double_3d


    !-------------------------------------------------------------------------------


    recursive subroutine readslice_logical_filename(filename, x1, x2, array, status)

        character(len=*), intent(in) :: filename
        integer*8, intent(in)        :: x1, x2
        logical*1, intent(inout)     :: array(x2 - x1 + 1)
        integer, intent(inout)       :: status

        integer                      :: unit

        if (status /= 0) return

        call readext_open(filename, unit, status)
        call readslice_logical_unit(unit, x1, x2, array, status)
        call ft_close(unit, status)

    end subroutine readslice_logical_filename


    !-------------------------------------------------------------------------------


    recursive subroutine readslice_logical_unit(unit, x1, x2, array, status)

        integer, intent(in)      :: unit
        integer*8, intent(in)    :: x1, x2
        logical*1, intent(inout) :: array(x2 - x1 + 1)
        integer, intent(inout)   :: status

        integer                  :: anynull

        if (status /= 0) return

        call ftgpvd(unit, GROUP, x1, x2 - x1 + 1, NULLVAL, array, anynull, status)

    end subroutine readslice_logical_unit


    !-------------------------------------------------------------------------------


    recursive subroutine readslice_int64_filename(filename, x1, x2, array, status)

        character(len=*), intent(in) :: filename
        integer*8, intent(in)        :: x1, x2
        integer*8, intent(inout)     :: array(x2 - x1 + 1)
        integer, intent(inout)       :: status

        integer                      :: unit

        if (status /= 0) return

        call readext_open(filename, unit, status)
        call readslice_int64_unit(unit, x1, x2, array, status)
        call ft_close(unit, status)

    end subroutine readslice_int64_filename


    !-------------------------------------------------------------------------------


    recursive subroutine readslice_int64_unit(unit, x1, x2, array, status)

        integer, intent(in)      :: unit
        integer*8, intent(in)    :: x1, x2
        integer*8, intent(inout) :: array(x2 - x1 + 1)
        integer, intent(inout)   :: status

        integer                  :: anynull

        if (status /= 0) return

        call ftgpvk(unit, GROUP, x1, x2 - x1 + 1, NULLVAL, array, anynull, status)

    end subroutine readslice_int64_unit


    !-------------------------------------------------------------------------------


    recursive subroutine readslice_double_filename(filename, x1, x2, array, status)

        character(len=*), intent(in) :: filename
        integer*8, intent(in)        :: x1, x2
        real*8, intent(inout)        :: array(x2 - x1 + 1)
        integer, intent(inout)       :: status

        integer                      :: unit

        if (status /= 0) return

        call readext_open(filename, unit, status)
        call readslice_double_unit(unit, x1, x2, array, status)
        call ft_close(unit, status)

    end subroutine readslice_double_filename


    !-------------------------------------------------------------------------------


    recursive subroutine readslice_double_unit(unit, x1, x2, array, status)

        integer, intent(in)    :: unit
        integer*8, intent(in)  :: x1, x2
        real*8, intent(inout)  :: array(x2 - x1 + 1)
        integer, intent(inout) :: status

        integer                :: anynull

        if (status /= 0) return

        call ftgpvd(unit, GROUP, x1, x2 - x1 + 1, NULLVAL, array, anynull, status)

    end subroutine readslice_double_unit


    !-------------------------------------------------------------------------------


    recursive subroutine readslice_logical_3d(unit, x1, x2, y, z, naxes, array, status)

        integer, intent(in)      :: unit
        integer*8, intent(in)    :: x1, x2
        integer*8, intent(in)    :: y, z
        integer, intent(in)      :: naxes(3)
        logical*1, intent(inout) :: array(x2 - x1 + 1)
        integer, intent(inout)   :: status

        integer                  :: anynull
        integer*8                :: firstpix

        if (status /= 0) return

        firstpix = x1 + naxes(1) * (y - 1 + naxes(2) * (z - 1))
        call ftgpvb(unit, GROUP, firstpix, x2 - x1 + 1, NULLVAL, array, anynull, status)

    end subroutine readslice_logical_3d


    !-------------------------------------------------------------------------------


    recursive subroutine readslice_double_3d(unit, x1, x2, y, z, naxes, array, status)

        integer, intent(in)    :: unit
        integer*8, intent(in)  :: x1, x2
        integer*8, intent(in)  :: y, z
        integer, intent(in)    :: naxes(3)
        real*8, intent(inout)  :: array(x2 - x1 + 1)
        integer, intent(inout) :: status

        integer                :: anynull
        integer*8              :: firstpix

        if (status /= 0) return

        firstpix = x1 + naxes(1) * (y - 1 + naxes(2) * (z - 1))
        call ftgpvd(unit, GROUP, firstpix, x2 - x1 + 1, NULLVAL, array, anynull, status)

    end subroutine readslice_double_3d


    !-------------------------------------------------------------------------------


    recursive subroutine readext_open(filename, unit, status)

        character(len=*), intent(in)  :: filename
        integer, intent(out)          :: unit
        integer, intent(inout)        :: status

        if (status /= 0) return
        call ftgiou(unit, status)

        ! open the fits file and move to the HDU specified in the filename
        call ftnopn(unit, filename, CFITSIO_READONLY, status)

    end subroutine readext_open


    !-------------------------------------------------------------------------------


    recursive subroutine ft_openimage(filename, unit, imagerank, imageshape, status, hdu)

        character(len=*), intent(in)        :: filename
        integer, intent(out)                :: unit
        integer, intent(in)                 :: imagerank
        integer, allocatable, intent(inout) :: imageshape(:)
        integer, intent(inout)              :: status
        integer, intent(in), optional       :: hdu

        integer                             :: nfound, hdutype, blocksize
        character(len=80)                   :: errmsg
        integer                             :: imageshape_(8)

        if (status /= 0) return

        call ftgiou(unit, status)

        if (present(hdu)) then

            ! open the fits file and point to the primary header
            call ftopen(unit, filename, CFITSIO_READONLY, blocksize, status)

            ! move to the specified HDU
            call ftmahd(unit, hdu, hdutype, status)
            if (status == 0 .and. hdutype /= 0) then
                write (errmsg,'(a,i5,a)') "HDU type is not an image: ", hdutype, "."
                call ftpmsg(errmsg)
                status = 1
            endif

        else

            ! open the fits file and move to the HDU specified in the filename
            ! otherwise, move to the first image in FITS file
            call ftiopn(unit, filename, CFITSIO_READONLY, status)

        end if

        if (status /= 0) return

        !  Determine the size of the image.
        call ftgknj(unit, 'NAXIS', 1, imagerank, imageshape_, nfound, status)

        !  Check that it found the NAXISn keywords.
        if (nfound /= imagerank) then
            call ftpmsg('Failed to read the NAXISn keywords.')
            status = 1
            return
        end if

        allocate(imageshape(imagerank))
        imageshape = imageshape_(1:imagerank)

    end subroutine ft_openimage


    !-------------------------------------------------------------------------------


    recursive subroutine ft_close(unit, status)

        integer, intent(in)    :: unit
        integer, intent(inout) :: status

        !  The FITS file must always be closed before exiting the program.
        !  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
        call ftclos(unit, status)
        call ftfiou(unit, status)
        call ft_printerror(status)

    end subroutine ft_close


    !-------------------------------------------------------------------------------


    recursive subroutine ft_printerror(status, filename, hdu)

        !  This subroutine prints out the descriptive text corresponding to the
        !  error status value and prints out the contents of the internal
        !  error message stack generated by FITSIO whenever an error occurs.

        integer, intent(in)                    :: status
        character(len=*), intent(in), optional :: filename
        integer, intent(in), optional          :: hdu
        character                              :: errtext*30, errmessage*80

        !  Check if status is OK (no error); if so, simply return
        if (status == 0) return

        ! flush the output unit
        flush(OUTPUT_UNIT)

        !  Print the filename
        if (present(filename)) write (ERROR_UNIT,*) 'In file ' // filename // ': '

        !  Print the HDU number
        if (present(hdu)) write (ERROR_UNIT,*) 'HDU: ', hdu

        !  The FTGERR subroutine returns a descriptive 30-character text string that
        !  corresponds to the integer error status number.  A complete list of all
        !  the error numbers can be found in the back of the FITSIO User's Guide.
        call ftgerr(status, errtext)
        write(ERROR_UNIT,*) 'FITSIO Error Status =',status,': ',errtext

        !  FITSIO usually generates an internal stack of error messages whenever
        !  an error occurs.  These messages provide much more information on the
        !  cause of the problem than can be provided by the single integer error
        !  status value.  The FTGMSG subroutine retrieves the oldest message from
        !  the stack and shifts any remaining messages on the stack down one
        !  position.  FTGMSG is called repeatedly until a blank message is
        !  returned, which indicates that the stack is empty.  Each error message
        !  may be up to 80 characters in length.  Another subroutine, called
        !  FTCMSG, is available to simply clear the whole error message stack in
        !  cases where one is not interested in the contents.
        call ftgmsg(errmessage)
        do while (errmessage /= ' ')
            write(ERROR_UNIT,*) errmessage
            call ftgmsg(errmessage)
        end do

    end subroutine ft_printerror


    !-------------------------------------------------------------------------------


    subroutine ft_create_header(naxis1, naxis2, cdelt1, cdelt2, crota2, crval1, crval2, crpix1, crpix2, header)

        integer, intent(in) :: naxis1, naxis2
        real*8, intent(in)  :: cdelt1, cdelt2, crval1, crval2, crpix1, crpix2, crota2
        character(len=2880), intent(out) :: header
        character(len=*), parameter   :: format_dbl = "(a8,'= ',f20.15)"
        character(len=*), parameter   :: format_int = "(a8,'= ',i20)"
        character(len=*), parameter   :: format_str = "(a8,'= ',a20)"

        header = ' '
        header(0*80+1:1*80) = 'SIMPLE  =                    T / Fits standard'
        header(1*80+1:2*80) = 'BITPIX  =                  -64 / Bits per pixel'
        header(2*80+1:3*80) = 'NAXIS   =                    2'
        write (header( 3*80+1: 4*80),format_int) 'NAXIS1  ', naxis1
        write (header( 4*80+1: 5*80),format_int) 'NAXIS2  ', naxis2
        write (header( 5*80+1: 6*80),format_dbl) 'CDELT1  ', cdelt1
        write (header( 6*80+1: 7*80),format_dbl) 'CDELT2  ', cdelt2
        write (header( 7*80+1: 8*80),format_dbl) 'CROTA2  ', crota2
        write (header( 8*80+1: 9*80),format_dbl) 'CRVAL1  ', crval1
        write (header( 9*80+1:10*80),format_dbl) 'CRVAL2  ', crval2
        write (header(10*80+1:11*80),format_dbl) 'CRPIX1  ', crpix1
        write (header(11*80+1:12*80),format_dbl) 'CRPIX2  ', crpix2
        write (header(12*80+1:13*80),format_str) 'CTYPE1  ', "'RA---TAN'"
        write (header(13*80+1:14*80),format_str) 'CTYPE2  ', "'DEC--TAN'"
        write (header(14*80+1:15*80),format_str) 'CUNIT1  ', "'deg'"
        write (header(15*80+1:16*80),format_str) 'CUNIT2  ', "'deg'"
        header(16*80+1:16*80+3) = 'END'

    end subroutine ft_create_header


    !-------------------------------------------------------------------------------


    subroutine writefits_double_1d(filename, data, wcs, status)

        character(len=*), intent(in)    :: filename
        real*8, intent(in)              :: data(:)
        integer, intent(in), optional   :: wcs(WCSLEN)
        integer, intent(inout)          :: status
        integer                         :: nkeyrec, irec, unit, blocksize, bitpix, naxis, naxes(1)
        logical                         :: simple, extend
        type(C_PTR)                     :: c_header
        character(len=1000*80), pointer :: header => null()

        if (status /= 0) return

        ! delete file if it exists
        open(10, file=filename, status='unknown')
        close(10, status='delete')

        ! open and initialise the fits file
        call ftgiou(unit, status)
        call ftinit(unit, filename, blocksize, status)
        simple = .true.
        bitpix = -64
        naxis  = 1
        naxes  = [ size(data) ]
        extend = .true.
        call FTPHPR(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

        ! write the astrometry keywords
        if (present(wcs)) then
            status = wcshdo(0, wcs, nkeyrec, c_header)
            call c_f_pointer(c_header, header)

            do irec=1, nkeyrec
                call FTPREC(unit,header((irec-1)*80+1:irec*80), status)
            end do
        end if

        ! write the image data
        call FTPPRD(unit, GROUP, 1, size(data), data, status)

        ! close the fits file
        call ft_close(unit, status)

    end subroutine writefits_double_1d


    !-------------------------------------------------------------------------------


    subroutine writefits_double_2d(filename, data, wcs, status)

        character(len=*), intent(in)    :: filename
        real*8, intent(in)              :: data(:,:)
        integer, intent(in), optional   :: wcs(WCSLEN)
        integer, intent(inout)          :: status
        integer                         :: nkeyrec, irec, unit, blocksize, bitpix, naxis, naxes(2)
        logical                         :: simple, extend
        type(C_PTR)                     :: c_header
        character(len=1000*80), pointer :: header => null()

        if (status /= 0) return

        ! delete file if it exists
        open(10, file=filename, status='unknown')
        close(10, status='delete')

        ! open and initialise the fits file
        call ftgiou(unit, status)
        call ftinit(unit, filename, blocksize, status)
        simple = .true.
        bitpix = -64
        naxis  = 2
        naxes  = [ size(data,1), size(data,2) ]
        extend = .true.
        call FTPHPR(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

        ! write the astrometry keywords
        if (present(wcs)) then
            status = wcshdo(0, wcs, nkeyrec, c_header)
            call c_f_pointer(c_header, header)

            do irec=1, nkeyrec
                call FTPREC(unit,header((irec-1)*80+1:irec*80), status)
            end do
        end if

        ! write the image data
        call FTPPRD(unit, GROUP, 1, size(data), data, status)

        ! close the fits file
        call ft_close(unit, status)

    end subroutine writefits_double_2d


    !-------------------------------------------------------------------------------


    ! translate the astrometry in a FITS header into a wcslib wcs
    ! XXX this routine currently leaks wcsp
    ! investigate how to keep wcs without wcsp
    subroutine ft_header2str(filename, header, status)

        use, intrinsic :: ISO_C_BINDING
        use module_cfitsio, only : CFITSIO_READONLY, fits_open_file, fits_hdr2str, fits_close_file, fits_report_error
        character(len=*), intent(in)     :: filename
        character(len=2880), intent(out) :: header
        integer, intent(inout)           :: status

        character(len=2880), pointer     :: f_header
        type(C_PTR)                      :: fptr, c_header
        integer                          :: nkeyrec

        if (status /= 0) return

        call fits_open_file(fptr, filename // C_NULL_CHAR, CFITSIO_READONLY, status)
        if (status /= 0) return

        call fits_hdr2str(fptr, 1, C_NULL_PTR, 0, c_header, nkeyrec, status)
        if (status /= 0) return

        call fits_close_file(fptr, status)
        if (status /= 0) return

        call c_f_pointer(c_header, f_header)
        header = f_header(1:nkeyrec*80)
        header(nkeyrec*80+1:) = ' '

    end subroutine ft_header2str


    !-------------------------------------------------------------------------------


    subroutine ft_header2wcs(header, wcs, nx, ny)

        use, intrinsic :: ISO_C_BINDING
        use module_wcslib
        character(len=*), intent(in)     :: header
        integer(kind=C_INT), intent(out) :: wcs(WCSLEN)
        integer, intent(out)             :: nx, ny
        integer(kind=C_INT)              :: wcs_(WCSLEN)

        integer(kind=C_INT) :: wcsp, status
        integer                          :: nreject, nwcs, statfix(WCSFIX_NWCS), iline
        character(len=80)                :: record

        status = wcspih(header // C_NULL_CHAR, len(header)/80, WCSHDR_all, 0, nreject, nwcs, wcsp)
        if (status /= 0) then
            write (*,*) 'WCSPIH error: ', status
            stop
        end if

        if (nwcs == 0) stop "Header has no astrometry."

        if (wcsp < 0) stop "wcsp negative"
        status = wcsvcopy(wcsp, 0, wcs_)
        status = wcsput(wcs, WCS_FLAG, -1_C_INT, 0_C_INT, 0_C_INT)
        status = wcscopy(wcs_, wcs)
        status = wcsfree(wcs_)
        status = wcsvfree(nwcs, wcsp)

        !Fix non-standard WCS keyvalues.
        status = wcsfix(ctrl=7, naxis=c_null_ptr, wcs=wcs, stat=statfix)
        if (status /= 0) then
            write (*,*) 'WCSFIX error: ', status
            stop
        end if

        status = wcsset(wcs)

        ! extract NAXIS1 and NAXIS2 keywords
        nx = -1
        ny = -1
        do iline = 1, len(header) / 80
           record = header((iline-1)*80+1:iline*80)
           if (record(1:8) == 'NAXIS1  ') then
              read(record, '(10x,i20)') nx
           end if
           if (record(1:8) == 'NAXIS2  ') then
              read(record, '(10x,i20)') ny
           end if
        end do
        if (nx < 0 .or. ny < 0) then
            write (*,'(a)') 'START header:'
            do iline = 1, len(header) / 80
                record = header((iline-1)*80+1:iline*80)
                write (*,'(a,a)') '>', record
            end do
            write (*,'(a)') 'END header.'
            stop "Invalid or missing NAXISn keywords in header."
        end if

    end subroutine ft_header2wcs


end module module_fitstools
