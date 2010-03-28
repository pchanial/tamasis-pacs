module module_fitstools

    use ISO_C_BINDING
    use ISO_FORTRAN_ENV, only : ERROR_UNIT, OUTPUT_UNIT
    use module_cfitsio
    use precision,       only : dp
    use string,          only : strlowcase, strsection, strupcase
    implicit none
    private

    public :: CFITSIO_IMAGE_HDU
    public :: CFITSIO_ASCII_TBL
    public :: CFITSIO_BINARY_TBL
    public :: CFITSIO_ANY_HDU
    public :: CFITSIO_READONLY
    public :: CFITSIO_READWRITE
    public :: ft_close
    public :: ft_create_header
    public :: ft_header2str
    public :: ft_open
    public :: ft_open_bintable
    public :: ft_open_image
    public :: ft_read_column
    public :: ft_readextension
    public :: ft_readparam
    public :: ft_readslice
    public :: ft_test_extension
    public :: ft_write

    integer, private, parameter :: GROUP = 1
    integer, private, parameter :: NULLVAL = 0
    integer, private, parameter :: BUFFERSIZE = 1024

    interface ft_read_column
        module procedure ft_read_column_character_filename,                    &
                         ft_read_column_character_unit,                        &
                         ft_read_column_int4_filename,                         &
                         ft_read_column_int4_unit,                             &
                         ft_read_column_int8_filename,                         &
                         ft_read_column_int8_unit,                             &
                         ft_read_column_double_filename,                       &
                         ft_read_column_double_unit
    end interface ft_read_column

    interface ft_readextension
        module procedure readext_logical_1d, readext_int8_1d,                  &
                         readext_double_1d,                                    &
                         readext_logical_2d, readext_double_2d,                &
                         readext_logical_3d, readext_double_3d
    end interface ft_readextension

    interface ft_readslice
        module procedure readslice_logical_filename, readslice_logical_unit,   &
                         readslice_int4_filename,   readslice_int4_unit,       &
                         readslice_int8_filename,   readslice_int8_unit,       &
                         readslice_double_filename,  readslice_double_unit,    &
                         readslice_logical_3d, readslice_int4_3d,              &
                         readslice_int8_3d, readslice_double_3d
    end interface ft_readslice

    interface ft_write
        module procedure writefits_double_1d, writefits_double_2d, writefits_double_3d
    end interface ft_write

    interface ft_readparam
        module procedure ft_readparam_logical, ft_readparam_int4, ft_readparam_int8, ft_readparam_double, ft_readparam_character
    end interface ft_readparam


contains


    subroutine ft_read_column_character_filename(filename, colname, data,status)
        character(len=*), intent(in)               :: filename
        character(len=*), intent(in)               :: colname
        character(len=*), allocatable, intent(out) :: data(:)
        integer, intent(out)                       :: status
        integer :: status_close, nrecords, unit
        
        call ft_open_bintable(filename, unit, nrecords, status)
        if (status /= 0) return
        
        allocate(data(nrecords))
        call ft_read_column_character_unit(unit, colname, 1, nrecords, data,   &
                                           status)
        
        call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine ft_read_column_character_filename


    !---------------------------------------------------------------------------


    subroutine ft_read_column_character_unit(unit, colname, first, last, data, &
                                             status)
        integer, intent(in)           :: unit
        character(len=*), intent(in)  :: colname
        integer, intent(in)           :: first, last
        character(len=*), intent(out) :: data(:)
        integer, intent(out)          :: status
        logical                       :: anyf   ! set to true if undef values
        integer                       :: colnum

        ! parameter check
        if (size(data) /= last-first+1) then
            status = 1
            write (ERROR_UNIT,'(a,i0,a)') "FT_READ_COLUMN: Output array has a s&
                &ize '", size(data), "' incompatible with the specified range '&
                &[" // strsection(first, last) // "]'."
            return
        end if

        status = 0

        ! get column number
        call ftgcno(unit, .true., colname, colnum, status)
        if (ft_checkerror_cfitsio(status)) return

        ! extract column
        call ftgcvs(unit, colnum, 1, first, last, '', data, anyf, status)
        if (ft_checkerror_cfitsio(status)) return

    end subroutine ft_read_column_character_unit


    !---------------------------------------------------------------------------


    subroutine ft_read_column_int4_filename(filename, colname, data, status)
        character(len=*), intent(in)        :: filename
        character(len=*), intent(in)        :: colname
        integer*4, allocatable, intent(out) :: data(:)
        integer, intent(out)                :: status
        integer                             :: status_close, nrecords, unit
        
        call ft_open_bintable(filename, unit, nrecords, status)
        if (status /= 0) return
        
        allocate(data(nrecords))
        call ft_read_column_int4_unit(unit, colname, 1, nrecords, data, status)
        
        call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine ft_read_column_int4_filename


    !---------------------------------------------------------------------------


    subroutine ft_read_column_int4_unit(unit, colname, first, last, data,status)
        integer, intent(in)          :: unit
        character(len=*), intent(in) :: colname
        integer, intent(in)          :: first, last
        integer*4, intent(out)       :: data(:)
        integer, intent(out)         :: status
        logical                      :: anyf   ! set to true if undef values
        integer                      :: colnum

        ! parameter check
        if (size(data) /= last-first+1) then
            status = 1
            write (ERROR_UNIT,'(a,i0,a)') "FT_READ_COLUMN: Output array has a s&
                &ize '", size(data), "' incompatible with the specified range '&
                &[" // strsection(first, last) // "]'."
            return
        end if

        status = 0

        ! get column number
        call ftgcno(unit, .true., colname, colnum, status)
        if (ft_checkerror_cfitsio(status)) return

        ! extract column
        call ftgcvj(unit, colnum, 1, first, last, nullval, data, anyf, status)
        if (ft_checkerror_cfitsio(status)) return

    end subroutine ft_read_column_int4_unit


    !---------------------------------------------------------------------------


    subroutine ft_read_column_int8_filename(filename, colname, data, status)
        character(len=*), intent(in)        :: filename
        character(len=*), intent(in)        :: colname
        integer*8, allocatable, intent(out) :: data(:)
        integer, intent(out)                :: status
        integer                             :: status_close, nrecords, unit
        
        call ft_open_bintable(filename, unit, nrecords, status)
        if (status /= 0) return
        
        allocate(data(nrecords))
        call ft_read_column_int8_unit(unit, colname, 1, nrecords, data, status)
        
        call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine ft_read_column_int8_filename


    !---------------------------------------------------------------------------


    subroutine ft_read_column_int8_unit(unit, colname, first, last, data,status)
        integer, intent(in)          :: unit
        character(len=*), intent(in) :: colname
        integer, intent(in)          :: first, last
        integer*8, intent(out)       :: data(:)
        integer, intent(out)         :: status
        logical                      :: anyf   ! set to true if undef values
        integer                      :: colnum

        ! parameter check
        if (size(data) /= last-first+1) then
            status = 1
            write (ERROR_UNIT,'(a,i0,a)') "FT_READ_COLUMN: Output array has a s&
                &ize '", size(data), "' incompatible with the specified range '&
                &[" // strsection(first, last) // "]'."
            return
        end if

        status = 0

        ! get column number
        call ftgcno(unit, .true., colname, colnum, status)
        if (ft_checkerror_cfitsio(status)) return

        ! extract column
        call ftgcvk(unit, colnum, 1, first, last, nullval, data, anyf, status)
        if (ft_checkerror_cfitsio(status)) return

    end subroutine ft_read_column_int8_unit


    !---------------------------------------------------------------------------


    subroutine ft_read_column_double_filename(filename, colname, data, status)
        character(len=*), intent(in)     :: filename
        character(len=*), intent(in)     :: colname
        real*8, allocatable, intent(out) :: data(:)
        integer, intent(out)             :: status
        integer                          :: status_close, nrecords, unit
        
        call ft_open_bintable(filename, unit, nrecords, status)
        if (status /= 0) return
        
        allocate(data(nrecords))
        call ft_read_column_double_unit(unit, colname, 1, nrecords, data,status)
        
        call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine ft_read_column_double_filename


    !---------------------------------------------------------------------------


    subroutine ft_read_column_double_unit(unit, colname, first, last, data,    &
                                          status)
        integer, intent(in)          :: unit
        character(len=*), intent(in) :: colname
        integer, intent(in)          :: first, last
        real(kind=dp), intent(out)   :: data(:)
        integer, intent(out)         :: status
        logical                      :: anyf   ! set to true if undef values
        integer                      :: colnum

        ! parameter check
        if (size(data) /= last-first+1) then
            status = 1
            write (ERROR_UNIT,'(a,i0,a)') "FT_READ_COLUMN: Output array has a s&
                &ize '", size(data), "' incompatible with the specified range '&
                &[" // strsection(first, last) // "]'."
            return
        end if

        status = 0

        ! get column number
        call ftgcno(unit, .true., colname, colnum, status)
        if (ft_checkerror_cfitsio(status)) return

        ! extract column
        call ftgcvd(unit, colnum, 1, first, last, nullval, data, anyf, status)
        if (ft_checkerror_cfitsio(status)) return

    end subroutine ft_read_column_double_unit


    !---------------------------------------------------------------------------


    subroutine readext_logical_1d(filename, output, status, hdu)

        character(len=*), intent(in)        :: filename
        logical*1, allocatable, intent(out) :: output(:)
        integer, intent(out)                :: status
        integer, optional, intent(in)       :: hdu

        integer                             :: unit, anynull, status_close
        integer, allocatable                :: imageshape(:)

        call ft_open_image(filename, unit, 1, imageshape, status, hdu=hdu)
        if (status /= 0) return

        !  Initialize variables
        allocate(output(imageshape(1)))

        call ftgpvb(unit, GROUP, 1, imageshape(1), NULLVAL, output, anynull,   &
                    status)
        if (ft_checkerror_cfitsio(status, filename)) continue

        call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine readext_logical_1d


    !---------------------------------------------------------------------------


    subroutine readext_int8_1d(filename, output, status, hdu)

        character(len=*), intent(in)        :: filename
        integer*8, allocatable, intent(out) :: output(:)
        integer, intent(out)                :: status
        integer, optional, intent(in)       :: hdu

        integer                             :: unit, anynull, status_close
        integer, allocatable                :: imageshape(:)

        call ft_open_image(filename, unit, 1, imageshape, status, hdu=hdu)
        if (status /= 0) return

        !  Initialize variables
        allocate(output(imageshape(1)))

        call ftgpvk(unit, GROUP, 1, imageshape(1), NULLVAL, output, anynull,   &
                    status)
        if (ft_checkerror_cfitsio(status, filename)) continue

        call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine readext_int8_1d


    !---------------------------------------------------------------------------


    subroutine readext_double_1d(filename, output, status, hdu)

        character(len=*), intent(in)     :: filename
        real*8, allocatable, intent(out) :: output(:)
        integer, intent(out)             :: status
        integer, optional, intent(in)    :: hdu

        integer                          :: unit, anynull, status_close
        integer, allocatable             :: imageshape(:)

        call ft_open_image(filename, unit, 1, imageshape, status, hdu=hdu)
        if (status /= 0) return

        !  Initialize variables
        allocate(output(imageshape(1)))

        call ftgpvd(unit, GROUP, 1, imageshape(1), NULLVAL, output, anynull,   &
                    status)
        if (ft_checkerror_cfitsio(status, filename)) continue

        call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine readext_double_1d


    !---------------------------------------------------------------------------


    subroutine readext_logical_2d(filename, output, status, hdu)

        character(len=*), intent(in)        :: filename
        integer, intent(out)                :: status
        logical*1, allocatable, intent(out) :: output(:,:)
        integer, optional, intent(in)       :: hdu

        integer                             :: unit, firstpix, anynull, j,     &
                                               status_close
        integer, allocatable                :: imageshape(:)

        call ft_open_image(filename, unit, 2, imageshape, status, hdu=hdu)
        if (status /= 0) return

        ! Initialize variables
        firstpix=1
        allocate(output(imageshape(1), imageshape(2)))

        do j = 1, imageshape(2)

           call ftgpvb(unit, GROUP, firstpix, imageshape(1), NULLVAL,          &
                       output(:, j), anynull, status)
           if (ft_checkerror_cfitsio(status, filename)) go to 999

           firstpix = firstpix + imageshape(1)

        end do

    999 call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine readext_logical_2d


    !---------------------------------------------------------------------------


    subroutine readext_double_2d(filename, output, status, hdu)
    
        character(len=*), intent(in)     :: filename
        real*8, allocatable, intent(out) :: output(:,:)
        integer, intent(out)             :: status
        integer, optional, intent(in)    :: hdu

        integer                          :: unit, firstpix, anynull, j,        &
                                            status_close
        integer, allocatable             :: imageshape(:)
    
        call ft_open_image(filename, unit, 2, imageshape, status, hdu=hdu)
        if (status /= 0) return

        !  Initialize variables
        firstpix=1
        allocate(output(imageshape(1),imageshape(2)))

        do j = 1, imageshape(2)
    
            call ftgpvd(unit, GROUP, firstpix, imageshape(1), NULLVAL,         &
                        output(:,j), anynull, status)
            if (ft_checkerror_cfitsio(status, filename)) go to 999

            firstpix = firstpix + imageshape(1)

        end do

    999 call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine readext_double_2d


    !---------------------------------------------------------------------------


    subroutine readext_logical_3d(filename, output, status)

        character(len=*), intent(in)        :: filename
        logical*1, allocatable, intent(out) :: output(:,:,:)
        integer, intent(out)                :: status

        integer                             :: unit, firstpix, anynull, j, k,  &
                                               status_close
        integer, allocatable                :: imageshape(:)

        call ft_open_image(filename, unit, 3, imageshape, status)
        if (status /= 0) return

        !  Initialize variables
        firstpix = 1
        allocate(output(imageshape(1), imageshape(2), imageshape(3)))

        outer: do k = 1, imageshape(3)
            do j = 1, imageshape(2)

                call ftgpvb(unit, GROUP, firstpix, imageshape(1), NULLVAL,     &
                            output(:,j,k), anynull, status)
                if (ft_checkerror_cfitsio(status, filename)) exit outer

                firstpix = firstpix + BUFFERSIZE

            end do
        end do outer

        call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine readext_logical_3d


    !---------------------------------------------------------------------------


    subroutine readext_double_3d(filename, output, status)

        character(len=*), intent(in)     :: filename
        real*8, allocatable, intent(out) :: output(:,:,:)
        integer, intent(out)             :: status

        integer                          :: unit, firstpix, anynull, j, k,     &
                                            status_close
        integer, allocatable             :: imageshape(:)

        call ft_open_image(filename, unit, 3, imageshape, status)
        if (status /= 0) return

        !  Initialize variables
        firstpix = 1
        allocate(output(imageshape(1), imageshape(2), imageshape(3)))

        outer: do k = 1, imageshape(3)

            do j = 1, imageshape(2)

                call ftgpvd(unit, GROUP, firstpix, imageshape(1), NULLVAL,     &
                            output(:,j,k), anynull, status)
                if (ft_checkerror_cfitsio(status, filename)) exit outer

                firstpix = firstpix + imageshape(1)

            end do

        end do outer

        call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine readext_double_3d


    !---------------------------------------------------------------------------


    subroutine readslice_logical_filename(filename, x1, x2, array, status)

        character(len=*), intent(in) :: filename
        integer*8, intent(in)        :: x1, x2
        logical*1, intent(inout)     :: array(x2 - x1 + 1)
        integer, intent(out)         :: status
        integer                      :: unit, status_close

        call ft_open(filename, unit, status)
        if (status /= 0) return
        call readslice_logical_unit(unit, x1, x2, array, status)
        call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine readslice_logical_filename


    !---------------------------------------------------------------------------


    subroutine readslice_logical_unit(unit, x1, x2, array, status)

        integer, intent(in)      :: unit
        integer*8, intent(in)    :: x1, x2
        logical*1, intent(inout) :: array(x2 - x1 + 1)
        integer, intent(out)     :: status
        integer                  :: anynull

        status = 0
        call ftgpvd(unit, GROUP, x1, x2 - x1 + 1, NULLVAL, array, anynull, status)
        if (ft_checkerror_cfitsio(status)) return

    end subroutine readslice_logical_unit


    !---------------------------------------------------------------------------


    subroutine readslice_int4_filename(filename, x1, x2, array, status)

        character(len=*), intent(in) :: filename
        integer*8, intent(in)        :: x1, x2
        integer*4, intent(inout)     :: array(x2 - x1 + 1)
        integer, intent(out)         :: status
        integer                      :: unit, status_close

        call ft_open(filename, unit, status)
        if (status /= 0) return
        call readslice_int4_unit(unit, x1, x2, array, status)
        call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine readslice_int4_filename


    !---------------------------------------------------------------------------


    subroutine readslice_int4_unit(unit, x1, x2, array, status)

        integer, intent(in)      :: unit
        integer*8, intent(in)    :: x1, x2
        integer*4, intent(inout) :: array(x2 - x1 + 1)
        integer, intent(out)     :: status
        integer                  :: anynull

        status = 0
        call ftgpvj(unit, GROUP, x1, x2-x1+1, NULLVAL, array, anynull, status)
        if (ft_checkerror_cfitsio(status)) return

    end subroutine readslice_int4_unit


    !---------------------------------------------------------------------------


    subroutine readslice_int8_filename(filename, x1, x2, array, status)

        character(len=*), intent(in) :: filename
        integer*8, intent(in)        :: x1, x2
        integer*8, intent(inout)     :: array(x2 - x1 + 1)
        integer, intent(out)         :: status
        integer                      :: unit, status_close

        call ft_open(filename, unit, status)
        if (status /= 0) return
        call readslice_int8_unit(unit, x1, x2, array, status)
        call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine readslice_int8_filename


    !---------------------------------------------------------------------------


    subroutine readslice_int8_unit(unit, x1, x2, array, status)

        integer, intent(in)      :: unit
        integer*8, intent(in)    :: x1, x2
        integer*8, intent(inout) :: array(x2 - x1 + 1)
        integer, intent(out)     :: status
        integer                  :: anynull

        status = 0
        call ftgpvk(unit, GROUP, x1, x2-x1+1, NULLVAL, array, anynull, status)
        if (ft_checkerror_cfitsio(status)) return

    end subroutine readslice_int8_unit


    !---------------------------------------------------------------------------


    subroutine readslice_double_filename(filename, x1, x2, array, status)

        character(len=*), intent(in) :: filename
        integer*8, intent(in)        :: x1, x2
        real*8, intent(inout)        :: array(x2 - x1 + 1)
        integer, intent(out)         :: status
        integer                      :: unit, status_close

        call ft_open(filename, unit, status)
        if (status /= 0) return
        call readslice_double_unit(unit, x1, x2, array, status)
        call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine readslice_double_filename


    !---------------------------------------------------------------------------


    subroutine readslice_double_unit(unit, x1, x2, array, status)

        integer, intent(in)    :: unit
        integer*8, intent(in)  :: x1, x2
        real*8, intent(inout)  :: array(x2 - x1 + 1)
        integer, intent(out)   :: status
        integer                :: anynull

        status = 0
        call ftgpvd(unit, GROUP, x1, x2-x1+1, NULLVAL, array, anynull, status)
        if (ft_checkerror_cfitsio(status)) return

    end subroutine readslice_double_unit


    !---------------------------------------------------------------------------


    subroutine readslice_logical_3d(unit, x1, x2, y, z, naxes, array, status)

        integer, intent(in)      :: unit
        integer*8, intent(in)    :: x1, x2
        integer*8, intent(in)    :: y, z
        integer, intent(in)      :: naxes(3)
        logical*1, intent(inout) :: array(x2 - x1 + 1)
        integer, intent(out)     :: status
        integer                  :: anynull
        integer*8                :: firstpix

        status = 0
        firstpix = x1 + naxes(1) * (y - 1 + naxes(2) * (z - 1))
        call ftgpvb(unit, GROUP, firstpix, x2 - x1 + 1, NULLVAL, array, anynull, status)
        if (ft_checkerror_cfitsio(status)) return

    end subroutine readslice_logical_3d


    !---------------------------------------------------------------------------


    subroutine readslice_int4_3d(unit, x1, x2, y, z, naxes, array, status)

        integer, intent(in)      :: unit
        integer*8, intent(in)    :: x1, x2
        integer*8, intent(in)    :: y, z
        integer, intent(in)      :: naxes(3)
        integer*4, intent(inout) :: array(x2 - x1 + 1)
        integer, intent(out)     :: status
        integer                  :: anynull
        integer*8                :: firstpix

        status = 0
        firstpix = x1 + naxes(1) * (y - 1 + naxes(2) * (z - 1))
        call ftgpvj(unit, GROUP, firstpix, x2 - x1 + 1, NULLVAL, array, anynull, status)
        if (ft_checkerror_cfitsio(status)) return

    end subroutine readslice_int4_3d


    !---------------------------------------------------------------------------


    subroutine readslice_int8_3d(unit, x1, x2, y, z, naxes, array, status)

        integer, intent(in)      :: unit
        integer*8, intent(in)    :: x1, x2
        integer*8, intent(in)    :: y, z
        integer, intent(in)      :: naxes(3)
        integer*8, intent(inout) :: array(x2 - x1 + 1)
        integer, intent(out)     :: status
        integer                  :: anynull
        integer*8                :: firstpix

        status = 0
        firstpix = x1 + naxes(1) * (y - 1 + naxes(2) * (z - 1))
        call ftgpvk(unit, GROUP, firstpix, x2 - x1 + 1, NULLVAL, array, anynull, status)
        if (ft_checkerror_cfitsio(status)) return

    end subroutine readslice_int8_3d


    !---------------------------------------------------------------------------


    subroutine readslice_double_3d(unit, x1, x2, y, z, naxes, array, status)

        integer, intent(in)   :: unit
        integer*8, intent(in) :: x1, x2
        integer*8, intent(in) :: y, z
        integer, intent(in)   :: naxes(3)
        real*8, intent(inout) :: array(x2 - x1 + 1)
        integer, intent(out)  :: status
        integer               :: anynull
        integer*8             :: firstpix

        status = 0
        firstpix = x1 + naxes(1) * (y - 1 + naxes(2) * (z - 1))
        call ftgpvd(unit, GROUP, firstpix, x2 - x1 + 1, NULLVAL, array, anynull, status)
        if (ft_checkerror_cfitsio(status)) return

    end subroutine readslice_double_3d


    !---------------------------------------------------------------------------


    ! return true if extension is false. status may be set to a non-zero value
    ! if a problem unrelated to the existence of the extension occurs.
    function ft_test_extension(filename, status)
        logical                      :: ft_test_extension
        character(len=*), intent(in) :: filename
        integer                      :: unit
        integer                      :: status, status_close

        ft_test_extension = .false.

        ! get an unused logical unit number for the FITS file
        status = 0
        call ftgiou(unit, status)
        if (ft_checkerror_cfitsio(status, filename)) return

        ! open the fits file and move to the specified extension
        call ftnopn(unit, filename, CFITSIO_READONLY, status)
        ft_test_extension = status == 0
        if (status == 301) then
            status = 0
        else
            if (ft_checkerror_cfitsio(status, filename)) continue
        end if

        ! close logical unit
        call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end function ft_test_extension


    !---------------------------------------------------------------------------


    subroutine ft_open(filename, unit, status)
        character(len=*), intent(in) :: filename
        integer, intent(out)         :: unit
        integer, intent(out)         :: status
        integer                      :: status_close

        status = 0

        ! get an unused logical unit number for the FITS file
        call ftgiou(unit, status)
        if (ft_checkerror_cfitsio(status, filename)) return

        ! open the fits file and move to the specified extension
        call ftnopn(unit, filename, CFITSIO_READONLY, status)
        if (ft_checkerror_cfitsio(status, filename)) then
            call ft_close(unit, status_close)
        end if

    end subroutine ft_open


    !---------------------------------------------------------------------------


    subroutine ft_open_image(filename, unit, imagerank, imageshape, status, hdu)
        character(len=*), intent(in)      :: filename
        integer, intent(out)              :: unit
        integer, intent(in)               :: imagerank
        integer, allocatable, intent(out) :: imageshape(:)
        integer, intent(out)              :: status
        integer, intent(in), optional     :: hdu

        integer                           :: nfound, hdutype, blocksize
        integer                           :: imageshape_(8), status_close

        status = 0

        ! get an unused logical unit number for the FITS file
        call ftgiou(unit, status)
        if (ft_checkerror_cfitsio(status, filename)) return

        if (present(hdu)) then

            ! open the fits file and point to the primary header
            call ftopen(unit, filename, CFITSIO_READONLY, blocksize, status)
            if (ft_checkerror_cfitsio(status, filename)) go to 999

            ! move to the specified HDU
            call ftmahd(unit, hdu, hdutype, status)
            if (status == 0 .and. hdutype /= CFITSIO_IMAGE_HDU) then
                write (ERROR_UNIT,'(a,i0,a)') "HDU type is not an image: ", hdutype, "."
                status = 1
                go to 999
            endif

        else

            ! open the fits file and move to the HDU specified in the filename
            ! otherwise, move to the first image in FITS file
            call ftiopn(unit, filename, CFITSIO_READONLY, status)
            if (ft_checkerror_cfitsio(status, filename)) go to 999

        end if

        !  Determine the size of the image.
        call ftgknj(unit, 'NAXIS', 1, size(imageshape_), imageshape_, nfound, status)
        if (ft_checkerror_cfitsio(status, filename)) go to 999

        !  Check that it found the NAXISn keywords.
        if (nfound /= imagerank) then
            write (ERROR_UNIT, '(a)') 'FT_OPEN_IMAGE: Incompatible NAXISn keywords.'
            status = 1
            go to 999
        end if

        if (allocated(imageshape)) deallocate(imageshape)
        allocate(imageshape(imagerank))
        imageshape = imageshape_(1:imagerank)
        return

    999 call ft_close(unit, status_close)

    end subroutine ft_open_image


    !---------------------------------------------------------------------------


    subroutine ft_open_bintable(filename, unit, nrecords, status)
        character(len=*), intent(in) :: filename
        integer, intent(out)         :: unit, nrecords
        integer, intent(out)         :: status

        integer                      :: naxes(8)
        integer                      :: nfound, ncolumns, status_close

        status = 0

        ! get an unused logical unit number for the FITS file
        call ftgiou(unit, status)
        if (ft_checkerror_cfitsio(status, filename)) return

        ! open the fits file and move to the specified extension
        call fttopn(unit, filename, CFITSIO_READONLY, status)
        if (ft_checkerror_cfitsio(status, filename)) go to 999

        !  Determine the size of the image.
        call ftgknj(unit, 'NAXIS', 1, size(naxes), naxes, nfound, status)
        if (ft_checkerror_cfitsio(status, filename)) go to 999

        !  Check that it found the NAXISn keywords.
        if (nfound /= 2) then
            write (ERROR_UNIT, '(a)') 'FT_OPEN_BINTABLE: Incompatible NAXISn keywords.'
            status = 1
            go to 999
        end if

        ncolumns = naxes(1)
        nrecords = naxes(2)
        return

    999 call ft_close(unit, status_close)

    end subroutine ft_open_bintable


    !---------------------------------------------------------------------------


    subroutine ft_close(unit, status)

        integer, intent(in)  :: unit
        integer, intent(out) :: status
        integer              :: status_release

        !  The FITS file must always be closed before exiting the program.
        !  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
        status = 0
        call ftclos(unit, status)
        if (ft_checkerror_cfitsio(status)) continue
        status_release = 0
        call ftfiou(unit, status_release)
        if (ft_checkerror_cfitsio(status_release) .and. status == 0)           &
            status = status_release

    end subroutine ft_close


    !---------------------------------------------------------------------------


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


    !---------------------------------------------------------------------------


    subroutine writefits_double_1d(filename, data, header, status)
        character(len=*), intent(in)    :: filename
        real*8, intent(in)              :: data(:)
        character(len=*), intent(in), optional :: header
        integer, intent(out)            :: status
        integer                         :: irec, unit, blocksize, bitpix,      &
                                           naxis, naxes(1), status_close
        logical                         :: simple, extend

        ! delete file if it exists
        open(10, file=filename, status='unknown')
        close(10, status='delete')

        ! open and initialise the fits file
        status = 0
        call ftgiou(unit, status)
        if (ft_checkerror_cfitsio(status, filename)) return

        call ftinit(unit, filename, blocksize, status)
        if (ft_checkerror_cfitsio(status, filename)) go to 999
        simple = .true.
        bitpix = -64
        naxis  = 1
        naxes  = [ size(data) ]
        extend = .true.
        call FTPHPR(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
        if (ft_checkerror_cfitsio(status, filename)) go to 999

        ! write the astrometry keywords
        if (present(header)) then
            do irec=1, len(header) / 80
                call FTPREC(unit,header((irec-1)*80+1:irec*80), status)
            end do
            if (ft_checkerror_cfitsio(status, filename)) go to 999
        end if

        ! write the image data
        call FTPPRD(unit, GROUP, 1, size(data), data, status)
        if (ft_checkerror_cfitsio(status, filename)) go to 999

        ! close the fits file
    999 call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine writefits_double_1d


    !---------------------------------------------------------------------------


    subroutine writefits_double_2d(filename, data, header, status)
        character(len=*), intent(in)    :: filename
        real*8, intent(in)              :: data(:,:)
        character(len=*), intent(in), optional :: header
        integer, intent(out)            :: status
        integer                         :: irec, unit, blocksize, bitpix,      &
                                           naxis, naxes(2), status_close
        logical                         :: simple, extend

        ! delete file if it exists
        open(10, file=filename, status='unknown')
        close(10, status='delete')

        ! open and initialise the fits file
        status = 0
        call ftgiou(unit, status)
        if (ft_checkerror_cfitsio(status, filename)) return

        call ftinit(unit, filename, blocksize, status)
        if (ft_checkerror_cfitsio(status, filename)) go to 999
        simple = .true.
        bitpix = -64
        naxis  = 2
        naxes  = [size(data,1), size(data,2)]
        extend = .true.
        call FTPHPR(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
        if (ft_checkerror_cfitsio(status, filename)) go to 999

        ! write the astrometry keywords
        if (present(header)) then
            do irec=1, len(header) / 80
                call FTPREC(unit,header((irec-1)*80+1:irec*80), status)
            end do
            if (ft_checkerror_cfitsio(status, filename)) go to 999
        end if

        ! write the image data
        call FTPPRD(unit, GROUP, 1, size(data), data, status)
        if (ft_checkerror_cfitsio(status, filename)) go to 999

        ! close the fits file
    999 call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine writefits_double_2d


    !---------------------------------------------------------------------------


    subroutine writefits_double_3d(filename, data, header, status)
        character(len=*), intent(in)    :: filename
        real*8, intent(in)              :: data(:,:,:)
        character(len=*), intent(in), optional :: header
        integer, intent(out)            :: status
        integer                         :: irec, unit, blocksize, bitpix,      &
                                           naxis, status_close
        logical                         :: simple, extend

        ! delete file if it exists
        open(10, file=filename, status='unknown')
        close(10, status='delete')

        ! open and initialise the fits file
        status = 0
        call ftgiou(unit, status)
        if (ft_checkerror_cfitsio(status, filename)) return

        call ftinit(unit, filename, blocksize, status)
        if (ft_checkerror_cfitsio(status, filename)) go to 999
        simple = .true.
        bitpix = -64
        naxis  = 3
        extend = .true.
        call FTPHPR(unit,simple,bitpix,naxis, shape(data),0,1,extend,status)
        if (ft_checkerror_cfitsio(status, filename)) go to 999

        ! write the astrometry keywords
        if (present(header)) then
            do irec=1, len(header) / 80
                call FTPREC(unit,header((irec-1)*80+1:irec*80), status)
            end do
            if (ft_checkerror_cfitsio(status, filename)) go to 999
        end if

        ! write the image data
        call FTPPRD(unit, GROUP, 1, size(data), data, status)
        if (ft_checkerror_cfitsio(status, filename)) go to 999

        ! close the fits file
    999 call ft_close(unit, status_close)
        if (status == 0) status = status_close

    end subroutine writefits_double_3d


    !---------------------------------------------------------------------------


    subroutine ft_header2str(filename, header, status)
        character(len=*), intent(in)      :: filename
        character(len=28800), intent(out) :: header
        integer, intent(out)              :: status

        integer(kind=C_INT)               :: fits_status
        character(len=28800), pointer     :: f_header
        type(C_PTR)                       :: fptr, c_header
        integer                           :: nkeyrec

        status = 0
        call fits_open_file(fptr, filename // C_NULL_CHAR, CFITSIO_READONLY, status)
        if (ft_checkerror_cfitsio(status, filename)) return

        call fits_hdr2str(fptr, 1, C_NULL_PTR, 0, c_header, nkeyrec, status)
        if (ft_checkerror_cfitsio(status, filename)) return

        call fits_close_file(fptr, status)
        if (ft_checkerror_cfitsio(status, filename)) return

        ! check we can hold the header in memory
        if (nkeyrec*80 > len(f_header)) then
            status = 1
            write (ERROR_UNIT,'(a)') "FT_HEADER2STR: The FITS file header in '"&
                // trim(filename) // "' is too large."
            return
        end if

        call c_f_pointer(c_header, f_header)
        header = f_header(1:nkeyrec*80)
        status = 0
        return

    end subroutine ft_header2str


    !---------------------------------------------------------------------------


    subroutine ft_getparam(header, param, count, value, comment, must_exist, status)
        character(len=*), intent(in)   :: header
        character(len=*), intent(in)   :: param
        integer, intent(out)           :: count
        character(len=70), intent(out) :: value
        character(len=70), intent(out), optional :: comment
        logical, intent(in), optional  :: must_exist
        integer, intent(out)           :: status
        character(len=70)              :: buffer
        character(len=len(param))      :: strlowparam
        integer                        :: ncards, i, iparam

        status = 1
        value = ' '

        ncards = (len(header)-1) / 80 + 1! allow for an additional \0
        iparam = 0
        count = 0
        if (present(comment)) comment = ' '
        strlowparam = strlowcase(param)

        ! find the record number associated to the parameter
        do i = 1, ncards
           if (strlowcase(header((i-1)*80+1:(i-1)*80+8)) /= strlowparam) cycle
           count = count + 1
           if (iparam > 0) cycle
           iparam = i
        end do

        ! the parameter is not found
        if (count == 0) then
            if (present(must_exist)) then
                if (must_exist) then
                    write (ERROR_UNIT, '(a)') "Missing keyword '" // strupcase(param) // "' in FITS header."
                    return
                end if
            end if
            goto 999
        end if

        buffer = adjustl(header((iparam-1)*80+11:iparam*80))

        ! simple case, the value is not enclosed in quotes
        if (buffer(1:1) /= "'") then
            i = index(buffer, '/')
            if (i == 0) i = 71
            goto 888 ! end of slash
        end if

        ! find ending quote
        i = 2
        do while (i <= 70)
            if (buffer(i:i) == "'") then
                i = i + 1
                if (i == 71) goto 888 ! end of slash
                if (buffer(i:i) /= "'") exit
            end if
            i = i + 1
        end do

        if (i == 71) then
            write (ERROR_UNIT,'(a)') 'FITS card value not properly terminated by a quote.' // header((iparam-1)*80:iparam*80)
            return
        end if

        ! i points right after ending quote, let's find '/'
        do while (i <= 70)
            if (buffer(i:i) == '/') exit
            if (buffer(i:i) /= ' ') then
                write (ERROR_UNIT,'(a)') 'FITS card value is malformed. ' // header((iparam-1)*80:iparam*80)
                return
            end if
            i = i + 1
        end do

        ! i points to '/' or is 71
    888 value   = buffer(1:i-1)
        if (present(comment)) comment = buffer(i+1:)

    999 status  = 0

    end subroutine ft_getparam


    !---------------------------------------------------------------------------


    subroutine ft_readparam_logical(header, param, count, value, comment, must_exist, status)
        character(len=*), intent(in)             :: header
        character(len=*), intent(in)             :: param
        integer, intent(out)                     :: count
        logical, intent(out)                     :: value
        character(len=70), optional, intent(out) :: comment
        logical, optional, intent(in)            :: must_exist
        integer, intent(out)                     :: status
        character(len=70)                        :: charvalue

        value = .false.

        call ft_getparam(header, param, count, charvalue, comment, must_exist, status)
        if (status /= 0 .or. count == 0) return

        if (charvalue /= 'F' .and. charvalue /= 'T') then
            status = 1
            write (ERROR_UNIT,'(a)') "ft_readparam_logical: invalid logical value '" // trim(charvalue) // "' for parameter '" // &
                                     param // "' in FITS header."
            return
        end if

        if (charvalue == 'T') value = .true.

    end subroutine ft_readparam_logical


    !---------------------------------------------------------------------------


    subroutine ft_readparam_int4(header, param, count, value, comment, must_exist, status)
        character(len=*), intent(in)             :: header
        character(len=*), intent(in)             :: param
        integer, intent(out)                     :: count
        integer*4, intent(out)                   :: value
        character(len=70), optional, intent(out) :: comment
        logical, optional, intent(in)            :: must_exist
        integer, intent(out)                     :: status
        character(len=70)                        :: charvalue

        value = 0

        call ft_getparam(header, param, count, charvalue, comment, must_exist, status)
        if (status /= 0 .or. count == 0) return

        read (charvalue,'(i20)',iostat=status) value
        if (status /= 0) then
            write (ERROR_UNIT,'(a)') "ft_readparam_int4: invalid integer value '" // trim(charvalue) // "' for parameter '" // &
                                     param // "' in FITS header."
            return
        end if

    end subroutine ft_readparam_int4


    !---------------------------------------------------------------------------


    subroutine ft_readparam_int8(header, param, count, value, comment, must_exist, status)
        character(len=*), intent(in)             :: header
        character(len=*), intent(in)             :: param
        integer, intent(out)                     :: count
        integer*8, intent(out)                   :: value
        character(len=70), optional, intent(out) :: comment
        logical, optional, intent(in)            :: must_exist
        integer, intent(out)                     :: status
        character(len=70)                        :: charvalue

        value = 0

        call ft_getparam(header, param, count, charvalue, comment, must_exist, status)
        if (status /= 0 .or. count == 0) return

        read (charvalue,'(i20)',iostat=status) value
        if (status /= 0) then
            write (ERROR_UNIT,'(a)') "ft_readparam_int8: invalid integer value '" // trim(charvalue) // "' for parameter '" // &
                                     param // "' in FITS header."
            return
        end if

    end subroutine ft_readparam_int8


    !---------------------------------------------------------------------------


    subroutine ft_readparam_double(header, param, count, value, comment, must_exist, status)
        character(len=*), intent(in)             :: header
        character(len=*), intent(in)             :: param
        integer, intent(out)                     :: count
        real*8, intent(out)                      :: value
        character(len=70), optional, intent(out) :: comment
        logical, optional, intent(in)            :: must_exist
        integer, intent(out)                     :: status
        character(len=70)                        :: charvalue

        value = 0.d0

        call ft_getparam(header, param, count, charvalue, comment, must_exist, status)
        if (status /= 0 .or. count == 0) return

        read (charvalue,'(bn,f20.0)',iostat=status) value
        if (status /= 0) then
            write (ERROR_UNIT,'(a)') "ft_readparam_double: invalid real value '" // trim(charvalue) // "' for parameter '" // &
                                     param // "' in FITS header."
            return
        end if

    end subroutine ft_readparam_double


    !---------------------------------------------------------------------------


    subroutine ft_readparam_character(header, param, count, value, comment, must_exist, status)
        character(len=*), intent(in)             :: header
        character(len=*), intent(in)             :: param
        integer, intent(out)                     :: count
        character(len=70), intent(out)           :: value
        character(len=70), optional, intent(out) :: comment
        logical, optional, intent(in)            :: must_exist
        integer, intent(out)                     :: status
        character(len=70)                        :: charvalue
        integer                                  :: ncharvalue, i, j

        value = " "

        call ft_getparam(header, param, count, charvalue, comment, must_exist, status)
        if (status /= 0 .or. count == 0) return

        ! remove quotes
        if (charvalue(1:1) == "'") then
           charvalue = charvalue(2:len_trim(charvalue)-1)
        end if

        ! make double quotes single
        ncharvalue = len_trim(charvalue)
        i = 1
        j = 1
        do while (j <= ncharvalue)
           value(i:i) = charvalue(j:j)
           i = i + 1
           if (j /= ncharvalue) then
               if (charvalue(j:j+1) == "''") j = j + 1
           end if
           j = j + 1
        end do

    end subroutine ft_readparam_character


    !---------------------------------------------------------------------------


    function ft_checkerror_cfitsio(status, filename, hdu)

        !  This subroutine prints out the descriptive text corresponding to the
        !  error status value and prints out the contents of the internal
        !  error message stack generated by FITSIO whenever an error occurs.
        logical                                :: ft_checkerror_cfitsio
        integer, intent(in)                    :: status
        character(len=*), intent(in), optional :: filename
        integer, intent(in), optional          :: hdu
        character                              :: errtext*30, errmessage*80

        ft_checkerror_cfitsio = status /= 0

        !  Check if status is OK (no error); if so, simply return
        if (status == 0) return

        ! flush the output unit
        flush(OUTPUT_UNIT)

        !  Print the filename
        if (present(filename)) then
            write (ERROR_UNIT,*) 'In file ' // trim(filename) // ':'
        end if

        !  Print the HDU number
        if (present(hdu)) write (ERROR_UNIT,*) 'HDU: ', hdu

        !  The FTGERR subroutine returns a descriptive 30-character text string that
        !  corresponds to the integer error status number.  A complete list of all
        !  the error numbers can be found in the back of the FITSIO User's Guide.
        call ftgerr(status, errtext)
        write(ERROR_UNIT,*) 'FITSIO Error Status =', status, ': ', errtext

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

    end function ft_checkerror_cfitsio


end module module_fitstools
