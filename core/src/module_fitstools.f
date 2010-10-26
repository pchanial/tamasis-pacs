module module_fitstools

    use iso_c_binding
    use iso_fortran_env,  only : ERROR_UNIT, OUTPUT_UNIT
    use module_cfitsio
    use module_precision, only : sp, dp, qp
    use module_string,    only : strinteger, strlowcase, strsection, strupcase
    use module_tamasis,   only : p
    implicit none
    private

    public :: CFITSIO_IMAGE_HDU
    public :: CFITSIO_ASCII_TBL
    public :: CFITSIO_BINARY_TBL
    public :: CFITSIO_ANY_HDU
    public :: CFITSIO_READONLY
    public :: CFITSIO_READWRITE
    public :: FLEN_KEYWORD
    public :: FLEN_CARD
    public :: FLEN_VALUE
    public :: FLEN_COMMENT
    public :: FLEN_ERRMSG
    public :: FLEN_STATUS
    public :: ft_check_error_cfitsio
    public :: ft_close
    public :: ft_create_header
    public :: ft_header2str
    public :: ft_open
    public :: ft_open_bintable
    public :: ft_open_image
    public :: ft_read_column
    public :: ft_read_image
    public :: ft_read_keyword
    public :: ft_read_keyword_hcss
    public :: ft_read_slice
    public :: ft_test_extension
    public :: ft_write_image

    integer, private, parameter :: GROUP = 1
    integer, private, parameter :: NULLVAL = 0
    integer, private, parameter :: BUFFERSIZE = 1024

    integer, parameter :: FLEN_KEYWORD = 71
    integer, parameter :: FLEN_CARD    = 80
    integer, parameter :: FLEN_VALUE   = 70
    integer, parameter :: FLEN_COMMENT = 72
    integer, parameter :: FLEN_ERRMSG  = 80
    integer, parameter :: FLEN_STATUS  = 30

    interface ft_open
        module procedure ft_open_must_exist, ft_open_must_not_exist
    end interface ft_open

    interface ft_read_column
        module procedure ft_read_column_filename_character, ft_read_column_character_unit,                                         &
                         ft_read_column_filename_int4, ft_read_column_unit_int4,                                                   &
                         ft_read_column_filename_int8, ft_read_column_unit_int8,                                                   &
                         ft_read_column_filename_real4, ft_read_column_unit_real4,                                                 &
                         ft_read_column_filename_real8, ft_read_column_unit_real8
#if PRECISION_REAL == 16
        module procedure ft_read_column_filename_real16, ft_read_column_unit_real16
#endif
    end interface ft_read_column

    interface ft_read_image
        module procedure ft_read_image_logical_1d, ft_read_image_logical_2d, ft_read_image_logical_3d,                             &
                         ft_read_image_int8_1d,                                                                                    &
                         ft_read_image_real4_1d, ft_read_image_real4_2d, ft_read_image_real4_3d,                                   &
                         ft_read_image_real8_1d, ft_read_image_real8_2d, ft_read_image_real8_3d
#if PRECISION_REAL == 16
        module procedure ft_read_image_real16_1d, ft_read_image_real16_2d, ft_read_image_real16_3d
#endif
    end interface ft_read_image

    interface ft_read_slice
        module procedure ft_read_slice_logical_3d, ft_read_slice_int4_3d, ft_read_slice_int8_3d,                                   &
                         ft_read_slice_real4_3d, ft_read_slice_real8_3d
#if PRECISION_REAL == 16
        module procedure ft_read_slice_real16_3d
#endif
    end interface ft_read_slice

    interface ft_write_image
        module procedure ft_write_image_int4_1d, ft_write_image_int4_2d, ft_write_image_int4_3d,                                   &
                         ft_write_image_real4_1d, ft_write_image_real4_2d, ft_write_image_real4_3d,                                &
                         ft_write_image_real8_1d, ft_write_image_real8_2d, ft_write_image_real8_3d
#if PRECISION_REAL == 16
        module procedure ft_write_image_real16_1d, ft_write_image_real16_2d, ft_write_image_real16_3d
#endif
    end interface ft_write_image

    interface ft_read_keyword
        module procedure ft_read_keyword_header_logical, ft_read_keyword_header_int4, ft_read_keyword_header_int8,                 &
                         ft_read_keyword_header_real4, ft_read_keyword_header_real8, ft_read_keyword_header_character,             &
                         ft_read_keyword_unit_logical, ft_read_keyword_unit_int4, ft_read_keyword_unit_int8,                       &
                         ft_read_keyword_unit_real4, ft_read_keyword_unit_real8, ft_read_keyword_unit_character
#if PRECISION_REAL == 16
        module procedure ft_read_keyword_header_real16, ft_read_keyword_unit_real16
#endif
    end interface ft_read_keyword

    interface ft_read_keyword_hcss
        module procedure ft_read_keyword_hcss_logical, ft_read_keyword_hcss_int4, ft_read_keyword_hcss_int8,                       &
                         ft_read_keyword_hcss_real4, ft_read_keyword_hcss_real8, ft_read_keyword_hcss_character
#if PRECISION_REAL == 16
        module procedure ft_read_keyword_hcss_real16
#endif
    end interface ft_read_keyword_hcss


contains


    subroutine ft_read_column_filename_character(filename, colname, data, status)

        character(len=*), intent(in)               :: filename
        character(len=*), intent(in)               :: colname
        character(len=*), allocatable, intent(out) :: data(:)
        integer, intent(out)                       :: status

        integer :: nrecords, unit
        
        call ft_open_bintable(filename, unit, nrecords, status)
        if (status /= 0) return
        
        allocate(data(nrecords))
        call ft_read_column(unit, colname, 1, nrecords, data, status)
        if (status /= 0) return

        call ft_close(unit, status)

    end subroutine ft_read_column_filename_character


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_column_character_unit(unit, colname, first, last, data, status)
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
        if (ft_check_error_cfitsio(status, unit)) return

        ! extract column
        call ftgcvs(unit, colnum, first, 1, last-first+1, '', data, anyf, status)
        if (ft_check_error_cfitsio(status, unit)) return

    end subroutine ft_read_column_character_unit


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_column_filename_int4(filename, colname, data, status)

        character(len=*), intent(in)        :: filename
        character(len=*), intent(in)        :: colname
        integer*4, allocatable, intent(out) :: data(:)
        integer, intent(out)                :: status

        integer :: nrecords, unit
        
        call ft_open_bintable(filename, unit, nrecords, status)
        if (status /= 0) return
        
        allocate(data(nrecords))
        call ft_read_column(unit, colname, 1, nrecords, data, status)
        if (status /= 0) return
        
        call ft_close(unit, status)

    end subroutine ft_read_column_filename_int4


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_column_unit_int4(unit, colname, first, last, data, status)

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
        if (ft_check_error_cfitsio(status, unit)) return

        ! extract column
        call ftgcvj(unit, colnum, first, 1, last-first+1, nullval, data, anyf, status)
        if (ft_check_error_cfitsio(status, unit)) return

    end subroutine ft_read_column_unit_int4


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_column_filename_int8(filename, colname, data, status)

        character(len=*), intent(in)        :: filename
        character(len=*), intent(in)        :: colname
        integer*8, allocatable, intent(out) :: data(:)
        integer, intent(out)                :: status

        integer                             :: nrecords, unit
        
        call ft_open_bintable(filename, unit, nrecords, status)
        if (status /= 0) return
        
        allocate(data(nrecords))
        call ft_read_column(unit, colname, 1, nrecords, data, status)
        if (status /= 0) return
        
        call ft_close(unit, status)

    end subroutine ft_read_column_filename_int8


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_column_unit_int8(unit, colname, first, last, data, status)

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
        if (ft_check_error_cfitsio(status, unit)) return

        ! extract column
        call ftgcvk(unit, colnum, first, 1, last-first+1, nullval, data, anyf, status)
        if (ft_check_error_cfitsio(status, unit)) return

    end subroutine ft_read_column_unit_int8


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_column_filename_real4(filename, colname, data, status)

        character(len=*), intent(in)       :: filename
        character(len=*), intent(in)       :: colname
        real(sp), allocatable, intent(out) :: data(:)
        integer, intent(out)               :: status

        integer :: nrecords, unit
        
        call ft_open_bintable(filename, unit, nrecords, status)
        if (status /= 0) return
        
        allocate(data(nrecords))
        call ft_read_column(unit, colname, 1, nrecords, data, status)
        if (status /= 0) return
        
        call ft_close(unit, status)

    end subroutine ft_read_column_filename_real4


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_column_filename_real8(filename, colname, data, status)

        character(len=*), intent(in)       :: filename
        character(len=*), intent(in)       :: colname
        real(dp), allocatable, intent(out) :: data(:)
        integer, intent(out)               :: status

        integer :: nrecords, unit
        
        call ft_open_bintable(filename, unit, nrecords, status)
        if (status /= 0) return
        
        allocate(data(nrecords))
        call ft_read_column(unit, colname, 1, nrecords, data, status)
        if (status /= 0) return
        
        call ft_close(unit, status)

    end subroutine ft_read_column_filename_real8


    !-------------------------------------------------------------------------------------------------------------------------------


#if PRECISION_REAL == 16
    subroutine ft_read_column_filename_real16(filename, colname, data, status)

        character(len=*), intent(in)       :: filename
        character(len=*), intent(in)       :: colname
        real(qp), allocatable, intent(out) :: data(:)
        integer, intent(out)               :: status

        real(dp), allocatable :: data_(:)
        
        call ft_read_column(filename, colname, data_, status)
        if (status /= 0) return
        
        allocate(data(size(data_)))
        data = data_

    end subroutine ft_read_column_filename_real16
#endif

    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_column_unit_real4(unit, colname, first, last, data, status)
        integer, intent(in)          :: unit
        character(len=*), intent(in) :: colname
        integer, intent(in)          :: first, last
        real(sp), intent(out)        :: data(:)
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
        if (ft_check_error_cfitsio(status, unit)) return

        ! extract column
        call ftgcve(unit, colnum, first, 1, last-first+1, nullval, data, anyf, status)
        if (ft_check_error_cfitsio(status, unit)) return

    end subroutine ft_read_column_unit_real4


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_column_unit_real8(unit, colname, first, last, data, status)
        integer, intent(in)          :: unit
        character(len=*), intent(in) :: colname
        integer, intent(in)          :: first, last
        real(dp), intent(out)        :: data(:)
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
        if (ft_check_error_cfitsio(status, unit)) return

        ! extract column
        call ftgcvd(unit, colnum, first, 1, last-first+1, nullval, data, anyf, status)
        if (ft_check_error_cfitsio(status, unit)) return

    end subroutine ft_read_column_unit_real8


    !-------------------------------------------------------------------------------------------------------------------------------


#if PRECISION_REAL == 16
    subroutine ft_read_column_unit_real16(unit, colname, first, last, data, status)
        integer, intent(in)          :: unit
        character(len=*), intent(in) :: colname
        integer, intent(in)          :: first, last
        real(qp), intent(out)        :: data(:)
        integer, intent(out)         :: status

        real(dp) :: data_(size(data))

        call ft_read_column(unit, colname, first, last, data_, status)
        if (status /= 0) return

        data = data_

    end subroutine ft_read_column_unit_real16
#endif

    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_image_logical_1d(filename, output, status, hdu)

        character(len=*), intent(in)        :: filename
        logical*1, allocatable, intent(out) :: output(:)
        integer, intent(out)                :: status
        integer, optional, intent(in)       :: hdu

        integer                             :: unit, anynull
        integer, allocatable                :: imageshape(:)

        call ft_open_image(filename, unit, 1, imageshape, status, hdu=hdu)
        if (status /= 0) return

        !  Initialize variables
        allocate(output(imageshape(1)))

        call ftgpvb(unit, GROUP, 1, imageshape(1), NULLVAL, output, anynull, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        call ft_close(unit, status)

    end subroutine ft_read_image_logical_1d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_image_int8_1d(filename, output, status, hdu)

        character(len=*), intent(in)        :: filename
        integer*8, allocatable, intent(out) :: output(:)
        integer, intent(out)                :: status
        integer, optional, intent(in)       :: hdu

        integer                             :: unit, anynull
        integer, allocatable                :: imageshape(:)

        call ft_open_image(filename, unit, 1, imageshape, status, hdu=hdu)
        if (status /= 0) return

        !  Initialize variables
        allocate(output(imageshape(1)))

        call ftgpvk(unit, GROUP, 1, imageshape(1), NULLVAL, output, anynull, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        call ft_close(unit, status)

    end subroutine ft_read_image_int8_1d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_image_real4_1d(filename, output, status, hdu)

        character(len=*), intent(in)       :: filename
        real(sp), allocatable, intent(out) :: output(:)
        integer, intent(out)               :: status
        integer, optional, intent(in)      :: hdu

        integer                            :: unit, anynull
        integer, allocatable               :: imageshape(:)

        call ft_open_image(filename, unit, 1, imageshape, status, hdu=hdu)
        if (status /= 0) return

        !  Initialize variables
        allocate(output(imageshape(1)))

        call ftgpve(unit, GROUP, 1, imageshape(1), NULLVAL, output, anynull, status)
        if (ft_check_error_cfitsio(status, unit, filename)) continue

        call ft_close(unit, status)

    end subroutine ft_read_image_real4_1d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_image_real8_1d(filename, output, status, hdu)

        character(len=*), intent(in)       :: filename
        real(dp), allocatable, intent(out) :: output(:)
        integer, intent(out)               :: status
        integer, optional, intent(in)      :: hdu

        integer                            :: unit, anynull
        integer, allocatable               :: imageshape(:)

        call ft_open_image(filename, unit, 1, imageshape, status, hdu=hdu)
        if (status /= 0) return

        !  Initialize variables
        allocate(output(imageshape(1)))

        call ftgpvd(unit, GROUP, 1, imageshape(1), NULLVAL, output, anynull, status)
        if (ft_check_error_cfitsio(status, unit, filename)) continue

        call ft_close(unit, status)

    end subroutine ft_read_image_real8_1d


    !-------------------------------------------------------------------------------------------------------------------------------


#if PRECISION_REAL == 16
    subroutine ft_read_image_real16_1d(filename, output, status, hdu)

        character(len=*), intent(in)       :: filename
        real(qp), allocatable, intent(out) :: output(:)
        integer, intent(out)               :: status
        integer, optional, intent(in)      :: hdu

        real(dp), allocatable              :: output_(:)

        call ft_read_image(filename, output_, status, hdu)
        if (status /= 0) return

        allocate(output(size(output_)))
        output = real(output_, qp)

    end subroutine ft_read_image_real16_1d
#endif


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_image_logical_2d(filename, output, status, hdu)

        character(len=*), intent(in)        :: filename
        integer, intent(out)                :: status
        logical*1, allocatable, intent(out) :: output(:,:)
        integer, optional, intent(in)       :: hdu

        integer                             :: unit, firstpix, anynull, j
        integer, allocatable                :: imageshape(:)

        call ft_open_image(filename, unit, 2, imageshape, status, hdu=hdu)
        if (status /= 0) return

        ! Initialize variables
        firstpix=1
        allocate(output(imageshape(1), imageshape(2)))

        do j = 1, imageshape(2)

           call ftgpvb(unit, GROUP, firstpix, imageshape(1), NULLVAL, output(:, j), anynull, status)
           if (ft_check_error_cfitsio(status, unit, filename)) return

           firstpix = firstpix + imageshape(1)

        end do

        call ft_close(unit, status)

    end subroutine ft_read_image_logical_2d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_image_real4_2d(filename, output, status, hdu)
    
        character(len=*), intent(in)       :: filename
        real(sp), allocatable, intent(out) :: output(:,:)
        integer, intent(out)               :: status
        integer, optional, intent(in)      :: hdu

        integer                          :: unit, firstpix, anynull, j
        integer, allocatable             :: imageshape(:)
    
        call ft_open_image(filename, unit, 2, imageshape, status, hdu=hdu)
        if (status /= 0) return

        !  Initialize variables
        firstpix=1
        allocate(output(imageshape(1),imageshape(2)))

        do j = 1, imageshape(2)
    
            call ftgpve(unit, GROUP, firstpix, imageshape(1), NULLVAL, output(:,j), anynull, status)
            if (ft_check_error_cfitsio(status, unit, filename)) return

            firstpix = firstpix + imageshape(1)

        end do

        call ft_close(unit, status)

    end subroutine ft_read_image_real4_2d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_image_real8_2d(filename, output, status, hdu)
    
        character(len=*), intent(in)       :: filename
        real(dp), allocatable, intent(out) :: output(:,:)
        integer, intent(out)               :: status
        integer, optional, intent(in)      :: hdu

        integer                          :: unit, firstpix, anynull, j
        integer, allocatable             :: imageshape(:)
    
        call ft_open_image(filename, unit, 2, imageshape, status, hdu=hdu)
        if (status /= 0) return

        !  Initialize variables
        firstpix=1
        allocate(output(imageshape(1),imageshape(2)))

        do j = 1, imageshape(2)
    
            call ftgpvd(unit, GROUP, firstpix, imageshape(1), NULLVAL, output(:,j), anynull, status)
            if (ft_check_error_cfitsio(status, unit, filename)) return

            firstpix = firstpix + imageshape(1)

        end do

        call ft_close(unit, status)

    end subroutine ft_read_image_real8_2d


    !-------------------------------------------------------------------------------------------------------------------------------


#if PRECISION_REAL == 16
    subroutine ft_read_image_real16_2d(filename, output, status, hdu)

        character(len=*), intent(in)       :: filename
        real(qp), allocatable, intent(out) :: output(:,:)
        integer, intent(out)               :: status
        integer, optional, intent(in)      :: hdu

        real(dp), allocatable              :: output_(:,:)

        call ft_read_image(filename, output_, status, hdu)
        if (status /= 0) return

        allocate(output(size(output_,1),size(output_,2)))
        output = real(output_, qp)

    end subroutine ft_read_image_real16_2d
#endif


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_image_logical_3d(filename, output, status)

        character(len=*), intent(in)        :: filename
        logical*1, allocatable, intent(out) :: output(:,:,:)
        integer, intent(out)                :: status

        integer                             :: unit, firstpix, anynull, j, k
        integer, allocatable                :: imageshape(:)

        call ft_open_image(filename, unit, 3, imageshape, status)
        if (status /= 0) return

        !  Initialize variables
        firstpix = 1
        allocate(output(imageshape(1), imageshape(2), imageshape(3)))

        do k = 1, imageshape(3)
            do j = 1, imageshape(2)

                call ftgpvb(unit, GROUP, firstpix, imageshape(1), NULLVAL, output(:,j,k), anynull, status)
                if (ft_check_error_cfitsio(status, unit, filename)) return

                firstpix = firstpix + BUFFERSIZE

            end do
        end do

        call ft_close(unit, status)

    end subroutine ft_read_image_logical_3d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_image_real4_3d(filename, output, status)

        character(len=*), intent(in)       :: filename
        real(sp), allocatable, intent(out) :: output(:,:,:)
        integer, intent(out)               :: status

        integer                          :: unit, firstpix, anynull, j, k
        integer, allocatable             :: imageshape(:)

        call ft_open_image(filename, unit, 3, imageshape, status)
        if (status /= 0) return

        !  Initialize variables
        firstpix = 1
        allocate(output(imageshape(1), imageshape(2), imageshape(3)))

        do k = 1, imageshape(3)

            do j = 1, imageshape(2)

                call ftgpve(unit, GROUP, firstpix, imageshape(1), NULLVAL, output(:,j,k), anynull, status)
                if (ft_check_error_cfitsio(status, unit, filename)) return

                firstpix = firstpix + imageshape(1)

            end do

        end do

        call ft_close(unit, status)

    end subroutine ft_read_image_real4_3d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_image_real8_3d(filename, output, status)

        character(len=*), intent(in)       :: filename
        real(dp), allocatable, intent(out) :: output(:,:,:)
        integer, intent(out)               :: status

        integer                          :: unit, firstpix, anynull, j, k
        integer, allocatable             :: imageshape(:)

        call ft_open_image(filename, unit, 3, imageshape, status)
        if (status /= 0) return

        !  Initialize variables
        firstpix = 1
        allocate(output(imageshape(1), imageshape(2), imageshape(3)))

        do k = 1, imageshape(3)

            do j = 1, imageshape(2)

                call ftgpvd(unit, GROUP, firstpix, imageshape(1), NULLVAL, output(:,j,k), anynull, status)
                if (ft_check_error_cfitsio(status, unit, filename)) return

                firstpix = firstpix + imageshape(1)

            end do

        end do

        call ft_close(unit, status)

    end subroutine ft_read_image_real8_3d


    !-------------------------------------------------------------------------------------------------------------------------------


#if PRECISION_REAL == 16
    subroutine ft_read_image_real16_3d(filename, output, status)

        character(len=*), intent(in)       :: filename
        real(qp), allocatable, intent(out) :: output(:,:,:)
        integer, intent(out)               :: status

        real(dp), allocatable              :: output_(:,:,:)

        call ft_read_image(filename, output_, status)
        if (status /= 0) return

        allocate (output(size(output_,1),size(output_,2),size(output_,3)))
        output = output_

    end subroutine ft_read_image_real16_3d
#endif


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_slice_logical_3d(unit, x1, x2, y, z, naxes, array, status)

        integer, intent(in)      :: unit
        integer, intent(in)      :: x1, x2
        integer, intent(in)      :: y, z
        integer, intent(in)      :: naxes(3)
        logical*1, intent(inout) :: array(x2 - x1 + 1)
        integer, intent(out)     :: status

        integer                  :: anynull
        integer                  :: firstpix

        status = 0
        firstpix = x1 + naxes(1) * (y - 1 + naxes(2) * (z - 1))
        call ftgpvb(unit, GROUP, firstpix, x2 - x1 + 1, NULLVAL, array, anynull, status)
        if (ft_check_error_cfitsio(status, unit)) return

    end subroutine ft_read_slice_logical_3d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_slice_int4_3d(unit, x1, x2, y, z, naxes, array, status)

        integer, intent(in)      :: unit
        integer, intent(in)      :: x1, x2
        integer, intent(in)      :: y, z
        integer, intent(in)      :: naxes(3)
        integer*4, intent(inout) :: array(x2 - x1 + 1)
        integer, intent(out)     :: status

        integer                  :: anynull
        integer                  :: firstpix

        status = 0
        firstpix = x1 + naxes(1) * (y - 1 + naxes(2) * (z - 1))
        call ftgpvj(unit, GROUP, firstpix, x2 - x1 + 1, NULLVAL, array, anynull, status)
        if (ft_check_error_cfitsio(status, unit)) return

    end subroutine ft_read_slice_int4_3d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_slice_int8_3d(unit, x1, x2, y, z, naxes, array, status)

        integer, intent(in)      :: unit
        integer, intent(in)      :: x1, x2
        integer, intent(in)      :: y, z
        integer, intent(in)      :: naxes(3)
        integer*8, intent(inout) :: array(x2 - x1 + 1)
        integer, intent(out)     :: status

        integer                  :: anynull
        integer                  :: firstpix

        status = 0
        firstpix = x1 + naxes(1) * (y - 1 + naxes(2) * (z - 1))
        call ftgpvk(unit, GROUP, firstpix, x2 - x1 + 1, NULLVAL, array, anynull, status)
        if (ft_check_error_cfitsio(status, unit)) return

    end subroutine ft_read_slice_int8_3d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_slice_real4_3d(unit, x1, x2, y, z, naxes, array, status)

        integer, intent(in)     :: unit
        integer, intent(in)     :: x1, x2
        integer, intent(in)     :: y, z
        integer, intent(in)     :: naxes(3)
        real(sp), intent(inout) :: array(x2 - x1 + 1)
        integer, intent(out)    :: status

        integer                 :: anynull
        integer                 :: firstpix

        status = 0
        firstpix = x1 + naxes(1) * (y - 1 + naxes(2) * (z - 1))
        call ftgpve(unit, GROUP, firstpix, x2 - x1 + 1, NULLVAL, array, anynull, status)
        if (ft_check_error_cfitsio(status, unit)) return

    end subroutine ft_read_slice_real4_3d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_slice_real8_3d(unit, x1, x2, y, z, naxes, array, status)

        integer, intent(in)     :: unit
        integer, intent(in)     :: x1, x2
        integer, intent(in)     :: y, z
        integer, intent(in)     :: naxes(3)
        real(dp), intent(inout) :: array(x2 - x1 + 1)
        integer, intent(out)    :: status

        integer                 :: anynull
        integer                 :: firstpix

        status = 0
        firstpix = x1 + naxes(1) * (y - 1 + naxes(2) * (z - 1))
        call ftgpvd(unit, GROUP, firstpix, x2 - x1 + 1, NULLVAL, array, anynull, status)
        if (ft_check_error_cfitsio(status, unit)) return

    end subroutine ft_read_slice_real8_3d


    !-------------------------------------------------------------------------------------------------------------------------------


#if PRECISION_REAL == 16
    subroutine ft_read_slice_real16_3d(unit, x1, x2, y, z, naxes, array, status)

        integer, intent(in)     :: unit
        integer, intent(in)     :: x1, x2
        integer, intent(in)     :: y, z
        integer, intent(in)     :: naxes(3)
        real(qp), intent(inout) :: array(x2 - x1 + 1)
        integer, intent(out)    :: status

        real(dp) :: array_(size(array))

        call ft_read_slice(unit, x1, x2, y, z, naxes, array_, status)
        if (status /= 0) return
        array = array_

    end subroutine ft_read_slice_real16_3d
#endif


    !-------------------------------------------------------------------------------------------------------------------------------


    ! return true if extension is false. status may be set to a non-zero value
    ! if a problem unrelated to the existence of the extension occurs.
    function ft_test_extension(filename, status)

        logical                      :: ft_test_extension
        character(len=*), intent(in) :: filename
        integer, intent(out)         :: status

        integer                      :: unit
 
        ft_test_extension = .false.

        ! get an unused logical unit number for the FITS file
        status = 0
        call ftgiou(unit, status)
        if (ft_check_error_cfitsio(status, filename=filename)) return

        ! open the fits file and move to the specified extension
        call ftnopn(unit, filename, CFITSIO_READONLY, status)
        ft_test_extension = status == 0
        if (status == 301) then
            status = 0
        else
            if (ft_check_error_cfitsio(status, unit, filename)) return
        end if

        ! close logical unit
        call ft_close(unit, status)

    end function ft_test_extension


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_open_must_exist(filename, unit, status)

        character(len=*), intent(in) :: filename
        integer, intent(out)         :: unit
        integer, intent(out)         :: status

        status = 0

        ! get an unused logical unit number for the FITS file
        call ftgiou(unit, status)
        if (ft_check_error_cfitsio(status, filename=filename)) return

        ! open the fits file and move to the specified extension
        call ftnopn(unit, filename, CFITSIO_READONLY, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

    end subroutine ft_open_must_exist


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_open_must_not_exist(filename, unit, found, status)

        character(len=*), intent(in) :: filename
        integer, intent(out)         :: unit
        logical, intent(out)         :: found
        integer, intent(out)         :: status

        found = .false.
        status = 0

        ! get an unused logical unit number for the FITS file
        call ftgiou(unit, status)
        if (ft_check_error_cfitsio(status, filename=filename)) return

        ! open the fits file and move to the specified extension
        call ftnopn(unit, filename, CFITSIO_READONLY, status)
        if (status == 301) then
            call ft_close(unit, status)
            return
        end if
        if (ft_check_error_cfitsio(status, unit, filename)) return

        found = .true.

    end subroutine ft_open_must_not_exist


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_open_image(filename, unit, imagerank, imageshape, status, hdu)

        character(len=*), intent(in)      :: filename
        integer, intent(out)              :: unit
        integer, intent(in)               :: imagerank
        integer, allocatable, intent(out) :: imageshape(:)
        integer, intent(out)              :: status
        integer, intent(in), optional     :: hdu

        integer                           :: nfound, hdutype, blocksize
        integer                           :: imageshape_(8)

        status = 0

        ! get an unused logical unit number for the FITS file
        call ftgiou(unit, status)
        if (ft_check_error_cfitsio(status, filename=filename)) return

        if (present(hdu)) then

            ! open the fits file and point to the primary header
            call ftopen(unit, filename, CFITSIO_READONLY, blocksize, status)
            if (ft_check_error_cfitsio(status, unit, filename)) return

            ! move to the specified HDU
            call ftmahd(unit, hdu, hdutype, status)
            if (status == 0 .and. hdutype /= CFITSIO_IMAGE_HDU) then
                write (ERROR_UNIT,'(a,i0,a)') "FT_OPEN_IMAGE: HDU type is not an image: ", hdutype, " in file '" // filename // "'."
                call ft_close(unit, status)
                status = 1
                return
            endif

        else

            ! open the fits file and move to the HDU specified in the filename
            ! otherwise, move to the first image in FITS file
            call ftiopn(unit, filename, CFITSIO_READONLY, status)
            if (ft_check_error_cfitsio(status, unit, filename)) return

        end if

        !  Determine the size of the image.
        call ftgknj(unit, 'NAXIS', 1, size(imageshape_), imageshape_, nfound, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        !  Check that it found the NAXISn keywords.
        if (nfound /= imagerank) then
            write (ERROR_UNIT, '(a,2(i0,a))') "FT_OPEN_IMAGE: Incompatible NAXIS value '", nfound, "' instead of '", imagerank,    &
                  "' in file '" // filename // "'."
            call ft_close(unit, status)
            status = 1
            return
        end if

        if (allocated(imageshape)) deallocate(imageshape)
        allocate(imageshape(imagerank))
        imageshape = imageshape_(1:imagerank)

    end subroutine ft_open_image


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_open_bintable(filename, unit, nrecords, status)

        character(len=*), intent(in) :: filename
        integer, intent(out)         :: unit, nrecords
        integer, intent(out)         :: status

        integer                      :: naxes(8)
        integer                      :: nfound, ncolumns

        status = 0

        ! get an unused logical unit number for the FITS file
        call ftgiou(unit, status)
        if (ft_check_error_cfitsio(status, filename=filename)) return

        ! open the fits file and move to the specified extension
        call fttopn(unit, filename, CFITSIO_READONLY, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        !  Determine the size of the image.
        call ftgknj(unit, 'NAXIS', 1, size(naxes), naxes, nfound, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        !  Check that it found the NAXISn keywords.
        if (nfound /= 2) then
            write (ERROR_UNIT, '(a)') 'FT_OPEN_BINTABLE: Incompatible NAXISn keywords.'
            call ft_close(unit, status)
            status = 1
        end if

        ncolumns = naxes(1)
        nrecords = naxes(2)

    end subroutine ft_open_bintable


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_close(unit, status)

        integer, intent(in)  :: unit
        integer, intent(out) :: status

        integer              :: status_release

        !  The FITS file must always be closed before exiting the program.
        !  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
        status = 0
        call ftclos(unit, status)
        if (ft_check_error_cfitsio(status)) continue
        status_release = 0
        call ftfiou(unit, status_release)
        if (ft_check_error_cfitsio(status_release) .and. status == 0) then
            status = status_release
        end if

    end subroutine ft_close


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_create_header(naxis1, naxis2, cd, crval1, crval2, crpix1, crpix2, header)

        integer, intent(in) :: naxis1, naxis2
        real(p), intent(in) :: cd(2,2), crval1, crval2, crpix1, crpix2
        character(len=2880), intent(out) :: header

        character(len=*), parameter :: format_dbl = "(a8,'= ',f20.15)"
        character(len=*), parameter :: format_int = "(a8,'= ',i20)"
        character(len=*), parameter :: format_str = "(a8,'= ',a20)"

        header = ' '
        header(0*80+1:1*80) = 'SIMPLE  =                    T / Fits standard'
        header(1*80+1:2*80) = 'BITPIX  =                  -64 / Bits per pixel'
        header(2*80+1:3*80) = 'NAXIS   =                    2'
        write (header( 3*80+1: 4*80),format_int) 'NAXIS1  ', naxis1
        write (header( 4*80+1: 5*80),format_int) 'NAXIS2  ', naxis2
        write (header( 5*80+1: 6*80),format_dbl) 'CRVAL1  ', crval1
        write (header( 7*80+1: 8*80),format_dbl) 'CRVAL2  ', crval2
        write (header( 6*80+1: 7*80),format_dbl) 'CRPIX1  ', crpix1
        write (header( 8*80+1: 9*80),format_dbl) 'CRPIX2  ', crpix2
        write (header( 9*80+1:10*80),format_dbl) 'CD1_1   ', cd(1,1)
        write (header(10*80+1:11*80),format_dbl) 'CD2_1   ', cd(2,1)
        write (header(11*80+1:12*80),format_dbl) 'CD1_2   ', cd(1,2)
        write (header(12*80+1:13*80),format_dbl) 'CD2_2   ', cd(2,2)
        write (header(13*80+1:14*80),format_str) 'CTYPE1  ', "'RA---TAN'"
        write (header(14*80+1:15*80),format_str) 'CTYPE2  ', "'DEC--TAN'"
        write (header(15*80+1:16*80),format_str) 'CUNIT1  ', "'deg'"
        write (header(16*80+1:17*80),format_str) 'CUNIT2  ', "'deg'"
        header(17*80+1:18*80+3) = 'END'

    end subroutine ft_create_header


    !-------------------------------------------------------------------------------------------------------------------------------


    function check_duplicate_keywords(keyword) result (is_dup)

        character(len=*), intent(in) :: keyword
        logical                                 :: is_dup

        is_dup = keyword == 'SIMPLE' .or. keyword == 'BITPIX' .or. keyword == 'NAXIS' .or. keyword == 'EXTEND'
        if (is_dup) return

        if (keyword(1:5) == 'NAXIS') then
            is_dup = .true.
        end if

    end function check_duplicate_keywords


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_write_image_int4_1d(filename, data, header, status)

        character(len=*), intent(in)    :: filename
        integer*4, intent(in)           :: data(:)
        character(len=*), intent(in), optional :: header
        integer, intent(out)            :: status

        integer                  :: irec, unit, blocksize, bitpix, naxis, naxes(1)
        logical                  :: simple, extend
        character(len=FLEN_CARD) :: record

        ! delete file if it exists
        open(10, file=filename, status='unknown')
        close(10, status='delete')

        ! open and initialise the fits file
        status = 0
        call ftgiou(unit, status)
        if (ft_check_error_cfitsio(status, filename=filename)) return

        call ftinit(unit, filename, blocksize, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        simple = .true.
        bitpix = 32
        naxis  = 1
        naxes  = [ size(data) ]
        extend = .true.
        call FTPHPR(unit, simple, bitpix, naxis, naxes, 0, 1, extend, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        ! write the astrometry keywords
        if (present(header)) then
            do irec=1, len(header) / 80
                record = header((irec-1)*80+1:irec*80)
                if (check_duplicate_keywords(record(1:8))) cycle
                call FTPREC(unit,header((irec-1)*80+1:irec*80), status)
            end do
            if (ft_check_error_cfitsio(status, unit, filename)) return
        end if

        ! write the image data
        call FTPPRJ(unit, GROUP, 1, size(data), data, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        ! close the fits file
        call ft_close(unit, status)

    end subroutine ft_write_image_int4_1d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_write_image_int4_2d(filename, data, header, status)

        character(len=*), intent(in)    :: filename
        integer*4, intent(in)           :: data(:,:)
        character(len=*), intent(in), optional :: header
        integer, intent(out)            :: status

        integer                  :: irec, unit, blocksize, bitpix, naxis, naxes(2)
        logical                  :: simple, extend
        character(len=FLEN_CARD) :: record

        ! delete file if it exists
        open(10, file=filename, status='unknown')
        close(10, status='delete')

        ! open and initialise the fits file
        status = 0
        call ftgiou(unit, status)
        if (ft_check_error_cfitsio(status, filename=filename)) return

        call ftinit(unit, filename, blocksize, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        simple = .true.
        bitpix = 32
        naxis  = 2
        naxes  = [size(data,1), size(data,2)]
        extend = .true.
        call FTPHPR(unit, simple, bitpix, naxis, naxes, 0, 1, extend, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        ! write the astrometry keywords
        if (present(header)) then
            do irec=1, len(header) / 80
                record = header((irec-1)*80+1:irec*80)
                if (check_duplicate_keywords(record(1:8))) cycle
                call FTPREC(unit,header((irec-1)*80+1:irec*80), status)
            end do
            if (ft_check_error_cfitsio(status, unit, filename)) return
        end if

        ! write the image data
        call FTPPRJ(unit, GROUP, 1, size(data), data, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        ! close the fits file
        call ft_close(unit, status)

    end subroutine ft_write_image_int4_2d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_write_image_int4_3d(filename, data, header, status)

        character(len=*), intent(in)    :: filename
        integer*4, intent(in)           :: data(:,:,:)
        character(len=*), intent(in), optional :: header
        integer, intent(out)            :: status

        integer                  :: irec, unit, blocksize, bitpix, naxis
        logical                  :: simple, extend
        character(len=FLEN_CARD) :: record

        ! delete file if it exists
        open(10, file=filename, status='unknown')
        close(10, status='delete')

        ! open and initialise the fits file
        status = 0
        call ftgiou(unit, status)
        if (ft_check_error_cfitsio(status, filename=filename)) return

        call ftinit(unit, filename, blocksize, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        simple = .true.
        bitpix = 32
        naxis  = 3
        extend = .true.
        call FTPHPR(unit, simple, bitpix, naxis, shape(data), 0, 1, extend, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        ! write the astrometry keywords
        if (present(header)) then
            do irec=1, len(header) / 80
                record = header((irec-1)*80+1:irec*80)
                if (check_duplicate_keywords(record(1:8))) cycle
                call FTPREC(unit,header((irec-1)*80+1:irec*80), status)
            end do
            if (ft_check_error_cfitsio(status, unit, filename)) return
        end if

        ! write the image data
        call FTPPRJ(unit, GROUP, 1, size(data), data, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        ! close the fits file
        call ft_close(unit, status)

    end subroutine ft_write_image_int4_3d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_write_image_real4_1d(filename, data, header, status)

        character(len=*), intent(in)           :: filename
        real(sp), intent(in)                   :: data(:)
        character(len=*), intent(in), optional :: header
        integer, intent(out)                   :: status

        integer                  :: irec, unit, blocksize, bitpix, naxis, naxes(1)
        logical                  :: simple, extend
        character(len=FLEN_CARD) :: record

        ! delete file if it exists
        open(10, file=filename, status='unknown')
        close(10, status='delete')

        ! open and initialise the fits file
        status = 0
        call ftgiou(unit, status)
        if (ft_check_error_cfitsio(status, filename=filename)) return

        call ftinit(unit, filename, blocksize, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        simple = .true.
        bitpix = -32
        naxis  = 1
        naxes  = [ size(data) ]
        extend = .true.
        call FTPHPR(unit, simple, bitpix, naxis, naxes, 0, 1, extend, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        ! write the astrometry keywords
        if (present(header)) then
            do irec=1, len(header) / 80
                record = header((irec-1)*80+1:irec*80)
                if (check_duplicate_keywords(record(1:8))) cycle
                call FTPREC(unit,header((irec-1)*80+1:irec*80), status)
            end do
            if (ft_check_error_cfitsio(status, unit, filename)) return
        end if

        ! write the image data
        call FTPPRE(unit, GROUP, 1, size(data), data, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        ! close the fits file
        call ft_close(unit, status)

    end subroutine ft_write_image_real4_1d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_write_image_real8_1d(filename, data, header, status)

        character(len=*), intent(in)           :: filename
        real(dp), intent(in)                   :: data(:)
        character(len=*), intent(in), optional :: header
        integer, intent(out)                   :: status

        integer                  :: irec, unit, blocksize, bitpix, naxis, naxes(1)
        logical                  :: simple, extend
        character(len=FLEN_CARD) :: record

        ! delete file if it exists
        open(10, file=filename, status='unknown')
        close(10, status='delete')

        ! open and initialise the fits file
        status = 0
        call ftgiou(unit, status)
        if (ft_check_error_cfitsio(status, filename=filename)) return

        call ftinit(unit, filename, blocksize, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        simple = .true.
        bitpix = -64
        naxis  = 1
        naxes  = [ size(data) ]
        extend = .true.
        call FTPHPR(unit, simple, bitpix, naxis, naxes, 0, 1, extend, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        ! write the astrometry keywords
        if (present(header)) then
            do irec=1, len(header) / 80
                record = header((irec-1)*80+1:irec*80)
                if (check_duplicate_keywords(record(1:8))) cycle
                call FTPREC(unit,header((irec-1)*80+1:irec*80), status)
            end do
            if (ft_check_error_cfitsio(status, unit, filename)) return
        end if

        ! write the image data
        call FTPPRD(unit, GROUP, 1, size(data), data, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        ! close the fits file
        call ft_close(unit, status)

    end subroutine ft_write_image_real8_1d


    !-------------------------------------------------------------------------------------------------------------------------------


#if PRECISION_REAL == 16
    subroutine ft_write_image_real16_1d(filename, data, header, status)

        character(len=*), intent(in)           :: filename
        real(qp), intent(in)                   :: data(:)
        character(len=*), intent(in), optional :: header
        integer, intent(out)                   :: status

        real(dp) :: data_(size(data))

        data_ = data
        call ft_write_image(filename, data_, header, status)

    end subroutine ft_write_image_real16_1d
#endif


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_write_image_real4_2d(filename, data, header, status)

        character(len=*), intent(in)           :: filename
        real(sp), intent(in)                   :: data(:,:)
        character(len=*), intent(in), optional :: header
        integer, intent(out)                   :: status

        integer                  :: irec, unit, blocksize, bitpix, naxis, naxes(2)
        logical                  :: simple, extend
        character(len=FLEN_CARD) :: record

        ! delete file if it exists
        open(10, file=filename, status='unknown')
        close(10, status='delete')

        ! open and initialise the fits file
        status = 0
        call ftgiou(unit, status)
        if (ft_check_error_cfitsio(status, filename=filename)) return

        call ftinit(unit, filename, blocksize, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        simple = .true.
        bitpix = -32
        naxis  = 2
        naxes  = [size(data,1), size(data,2)]
        extend = .true.
        call FTPHPR(unit, simple, bitpix, naxis, naxes, 0, 1, extend, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        ! write the astrometry keywords
        if (present(header)) then
            do irec=1, len(header) / 80
                record = header((irec-1)*80+1:irec*80)
                if (check_duplicate_keywords(record(1:8))) cycle
                call FTPREC(unit,header((irec-1)*80+1:irec*80), status)
            end do
            if (ft_check_error_cfitsio(status, unit, filename)) return
        end if

        ! write the image data
        call FTPPRE(unit, GROUP, 1, size(data), data, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        ! close the fits file
        call ft_close(unit, status)

    end subroutine ft_write_image_real4_2d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_write_image_real8_2d(filename, data, header, status)

        character(len=*), intent(in)           :: filename
        real(dp), intent(in)                   :: data(:,:)
        character(len=*), intent(in), optional :: header
        integer, intent(out)                   :: status

        integer                  :: irec, unit, blocksize, bitpix, naxis, naxes(2)
        logical                  :: simple, extend
        character(len=FLEN_CARD) :: record

        ! delete file if it exists
        open(10, file=filename, status='unknown')
        close(10, status='delete')

        ! open and initialise the fits file
        status = 0
        call ftgiou(unit, status)
        if (ft_check_error_cfitsio(status, filename=filename)) return

        call ftinit(unit, filename, blocksize, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        simple = .true.
        bitpix = -64
        naxis  = 2
        naxes  = [size(data,1), size(data,2)]
        extend = .true.
        call FTPHPR(unit, simple, bitpix, naxis, naxes, 0, 1, extend, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        ! write the astrometry keywords
        if (present(header)) then
            do irec=1, len(header) / 80
                record = header((irec-1)*80+1:irec*80)
                if (check_duplicate_keywords(record(1:8))) cycle
                call FTPREC(unit,header((irec-1)*80+1:irec*80), status)
            end do
            if (ft_check_error_cfitsio(status, unit, filename)) return
        end if

        ! write the image data
        call FTPPRD(unit, GROUP, 1, size(data), data, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        ! close the fits file
        call ft_close(unit, status)

    end subroutine ft_write_image_real8_2d


    !-------------------------------------------------------------------------------------------------------------------------------


#if PRECISION_REAL == 16
    subroutine ft_write_image_real16_2d(filename, data, header, status)

        character(len=*), intent(in)           :: filename
        real(qp), intent(in)                   :: data(:,:)
        character(len=*), intent(in), optional :: header
        integer, intent(out)                   :: status

        real(dp) :: data_(size(data,1),size(data,2))

        data_ = data
        call ft_write_image(filename, data_, header, status)

    end subroutine ft_write_image_real16_2d
#endif


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_write_image_real4_3d(filename, data, header, status)

        character(len=*), intent(in)           :: filename
        real(sp), intent(in)                   :: data(:,:,:)
        character(len=*), intent(in), optional :: header
        integer, intent(out)                   :: status

        integer                  :: irec, unit, blocksize, bitpix, naxis
        logical                  :: simple, extend
        character(len=FLEN_CARD) :: record

        ! delete file if it exists
        open(10, file=filename, status='unknown')
        close(10, status='delete')

        ! open and initialise the fits file
        status = 0
        call ftgiou(unit, status)
        if (ft_check_error_cfitsio(status, filename=filename)) return

        call ftinit(unit, filename, blocksize, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        simple = .true.
        bitpix = -32
        naxis  = 3
        extend = .true.
        call FTPHPR(unit, simple, bitpix, naxis, shape(data), 0, 1, extend, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        ! write the astrometry keywords
        if (present(header)) then
            do irec=1, len(header) / 80
                record = header((irec-1)*80+1:irec*80)
                if (check_duplicate_keywords(record(1:8))) cycle
                call FTPREC(unit,header((irec-1)*80+1:irec*80), status)
            end do
            if (ft_check_error_cfitsio(status, unit, filename)) return
        end if

        ! write the image data
        call FTPPRE(unit, GROUP, 1, size(data), data, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        ! close the fits file
        call ft_close(unit, status)

    end subroutine ft_write_image_real4_3d


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_write_image_real8_3d(filename, data, header, status)

        character(len=*), intent(in)           :: filename
        real(dp), intent(in)                   :: data(:,:,:)
        character(len=*), intent(in), optional :: header
        integer, intent(out)                   :: status

        integer                  :: irec, unit, blocksize, bitpix, naxis
        logical                  :: simple, extend
        character(len=FLEN_CARD) :: record

        ! delete file if it exists
        open(10, file=filename, status='unknown')
        close(10, status='delete')

        ! open and initialise the fits file
        status = 0
        call ftgiou(unit, status)
        if (ft_check_error_cfitsio(status, filename=filename)) return

        call ftinit(unit, filename, blocksize, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        simple = .true.
        bitpix = -64
        naxis  = 3
        extend = .true.
        call FTPHPR(unit, simple, bitpix, naxis, shape(data), 0, 1, extend, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        ! write the astrometry keywords
        if (present(header)) then
            do irec=1, len(header) / 80
                record = header((irec-1)*80+1:irec*80)
                if (check_duplicate_keywords(record(1:8))) cycle
                call FTPREC(unit,header((irec-1)*80+1:irec*80), status)
            end do
            if (ft_check_error_cfitsio(status, unit, filename)) return
        end if

        ! write the image data
        call FTPPRD(unit, GROUP, 1, size(data), data, status)
        if (ft_check_error_cfitsio(status, unit, filename)) return

        ! close the fits file
        call ft_close(unit, status)

    end subroutine ft_write_image_real8_3d


    !-------------------------------------------------------------------------------------------------------------------------------


#if PRECISION_REAL == 16
    subroutine ft_write_image_real16_3d(filename, data, header, status)

        character(len=*), intent(in)           :: filename
        real(qp), intent(in)                   :: data(:,:,:)
        character(len=*), intent(in), optional :: header
        integer, intent(out)                   :: status

        real(dp) :: data_(size(data,1),size(data,2),size(data,3))

        data_ = data
        call ft_write_image(filename, data_, header, status)

    end subroutine ft_write_image_real16_3d
#endif


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_header2str(filename, header, status)

        character(len=*), intent(in)        :: filename
        character(len=*), intent(out)       :: header
        integer, intent(out)                :: status

        character(len=len(header)), pointer :: f_header
        type(C_PTR)                         :: fptr, c_header
        integer                             :: nkeyrec

        status = 0
        call fits_open_file(fptr, filename // C_NULL_CHAR, CFITSIO_READONLY, status)
        if (ft_check_error_cfitsio(status, filename=filename)) return

        call fits_hdr2str(fptr, 1, C_NULL_PTR, 0, c_header, nkeyrec, status)
        if (ft_check_error_cfitsio(status, filename=filename)) return

        call fits_close_file(fptr, status)
        if (ft_check_error_cfitsio(status, filename=filename)) return

        ! check we can hold the header in memory
        if (nkeyrec*80 > len(f_header)) then
            status = 1
            write (ERROR_UNIT,'(a)') "FT_HEADER2STR: The argument cannot store &
                &the FITS header of file '" // trim(filename) // "'. The length&
                & should be at least equal to '"//strinteger(nkeyrec*80)//"'."
            return
        end if

        call c_f_pointer(c_header, f_header)
        header = f_header(1:nkeyrec*80)

    end subroutine ft_header2str


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine get_keyword(header, keyword, value, found, comment, must_exist, status)

        character(len=*), intent(in)                       :: header
        character(len=*), intent(in)                       :: keyword
        character(len=FLEN_VALUE), intent(out)             :: value
        logical, intent(out)                               :: found
        character(len=FLEN_COMMENT), intent(out), optional :: comment
        logical, intent(in), optional                      :: must_exist
        integer, intent(out)                               :: status

        character(len=FLEN_VALUE)   :: buffer
        character(len=FLEN_KEYWORD) :: strlowkeyword
        integer                     :: i, ikeyword, ncards

        status = 0
        value = ' '

        ncards = (len(header)-1) / 80 + 1! allow for an additional \0
        found = .false.
        if (present(comment)) comment = ' '
        strlowkeyword = strlowcase(keyword)

        ! find the record number associated to the keyword
        do ikeyword = 1, ncards
           if (strlowcase(header((ikeyword-1)*80+1:(ikeyword-1)*80+8)) == strlowkeyword) then
              found = .true.
              exit
           end if
        end do

        ! the keyword is not found
        if (.not. found) then
            if (present(must_exist)) then
                if (must_exist) then
                    status = 202
                    write (ERROR_UNIT, '(a)') "Missing keyword '" // strupcase(keyword) // "' in FITS header."
                    return
                end if
            end if
            return
        end if

        buffer = adjustl(header((ikeyword-1)*80+11:ikeyword*80))

        ! simple case, the value is not enclosed in quotes
        if (buffer(1:1) /= "'") then
            i = index(buffer, '/')
            if (i == 0) i = 71
            goto 999 ! end of slash
        end if

        ! find ending quote
        i = 2
        do while (i <= 70)
            if (buffer(i:i) == "'") then
                i = i + 1
                if (i == 71) goto 999 ! end of slash
                if (buffer(i:i) /= "'") exit
            end if
            i = i + 1
        end do

        ! i points right after ending quote, let's find '/' ignoring what's before it
        do while (i <= 70)
            if (buffer(i:i) == '/') exit
            i = i + 1
        end do

        ! i points to '/' or is 71
    999 value = buffer(1:i-1)
        if (present(comment)) comment = buffer(i+1:)

    end subroutine get_keyword


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_keyword_header_logical(header, keyword, value, found, status, comment)

        character(len=*), intent(in)                       :: header
        character(len=*), intent(in)                       :: keyword
        logical, intent(out)                               :: value
        logical, intent(out), optional                     :: found
        integer, intent(out)                               :: status
        character(len=FLEN_COMMENT), optional, intent(out) :: comment

        character(len=FLEN_VALUE) :: charvalue
        logical                   :: found_

        value = .false.

        call get_keyword(header, keyword, charvalue, found_, comment, .not. present(found), status)
        if (present(found)) found = found_
        if (status /= 0 .or. .not. found_) return

        charvalue = strupcase(charvalue)
        if (charvalue /= 'FALSE'(1:min(5,len_trim(charvalue))) .and. charvalue /= 'TRUE'(1:min(4,len_trim(charvalue)))) then
            status = 404
            write (ERROR_UNIT,'(a)') "ft_read_keyword_header_logical: invalid logical value '" // trim(charvalue) //             &
                  "' for keyword '" // keyword // "' in FITS header."
            return
        end if

        if (charvalue(1:1) == 'T') value = .true.

    end subroutine ft_read_keyword_header_logical


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_keyword_header_int4(header, keyword, value, found, status, comment)

        character(len=*), intent(in)                       :: header
        character(len=*), intent(in)                       :: keyword
        integer*4, intent(out)                             :: value
        logical, intent(out), optional                     :: found
        integer, intent(out)                               :: status
        character(len=FLEN_COMMENT), optional, intent(out) :: comment

        character(len=FLEN_VALUE) :: charvalue
        logical                   :: found_

        value = 0

        call get_keyword(header, keyword, charvalue, found_, comment, .not. present(found), status)
        if (present(found)) found = found_
        if (status /= 0 .or. .not. found_) return

        read (charvalue,'(i20)',iostat=status) value
        if (status /= 0) then
            status = 407
            write (ERROR_UNIT,'(a)') "ft_read_keyword_header_int4: invalid integer value '" // trim(charvalue) //                &
                  "' for keyword '" // keyword // "' in FITS header."
            return
        end if

    end subroutine ft_read_keyword_header_int4


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_keyword_header_int8(header, keyword, value, found, status, comment)

        character(len=*), intent(in)                       :: header
        character(len=*), intent(in)                       :: keyword
        integer*8, intent(out)                             :: value
        logical, intent(out), optional                     :: found
        integer, intent(out)                               :: status
        character(len=FLEN_COMMENT), optional, intent(out) :: comment

        character(len=FLEN_VALUE) :: charvalue
        logical                   :: found_

        value = 0

        call get_keyword(header, keyword, charvalue, found_, comment, .not. present(found), status)
        if (present(found)) found = found_
        if (status /= 0 .or. .not. found_) return

        read (charvalue,'(i20)',iostat=status) value
        if (status /= 0) then
            status = 407
            write (ERROR_UNIT,'(a)') "ft_read_keyword_header_int8: invalid integer value '" // trim(charvalue) //                &
                  "' for keyword '" // keyword // "' in FITS header."
            return
        end if

    end subroutine ft_read_keyword_header_int8


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_keyword_header_real4(header, keyword, value, found, status, comment)

        character(len=*), intent(in)                       :: header
        character(len=*), intent(in)                       :: keyword
        real(sp), intent(out)                              :: value
        logical, intent(out), optional                     :: found
        integer, intent(out)                               :: status
        character(len=FLEN_COMMENT), optional, intent(out) :: comment

        real(dp) :: value_

        call ft_read_keyword(header, keyword, value_, found, status, comment)
        value = real(value_, sp)

    end subroutine ft_read_keyword_header_real4


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_keyword_header_real8(header, keyword, value, found, status, comment)

        character(len=*), intent(in)                       :: header
        character(len=*), intent(in)                       :: keyword
        real(dp), intent(out)                              :: value
        logical, intent(out), optional                     :: found
        integer, intent(out)                               :: status
        character(len=FLEN_COMMENT), optional, intent(out) :: comment

        character(len=FLEN_VALUE) :: charvalue
        logical                   :: found_

        value = 0

        call get_keyword(header, keyword, charvalue, found_, comment, .not. present(found), status)
        if (present(found)) found = found_
        if (status /= 0 .or. .not. found_) return

        read (charvalue,'(bn,f'//strinteger(FLEN_VALUE)//'.0)',iostat=status) value
        if (status /= 0) then
            status = 409
            write (ERROR_UNIT,'(a)') "ft_read_keyword: invalid real value '" // trim(charvalue) //                 &
                  "' for keyword '" // keyword // "' in FITS header."
            return
        end if

    end subroutine ft_read_keyword_header_real8


    !-------------------------------------------------------------------------------------------------------------------------------


#if PRECISION_REAL == 16
    subroutine ft_read_keyword_header_real16(header, keyword, value, found, status, comment)

        character(len=*), intent(in)                       :: header
        character(len=*), intent(in)                       :: keyword
        real(qp), intent(out)                              :: value
        logical, intent(out), optional                     :: found
        integer, intent(out)                               :: status
        character(len=FLEN_COMMENT), optional, intent(out) :: comment

        real(dp) :: value_

        call ft_read_keyword(header, keyword, value_, found, status, comment)
        value = real(value_, qp)

    end subroutine ft_read_keyword_header_real16
#endif


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_keyword_header_character(header, keyword, value, found, status, comment)

        character(len=*), intent(in)                       :: header
        character(len=*), intent(in)                       :: keyword
        character(len=*), intent(out)                      :: value
        logical, intent(out), optional                     :: found
        integer, intent(out)                               :: status
        character(len=FLEN_COMMENT), optional, intent(out) :: comment

        character(len=FLEN_VALUE) :: charvalue
        integer                   :: ncharvalue, i, j
        logical                   :: found_

        call get_keyword(header, keyword, charvalue, found_, comment, .not. present(found), status)
        if (present(found)) found = found_
        if (status /= 0 .or. .not. found_) return

        if (charvalue(1:1) /= "'") then
           value = charvalue
           return
        end if

        value = ' '

        ! remove leading quote
        charvalue = charvalue(2:len_trim(charvalue))

        ! scan until next '
        ! make double quotes single
        ncharvalue = len_trim(charvalue)
        i = 1
        j = 1
        do while (j <= ncharvalue)
           if (charvalue(j:j) == "'") then
              if (j == ncharvalue .or. charvalue(j+1:j+1) /= "'") exit
           end if
           value(i:i) = charvalue(j:j)
           i = i + 1
           if (j /= ncharvalue) then
               if (charvalue(j:j+1) == "''") j = j + 1
           end if
           j = j + 1
        end do

    end subroutine ft_read_keyword_header_character


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_keyword_unit_logical(unit, keyword, value, found, status, comment)

        integer, intent(in)                                :: unit
        character(len=*), intent(in)                       :: keyword
        logical*4, intent(out)                             :: value
        logical, intent(out), optional                     :: found
        integer, intent(out)                               :: status
        character(len=FLEN_COMMENT), optional, intent(out) :: comment

        logical                     :: found_, junk
        character(len=FLEN_COMMENT) :: comment_
        
        status = 0
        call ftgkyl(unit, keyword, value, comment_, status)
        found_ = status == 0
        if (.not. found_) then
            value = .false.
            if (status == 202 .and. present(found)) status = 0
            junk = ft_check_error_cfitsio(status, unit)
        end if

        if (present(found))   found = found_
        if (present(comment)) comment = comment_

    end subroutine ft_read_keyword_unit_logical


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_keyword_unit_int4(unit, keyword, value, found, status, comment)

        integer, intent(in)                                :: unit
        character(len=*), intent(in)                       :: keyword
        integer*4, intent(out)                             :: value
        logical, intent(out), optional                     :: found
        integer, intent(out)                               :: status
        character(len=FLEN_COMMENT), optional, intent(out) :: comment

        logical                     :: found_, junk
        character(len=FLEN_COMMENT) :: comment_
        
        status = 0
        call ftgkyj(unit, keyword, value, comment_, status)
        found_ = status == 0
        if (.not. found_) then
            value = 0
            if (status == 202 .and. present(found)) status = 0
            junk = ft_check_error_cfitsio(status, unit)
        end if

        if (present(found))   found = found_
        if (present(comment)) comment = comment_

    end subroutine ft_read_keyword_unit_int4


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_keyword_unit_int8(unit, keyword, value, found, status, comment)

        integer, intent(in)                                :: unit
        character(len=*), intent(in)                       :: keyword
        integer*8, intent(out)                             :: value
        logical, intent(out), optional                     :: found
        integer, intent(out)                               :: status
        character(len=FLEN_COMMENT), optional, intent(out) :: comment

        logical                     :: found_, junk
        character(len=fLEN_COMMENT) :: comment_
                
        status = 0
        call ftgkyk(unit, keyword, value, comment_, status)
        found_ = status == 0
        if (.not. found_) then
            value = 0
            if (status == 202 .and. present(found)) status = 0
            junk = ft_check_error_cfitsio(status, unit)
        end if

        if (present(found))   found = found_
        if (present(comment)) comment = comment_

    end subroutine ft_read_keyword_unit_int8


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_keyword_unit_real4(unit, keyword, value, found, status, comment)

        integer, intent(in)                                :: unit
        character(len=*), intent(in)                       :: keyword
        real(sp), intent(out)                              :: value
        logical, intent(out), optional                     :: found
        integer, intent(out)                               :: status
        character(len=FLEN_COMMENT), optional, intent(out) :: comment

        real(dp) :: value_

        call ft_read_keyword(unit, keyword, value_, found, status, comment)
        value = real(value_, sp)

    end subroutine ft_read_keyword_unit_real4


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_keyword_unit_real8(unit, keyword, value, found, status, comment)

        integer, intent(in)                                :: unit
        character(len=*), intent(in)                       :: keyword
        real(dp) , intent(out)                             :: value
        logical, intent(out), optional                     :: found
        integer, intent(out)                               :: status
        character(len=FLEN_COMMENT), optional, intent(out) :: comment

        logical                     :: found_, junk
        character(len=FLEN_COMMENT) :: comment_
        
        status = 0
        call ftgkyd(unit, keyword, value, comment_, status)
        found_ = status == 0
        if (.not. found_) then
            value = 0._dp
            if (status == 202 .and. present(found)) status = 0
            junk = ft_check_error_cfitsio(status, unit)
        end if

        if (present(found))   found = found_
        if (present(comment)) comment = comment_

    end subroutine ft_read_keyword_unit_real8


    !-------------------------------------------------------------------------------------------------------------------------------


#if PRECISION_REAL == 16
    subroutine ft_read_keyword_unit_real16(unit, keyword, value, found, status, comment)

        integer, intent(in)                                :: unit
        character(len=*), intent(in)                       :: keyword
        real(qp), intent(out)                              :: value
        logical, intent(out), optional                     :: found
        integer, intent(out)                               :: status
        character(len=FLEN_COMMENT), optional, intent(out) :: comment

        real(dp) :: value_

        call ft_read_keyword(unit, keyword, value_, found, status, comment)
        value = real(value_, qp)

    end subroutine ft_read_keyword_unit_real16
#endif


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_keyword_unit_character(unit, keyword, value, found, status, comment)

        integer, intent(in)                                :: unit
        character(len=*), intent(in)                       :: keyword
        character(len=*), intent(out)                      :: value
        logical, intent(out), optional                     :: found
        integer, intent(out)                               :: status
        character(len=FLEN_COMMENT), optional, intent(out) :: comment

        logical                     :: found_, junk
        character(len=FLEN_COMMENT) :: comment_
        
        status = 0
        call ftgkys(unit, keyword, value, comment_, status)
        found_ = status == 0
        if (.not. found_) then
            if (status == 202 .and. present(found)) status = 0
            junk = ft_check_error_cfitsio(status, unit)
        end if

        if (present(found))   found = found_
        if (present(comment)) comment = comment_

    end subroutine ft_read_keyword_unit_character


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_keyword_hcss_logical(unit, keyword, value, found, status, comment)

        integer, intent(in)                                :: unit
        character(len=*), intent(in)                       :: keyword
        logical, intent(out)                               :: value
        logical, intent(out), optional                     :: found
        integer, intent(out)                               :: status
        character(len=FLEN_COMMENT), optional, intent(out) :: comment

        logical                     :: found_
        integer                     :: ikey
        character(len=FLEN_COMMENT) :: comment_
        character(len=FLEN_KEYWORD) :: kwd
        
        ikey = 0
        do
            call ft_read_keyword(unit, 'key.META_' // strinteger(ikey), kwd, found_, status)
            if (status /= 0) return
            if (.not. found_) go to 999
            if (kwd == keyword) exit
            ikey = ikey + 1
        end do

        call ft_read_keyword(unit, 'META_' // strinteger(ikey), value, found_, status, comment_)
        if (status /= 0) return
        if (.not. found_) go to 999

        if (present(found)) found = .true.
        if (present(comment)) comment = comment_
        return

    999 value = .false.
        if (present(found)) then
            found = .false.
            return
        end if
        status = 1
        write (ERROR_UNIT,'(a)') "FT_READ_KEYWORD_HCSS: FITS keyword '" // keyword // "' is not found."

    end subroutine ft_read_keyword_hcss_logical


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_keyword_hcss_int4(unit, keyword, value, found, status, comment)

        integer, intent(in)                                :: unit
        character(len=*), intent(in)                       :: keyword
        integer*4, intent(out)                             :: value
        logical, intent(out), optional                     :: found
        integer, intent(out)                               :: status
        character(len=FLEN_COMMENT), optional, intent(out) :: comment

        logical                     :: found_
        integer                     :: ikey
        character(len=FLEN_COMMENT) :: comment_
        character(len=FLEN_KEYWORD) :: kwd
        
        ikey = 0
        do
            call ft_read_keyword(unit, 'key.META_' // strinteger(ikey), kwd, found_, status)
            if (status /= 0) return
            if (.not. found_) go to 999
            if (kwd == keyword) exit
            ikey = ikey + 1
        end do

        call ft_read_keyword(unit, 'META_' // strinteger(ikey), value, found_, status, comment_)
        if (status /= 0) return
        if (.not. found_) go to 999

        if (present(found)) found = .true.
        if (present(comment)) comment = comment_
        return

    999 value = 0
        if (present(found)) then
            found = .false.
            return
        end if
        status = 1
        write (ERROR_UNIT,'(a)') "FT_READ_KEYWORD_HCSS: FITS keyword '" // keyword // "' is not found."

    end subroutine ft_read_keyword_hcss_int4


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_keyword_hcss_int8(unit, keyword, value, found, status, comment)

        integer, intent(in)                                :: unit
        character(len=*), intent(in)                       :: keyword
        integer*8, intent(out)                             :: value
        logical, intent(out), optional                     :: found
        integer, intent(out)                               :: status
        character(len=FLEN_COMMENT), optional, intent(out) :: comment

        logical                     :: found_
        integer                     :: ikey
        character(len=FLEN_COMMENT) :: comment_
        character(len=FLEN_KEYWORD) :: kwd
        
        ikey = 0
        do
            call ft_read_keyword(unit, 'key.META_' // strinteger(ikey), kwd, found_, status)
            if (status /= 0) return
            if (.not. found_) go to 999
            if (kwd == keyword) exit
            ikey = ikey + 1
        end do

        call ft_read_keyword(unit, 'META_' // strinteger(ikey), value, found_, status, comment_)
        if (status /= 0) return
        if (.not. found_) go to 999

        if (present(found)) found = .true.
        if (present(comment)) comment = comment_
        return

    999 value = 0
        if (present(found)) then
            found = .false.
            return
        end if
        status = 1
        write (ERROR_UNIT,'(a)') "FT_READ_KEYWORD_HCSS: FITS keyword '" // keyword // "' is not found."

    end subroutine ft_read_keyword_hcss_int8


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_keyword_hcss_real4(unit, keyword, value, found, status, comment)

        integer, intent(in)                                :: unit
        character(len=*), intent(in)                       :: keyword
        real(sp), intent(out)                              :: value
        logical, intent(out), optional                     :: found
        integer, intent(out)                               :: status
        character(len=FLEN_COMMENT), optional, intent(out) :: comment

        real(dp) :: value_

        call ft_read_keyword_hcss(unit, keyword, value_, found, status, comment)
        value = real(value_, sp)

    end subroutine ft_read_keyword_hcss_real4


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_keyword_hcss_real8(unit, keyword, value, found, status, comment)

        integer, intent(in)                                :: unit
        character(len=*), intent(in)                       :: keyword
        real(dp), intent(out)                              :: value
        logical, intent(out), optional                     :: found
        integer, intent(out)                               :: status
        character(len=FLEN_COMMENT), optional, intent(out) :: comment

        logical                     :: found_
        integer                     :: ikey
        character(len=FLEN_COMMENT) :: comment_
        character(len=FLEN_KEYWORD) :: kwd
        
        ikey = 0
        do
            call ft_read_keyword(unit, 'key.META_' // strinteger(ikey), kwd, found_, status)
            if (status /= 0) return
            if (.not. found_) go to 999
            if (kwd == keyword) exit
            ikey = ikey + 1
        end do

        call ft_read_keyword(unit, 'META_' // strinteger(ikey), value, found_, status, comment_)
        if (status /= 0) return
        if (.not. found_) go to 999

        if (present(found)) found = .true.
        if (present(comment)) comment = comment_
        return

    999 value = 0._dp
        if (present(found)) then
            found = .false.
            return
        end if
        status = 1
        write (ERROR_UNIT,'(a)') "FT_READ_KEYWORD_HCSS: FITS keyword '" // keyword // "' is not found."

    end subroutine ft_read_keyword_hcss_real8


    !-------------------------------------------------------------------------------------------------------------------------------


#if PRECISION_REAL == 16
    subroutine ft_read_keyword_hcss_real16(unit, keyword, value, found, status, comment)

        integer, intent(in)                                :: unit
        character(len=*), intent(in)                       :: keyword
        real(qp), intent(out)                              :: value
        logical, intent(out), optional                     :: found
        integer, intent(out)                               :: status
        character(len=FLEN_COMMENT), optional, intent(out) :: comment

        real(dp) :: value_

        call ft_read_keyword_hcss(unit, keyword, value_, found, status, comment)
        value = real(value_, qp)

    end subroutine ft_read_keyword_hcss_real16
#endif


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine ft_read_keyword_hcss_character(unit, keyword, value, found, status, comment)

        integer, intent(in)                                :: unit
        character(len=*), intent(in)                       :: keyword
        character(len=*), intent(out)                      :: value
        logical, intent(out), optional                     :: found
        integer, intent(out)                               :: status
        character(len=FLEN_COMMENT), optional, intent(out) :: comment

        logical           :: found_
        integer           :: ikey
        character(len=FLEN_COMMENT) :: comment_
        character(len=FLEN_KEYWORD) :: kwd
        
        ikey = 0
        do
            call ft_read_keyword(unit, 'key.META_' // strinteger(ikey), kwd, found_, status)
            if (status /= 0) return
            if (.not. found_) go to 999
            if (kwd == keyword) exit
            ikey = ikey + 1
        end do

        call ft_read_keyword(unit, 'META_' // strinteger(ikey), value, found_, status, comment_)
        if (status /= 0) return
        if (.not. found_) go to 999

        if (present(found)) found = .true.
        if (present(comment)) comment = comment_
        return

    999 value = ''
        if (present(found)) then
            found = .false.
            return
        end if
        status = 1
        write (ERROR_UNIT,'(a)') "FT_READ_KEYWORD_HCSS: FITS keyword '" // keyword // "' is not found."

    end subroutine ft_read_keyword_hcss_character


    !-------------------------------------------------------------------------------------------------------------------------------


    function ft_check_error_cfitsio(status, unit, filename, hdu)

        !  This subroutine prints out the descriptive text corresponding to the
        !  error status value and prints out the contents of the internal
        !  error message stack generated by FITSIO whenever an error occurs.
        logical                                :: ft_check_error_cfitsio
        integer, intent(in)                    :: status
        integer, intent(in), optional          :: unit
        character(len=*), intent(in), optional :: filename
        integer, intent(in), optional          :: hdu

        character(len=FLEN_STATUS)             :: errtext
        character(len=FLEN_ERRMSG)             :: errmessage
        integer                                :: status_close

        ft_check_error_cfitsio = status /= 0

        !  Check if status is OK (no error); if so, simply return
        if (status == 0) return

        ! file not found
        if (status == 104 .and. present(filename)) then
            write (ERROR_UNIT,'(a)') "ERROR: Failure to open or find file '" // trim(filename) // "'."
            go to 999
        end if

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

    999 if (present(unit)) then
            call ft_close(unit, status_close)
        end if

    end function ft_check_error_cfitsio


end module module_fitstools
