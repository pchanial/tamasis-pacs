module module_madcap

    use iso_fortran_env,       only : ERROR_UNIT
    use module_filtering,      only : filterset
    use module_pointingmatrix, only : pointingelement
    use module_precision,      only : p
    use module_string,         only : strinteger
    implicit none
    private

    public :: read_filter
    public :: read_tod
    public :: read_tod_header

    interface read_tod
        module procedure read_tod_pix_per_sample_is_one, read_tod_madmap1
    end interface read_tod

    interface read_tod_header
        module procedure read_tod_header_pix_per_sample_is_one, read_tod_header_madmap1
    end interface read_tod_header


contains


    subroutine read_tod_header_pix_per_sample_is_one(signalfile, weightfile, pixelfile, nsamples, npixels_per_sample,     &
                                                     status)

        character(len=*), intent(in)       :: signalfile, weightfile, pixelfile
        integer*8, intent(out)             :: nsamples
        integer, intent(out)               :: npixels_per_sample, status

        integer*8 :: signalsize, othersize
        logical   :: exist
!XXX gfortran 4.5 inquire size does not work
#ifdef GFORTRAN
        integer*4 :: sarray(13)
#endif

        status = 1
        inquire (file=signalfile, exist=exist)
        if (.not. exist) then
            write (ERROR_UNIT,'(a)') "Error: Failed to open signal file '" // signalfile // "'."
            return
        end if
        inquire (file=weightfile, exist=exist)
        if (.not. exist) then
            write (ERROR_UNIT,'(a)') "Error: Failed to open weight file '" // weightfile // "'."
            return
        end if
        inquire (file=pixelfile, exist=exist)
        if (.not. exist) then
            write (ERROR_UNIT,'(a)') "Error: Failed to open pixel file '" // pixelfile // "'."
            return
        end if

        inquire (file=signalfile, size=signalsize)
#ifdef GFORTRAN
        call stat(signalfile, sarray, status)
        signalsize = sarray(8)
        status = 1
#endif
        if (mod(signalsize,8) /= 0) then
            write (ERROR_UNIT,'(a)') "Error: Size of file '" // signalfile // "' is not a multiple of 8 (" //                      &
                  strinteger(signalsize) // ")."
            return
        end if

        inquire (file=weightfile, size=othersize)
#ifdef GFORTRAN
        call stat(weightfile, sarray, status)
        othersize = sarray(8)
        status = 1
#endif
        if (othersize*2 /= signalsize) then
            write (ERROR_UNIT,'(a,2(i0,a))') "Error: Size of weight file '", othersize, "' is not compatible with that of the signa&
                  &l file (", signalsize, ")."
            return
        end if

        inquire (file=pixelfile, size=othersize)
#ifdef GFORTRAN
        call stat(pixelfile, sarray, status)
        othersize = sarray(8)
        status = 1
#endif
        if (othersize*2 /= signalsize) then
            write (ERROR_UNIT,'(a,2(i0,a))') "Error: Size of pixel file '",  othersize, "' is not compatible with that of the signa&
                  &l file (", signalsize, ")."
            return
        end if

        nsamples = signalsize / 8
        npixels_per_sample = 1
        status = 0

    end subroutine read_tod_header_pix_per_sample_is_one


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine read_tod_header_madmap1(filename, convert, nsamples, npixels_per_sample, status)

        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: convert
        integer*8, intent(out)       :: nsamples
        integer, intent(out)         :: npixels_per_sample, status

        integer*8 :: npixels_map, npixels_per_sample_, filesize, first, last
        logical   :: exist
#ifdef GFORTRAN
        integer*4 :: sarray(13)
#endif

        status = 1
        inquire (file=filename, exist=exist)
        if (.not. exist) then
            write (ERROR_UNIT,'(a)') "Error: Failed to open file '" // filename // "'."
            return
        end if

        inquire (file=filename, size=filesize)
#ifdef GFORTRAN
        call stat(filename, sarray, status)
        filesize = sarray(8)
        status = 1
#endif
        if (filesize < 4*8) then
            write (ERROR_UNIT,'(a)') "Error: Size of file '" // filename // "' is too small (" // strinteger(filesize) // ")."
            return
        end if

        open  (11, file=filename, access='stream', form='unformatted', convert=convert, action='read')
        read  (11) first, last, npixels_per_sample_, npixels_map
        close (11)

        first = first + 1
        last  = last + 1
        nsamples = last - first + 1
        npixels_per_sample = npixels_per_sample_

        if (first < 1 .or. last < 1 .or. npixels_per_sample < 1 .or. npixels_per_sample > 1000 .or. npixels_map < 1) then
            write (ERROR_UNIT,'(a,4(i0,a))') 'Error: Invalid tod header (first=', first, ', last=', last, ', npixels_per_sample=', &
                  npixels_per_sample, ', npixels in map=', npixels_map, '). Check endianness.'
            return
        end if

        status = 0

    end subroutine read_tod_header_madmap1


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine read_tod_pix_per_sample_is_one(signalfile, weightfile, pixelfile, convert, first, last, tod, pmatrix, status)

        character(len=*), intent(in)       :: signalfile, weightfile, pixelfile
        character(len=*), intent(in)       :: convert
        integer*8, intent(in)              :: first(:), last(:)
        real(p), intent(out)               :: tod(:,:)
        type(pointingelement), intent(out) :: pmatrix(:,:,:)
        integer, intent(out)               :: status

        integer   :: ndetectors, nslices, idetector, islice
#ifdef GFORTRAN
        integer*4 :: sarray(13)
#endif

        nslices = size(first)
        ndetectors = size(tod,2)
        
        status = 1

        open (11, file=signalfile, access='stream', form='unformatted', convert=convert, action='read', iostat=status)
        if (status /= 0) then
            write (ERROR_UNIT,'(a)') "Error: Failed to open signal file '" // signalfile // "'."
            return
        end if

        open (12, file=weightfile, access='stream', form='unformatted', convert=convert, action='read', iostat=status)
        if (status /= 0) then
            write (ERROR_UNIT,'(a)') "Error: Failed to open weight file '" // weightfile // "'."
            close (11)
            return
        end if

        open (13, file=pixelfile, access='stream', form='unformatted', convert=convert, action='read', iostat=status)
        if (status /= 0) then
            write (ERROR_UNIT,'(a)') "Error: Failed to open pixel file '" // pixelfile // "'."
            close (11)
            close (12)
            return
        end if

        do islice = 1, nslices

            do idetector = 1, ndetectors

                read (11, iostat=status) tod(first(islice):last(islice),idetector)
                if (status /= 0) then
                    write (ERROR_UNIT,'(a)') "Error: Failed to read signal file '" // signalfile // "'."
                    go to 999
                end if
                read (12, iostat=status) pmatrix(1,first(islice):last(islice),idetector)%weight
                if (status /= 0) then
                    write (ERROR_UNIT,'(a)') "Error: Failed to read weight file '" // weightfile // "'."
                    go to 999
                end if
                read (13, iostat=status) pmatrix(1,first(islice):last(islice),idetector)%pixel
                if (status /= 0) then
                    write (ERROR_UNIT,'(a)') "Error: Failed to read pixel file '" // pixelfile // "'."
                    go to 999
                end if

            end do

        end do

    999 close (11)
        close (12)
        close (13)

    end subroutine read_tod_pix_per_sample_is_one


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine read_tod_madmap1(filename, convert, first, last, tod, pmatrix, status)

        character(len=*), intent(in)       :: filename
        character(len=*), intent(in)       :: convert
        integer*8, intent(in)              :: first(:), last(:)
        real(p), intent(out)               :: tod(:,:)
        type(pointingelement), intent(out) :: pmatrix(:,:,:)
        integer, intent(out)               :: status

        integer*8 :: filesize, nsamples, isample
        integer   :: npixels_per_sample, ndetectors, nslices, idetector, islice
#ifdef GFORTRAN
        integer*4 :: sarray(13)
#endif

        ndetectors = size(tod,2)
        npixels_per_sample = size(pmatrix,1)
        nslices = size(first)
        nsamples = sum(last-first+1)

        inquire (file=filename, size=filesize)
#ifdef GFORTRAN
        call stat(filename, sarray, status)
        filesize = sarray(8)
#endif

print *, 'f', first
print *, last
print *, 'nsamples:',nsamples, 'npixels_per_sample', npixels_per_sample, 'ndetector', ndetectors
        
        if (filesize /= 4*8 + 8 * nsamples * (1 + npixels_per_sample) * ndetectors) then
            status = 1
            write (ERROR_UNIT,'(a,2(i0,a))') "Error: Invalid file size ('", filesize, "' instead of '",                            &
                  4*8 + 8 * nsamples * (1 + npixels_per_sample) * ndetectors, "')."
            return
        end if

        open (11, file=filename, access='stream', form='unformatted', convert=convert, action='read', iostat=status)
        if (status /= 0) then
            write (ERROR_UNIT,'(a)') "Error: Failed to open file '" // filename // "'."
            return
        end if

        if (sum(last-first+1) /= size(tod,1)) then
            status = 1
            write (ERROR_UNIT,'(a,2(i0,a))') "Error: Invalid number of samples in input tod ('", size(tod,1), "' instead of '",    &
                  sum(last-first+1), "')."
            return
        end if

        ! read the tod
        do islice = 1, nslices

            do idetector = 1, ndetectors

                do isample = first(islice), last(islice)

                    read (11, iostat=status) tod(isample,idetector)
                    if (status /= 0) go to 999
                    read (11, iostat=status) pmatrix(:,isample,idetector)
                    if (status /= 0) go to 999

                end do

            end do

        end do

    999 close (11)
        
    end subroutine read_tod_madmap1


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine read_filter_header(filename, convert, first, last, ncorrelations, status)

        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: convert
        integer*8, intent(out)       :: first, last
        integer, intent(out)         :: ncorrelations
        integer, intent(out)         :: status

        integer   :: filesize
        integer*8 :: ncorrelations_
#ifdef GFORTRAN
        integer*4 :: sarray(13)
#endif
    
        status = 1

        inquire (file=filename, size=filesize)
!XXX gfortran 4.5 inquire size does not work
#ifdef GFORTRAN
        call stat(filename, sarray, status)
        filesize = sarray(8)
#endif
        if (filesize < 4*8) then
            write (ERROR_UNIT,'(a)') "Error: Size of file '" // filename // "' is too small (" // strinteger(filesize) // ")."
            return
        end if

        if (mod(filesize,8) /= 0) then
            write (ERROR_UNIT,'(a)') "Error: Size of file '" // filename // "' is not a multiple of 8."
            return
        end if

        open  (11, file=filename, form='unformatted', access='stream', convert=convert, action='read')
        read  (11, pos=1) first, last, ncorrelations_
        close (11)

        first = first + 1
        last  = last  + 1
        ncorrelations = ncorrelations_

        if (first < 1 .or. last < 1 .or. ncorrelations < 0) then
            write (ERROR_UNIT,'(a,3(i0,a))') 'Error: Invalid filter header (first=', first, ', last=', last, ', ncorrelations=',   &
                  ncorrelations, '). Check endianness.'
            return
        end if

        if (filesize / 8 - 3 /= ncorrelations + 1) then
            write (ERROR_UNIT,'(a,2(i0,a))') "Error: Invalid number of values in file '" // filename // "' ('", filesize / 8 - 3,  &
                  "' instead of '", ncorrelations+1, "'). Check endianness."
            return
        end if

        status = 0

    end subroutine read_filter_header


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine read_filter_headers(filename, convert, ndetectors, first, last, ncorrelations, status)

        character(len=*), intent(in)        :: filename
        character(len=*), intent(in)        :: convert
        integer, intent(in)                 :: ndetectors
        integer*8, allocatable, intent(out) :: first(:), last(:)
        integer, allocatable, intent(out)   :: ncorrelations(:)
        integer, intent(out)                :: status

        logical                             :: exist
        integer                             :: nfiles, ifile

        status = 1

        ! determine the number of filter files
        nfiles = 0
        do

            inquire (file=filename // '.' // strinteger(nfiles), exist=exist)
            if (.not. exist) exit
            nfiles = nfiles + 1

        end do

        ! some checks
        if (nfiles == 0) then
            write (ERROR_UNIT,'(a)') "Error: Failed to open filter '" // filename // "'."
            return
        end if

        if (mod(nfiles, ndetectors) /= 0) then
           write (ERROR_UNIT,'(a,2(i0,a))') "Error: In file '" // filename // "', the number of filters '", nfiles, "' is not an in&
                 &teger time the number of detectors '", ndetectors, "'."
           return
        end if

        allocate (first(nfiles))
        allocate (last(nfiles))
        allocate (ncorrelations(nfiles))

        ! loop over the filters
        do ifile = 1, nfiles
            call read_filter_header(filename // '.' // strinteger(ifile-1), convert, first(ifile), last(ifile),                    &
                 ncorrelations(ifile), status)
            if (status /= 0) return
        end do
        
        status = 0
        
    end subroutine read_filter_headers


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine read_filter(filename, convert, ndetectors, filter, status)

        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: convert
        integer, intent(in)          :: ndetectors
        type(filterset), intent(out) :: filter
        integer, intent(out)         :: status

        integer*8, allocatable :: first(:), last(:)
        integer, allocatable   :: ncorrelations(:)
        integer*8              :: nsamples, isample
        integer                :: islice, idetector, ifilter

        call read_filter_headers(filename, convert, ndetectors, first, last, ncorrelations, status)
        if (status /= 0) return

        ! set number of correlations and bandwidth
        if (any(ncorrelations /= ncorrelations(1))) then
            write (ERROR_UNIT,'(a)') 'Error: Filter files do not have the same correlation length.'
            return
        end if        
        filter%ncorrelations = ncorrelations(1)
        filter%bandwidth = 2*ncorrelations(1) + 1

        ! set number of detectors
        filter%ndetectors = ndetectors

        ! set number of slices, update ranges and read invntt values
        filter%nslices = size(first) / ndetectors

        allocate(filter%first(filter%nslices))
        allocate(filter%last (filter%nslices))
        allocate(filter%data (filter%ncorrelations,filter%ndetectors,filter%nslices))
        isample = 1

        do islice=1, filter%nslices

           nsamples = last((islice-1)*ndetectors+1) - first((islice-1)*ndetectors+1) + 1

           do idetector = 1, ndetectors

               ifilter = (islice-1)*ndetectors + idetector
               if (last(ifilter) - first(ifilter) + 1 /= nsamples) then
                   write (ERROR_UNIT,'(a)') "Error: Filter '" // filename // '.' // strinteger(ifilter-1)// "' does not apply to th&
                         &e same number of samples than filter '" // filename // '.' // strinteger((islice-1)*ndetectors) // "'."
                   status = 1
                   deallocate (filter%first)
                   deallocate (filter%last)
                   deallocate (filter%data)
                   return
               end if

               open (11, file=filename // '.' // strinteger(ifilter-1), form='unformatted', access='stream', convert=convert,      &
                     action='read')
               read (11, pos=3*8+1) filter%data(:,idetector,islice)
               close (11)

           end do

           filter%first(islice) = isample
           filter%last (islice) = isample + nsamples - 1
           isample = isample + nsamples

        end do

    end subroutine read_filter


end module module_madcap