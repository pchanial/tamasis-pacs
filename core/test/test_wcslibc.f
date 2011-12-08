program test_wcslibc

    use iso_c_binding
    use module_cfitsio
    use module_stdio
    use module_string, only : strjoin
    use module_wcslibc
    implicit none

    character(len=*), parameter       :: filename_header = 'core/test/data/pih.fits'
    integer                           :: i, iwcs, ifix
    integer(kind=C_INT)               :: alts(27), stat(NWCSFIX)
    integer(kind=C_INT)               :: nkeyrec, relax, ctrl, nreject, nwcs, status
    type(C_PTR)                       :: c_wcs, c_header, fptr
    character(len=80*1000+1), pointer :: header => null()
    type(wcsprm), pointer             :: wcs(:) => null()

    ! not working...
    stop

    status = 0

    call init_stdio()
    status = 0

    call fits_open_file(fptr, filename_header // C_NULL_CHAR, CFITSIO_READONLY, status)
    if (status /= 0) then
        call fits_report_error(stderr, status)
        stop
    end if

    call fits_hdr2str(fptr, 1_C_INT, c_null_ptr, 0_C_INT, c_header, nkeyrec, status)
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

    relax = WCSHDR_all
    ctrl = -2_C_INT

    status = wcspih(header, nkeyrec, relax, ctrl, nreject, nwcs, c_wcs)

    call c_f_pointer(c_wcs, wcs, shape = [nwcs])

    !/* Summarize what was found. */
    !status = wcsidx(nwcs, &wcs, alts);
    status = wcsidx(nwcs, c_wcs, alts);
    !printf("\n\nFound %d alternate coordinate descriptions with indices:\n  ",
    !       nwcs);
    write (*,'(a,i2,a)') 'Found ', nwcs, ' alternate coordinate descriptions with indices:'
    !for (a = 'A'; a <= 'Z'; a++) {
    !  printf("%2c", a);
    !}
    write (*,'(a)') '   A B C D E F G H I J K L M N O P Q R S T U V W X Y Z'

    !for (ialt = 0; ialt < 27; ialt++) {
    !  if (alts[ialt] < 0) {
    !    printf(" -");
    !  } else {
    !    printf("%2d", alts[ialt]);
    !  }
    !}
    do i=1, 27
       if (alts(i) < 0) then
           write (*,'(a,$)') ' -'
       else
           write (*,'(i2,$)') alts(i)
       end if
    end do
    write (*,*)

    !/* Fix non-standard usage and print each of the wcsprm structs. */
    !for (iwcs = 0; iwcs < nwcs; iwcs++) {
    !  printf("\n------------------------------------"
    !         "------------------------------------\n");
    !
    !  /* Fix non-standard WCS keyvalues. */
    !  if ((status = wcsfix(7, 0, wcs+iwcs, stat))) {
    !    printf("wcsfix ERROR, status returns: (");
    !    for (ifix = 0; ifix < NWCSFIX; ifix++) {
    !      printf(ifix ? ", %d" : "%d", stat[ifix]);
    !    }
    !    printf(")\n\n");
    !  }
    !

  do iwcs = 1, nwcs
      write (*,*)
      write (*,'(a)') '------------------------------------------------------------------------'
      ! fix non-standard WCS keyvalues
      status = wcsfix(7_C_INT, c_null_ptr, c_loc(wcs(iwcs)), stat)
      if (status /= 0) then
          write (*,'(a,$)') "wcsfix ERROR, status returns: ("
          do ifix = 1, NWCSFIX
              if (ifix > 1) write (*,'(a,$)') ', '
              write (*,'(i3,$)') stat(ifix)
          end do
          write (*,*) ')'
      end if

      !  if ((status = wcsset(wcs+iwcs))) {
      !    fprintf(stderr, "wcsset ERROR %d: %s.\n", status, wcs_errmsg[status]);
      !    continue;
      !  }
      !
      status = wcsset(c_loc(wcs(iwcs)))
      if (status /= 0) then
          write (*,'(a,i3)') 'wcsset ERROR ', status
          cycle
      end if

      !  if ((status = wcsprt(wcs+iwcs))) {
      !    fprintf(stderr, "wcsprt ERROR %d: %s.\n", status, wcs_errmsg[status]);
      !  }
      status = wcsprt(c_loc(wcs(iwcs)))
      if (status /= 0) then
          write (*,'(a,i3)') 'wcsprt ERROR ', status
      end if

  end do

  !status = wcsvfree(&nwcs, &wcs);

  stop 'OK.'

end program test_wcslibc
