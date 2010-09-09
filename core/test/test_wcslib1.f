program test_wcslib1

    use, intrinsic :: ISO_C_BINDING
    use module_wcslib
    implicit none

    INTEGER   ALTS(0:26), CTRL, I, IERR, J, K, NKEYRC, NREJECT, NWCS, RELAX, WCSP
    CHARACTER CALTS(0:26)*2, KEYREC*80, HEADER*288001
    character(len=*), parameter :: infile = 'core/test/data/pih.fits'

    INTEGER WCS(WCSLEN), STAT(WCSFIX_NWCS)

!-----------------------------------------------------------------------
      WRITE (*, 10)
 10   FORMAT (&
     &  'Testing WCSLIB parser for FITS image headers (tpih1.f)',/,&
     &  '------------------------------------------------------',/)

!     Open the FITS WCS test header for formatted, direct I/O.
      OPEN (UNIT=1, FILE=INFILE, FORM='FORMATTED', ACCESS='DIRECT',&
     &      RECL=80, IOSTAT=IERR)
      IF (IERR.NE.0) THEN
        WRITE (*, 20) IERR, INFILE
 20     FORMAT ('ERROR',I3,' opening ',A)
        stop
      END IF

!     Read in the FITS header, excluding COMMENT and HISTORY keyrecords.
      K = 1
      NKEYRC = 0
      outer: DO J = 0, 100
        DO I = 1, 36
          READ (1, '(A80)', REC=36*J+I, IOSTAT=IERR) KEYREC
          IF (IERR.NE.0) THEN
            WRITE (*, 30) IERR
 30         FORMAT ('ERROR',I3,' reading header.')
            stop
          END IF

          IF (KEYREC(:8).EQ.'        ') cycle
          IF (KEYREC(:8).EQ.'COMMENT ') cycle
          IF (KEYREC(:8).EQ.'HISTORY ') cycle

          HEADER(K:) = KEYREC
          K = K + 80
          NKEYRC = NKEYRC + 1

!         An END keyrecord was read, read the rest of the block.
          IF (KEYREC(:8).EQ.'END     ') exit outer
         end do

      end do outer

      CLOSE (UNIT=1)

      HEADER(K:K) = CHAR (0)
      WRITE (*, 70) NKEYRC
 70   FORMAT ('Found',I4,' non-comment header keyrecords.',/)


!     Cull all WCS keyrecords from the header but report illegal ones.
      WRITE (*, 80)
 80   FORMAT (/,'Illegal-WCS header keyrecords rejected by wcspih():')
      RELAX = WCSHDR_all
      CTRL = -2
!     WCSPIH will allocate memory for NWCS intialized WCSPRM structs.
      IERR = WCSPIH(HEADER, NKEYRC, RELAX, CTRL, NREJECT, NWCS, WCSP)
      IF (IERR.NE.0) THEN
        WRITE (*, 90) IERR
 90     FORMAT ('WCSPIH ERROR',I2,'.')
        stop
      END IF

!     List keyrecords that were not consumed by WCSPIH.
      WRITE (*, 100)
 100  FORMAT (//,'Non-WCS header keyrecords not used by WCSPIH:')
      DO I = 1, 288001, 80
        IF (HEADER(I:I).EQ.CHAR(0)) exit
        WRITE (*, '(A)') HEADER(I:I+79)
      end do

      IERR = WCSIDX (NWCS, WCSP, ALTS)
      WRITE (*, 130)
 130  FORMAT (//,'Index of alternate coordinate descriptions found:',/,&
     &        '   A B C D E F G H I J K L M N O P Q R S T U V W X Y Z')
      DO I = 0, 26
        IF (ALTS(I).LT.0) THEN
          CALTS(I) = ' -'
        ELSE
          WRITE (CALTS(I), '(I2)') ALTS(I)
        END IF
      end do
      WRITE (*, '(27A)') CALTS

      DO I = 0, NWCS-1
        WRITE (*, 150)
 150    FORMAT (/,'------------------------------------', '------------------------------------')

!       Copy into our WCSPRM struct.
        IERR = WCSVCOPY (WCSP, I, WCS)

!       Fix non-standard WCS keyvalues.
        IERR = WCSFIX (ctrl=7, naxis=c_null_ptr, wcs=WCS, stat=STAT)
        IF (IERR.NE.0) THEN
            WRITE (*, 160) (STAT(J), J=1,WCSFIX_NWCS)
 160        FORMAT ('WCSFIX ERROR, status returns: (',(I2,:,','),')',/)
        END IF

        IERR = WCSSET (WCS)
        IF (IERR.NE.0) THEN
            WRITE (*, 170) IERR
 170        FORMAT ('WCSSET ERROR',I2,'.')
            cycle
        END IF

        IERR = WCSPRT (WCS)
        IF (IERR.NE.0) THEN
            WRITE (*, 180) IERR
 180        FORMAT ('WCSPRT ERROR',I2,'.')
            cycle
        END IF

!       Free memory (doesn't free memory allocated by WCSPIH).
        IERR = WCSFREE (WCS)
      end do

!     Free the memory allocated by WCSPIH.
      IERR = WCSVFREE (NWCS, WCSP)

      stop 'OK.'

end program test_wcslib1
