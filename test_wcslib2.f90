program test_wcslib2

    use, intrinsic :: ISO_C_BINDING
    use module_wcslib
    use module_fitstools

    implicit none

    integer, parameter          :: ncoords = 3
    character(len=*), parameter :: filename = '/home/pchanial/work/pacs/tamasis/csh_header_ngc6946.fits'
    character(len=2880)         :: header
    integer                     :: wcs(WCSLEN), statfix(WCSFIX_NWCS), statproj(ncoords)
    integer                     :: alts(0:26)
    character(len=2)            :: calts(0:26)
    integer                     :: i, status, relax, ctrl, nreject, nwcs, wcsp

    real*8                      :: ad(2*ncoords), xy(2*ncoords), imgcrd(2*ncoords), ra_
    real*8                      :: phi(ncoords), theta(ncoords)

    status = 0
    call ft_header2str(filename, header, status)
    call ft_printerror(status, filename)

    relax = WCSHDR_all
    ctrl  = 0
    write(*,*) header
    status = wcspih(header, len(header)/80, relax, ctrl, nreject, nwcs, wcsp)

    write (*,*) 'nwcs:', nwcs

    ! print the alternate representation, if any
    status = WCSIDX (NWCS, WCSP, ALTS)
    write (*, '(/,a,/,a)') 'Index of alternate coordinate descriptions found:', &
                           '   A B C D E F G H I J K L M N O P Q R S T U V W X Y Z'
    do i = 0, 26
        if (alts(i) < 0) then
            calts(i) = ' -'
        else
            write (calts(i), '(i2)') alts(i)
        end if
    end do
    write (*, '(27a)') calts

    status = wcsvcopy(wcsp, 0, wcs)

    !Fix non-standard WCS keyvalues.
    status = WCSFIX(ctrl=7, naxis=c_null_ptr, wcs=WCS, stat=statfix)
    IF (status /= 0) THEN
        WRITE (*, 160) (statfix(i), i=1,WCSFIX_NWCS)
160     FORMAT ('WCSFIX ERROR, status returns: (',(I2,:,','),')',/)
        stop
    END IF

    status = WCSSET(WCS)
    IF (status /= 0) THEN
        WRITE (*, 170) status
170     FORMAT('WCSSET ERROR',I2,'.')
        stop
    END IF

    ra_ = 300d0
    ad = [ra_, 50d0, ra_, 60d0, ra_, 62d0]

    status = wcss2p(wcs, ncoords, 2, ad, phi, theta, imgcrd, xy, statproj)
    if (status /= 0) then
        write (*,*) 'Error in wcss2p: ', status
    else
        do i=1, ncoords
            write (*,*) ad(i*2-1), ad(i*2), ':', xy(i*2-1), xy(i*2)
        end do
    end if

    !IDL> forprint, ra, dec, x+1, y+1
    !   300.00000       50.000000      -9113.9278       10412.664
    !   300.00000       60.000000      -4896.1894      -1012.9964
    !   300.00000       62.000000      -4062.7907      -3270.6355

    ! some cleanup
    status = wcsfree(wcs)
    status = wcsvfree(nwcs, wcsp)

    stop "OK."

end program test_wcslib2
