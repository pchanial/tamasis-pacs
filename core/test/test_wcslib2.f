program test_wcslib2

    use iso_c_binding
    use module_math,      only : neq_real
    use module_fitstools, only : ft_header2str
    use module_wcslib

    implicit none

    integer, parameter          :: ncoords = 3
    character(len=*), parameter :: filename = 'core/test/data/header.fits'
    character(len=2880)         :: header
    integer                     :: wcs(WCSLEN), statfix(WCSFIX_NWCS), statproj(ncoords)
    integer                     :: alts(0:26)
    character(len=2)            :: calts(0:26)
    integer                     :: i, status, relax, ctrl, nreject, nwcs, wcsp

    real*8                      :: ad(2*ncoords), xy(2*ncoords), imgcrd(2*ncoords), ra_
    real*8                      :: phi(ncoords), theta(ncoords)

    call ft_header2str(filename, header, status)
    if (status /= 0) stop 'ft_header2str: FAILED.'

    relax = WCSHDR_all
    ctrl  = 0
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

    if (any(neq_real(xy, [-6538.19790689358d0,12196.226642292400d0, -4998.77434749560d0,114.620164770842d0,-4694.59388530175d0,    &
        -2272.629672963080d0], 1d-12))) then
        stop 'wcss2p: Wrong result.'
    end if

    ! some cleanup
    status = wcsfree(wcs)
    status = wcsvfree(nwcs, wcsp)

    stop "OK."

end program test_wcslib2
