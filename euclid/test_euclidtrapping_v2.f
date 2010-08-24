program test_euclidtrapping_v2

    use module_euclidtrapping
    use module_euclidtrapping_v2
    use module_fitstools, only : ft_write_image
    implicit none

    integer               :: imain(nlines)
    real(p),  allocatable :: time(:)
    integer,  allocatable :: tirage_photons(:,:)
    integer, allocatable  :: imaout(:,:)
    integer               :: status

    imain = 94
    call simul_integration_v2(imain, time, tirage_photons, imaout)

    print *, 'IMAOUT FINAL:', imaout(:,size(imaout,2))

    !call ft_write_image('/home/pchanial/work/euclid/simul_trapping/simul_trapping_proba_integration_full_nv.fits', &
    !     imaout, status=status)

end program test_euclidtrapping_v2
