program test_euclidtrapping2

    use module_euclidtrapping
    use module_fitstools, only : ft_write_image
    use module_precision, only : p => dp

    integer  :: imain(nlines)
    real(p),  allocatable :: time(:)
    integer,  allocatable :: tirage_photons(:,:)
    integer, allocatable :: imaout(:,:)
    integer :: status

    !call run_tirage_proba(ntirages)
    
    imain = 94
    call simul_integration(imain, time, tirage_photons, imaout)

    print *, 'IMAOUT FINAL:', imaout(:,size(imaout,2))

    !print *, 'coord:', trapsin(1,1)%coord(1:6,1)
    !print *, 'species:', trapsin(1,1)%species(1:20)
    !print *, 'state:', trapsin(1,1)%state(1:20)
    !print *, 'sig:', trapsin(1,1)%sig
    !print *, 'tau_r:', trapsin(1,1)%tau_r
    !print *, 'tau_c:', trapsin(1,1)%tau_c(:)
    !print *, 'proba_r:', trapsin(1,1)%proba_r
    !print *, 'proba_c:', trapsin(1,1)%proba_c
    !print *, 'nv:', trapsin(1,1)%nv(:)
    !print *, 'nc:', trapsin(1,1)%nc(:)

    !call ft_write_image('/home/pchanial/work/euclid/simul_trapping/simul_trapping_proba_integration_full_nv.fits', &
    !     trapsout(1,1)%nv(3,:), status=status)

end program test_euclidtrapping2
