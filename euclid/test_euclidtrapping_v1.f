program test_euclidtrapping

    use module_euclidtrapping
    use module_euclidtrapping_v1
    use module_fitstools, only : ft_write_image
    implicit none

    integer, parameter         :: ntirages = 100000
    integer                    :: imain(nlines, ncolumns)
    real(p),  allocatable      :: time(:)
    integer,  allocatable      :: tirage_photons(:,:,:)
    real(sp), allocatable      :: coord(:,:)
    type(TrapIn), allocatable  :: trapsin(:,:)
    type(TrapOut), allocatable :: trapsout(:,:)
    integer, allocatable       :: imaout(:,:,:)
    real(p), allocatable       :: conf_z(:,:,:)
    integer :: status

    !call run_tirage_proba(ntirages)
    
    imain = 94
    call simul_integration_v1(imain, time, tirage_photons, coord, trapsin, imaout, trapsout, conf_z, status)

    print *, 'IMAOUT FINAL:', imaout(:,1,size(imaout,3))

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

end program test_euclidtrapping
