program test_euclidtrapping

    use module_euclidtrapping
    use module_fitstools, only : ft_write_image
    use module_precision, only : p => dp

    integer, parameter :: ntirages = 100000
    integer  :: imain(nlines, ncolumns)
    real(p),  allocatable :: time(:)
    integer,  allocatable :: tirage_photons(:,:,:)
    real(sp), allocatable :: coord(:,:)
    type(TrapIn), allocatable :: trapsin(:,:)
    integer, allocatable :: imaout(:,:,:)
    type(Trap), allocatable :: trapsout(:,:)
    real(p), allocatable :: conf_z(:,:,:)
    integer :: status

    !call run_tirage_proba(ntirages)
    
    imain = 94
    call simul_trapping_proba_integration_full(imain, time, tirage_photons, coord, trapsin, imaout, trapsout, conf_z, status)

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

contains

    subroutine run_tirage_proba(ntirages)

        integer, intent(in) :: ntirages
        
        integer :: output(ntirages, nlines, ncolumns)
        real(p) :: proba(nlines, ncolumns)

        print *, 'starting with ', shape(output)
        proba = 0.5_p
        call tirage_proba(proba, ntirages, output)

    end subroutine run_tirage_proba

end program test_euclidtrapping
