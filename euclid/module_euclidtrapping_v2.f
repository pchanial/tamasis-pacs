!d'apres les equations de alex short
!on calcule le trapping-detrapping qui intervient dans un pixel pendant l'integration puis pendant la lecture
!approche probabiliste
!imain: image a processer
!tirage_photons: tirage aleatoire de l'arrivee des photons
!traps: structure contenant la description de chaque piege
!t_ccd: temperature de fonctionnement
!int_time: temps d'integration (s)
!time_int: intervalle de temps (s)
!imaout: image de sortie apres calcul de la dynamique des pieges
!tau_c: constantes de temp de capture des pieges
!proba_c: probabilites de capture
!proba_r: probabilites de relache

!nv: nombre de pieges vacants
!nc: nombre de pieges occupes
!nt nombre total de pieges
!nt = nv+nc
!dynamique de capture et relache d'un electron
!dnc/dt = +-nc/tau_r
!dnv/dt = +-nv/tau_c
!tau_r: constante de temps de relache d'un electron piege
!tau_c: constante de temps de capture d'un electron libre
!Ne:nombre d'electrons libres
!dynamique des electrons libres
!dNe/dt = nc/tau_r-nv/tau_c
!Ne+Nc = cte

module module_euclidtrapping_v2

    use module_euclidtrapping
    use module_fitstools, only : ft_read_image
    use module_math,      only : PI, linspace
    use module_sort,      only : qsorti
    implicit none

    character, parameter :: mode = 'd'

    type Detector
        integer :: nv
        integer :: nc
        logical :: filled(ntraps)
    end type Detector

    type DetectorColumn
        integer        :: ndetectors
        integer        :: ntraps
        real(p)        :: vth
        real(p)        :: z(ntraps)
        real(p)        :: sig(ntraps)
        real(p)        :: proba_r(ntraps)
        type(Detector) :: detector(nlines)
    end type DetectorColumn

    real(p) :: vth, nstates
    real(p), dimension(ntraptypes) :: sig, tau_r, e_traps
    integer :: nt(ntraptypes) !nt: nombre de pieges par volume de confinement (reel)
    integer :: tottraps, isig(ntraptypes)
    type(DetectorColumn) :: column
    

contains


    subroutine init_v2

        real(p), allocatable :: coord_(:,:)
        integer              :: itrap, idetector, status
        integer              :: nt_debut(ntraptypes), nt_fin(ntraptypes)

        !section efficace de capture
        sig = [1.092e-18_p, 1.179e-17_p, 1.029e-19_p]*1.e-4_p !m2
        !constante de temps d'emission des electrons pieges
        tau_r = [5.913e-4_p, 7.623_p, 15.282_p] !s

        !nt: nombre de pieges par volume de confinement (reel)
        nt = nint(nt_gaia*(fwc_euclid/fwc_gaia))

        !scaling des constantes de temps avec la temperature
        !d'apres a. short varie en exp(Et/kT)/T^2 (eq. 4 et 5)
        !vth: velocite thermique des electron (cm/s), d'apres alex short
        vth     = sqrt(3_p*k*t_ccd/me) !m/s
        nstates = 2*(2*PI*me*k*t_ccd/(h**2))**1.5_p !m-3
        e_traps = k*t_ccd*log(tau_r*sig*vth*nstates) !J
        tau_r   = exp(e_traps/(k*t_ccd))/(sig*vth*nstates)
        tottraps = sum(nt)

        !il faut calculer la dynamique des pieges en commencant par les pieges de constantes de temps de capture la plus courte
        !donc la section efficace de capture la plus grande donc il faut ordonner les pieges par sig decroissant
        call qsorti(sig, isig)
        isig = isig(ntraptypes:1:-1)
        sig = sig(isig)
        tau_r = tau_r(isig)
        nt = nt(isig)
        
        !indices des pieges des differentes familles
        nt_debut(1) = 1
        nt_fin(1) = nt(1)
        do itrap = 2, ntraptypes
            nt_debut(itrap) = nt_fin(itrap-1) + 1
            nt_fin(itrap) = nt_fin(itrap-1) + nt(itrap)
        end do

        column%ndetectors = nlines
        column%ntraps = sum(nt)
        column%vth = vth

        !generation des pieges dans le volume du buried channel
        call ft_read_image('~/work/tamasis/tamasis-unstable/euclid/data/coord_traps_integ_proba_comp_94.fits', coord_,status=status)
        if (status /= 0) return
        column%z = real(coord_(:,3), kind=p)

        call fill_trap(column%sig, sig, nt_debut, nt_fin)
        call fill_trap(column%proba_r, 1._p - exp(-time_int / tau_r), nt_debut, nt_fin)

        ! we start with empty traps
        column%detector%nv = sum(nt)
        column%detector%nc = 0
        do idetector = 1, column%ndetectors
            column%detector(idetector)%filled = .false.
        end do

    end subroutine init_v2


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine fill_trap(dest, source, start, end)
        real(p), intent(out) :: dest(:)
        real(p), intent(in)  :: source(:)
        integer, intent(in)  :: start(:), end(:)

        integer :: itype

        do itype = 1, size(source)
            dest(start(itype):end(itype)) = source(itype)
        end do

    end subroutine


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine simul_integration(imain, time, tirage_photons, imaout)

        integer, intent(in)  :: imain(nlines)
        real(p), intent(out), allocatable :: time(:)
        integer, intent(out), allocatable :: tirage_photons(:,:)
        integer, intent(out), allocatable :: imaout(:,:)

        real(p) :: model
        integer :: nsteps, i

        integer :: flux(nlines)
        real(p), dimension(ntraps) :: alp, tau_c, proba_c
        real(p), dimension(nlines) :: conf_vol, conf_z
        real(p) :: rnd0
        integer :: idetector, itrap, ncaptures, nreleases

        integer :: count1, count2, count_rate, count_max
        
        call init_v2

        !calcul du nombre d'e- rajoutes par effet photo-electrique a chaque pas de temps
        !generation de la sequence temporelle d'arrivee des photons pendant l'integration
        !nombre de pas d'integration
        nsteps = nint(int_time/time_int)
        print *, 'nsteps:', nsteps
        allocate (time(nsteps+1))
        time = linspace(0._p, nsteps * time_int, nsteps+1)
        flux = imain / nsteps !XXX check me

        allocate (tirage_photons(nlines, nsteps))
        call tirage_proba_slow2(imain/real(nsteps, kind=p)-flux, nsteps, tirage_photons)

        do i=1, nsteps
            tirage_photons(:,i) = tirage_photons(:,i) + flux
        end do
        
        allocate (imaout(column%ndetectors, nsteps+1))
        imaout = 0
        
        !choix du modele density driven ou volume driven
        if (mode == 'v') then
            model = 1._p/3._p
        else if (mode == 'd') then
            model = 0_p
        else
            stop 'mode should be either d or v'
        end if
        
        !integration de l'image
        call system_clock(count1, count_rate, count_max)
        
        do i=1, nsteps

            if (i == (i/1000)*1000) print *, 'step:', i

            !calcul du nombre d'electrons rajoutes a chaque pas d'integration
            imaout(:,i) = imaout(:,i) + tirage_photons(:,i)

            !calcul de la profondeur du volume de confinement des electrons libres
            !origine en profondeur a z=0
            conf_z = (imaout(:,i) / fwc_euclid)**model

            !calcul du volume de confinement des electrons libres
            conf_vol = conf_z * bc_vol

            do idetector = 1, column%ndetectors

                ! Capture
                if (imaout(idetector,i) /= 0 .and. column%detector(idetector)%nv > 0) then
                    
                    !constante de temps de capture
                    alp = (column%sig*vth) / conf_vol(idetector)
                    tau_c = 1_p / (alp * imaout(idetector,i))

                    !probabilite de capture
                    proba_c = 1_p - exp(-time_int / tau_c)

                    !on identifie les pieges vides de chaque pixel
                    ncaptures = 0
                    do itrap = 1, column%ntraps
                        if (column%z(itrap) <= conf_z(idetector) .and. .not. column%detector(idetector)%filled(itrap)) then
                            call random_number(rnd0)
                            if (floor(rnd0 + proba_c(itrap)) == 1) then
                                column%detector(idetector)%filled(itrap) = .true.
                                ncaptures = ncaptures + 1
                                if (ncaptures == imaout(idetector,i)) exit
                            end if
                        end if
                    end do
                    imaout(idetector,i) = imaout(idetector,i) - ncaptures

                end if

                ! Relache
                nreleases = 0
                do itrap = 1, column%ntraps
                    if (column%detector(idetector)%filled(itrap)) then
                        call random_number(rnd0)
                        column%detector(idetector)%filled(itrap) = floor(rnd0+column%proba_r(itrap)) == 0
                        if (.not. column%detector(idetector)%filled(itrap)) then
                            nreleases = nreleases + 1
                        end if
                    end if
                end do
                imaout(idetector,i) = imaout(idetector,i) + nreleases
                
                column%detector(idetector)%nc = column%detector(idetector)%nc - nreleases
                column%detector(idetector)%nv = column%ntraps - column%detector(idetector)%nc

            end do

            imaout(:,i+1)=imaout(:,i)
            
        end do

        call system_clock(count2, count_rate, count_max)
        write (*,'(a,f7.2,a)') 'proba_integration_full: ', real(count2-count1)/count_rate, 's'

    end subroutine simul_integration


end module module_euclidtrapping_v2
