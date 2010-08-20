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

module module_euclidtrapping

    use module_math, only : PI, linspace
    use module_sort, only : qsorti, where
    use module_precision, only : p => dp, sp
    implicit none

    real(p), parameter   :: t_ccd = 163.d0
    real(p), parameter   :: time_par = 1.d-2
    character, parameter :: mode = 'd'

    integer, parameter   :: npixal = 1024
    integer, parameter   :: npixac = 1
    integer, parameter   :: nlines = npixal
    integer, parameter   :: ncolumns = npixac
    integer, parameter   :: ntraps = 3 !nombre de familles de pieges considerees
    real(p), parameter   :: k = 1.381e-23_p !J/K
    real(p), parameter   :: me = 0.5_p * 9e-31_p !kg ! WHAT'S THAT?
    real(p), parameter   :: h = 6.626e-34_p !J.s

    !parametres des traps deduits de gaia
    !sig: section efficace de capture (cm2)
    !tau_r: constante de temps de relache des pieges (s)
    !e_vol: volume de confinement des electrons dans le pixel
    !1.26e-9 cm3 d'apres short, 1.296e-10 cm3 d'apres seabroke (volume du BC 24x10x0.54 µm3)
    !scaling BC euclid: FWC GAIA 200ke-, FWC Euclid 250ke-
    !d'ou BC euclid=BC gaia*1.25 d'ou volume BC euclid=volume BC gaia*1.25
    !or surface BC euclid=144 µm2 compare a surface BC gaia 240 µm2
    !d'ou profondeur BC euclid=1.25*240/144 profondeur BC gaia = 1.125 µm?
    real(p), parameter :: fwc_gaia   = 2.0e5_p
    real(p), parameter :: fwc_euclid = 2.5e5_p
    real(p), parameter :: bc_vol     = 1.26e-9_p * (fwc_euclid/fwc_gaia) * 1.e-6_p !m3

    !familles de pieges: 0=divacancy, 1=?, 2=E-centre
    !d'apres alex short, best fit to gaia test data at 163K
    !nombre de pieges par pixel gaia par famille
    real(p), parameter :: nt_gaia(3) = [23.0_p, 1.6421_p, 110.332_p]
    !nombre de pieges par volume buried channel euclid
    integer, parameter :: nt_(3) = nint(nt_gaia*(fwc_euclid/fwc_gaia))
#ifdef IFORT
    integer, parameter :: tottraps_ = nt_(1) + nt_(2) + nt_(3)
#else
    integer, parameter :: tottraps_ = sum(nt_)
#endif
    integer, parameter :: npar = 10
    
    ! pieges mode stochastique
    type TrapIn
        real(sp) :: coord(tottraps_, 3)
        integer  :: species(tottraps_)
        integer  :: state(tottraps_)
        real(p)  :: sig(ntraps)
        real(p)  :: tau_r(ntraps)
        real(p)  :: tau_c(ntraps)
        real(p)  :: proba_r(ntraps)
        real(p)  :: proba_c(ntraps)
        integer  :: nv(ntraps)
        integer  :: nc(ntraps)
    end type TrapIn
        
    type Trap
        real(sp) :: coord(tottraps_, 3)
        integer  :: species(tottraps_)
        integer  :: state(tottraps_, npar+1)
        real(p)  :: sig(ntraps)
        real(p)  :: tau_r(ntraps)
        real(p)  :: tau_c(ntraps, npar)
        real(p)  :: proba_r(ntraps)
        real(p)  :: proba_c(ntraps, npar)
        integer  :: nv(ntraps, npar+1)
        integer  :: nc(ntraps, npar+1)
    end type Trap

    real(p) :: vth, nstates
    real(p), dimension(ntraps) :: sig, tau_r, e_traps
    integer :: nt(ntraps) !nt: nombre de pieges par volume de confinement (reel)
    integer :: tottraps
    integer :: nt_debut(ntraps), nt_fin(ntraps)


contains


    subroutine init

        integer :: itrap

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
        sig = sig(ntraps:1:-1)
        tau_r = tau_r(ntraps:1:-1)
        nt = nt(ntraps:1:-1)
        
        !indices des pieges des differentes familles
        nt_debut(1) = 1
        nt_fin(1) = nt(1)
        do itrap = 2, ntraps
            nt_debut(itrap) = nt_fin(itrap-1) + 1
            nt_fin(itrap) = nt_fin(itrap-1) + nt(itrap)
        end do

    end subroutine init


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine tirage_proba(proba, ntirages, tirage)
        !procedure de tirage d'un evenement ayant une certaine probabilite d'arriver
        !proba: probabilite de l'evenement
        !ntirages: nombre de tirages
        !on cree un vecteur de "ntirages" tirages compose de 1 et 0 avec un rapport 1/0 egal a "proba"
        !on reordonne aleatoirement ce vecteur
        
        implicit none
        
        real(p), intent(in) :: proba(:,:)
        integer, intent(in) :: ntirages
        integer, intent(out) :: tirage(ntirages,size(proba,1),size(proba,2))

        integer :: n_un(size(proba,1),size(proba,2)), ordre(ntirages)
        integer :: icolumn, iline
        real(p) :: rnd(ntirages)
        integer :: count1, count2, count_rate, count_max

        tirage = 0
        n_un = nint(proba * ntirages)
 
        call system_clock(count1, count_rate, count_max)
        call random_number(rnd)

        call qsorti(rnd, ordre)

        do icolumn = 1, ncolumns
            do iline = 1, nlines
                tirage(ordre(1:n_un(iline,icolumn)),iline,icolumn) = 1
            end do
        end do

        call system_clock(count2, count_rate, count_max)
        write (*,'(a,f7.2,a)') 'tirage_proba: ', real(count2-count1)/count_rate, 's'

    end subroutine tirage_proba

    subroutine tirage_proba_slow(proba, ntirages, tirage)
        !procedure de tirage d'un evenement ayant une certaine probabilite d'arriver
        !proba: probabilite de l'evenement
        !ntirages: nombre de tirages
        !on cree un vecteur de "ntirages" tirages compose de 1 et 0 avec un rapport 1/0 egal a "proba"
        !on reordonne aleatoirement ce vecteur
        
        implicit none
        
        real(p), intent(in) :: proba(:,:)
        integer, intent(in) :: ntirages
        integer, intent(out) :: tirage(size(proba,1),size(proba,2),ntirages)

        integer :: n_un(size(proba,1),size(proba,2)), ordre(ntirages)
        integer :: icolumn, iline
        real(p) :: rnd(ntirages)
        integer :: count1, count2, count_rate, count_max

        tirage = 0
        n_un = nint(proba * ntirages)
 
        call system_clock(count1, count_rate, count_max)
        call random_number(rnd)

        call qsorti(rnd, ordre)

        do icolumn = 1, ncolumns
            do iline = 1, nlines
                tirage(iline,icolumn,ordre(1:n_un(iline,icolumn))) = 1
            end do
        end do

        call system_clock(count2, count_rate, count_max)
        write (*,'(a,f7.2,a)') 'tirage_proba: ', real(count2-count1)/count_rate, 's'

    end subroutine tirage_proba_slow


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine simul_trapping_proba_integration_full(imain,time,tirage_photons,coord,trapsin,imaout,trapsout,conf_z,t_ccd=t_ccd,int_time=int_time,time_int=time_int, status)

        implicit none

        integer, intent(in)  :: imain(nlines, ncolumns)
        type(Trap), intent(out), allocatable :: trapsout(:,:)
        integer, intent(out) :: status

        type(TrapIn), allocatable :: trapsin(:,:)
        type(TrapIn)              :: trapin_
        real(sp), allocatable :: coord(:,:)
        real(p), allocatable :: time(:), conf_vol(:), conf_z(:)
        integer, allocatable :: imaout(:,:,:), tirage_photon(:,:,:)
        real(p) :: int_time, time_int
        integer :: nstep
        integer :: count1, count2, count_rate, count_max

        
        !parametres par defaut
        !if keyword_set(int_time) then int_time=int_time else int_time=450.
        !if keyword_set(time_int) then time_int=time_int else time_int=1.e-1
        !if keyword_set(init) then init=init else init=0
        print *, 'int_time: ',int_time
        print *, 'time_int: ',time_int

        call init

        !calcul du nombre d'e- rajoutes par effet photo-electrique a chaque pas de temps
        !generation de la sequence temporelle d'arrivee des photons pendant l'integration
        !nombre de pas d'integration
        nstep = nint(int_time/time_int)
        allocate (time(nstep+1))
        time = linspace(0._p, nstep * time_int, nstep+1)
        flux = imain / nstep

        allocate (tirage_photons(nlines, ncolumns, nstep))
        call tirage_proba_slow(imain/float(nstep)-flux, nstep, tirage_photons)
        do i=1, nstep
            tirage_photons(:,:,i) = tirage_photons(:,:,i)+flux
        end do
        
        allocate (imaout(nlines, ncolumns, nstep+1))
        
        !generation des pieges dans le volume du buried channel
        call ft_read_image('coord_traps_integ_proba_comp_94.fits', coord, status=status)
        if (status /= 0) return

        !coordonnees normalisees par rapport aux dimensions du buried channel, grille de 1000x1000x1000
        trapin_%coord=coord
        trapin_%sig=sig
        trapin_%tau_r=tau_r
        trapin_%proba_r=1.-exp((-1.d)*time_int/tau_r)
        do i=1, ntraps
            trapin_%species(nt_debut(i):nt_fin(i)) = i
            trapin_%nv(i) = nt(i)
        end do
        !etat (0 vide 1 occupe)
        trapin_%state=0
        
        allocate (trapsin(nlines, ncolumns))
        trapsin = trapin_

        !generation des pieges mode stochastique
        allocate (trapsout(nlines, ncolumns))
        trapsout%coord        = trapsin%coord
        trapsout%species      = trapsin%species
        trapsout%state(:,1)   = trapsin%state
        trapsout%sig          = trapsin%sig
        trapsout%tau_r        = trapsin%tau_r
        trapsout%proba_r      = trapsin%proba_r
        trapsout%tau_c(:,1)   = trapsin%tau_c
        trapsout%proba_c(:,1) = trapsin%proba_c
        trapsout%nv(:,1)      = trapsin%nv
        trapsout%nc(:,1)      = trapsin%nc
        
        !volume occupe par les electrons libres
        allocate (conf_vol(nlines, ncolumns, nstep))

        !profondeur du volume occupe par les electrons libres
        allocate (conf_z(nlines, ncolumns, nstep))
        
        !choix du modele density driven ou volume driven
        if (mode == 'v') then
            model = 1.d/3.d
        else if (mode == 'd') then
            model = 0
        else
            stop 'mode should be either d or v'
        end if
        
        !integration de l'image
        call system_clock(count1, count_rate, count_max)
        
        do i=1, nstep
            if (i == (i/1000)*1000) print *,'step:',i
            !calcul du nombre d'electrons rajoutes a chaque pas d'integration
            imaout(:,:,i) = imaout(:,:,i)+tirage_photons(:,:,i)

            !calcul du trapping et detrapping
            do itrap=1, ntraps

                !calcul du volume de confinement des electrons
                !  conf_vol(:,:,i)=(bc_vol/3.)*(imaout(:,:,i)/fwc_euclid)^(double(1.)/double(3.))
                conf_vol(:,:,i) = bc_vol*(imaout(:,:,i)/fwc_euclid)**model
                !calcul de la profondeur du volume de confinement
                !origine en profondeur a z=0
                conf_z(:,:,i) = conf_vol(:,:,i) / bc_vol  
                !capture par les pieges vides
                temp_itrap=trapsout%nv(itrap,i)

                wv=where(imaout(:,:,i) /= 0 .and. temp_itrap > 0, count_wv)
                !capture par les pieges vides
                if (count_wv > 0) then
                    !vectorisations
                    temp_ima = imaout(:,:,i)
                    temp_ima = temp_ima(wv)
                    temp_conf_vol = conf_vol(:,:,i)
                    temp_conf_vol = temp_conf_vol(wv)
                    temp_conf_z   = conf_z(:,:,i)
                    temp_conf_z   = temp_conf_z(wv)
                    temp_coord    = trapsout(wv)%coord(:,2)
                    temp_species  = trapsout(wv)%species
                    temp_state    = trapsout(wv)%state(:,i)
                    !constante de temps de capture
                    alp = (sig(itrap)*vth) / temp_conf_vol
                    temp_tau_c = 1_p / (alp*temp_ima)
                    trapsout(wv)%tau_c(itrap,i) = temp_tau_c
                    !probabilite de capture
                    temp_proba_c=1.-exp((-1.d)*time_int/temp_tau_c)
                    trapsout(wv)%proba_c(itrap,i)=temp_proba_c
                    !tirage du resultat
                    do iv=0,count_wv-1
                        !capture?
                        !on identifie les pieges vides de chaque pixel
                        wt=where(temp_species(:,iv) == itrap .and. temp_coord(:,iv) <= temp_conf_z(iv) .and. temp_state(:,iv) == 0,count_wt)
                        if (count_wt > 0) then
                            tir=fix(randomu(seed,count_wt)+temp_proba_c(iv))
                            !on cherche les pieges dont le tir est "positif"
                            wun=where(tir eq 1,count_wun)
                            !les pieges capturent au maximum le nombre d'e- libres du pixel
                            if (count_wun > temp_ima(iv)) then
                                wun=wun(0:(temp_ima(iv)-1))
                                count_wun=temp_ima(iv)
                            end if
                            !changement d'etat et mise a jour du contenu en electrons libres
                            if (count_wun > 0) then
                                trapsout(wv(iv))%state(wt(wun),i)=1
                                temp_ima(iv)=temp_ima(iv)-count_wun
                            end if
                        end if
                    end do
                    !reformatage de l'image
                    nc = imaout(:,:,i)
                    nc(wv) = temp_ima
                    imaout(:,:,i) = nc
                end if
                trapsout%nc(itrap,i)=sum(trapsout%state(nt_debut(itrap):nt_fin(itrap),i),1)
                trapsout%nv(itrap,i)=nt(itrap)-trapsout%nc(itrap,i)
                !relache par les pieges occupes
                !on identifie les pieges pleins de chaque pixel
                temp_ima=imaout(:,:,i)
                temp_species=trapsout%species
                temp_state=trapsout%state(:,i)
                !probabilite de relache
                temp_proba_r=trapsout(0)%proba_r(itrap)
                wt=where(temp_species == itrap .and. temp_state == 1,count_wt)
                if (count_wt > 0) then
                    !relache?
                    tir=fix(randomu(seed,count_wt)+temp_proba_r)
                    !on cherche les pieges dont le tir est "positif"
                    wun=where(tir == 1,count_wun)
                    !changement d'etat et mise a jour du contenu en electrons libres
                    if (count_wun > 0) then
                        !pieges changeant d'etat
                        temp_state(wt(wun))=0
                        trapsout%state(:,i)=temp_state
                        !pixels comprenant des pieges changeant d'etat
                        pixt = wt(wun)/tottraps
                        !indices des pixels
                        pixind = uniq(pixt)
                        !numeros des pixels
                        pixn = pixt(pixind)
                        !nombre de pieges changeant d'etat par pixel
                        totpix = n_elements(pixind)
                        if (totpix == 1) then
                            count=pixind(1)+1
                        else
                            count=[pixind(1)+1,pixind(2:totpix)-pixind(1:totpix-1)]
                        end if
                        temp_ima(pixn)=temp_ima(pixn)+count
                    end if
                end if

                !reformatage de l'image
                imaout(:,:,i) = temp_ima
                trapsout%nc(itrap,i) = sum(trapsout%state(nt_debut(itrap):nt_fin(itrap),i),1)
                trapsout%nv(itrap,i) = nt(itrap)-trapsout%nc(itrap,i)

            end do
            
            trapsout%state(:,i+1)=trapsout%state(:,i)
            imaout(:,:,i+1)=imaout(:,:,i)
            trapsout%nv(:,i+1)=trapsout%nv(:,i)
            trapsout%nc(:,i+1)=trapsout%nc(:,i)
            
        end do
        
        call system_clock(count2, count_rate, count_max)
        write (*,'(a,f7.2,a)') 'proba_integration_full: ', real(count2-count1)/count_rate, 's'

    end subroutine simul_trapping_proba_integration_full


    !-------------------------------------------------------------------------------------------------------------------------------


!!$    subroutine simul_trapping_proba_parallel_full_bis, imaout, time, trapsout, imapar, trapspar, conf_z, imarem
!!$
!!$        integer, intent(out) :: imaout(:,:)
!!$        integer, intent(out) :: imapar(size(imaout,1), size(imaout,2))
!!$        integer, intent(out), allocatable :: imarem(:,:,:)
!!$
!!$        real(p) :: model
!!$        integer :: i
!!$        integer, allocatable :: tir_r(:,:,:)
!!$        real(p), allocatable :: conf_vol(:,:,:), conf_z(:,:,:), temp_conf_vol(:,:)
!!$        integer :: count1, count2, count_max, count_ratex
!!$        type(trap(:,:)), allocatable :: trapspar(:,:)
!!$
!!$        integer, allocatable :: wt(:), wun(:)
!!$        integer :: count_wt, count_wv, count_wun
!!$        integer :: imapartmp(size(imaout,1), size(imaout,2))
!!$
!!$        !parametres par defaut
!!$        print *, 'time_par: ',time_par
!!$        print *, 'mode: ',mode
!!$
!!$        npixal = size(imaout, 1)
!!$        npixac = size(imaout, 2)
!!$  
!!$        !nombre de pas de transfert parallel
!!$        npar = npixal
!!$        time = findgen(npar+1) * time_par
!!$
!!$        !image en debut de transfert parallel
!!$        imapar = 0
!!$
!!$        !pieges en debut de transfert parallel
!!$        allocate (trapspar(npixal, npixac))
!!$        trapspar%coord = trapsout%coord
!!$        trapspar%species = trapsout%species
!!$        trapspar%state(:,1) = trapsout%state
!!$        trapspar%sig = trapsout%sig
!!$        trapspar%tau_r = trapsout%tau_r
!!$        do itrap = 1, ntraps do
!!$            trapspar%proba_r(itrap) = 1 - exp(-time_par/tau_r(itrap))
!!$        end do
!!$        trapspar%tau_c(:,1) = trapsout%tau_c
!!$        trapspar%proba_c(:,1) = trapsout%proba_c
!!$        trapspar%nv(:,1) = trapsout%nv
!!$        trapspar%nc(:,1) = trapsout%nc
!!$
!!$        !image remanente
!!$        allocate (imarem(npixal, npixac, npar))
!!$        imarem = 0
!!$
!!$        !tirages de probabilites
!!$        allocate (tir_r(tottraps, npixal, npixac))
!!$
!!$        !volume occupe par les electrons libres
!!$        allocate (conf_vol(npixal, npixac, npar))
!!$        allocate (temp_conf_vol(npixal, npixac))
!!$
!!$        !profondeur du volume occupe par les electrons libres
!!$        allocate (conf_z(npixal, npixac, npar))
!!$
!!$        !choix du modele density driven ou volume driven
!!$        if (mode ==  'v') then
!!$            model = 1_p/3_p
!!$        else if (mode ==  'd') then
!!$            model = 0
!!$        else
!!$            stop 'mode should be either d or v'
!!$        end if
!!$
!!$        !integration de l'image
!!$        call system_clock(count1, count_rate, count_max)
!!$
!!$        !un pas = un transfert ligne
!!$        imapartmp = imaout
!!$        do i = 1, npar
!!$
!!$            if (i == (i/100)*100) print *, 'step:', i
!!$            
!!$            !transfert de la 1ere ligne de la zone image dans le registre serie
!!$            imapar(i,*) = imapartmp(1,:)
!!$            !decalage d'1 ligne de la zone image
!!$            imapartmp = [imapartmp(1:(npixal-1),*),lonarr(1,npixac)]
!!$            !calcul du trapping et detrapping
!!$            do itrap = 1, ntraps
!!$                !calcul du volume de confinement des electrons
!!$                conf_vol(:,:,i) = bc_vol*(imapartmp/fwc_euclid)**model
!!$                !calcul de la profondeur du volume de confinement
!!$                !origine en profondeur a z = 0
!!$                conf_z(:,:,i) = conf_vol(:,:,i) / bc_vol
!!$                !capture par les pieges vides
!!$                temp_itrap = trapspar%nv(itrap,i)
!!$                wv = where(imapartmp .ne. 0 .and. temp_itrap > 0, count_wv)
!!$                !capture par les pieges vides
!!$                if (count_wv > 0) then
!!$                    !vectorisations
!!$                    temp_ima = imapartmp
!!$                    temp_ima = temp_ima(wv)
!!$                    temp_conf_vol = conf_vol(:,:,i)
!!$                    temp_conf_vol = temp_conf_vol(wv)
!!$                    temp_conf_z = conf_z(:,:,i)
!!$                    temp_conf_z = temp_conf_z(wv)
!!$                    temp_coord = trapspar(wv)%coord(:,2)
!!$                    temp_species = trapspar(wv)%species
!!$                    temp_state = trapspar(wv)%state(:,i)
!!$                    !constante de temps de capture
!!$                    alp = sig(itrap) * vth / temp_conf_vol
!!$                    temp_tau_c = 1 / (alp * temp_ima)
!!$                    trapspar(wv)%tau_c(itrap,i) = temp_tau_c
!!$                    !probabilite de capture
!!$                    temp_proba_c = 1 - exp(-time_par/temp_tau_c)
!!$                    trapspar(wv)%proba_c(itrap,i) = temp_proba_c
!!$                    !tirage du resultat
!!$                    do iv = 1, n_elements(wv)
!!$                        !capture?
!!$                        !on identifie les pieges vides de chaque pixel
!!$                        call where(temp_species(:,iv) == itrap .and. temp_coord(:,iv) <= temp_conf_z(iv) .and. temp_state(:,iv) == 0, wt, count_wt)
!!$                        if (count_wt > 0) then
!!$                            tir = fix(randomu(seed, count_wt)+temp_proba_c(iv))
!!$                            !on cherche les pieges dont le tir est "positif"
!!$                            call where(tir == 1, wun, count_wun)
!!$                            !les pieges capturent au maximum le nombre d'e- libres du pixel
!!$                            if (count_wun > temp_ima(iv)) then
!!$                                wun = wun(1:temp_ima(iv))
!!$                                count_wun = temp_ima(iv)
!!$                            end if
!!$                            !changement d'etat et mise a jour du contenu en electrons libres
!!$                            if (count_wun > 0) then
!!$                                trapspar(wv(iv))%state(wt(wun),i) = 1
!!$                                temp_ima(iv) = temp_ima(iv)-count_wun
!!$                            end if
!!$                        end if
!!$                    end do
!!$
!!$                    !reformatage de l'image
!!$                    nc = imapartmp
!!$                    nc(wv) = temp_ima
!!$                    imapartmp = nc
!!$
!!$                end if
!!$
!!$                trapspar%nc(itrap,i) = sum(trapspar%state(nt_debut(itrap):nt_fin(itrap),i),1)
!!$                trapspar%nv(itrap,i) = nt(itrap) - trapspar%nc(itrap,i)
!!$
!!$            end do
!!$
!!$            !relache par les pieges occupes
!!$            !on identifie les pieges pleins de chaque pixel
!!$            temp_ima = imapartmp
!!$            !    temp_species = trapsout%species
!!$            temp_state = trapspar%state(:,i)
!!$            temp_proba_r = trapspar(1)%proba_r
!!$
!!$            do itrap = 1,ntraps
!!$                tir_r(nt_debut(itrap):nt_fin(itrap),:,:) = fix(randomu(seed,nt(itrap),npixal,npixac)+temp_proba_r(itrap))
!!$            end do
!!$
!!$            wt = where(tir_r ==  1 .and. temp_state ==  1, count_wt)
!!$            if (count_wt > 0) then
!!$                !relache?
!!$                !      tir = fix(randomu(seed,count_wt)+temp_proba_r)
!!$                !on cherche les pieges dont le tir est "positif"
!!$                !      wun = where(tir eq 1,count_wun)
!!$                !changement d'etat et mise a jour du contenu en electrons libres
!!$                !      if count_wun gt 0 then begin
!!$                !pieges changeant d'etat
!!$                !        temp_state(wt(wun)) = 0
!!$                temp_state(wt) = 0
!!$                trapspar%state(:,i) = temp_state
!!$                !pixels comprenant des pieges changeant d'etat
!!$                !        pixt = wt(wun)/tottraps
!!$                pixt = wt/tottraps
!!$                !indices des pixels
!!$                pixind = uniq(pixt)
!!$                !numeros des pixels
!!$                pixn = pixt(pixind)
!!$                !nombre de pieges changeant d'etat par pixel
!!$                totpix = n_elements(pixind)
!!$                if (totpix ==  1) then
!!$                    count = pixind(0)+1
!!$                else
!!$                    count = [pixind(0)+1,pixind(1:(totpix-1))-pixind(0:(totpix-2))]
!!$                end if
!!$                temp_ima(pixn) = temp_ima(pixn)+count
!!$            end if
!!$
!!$            !reformatage de l'image
!!$            imapartmp = temp_ima
!!$
!!$            !fabrication de l'image remanente
!!$            imarem(:,:,i) = imapartmp
!!$
!!$            do itr = 1, ntraps
!!$                trapspar%nc(itr,i) = sum(trapspar%state(nt_debut(itr):nt_fin(itr),i),1)
!!$                trapspar%nv(itr,i) = nt(itr)-trapspar%nc(itr,i)
!!$            end do
!!$
!!$            trapspar%state(:,i+1) = trapspar%state(:,i)
!!$            trapspar%nv(:,i+1) = trapspar%nv(:,i)
!!$            trapspar%nc(:,i+1) = trapspar%nc(:,i)
!!$
!!$        end do
!!$
!!$        call system_clock(count2, count_rate, count_max)
!!$        write (*,'(a,f7.2,a)') 'Elapsed time: ', real(count2-count1)/count_rate, 's'
!!$        
!!$    end subroutine simul_trapping_proba_parallel_full_bis

end module module_euclidtrapping
