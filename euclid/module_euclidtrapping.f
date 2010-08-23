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

#ifdef IFORT
    use ifport, only : rand
#endif
    use module_fitstools, only : ft_read_image, ft_write_image
    use module_math,      only : PI, linspace
    use module_sort,      only : qsorti, uniq, where
    use module_precision, only : p => dp, sp
    implicit none

    real(p), parameter   :: t_ccd = 163._p
    real(p), parameter   :: int_time = 450._p
    real(p), parameter   :: time_int = 0.1_p
    real(p), parameter   :: time_par = 1.e-2_p
    character, parameter :: mode = 'd'

    integer, parameter   :: npixal = 256
    integer, parameter   :: npixac = 1
    integer, parameter   :: nlines = npixal
    integer, parameter   :: ncolumns = npixac
    integer, parameter   :: ntraps = 3 !nombre de familles de pieges considerees
    real(p), parameter   :: k = 1.381e-23_p !J/K
    real(p), parameter   :: me = 0.5_p * 9e-31_p !kg ! XXXPCH WHAT'S THAT?
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
    integer, parameter :: nstep = nint(int_time / time_int)
    
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
        integer  :: state(tottraps_, nstep+1)
        real(p)  :: sig(ntraps)
        real(p)  :: tau_r(ntraps)
        real(p)  :: tau_c(ntraps, nstep)
        real(p)  :: proba_r(ntraps)
        real(p)  :: proba_c(ntraps, nstep)
        integer  :: nv(ntraps, nstep+1)
        integer  :: nc(ntraps, nstep+1)
    end type Trap

    type Detector
        integer :: nv
        integer :: nc
        logical :: filled(tottraps_)
    end type Detector

    type DetectorColumn
        integer        :: ndetectors
        integer        :: ntraps
        real(p)        :: vth
        real(p)        :: z(tottraps_)
        real(p)        :: sig(tottraps_)
        real(p)        :: proba_r(tottraps_)
        type(Detector) :: detector(nlines)
    end type DetectorColumn

    real(p) :: vth, nstates
    real(p), dimension(ntraps) :: sig, tau_r, e_traps
    integer :: nt(ntraps) !nt: nombre de pieges par volume de confinement (reel)
    integer :: tottraps, isig(ntraps)
    integer :: nt_debut(ntraps), nt_fin(ntraps)
    type(DetectorColumn) :: column
    

contains


    subroutine init

        integer :: itrap, idetector, status
        real(p), allocatable :: coord_(:,:)

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
        isig = isig(ntraps:1:-1)
        sig = sig(isig)
        tau_r = tau_r(isig)
        nt = nt(isig)
        
        !indices des pieges des differentes familles
        nt_debut(1) = 1
        nt_fin(1) = nt(1)
        do itrap = 2, ntraps
            nt_debut(itrap) = nt_fin(itrap-1) + 1
            nt_fin(itrap) = nt_fin(itrap-1) + nt(itrap)
        end do

        column%ndetectors = nlines
        column%ntraps = sum(nt)
        column%vth = vth

        !generation des pieges dans le volume du buried channel
        call ft_read_image('~/work/tamasis/tamasis-unstable/euclid/data/coord_traps_integ_proba_comp_94.fits', coord_, status=status)
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

    end subroutine init


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


    subroutine tirage_proba(proba, ntirages, tirage)
        !procedure de tirage d'un evenement ayant une certaine probabilite d'arriver
        !proba: probabilite de l'evenement
        !ntirages: nombre de tirages
        !on cree un vecteur de "ntirages" tirages compose de 1 et 0 avec un rapport 1/0 egal a "proba"
        !on reordonne aleatoirement ce vecteur
        
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
        call notrandom_number(rnd)

        call qsorti(rnd, ordre)

        do icolumn = 1, ncolumns
            do iline = 1, nlines
                tirage(iline,icolumn,ordre(1:n_un(iline,icolumn))) = 1
            end do
        end do

        call system_clock(count2, count_rate, count_max)
        write (*,'(a,f7.2,a)') 'tirage_proba: ', real(count2-count1)/count_rate, 's'

    end subroutine tirage_proba_slow
    subroutine tirage_proba_slow2(proba, ntirages, tirage)
        !procedure de tirage d'un evenement ayant une certaine probabilite d'arriver
        !proba: probabilite de l'evenement
        !ntirages: nombre de tirages
        !on cree un vecteur de "ntirages" tirages compose de 1 et 0 avec un rapport 1/0 egal a "proba"
        !on reordonne aleatoirement ce vecteur
        
        real(p), intent(in) :: proba(:)
        integer, intent(in) :: ntirages
        integer, intent(out) :: tirage(size(proba),ntirages)

        integer :: n_un(size(proba)), ordre(ntirages)
        integer :: iline
        real(p) :: rnd(ntirages)
        integer :: count1, count2, count_rate, count_max

        tirage = 0
        n_un = nint(proba * ntirages)
 
        call system_clock(count1, count_rate, count_max)
        call notrandom_number(rnd)

        call qsorti(rnd, ordre)

        do iline = 1, nlines
            tirage(iline,ordre(1:n_un(iline))) = 1
        end do

        call system_clock(count2, count_rate, count_max)
        write (*,'(a,f7.2,a)') 'tirage_proba: ', real(count2-count1)/count_rate, 's'

    end subroutine tirage_proba_slow2


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine simul_trapping_proba_integration_full(imain, time, tirage_photons, coord, trapsin, imaout, trapsout, conf_z, status)

        integer, intent(in)  :: imain(nlines, ncolumns)
        real(p), intent(out), allocatable :: time(:)
        integer, intent(out), allocatable :: tirage_photons(:,:,:)
        real(sp), intent(out), allocatable :: coord(:,:)
        type(TrapIn), intent(out), allocatable :: trapsin(:,:)
        integer, intent(out), allocatable :: imaout(:,:,:)
        type(Trap), intent(out), allocatable :: trapsout(:,:)
        real(p), intent(out), allocatable :: conf_z(:,:,:)
        integer, intent(out) :: status

        type(TrapIn)         :: trapin_
        real(p), allocatable :: conf_vol(:,:,:), coord_(:,:)
        real(p) :: model, temp_proba_r
        integer :: nstep, i, iv, itrap, itraptot, count_wt, count_wv, count_wun, iline, icolumn, totpix, iun
        integer, allocatable :: wv(:), wt(:), wun(:), wun_(:), temp_ima(:), temp_species(:,:), temp_state(:,:), tir(:)
        integer, allocatable :: pixt(:), pixn(:), pixind(:), count(:)
        real(p), allocatable :: temp_conf_vol(:), temp_conf_z(:), temp_coord(:,:), alp(:), temp_tau_c(:), temp_proba_c(:), rnd(:)
        integer :: temp_itrap(nlines), flux(nlines,ncolumns)

        integer :: count1, count2, count_rate, count_max
        
        call init

        !calcul du nombre d'e- rajoutes par effet photo-electrique a chaque pas de temps
        !generation de la sequence temporelle d'arrivee des photons pendant l'integration
        !nombre de pas d'integration
        nstep = nint(int_time/time_int)
        print *, 'nstep:', nstep
        allocate (time(nstep+1))
        time = linspace(0._p, nstep * time_int, nstep+1)
        flux = imain / nstep

        allocate (tirage_photons(nlines, ncolumns, nstep))
        call tirage_proba_slow(imain/real(nstep, kind=p)-flux, nstep, tirage_photons)

        do i=1, nstep
            tirage_photons(:,:,i) = tirage_photons(:,:,i) + flux
        end do
        
        allocate (imaout(nlines, ncolumns, nstep+1))
        imaout = 0
        
        !generation des pieges dans le volume du buried channel
        call ft_read_image('~/work/tamasis/tamasis-latest/euclid/data/coord_traps_integ_proba_comp_94.fits', coord_, status=status)
        if (status /= 0) return
        allocate (coord(size(coord_,1),size(coord_,2)))
        coord = real(coord_, kind=sp)

        !coordonnees normalisees par rapport aux dimensions du buried channel, grille de 1000x1000x1000
        trapin_%coord = coord
        trapin_%sig = sig
        trapin_%tau_r = tau_r
        trapin_%proba_r = 1._p-exp(-time_int/tau_r)
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
        do icolumn = 1, ncolumns
            do iline = 1, nlines
            
                trapsout(iline,icolumn)%coord        = trapsin(iline,icolumn)%coord
                trapsout(iline,icolumn)%species      = trapsin(iline,icolumn)%species
                trapsout(iline,icolumn)%state(:,1)   = trapsin(iline,icolumn)%state
                trapsout(iline,icolumn)%sig          = trapsin(iline,icolumn)%sig
                trapsout(iline,icolumn)%tau_r        = trapsin(iline,icolumn)%tau_r
                trapsout(iline,icolumn)%proba_r      = trapsin(iline,icolumn)%proba_r
                trapsout(iline,icolumn)%tau_c(:,1)   = trapsin(iline,icolumn)%tau_c
                trapsout(iline,icolumn)%proba_c(:,1) = trapsin(iline,icolumn)%proba_c
                trapsout(iline,icolumn)%nv(:,1)      = trapsin(iline,icolumn)%nv
                trapsout(iline,icolumn)%nc(:,1)      = trapsin(iline,icolumn)%nc

            end do
        end do
        
        !volume occupe par les electrons libres
        allocate (conf_vol(nlines, ncolumns, nstep))

        !profondeur du volume occupe par les electrons libres
        allocate (conf_z(nlines, ncolumns, nstep))
        
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
        
        do i=1, nstep

            if (i == (i/1000)*1000) print *,'step:',i

            !calcul du nombre d'electrons rajoutes a chaque pas d'integration
            imaout(:,:,i) = imaout(:,:,i) + tirage_photons(:,:,i)

            !calcul du trapping et detrapping
            do itrap=1, ntraps

                !calcul du volume de confinement des electrons
                conf_vol(:,:,i) = bc_vol * (imaout(:,:,i)/fwc_euclid)**model
                !calcul de la profondeur du volume de confinement
                !origine en profondeur a z=0
                conf_z(:,:,i) = conf_vol(:,:,i) / bc_vol  
                !capture par les pieges vides
                temp_itrap = trapsout(:,1)%nv(itrap,i)

                call where(imaout(:,1,i) /= 0 .and. temp_itrap > 0, wv, count_wv)

                !capture par les pieges vides
                if (count_wv > 0) then

                    !vectorisations
                    if (allocated(temp_ima))      deallocate (temp_ima)
                    if (allocated(temp_conf_vol)) deallocate (temp_conf_vol)
                    if (allocated(temp_conf_z))   deallocate (temp_conf_z)
                    if (allocated(temp_coord))    deallocate (temp_coord)
                    if (allocated(temp_species))  deallocate (temp_species)
                    if (allocated(temp_state))    deallocate (temp_state)
                    if (allocated(alp))           deallocate (alp)
                    if (allocated(temp_tau_c))    deallocate (temp_tau_c)
                    if (allocated(temp_proba_c))  deallocate (temp_proba_c)
                    
                    allocate (temp_ima(count_wv))
                    allocate (temp_conf_vol(count_wv))
                    allocate (temp_conf_z(count_wv))
                    allocate (temp_coord(tottraps,count_wv))
                    allocate (temp_species(tottraps,count_wv))
                    allocate (temp_state(tottraps,count_wv))
                    allocate (alp(count_wv))
                    allocate (temp_tau_c(count_wv))
                    allocate (temp_proba_c(count_wv))

                    temp_ima      = imaout(wv,1,i)
                    temp_conf_vol = conf_vol(wv,1,i)
                    temp_conf_z   = conf_z(wv,1,i)                    
                    do itraptot=1, tottraps
                        temp_coord(itraptot,:)   = trapsout(wv,1)%coord(itraptot,3)
                        temp_species(itraptot,:) = trapsout(wv,1)%species(itraptot)
                        temp_state(itraptot,:)   = trapsout(wv,1)%state(itraptot,i)
                    end do

                    !constante de temps de capture
                    alp = (sig(itrap)*vth) / temp_conf_vol
                    temp_tau_c = 1_p / (alp*temp_ima)
                    trapsout(wv,1)%tau_c(itrap,i) = temp_tau_c

                    !probabilite de capture
                    temp_proba_c = 1-exp(-time_int/temp_tau_c)
                    trapsout(wv,1)%proba_c(itrap,i) = temp_proba_c

                    !tirage du resultat
                    do iv=1, count_wv

                        !capture?
                        !on identifie les pieges vides de chaque pixel
                        call where(temp_species(:,iv) == itrap .and. temp_coord(:,iv) <= temp_conf_z(iv) .and. &
                             temp_state(:,iv) == 0, wt, count_wt)
                        if (count_wt > 0) then
                            if (allocated(rnd)) deallocate (rnd)
                            if (allocated(tir)) deallocate (tir)
                            allocate (rnd(count_wt))
                            allocate (tir(count_wt))
                            call notrandom_number(rnd)
                            tir = floor(rnd + temp_proba_c(iv))
                            !on cherche les pieges dont le tir est "positif"
                            call where(tir == 1, wun, count_wun)
                            !les pieges capturent au maximum le nombre d'e- libres du pixel
                            if (count_wun > temp_ima(iv)) then
                                if (allocated(wun_)) deallocate(wun_)
                                allocate (wun_(temp_ima(iv)))
                                wun_ = wun(1:temp_ima(iv))
                                call move_alloc(wun_, wun)
                                count_wun = temp_ima(iv)
                            end if
                            !changement d'etat et mise a jour du contenu en electrons libres
                            if (count_wun > 0) then
                                trapsout(wv(iv),1)%state(wt(wun),i) = 1
                                temp_ima(iv) = temp_ima(iv) - count_wun
                            end if
                        end if
                    end do

                    !reformatage de l'image
                    imaout(wv,1,i) = temp_ima

                end if

                do iline=1, nlines
                    trapsout(iline,1)%nc(itrap,i) = sum(trapsout(iline,1)%state(nt_debut(itrap):nt_fin(itrap),i))
                    trapsout(iline,1)%nv(itrap,i) = nt(itrap)-trapsout(iline,1)%nc(itrap,i)
                end do

                !relache par les pieges occupes
                !on identifie les pieges pleins de chaque pixel
                deallocate (temp_ima)
                deallocate (temp_species)
                deallocate (temp_state)
                allocate (temp_ima(nlines))
                allocate (temp_species(tottraps,nlines))
                allocate (temp_state(tottraps,nlines))
                temp_ima=imaout(:,1,i)
                do iline=1, nlines
                    temp_species(:,iline) = trapsout(iline,1)%species(:)
                    temp_state  (:,iline) = trapsout(iline,1)%state(:,i)
                end do

                !probabilite de relache
                temp_proba_r=trapsout(1,1)%proba_r(itrap)
                call where(temp_species == itrap .and. temp_state == 1, wt, count_wt)
                if (count_wt > 0) then
                    !relache?
                    if (allocated(rnd)) deallocate (rnd)
                    if (allocated(tir)) deallocate (tir)
                    allocate (rnd(count_wt))
                    allocate (tir(count_wt))
                    call notrandom_number(rnd)
                    tir = floor(rnd+temp_proba_r)

                    !on cherche les pieges dont le tir est "positif"
                    call where(tir == 1, wun, count_wun)
                    !changement d'etat et mise a jour du contenu en electrons libres
                    if (count_wun > 0) then

                        !pieges changeant d'etat
                        do iun=1, count_wun
                            iline = (wt(wun(iun))-1) / tottraps + 1
                            itraptot = modulo(wt(wun(iun))-1, tottraps) + 1
                            temp_state(itraptot,iline) = 0
                        end do
                        do iline = 1, nlines
                            trapsout(iline,1)%state(:,i)=temp_state(:,iline)
                        end do
                        if (allocated(pixt)) deallocate (pixt)
                        allocate (pixt(count_wun))
                        !pixels comprenant des pieges changeant d'etat
                        pixt = (wt(wun)-1)/tottraps + 1

                        !indices des pixels
                        call uniq(pixt, pixind)
                        totpix = size(pixind)
                        if (allocated(pixn)) deallocate (pixn)
                        allocate (pixn(totpix))
                        !numeros des pixels
                        pixn = pixt(pixind)

                        if (allocated(count)) deallocate (count)
                        allocate (count(totpix))
                        !nombre de pieges changeant d'etat par pixel
                        count = [ pixind(1), pixind(2:totpix)-pixind(1:totpix-1) ]
                        temp_ima(pixn) = temp_ima(pixn) + count
                    end if

                end if

                !reformatage de l'image
                imaout(:,1,i) = temp_ima
                
                do iline=1,nlines
                    trapsout(iline,1)%nc(itrap,i) = sum(trapsout(iline,1)%state(nt_debut(itrap):nt_fin(itrap),i))
                end do
                trapsout(:,1)%nv(itrap,i) = nt(itrap)-trapsout(:,1)%nc(itrap,i)

            end do
            
            do iline = 1, nlines
                trapsout(iline,1)%state(:,i+1) = trapsout(iline,1)%state(:,i)
                trapsout(iline,1)%nv   (:,i+1) = trapsout(iline,1)%nv   (:,i)
                trapsout(iline,1)%nc   (:,i+1) = trapsout(iline,1)%nc   (:,i)
            end do
            imaout(:,:,i+1)=imaout(:,:,i)
            
        end do

        call system_clock(count2, count_rate, count_max)
        write (*,'(a,f7.2,a)') 'proba_integration_full: ', real(count2-count1)/count_rate, 's'

    end subroutine simul_trapping_proba_integration_full


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine simul_integration(imain, time, tirage_photons, imaout)

        integer, intent(in)  :: imain(nlines)
        real(p), intent(out), allocatable :: time(:)
        integer, intent(out), allocatable :: tirage_photons(:,:)
        integer, intent(out), allocatable :: imaout(:,:)

        real(p) :: model
        integer :: nstep, i

        integer :: flux(nlines)
        real(p), dimension(tottraps_) :: alp, tau_c, proba_c, test
        real(p), dimension(nlines) :: conf_vol, conf_z
        real(p) :: rnd0
        integer :: idetector, itrap, nc, ncaptures, nreleases

        integer :: count1, count2, count_rate, count_max
        
        call init

        !calcul du nombre d'e- rajoutes par effet photo-electrique a chaque pas de temps
        !generation de la sequence temporelle d'arrivee des photons pendant l'integration
        !nombre de pas d'integration
        nstep = nint(int_time/time_int)
        print *, 'nstep:', nstep
        allocate (time(nstep+1))
        time = linspace(0._p, nstep * time_int, nstep+1)
        flux = imain / nstep !XXX check me

        allocate (tirage_photons(nlines, nstep))
        call tirage_proba_slow2(imain/real(nstep, kind=p)-flux, nstep, tirage_photons)

        do i=1, nstep
            tirage_photons(:,i) = tirage_photons(:,i) + flux
        end do
        
        allocate (imaout(column%ndetectors, nstep+1))
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
        
        do i=1, nstep

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
                        column%detector(idetector)%filled(itrap) = floor(rand()+column%proba_r(itrap)) == 0
                        if (.not. column%detector(idetector)%filled(itrap)) then
                            nreleases = nreleases + 1
                        end if
                    end if
                end do

!!$                nc = count(column%detector(idetector)%filled)
!!$                where (column%detector(idetector)%filled)
!!$                    test = rand()
!!$                    column%detector(idetector)%filled = floor(test+column%proba_r) == 0
!!$                end where
                imaout(idetector,i) = imaout(idetector,i) + nreleases
                
                column%detector(idetector)%nc = column%detector(idetector)%nc - nreleases
                column%detector(idetector)%nv = column%ntraps - column%detector(idetector)%nc

            end do

            imaout(:,i+1)=imaout(:,i)
            
        end do

        call system_clock(count2, count_rate, count_max)
        write (*,'(a,f7.2,a)') 'proba_integration_full: ', real(count2-count1)/count_rate, 's'

    end subroutine simul_integration


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


    subroutine notrandom_number(array)

        real(p), intent(out) :: array(:)

        real(p), allocatable, save :: data(:)
        integer, save              :: index
        integer                    :: count, status

        count = size(array)
        if (.not. allocated(data)) then
            call ft_read_image('~/work/tamasis/tests/test_random/random.fits', data, status=status)
            if (status /= 0) then
                stop 'NOTRANDOM_NUMBER: Cannot read data file.'
            end if
            index = 0
        end if

        if (count > size(data)) then
            stop 'NOTRANDOM_NUMBER: buffer is not large enough.'
        endif
 
        if (index + count > size(data)) then
            data = cshift(data, index)
            index = 0
        endif
 
        array = data(index+1:index+count)
        index = index + count

    end subroutine notrandom_number

end module module_euclidtrapping
