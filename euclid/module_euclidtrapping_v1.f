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

module module_euclidtrapping_v1

    use module_euclidtrapping
    use module_fitstools, only : ft_read_image, ft_write_image
    use module_math,      only : PI, linspace
    use module_sort,      only : qsorti, uniq, where
    use module_precision, only : p => dp, sp
    implicit none

    character, parameter :: mode = 'd'
    integer, parameter   :: ncolumns = 1
    
    ! pieges mode stochastique
    type TrapIn
        real(sp) :: coord(ntraps, 3)
        integer  :: species(ntraps)
        integer  :: state(ntraps)
        real(p)  :: sig(ntraptypes)
        real(p)  :: tau_r(ntraptypes)
        real(p)  :: tau_c(ntraptypes)
        real(p)  :: proba_r(ntraptypes)
        real(p)  :: proba_c(ntraptypes)
        integer  :: nv(ntraptypes)
        integer  :: nc(ntraptypes)
    end type TrapIn
        
    type TrapOut
        real(sp) :: coord(ntraps, 3)
        integer  :: species(ntraps)
        integer  :: state(ntraps, nsteps+1)
        real(p)  :: sig(ntraptypes)
        real(p)  :: tau_r(ntraptypes)
        real(p)  :: tau_c(ntraptypes, nsteps)
        real(p)  :: proba_r(ntraptypes)
        real(p)  :: proba_c(ntraptypes, nsteps)
        integer  :: nv(ntraptypes, nsteps+1)
        integer  :: nc(ntraptypes, nsteps+1)
    end type TrapOut

    real(p) :: vth, nstates
    real(p), dimension(ntraptypes) :: sig, tau_r, e_traps
    integer :: nt(ntraptypes) !nt: nombre de pieges par volume de confinement (reel)
    integer :: tottraps, isig(ntraptypes)
    integer :: nt_debut(ntraptypes), nt_fin(ntraptypes)
    

contains


    subroutine init_v1

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

    end subroutine init_v1


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine simul_integration_v1(imain, time, tirage_photons, coord, trapsin, imaout, trapsout, conf_z, status)

        integer, intent(in)                     :: imain(nlines,ncolumns)
        real(p), intent(out), allocatable       :: time(:)
        integer, intent(out), allocatable       :: tirage_photons(:,:,:)
        real(sp), intent(out), allocatable      :: coord(:,:)
        integer, intent(out), allocatable       :: imaout(:,:,:)
        type(TrapIn), intent(out), allocatable  :: trapsin(:,:)
        type(TrapOut), intent(out), allocatable :: trapsout(:,:)
        real(p), intent(out), allocatable       :: conf_z(:,:,:)
        integer, intent(out)                    :: status

        type(TrapIn)         :: trapin_
        real(p), allocatable :: conf_vol(:,:,:), coord_(:,:)
        integer, allocatable :: wv(:), wt(:), wun(:), wun_(:), temp_ima(:), temp_species(:,:), temp_state(:,:), tir(:)
        integer, allocatable :: pixt(:), pixn(:), pixind(:), count(:)
        real(p), allocatable :: temp_conf_vol(:), temp_conf_z(:), temp_coord(:,:), alp(:), temp_tau_c(:), temp_proba_c(:), rnd(:)
        real(p)              :: model, temp_proba_r
        integer              :: nsteps, i, iv, itrap, itraptot, count_wt, count_wv, count_wun, iline, icolumn, totpix, iun
        integer              :: temp_itrap(nlines), flux(nlines,ncolumns)
        integer              :: count1, count2, count_rate, count_max
        
        call init_v1

        !calcul du nombre d'e- rajoutes par effet photo-electrique a chaque pas de temps
        !generation de la sequence temporelle d'arrivee des photons pendant l'integration
        !nombre de pas d'integration
        nsteps = nint(int_time/time_int)
        print *, 'nsteps:', nsteps
        allocate (time(nsteps+1))
        time = linspace(0._p, nsteps * time_int, nsteps+1)
        flux = imain / nsteps

        allocate (tirage_photons(nlines, ncolumns, nsteps))
        call tirage_proba(imain(:,1)/real(nsteps, kind=p)-flux(:,1), nsteps, tirage_photons(:,1,:))

        do i=1, nsteps
            tirage_photons(:,:,i) = tirage_photons(:,:,i) + flux
        end do
        
        allocate (imaout(nlines, ncolumns, nsteps+1))
        imaout = 0
        
        !generation des pieges dans le volume du buried channel
        call ft_read_image('~/work/tamasis/tamasis-unstable/euclid/data/coord_traps_integ_proba_comp_94.fits', coord_,status=status)
        if (status /= 0) return
        allocate (coord(size(coord_,1),size(coord_,2)))
        coord = real(coord_, kind=sp)

        !coordonnees normalisees par rapport aux dimensions du buried channel, grille de 1000x1000x1000
        trapin_%coord = coord
        trapin_%sig = sig
        trapin_%tau_r = tau_r
        trapin_%proba_r = 1._p-exp(-time_int/tau_r)
        do i=1, ntraptypes
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
        allocate (conf_vol(nlines, ncolumns, nsteps))

        !profondeur du volume occupe par les electrons libres
        allocate (conf_z(nlines, ncolumns, nsteps))
        
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

            if (i == (i/1000)*1000) print *,'step:',i

            !calcul du nombre d'electrons rajoutes a chaque pas d'integration
            imaout(:,:,i) = imaout(:,:,i) + tirage_photons(:,:,i)

            !calcul du trapping et detrapping
            do itrap=1, ntraptypes

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
                if (allocated(temp_ima))     deallocate (temp_ima)
                if (allocated(temp_species)) deallocate (temp_species)
                if (allocated(temp_state))   deallocate (temp_state)
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

    end subroutine simul_integration_v1

end module module_euclidtrapping_v1
