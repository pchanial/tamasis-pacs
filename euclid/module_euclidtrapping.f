module module_euclidtrapping

    use module_fitstools, only : ft_read_image
    use module_precision, only : p => dp
    use module_sort,      only : qsorti

    real(p), parameter   :: t_ccd = 163._p
    real(p), parameter   :: int_time = 450._p
    real(p), parameter   :: time_int = 0.1_p
    real(p), parameter   :: time_par = 1.e-2_p

    integer, parameter   :: nlines = 256
    integer, parameter   :: ntraptypes = 3

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
    integer, parameter :: ntraps = nt_(1) + nt_(2) + nt_(3)
#else
    integer, parameter :: ntraps = sum(nt_)
#endif
    integer, parameter :: nsteps = nint(int_time / time_int)

contains


    subroutine tirage_proba(proba, ntirages, tirage)
        !procedure de tirage d'un evenement ayant une certaine probabilite d'arriver
        !proba: probabilite de l'evenement
        !ntirages: nombre de tirages
        !on cree un vecteur de "ntirages" tirages compose de 1 et 0 avec un rapport 1/0 egal a "proba"
        !on reordonne aleatoirement ce vecteur
        
        real(p), intent(in)  :: proba(:)
        integer, intent(in)  :: ntirages
        integer, intent(out) :: tirage(ntirages,size(proba))

        integer :: tirage_(size(proba),ntirages), n_un(size(proba)), ordre(ntirages)
        integer :: iline
        real(p) :: rnd(ntirages)
        integer :: count1, count2, count_rate, count_max

        tirage_ = 0
        n_un = nint(proba * ntirages)
 
        call system_clock(count1, count_rate, count_max)
        call random_number(rnd)

        call qsorti(rnd, ordre)

        do iline = 1, nlines
            tirage(ordre(1:n_un(iline)),iline) = 1
        end do
        
        tirage = transpose(tirage_)

        call system_clock(count2, count_rate, count_max)
        write (*,'(a,f7.2,a)') 'tirage_proba: ', real(count2-count1)/count_rate, 's'

    end subroutine tirage_proba


    !-------------------------------------------------------------------------------------------------------------------------------


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
