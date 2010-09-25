program test_wcs

    use module_math,      only : neq_real
    use module_fitstools, only : ft_header2str
    use module_tamasis,   only : p
    use module_wcs
    implicit none

    character(len=*), parameter :: filename_header = 'core/test/data/header.fits'
    integer                     :: n

    type(astrometry)            :: astr
    character(len=2880*2)       :: header
    integer                     :: status, i
    integer                     :: count1, count2, count_rate, count_max
    real(p), dimension(5)     :: a, d, x, y, xref, yref
    real(p), dimension(:), allocatable   :: a_vect, d_vect, x_vect, y_vect
    real(p), dimension(:,:), allocatable :: ad_vect, xy_vect

    a = [308.3404249521426_p, 308.3404249521426_p-2, 308.3404249521426_p+10, 308.3404249521426_p-2, 308.3404249521426_p+10]
    d = [60.357773247484651_p, 60.357773247484651_p-10, 60.357773247484651_p-10, 60.3577732474846510_p+10, 60.357773247484651_p+10]

    xref = [ 1.0_p, 1555.78869806518_p, -7771.43047054905_p, 820.11316879919_p, -4084.69147364694_p ]
    yref = [ 1.0_p, -12101.1296270193_p, -11590.6432782866_p, 12138.0230720130_p,  12466.1951260863_p ]

    call ft_header2str(filename_header, header, status=status)
    if (status /= 0) stop 'FAILED: ft_header2str.'

    call init_astrometry(header, astr, status)
    if (status /= 0) stop 'FAILED: init_astrometry.'

    call print_astrometry(astr)

    do i=1,5
        call ad2xy_gnomonic_vect(a(i), d(i), x(i), y(i))
        if (neq_real(x(i), xref(i), 1e-7_p)) stop 'FAILED: ad2xy_gnomonic'
        if (neq_real(y(i), yref(i), 1e-7_p)) stop 'FAILED: ad2xy_gnomonic'
    end do

    n = 10000000
    allocate(a_vect(n))
    allocate(d_vect(n))
    allocate(ad_vect(2,n))

    a_vect = [(a(1)+(10._p*i)/n, i=1,n)]
    d_vect = [(d(1)+(10._p*i)/n, i=1,n)]
    ad_vect(1,:) = a_vect
    ad_vect(2,:) = d_vect

    call system_clock(count1, count_rate, count_max)
    call ad2xy_gnomonic_vect(a_vect, d_vect, a_vect, d_vect)
    call system_clock(count2, count_rate, count_max)
    write(*,'(a,f6.2,a)') 'ad2xy_vect: ',real(count2-count1)/count_rate, 's'

    call system_clock(count1, count_rate, count_max)
    ad_vect = ad2xy_gnomonic(ad_vect)
    call system_clock(count2, count_rate, count_max)
    write(*,'(a,f6.2,a)') 'ad2xy: ',real(count2-count1)/count_rate, 's'

!!$    call init_wcslib(header, status)
!!$    call system_clock(count1, count_rate, count_max)
!!$    do i=1,1000
!!$        xy_vect = ad2xy_wcslib(ad_vect)
!!$    end do
!!$    call system_clock(count2, count_rate, count_max)
!!$    write(*,'(a,f6.2,a)') 'wcslib: ',real(count2-count1)/count_rate, 's'

    stop 'OK.'
   
end program test_wcs
