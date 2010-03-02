program test_wcs

    use precision, only : test_real_eq
    use module_fitstools, only : ft_header2str
    use module_wcs
    implicit none

    character(len=*), parameter :: filename_header = 'tests/csh_header_ngc6946.fits'
    integer                     :: n

    type(astrometry)            :: astr
    character(len=2880)         :: header
    integer                     :: status, i, j
    real*8, dimension(5)        :: a, d, x, y, xref, yref
    real*8, dimension(:), allocatable   :: a_vect, d_vect, x_vect, y_vect
    real*8, dimension(:,:), allocatable :: ad_vect, xy_vect
    integer                     :: count1, count2, count_rate, count_max

    a = [308.3404249521426d0, 308.3404249521426d0-2d0,        &
         308.3404249521426d0+10d0, 308.3404249521426d0-2d0,   &
         308.3404249521426d0+10d0]
    d = [60.357773247484651d0, 60.357773247484651d0-10d0,     &
         60.357773247484651d0-10d0, 60.3577732474846510+10d0, &
         60.357773247484651d0+10d0]

    ! header =headfits('tests/csh_header_ngc6946.fits')
    ! adxy, header, a, d, x, y
    ! print_text, x+1, y+1

    xref = [ 1.00000000000, 1555.78869806518, -7771.43047054905, &
             820.11316879919, -4084.69147364694 ]
    yref = [ 1.0000000000, -12101.1296270193, -11590.6432782866, &
              12138.0230720130,  12466.1951260863]

    call ft_header2str(filename_header, header, status=status)
    if (status /= 0) stop 'FAILED: ft_header2str.'

    call init_astrometry(header, astr, status)
    if (status /= 0) stop 'FAILED: init_astrometry.'

    call print_astrometry(astr)

    do i=1,5
        call ad2xy_gnomonic_vect(a(i), d(i), x(i), y(i))
        if (.not. test_real_eq(x(i), xref(i), 7)) stop 'FAILED: ad2xy_gnomonic'
        if (.not. test_real_eq(y(i), yref(i), 7)) stop 'FAILED: ad2xy_gnomonic'
    end do

    n = 10000000
    allocate(a_vect(n))
    allocate(d_vect(n))
    allocate(ad_vect(2,n))

    a_vect = [(a(1)+(10.d0*i)/n, i=1,n)]
    d_vect = [(d(1)+(10.d0*i)/n, i=1,n)]
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
