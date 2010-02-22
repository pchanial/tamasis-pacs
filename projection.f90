module projection

 use precision, only : p
 implicit none

contains

!===============================================================================

subroutine radec2yz(ra, dec, y, z, ra0, dec0, pa0)
 
 use precision, only : p

 implicit none
 real(kind=p), intent(in)  :: y(*), z(*)     ! in degrees
 real(kind=p), intent(in)  :: ra0, dec0, pa0 ! PA in degrees
 real(kind=p), intent(out) :: ra(*), dec(*)
 real(kind=p)              :: cospa, sinpa
 real(kind=p), parameter   :: pi = 4.0_p * atan(1.0_p)
 real(kind=p), parameter   :: radeg = pi / 180_p

 cospa = cos(pa0*radeg)
 sinpa = sin(pa0*radeg)

 !y * sinpa + z * cospa = dec - dec0
 !y * cospa - z * sinpa = (ra - ra0) * cos(dec * radeg)

 y = (dec - dec0) * sinpa + (ra - ra0) * cos(dec * radeg) * cospa
 z = (dec - dec0) * cospa - (ra - ra0) * cos(dec * radeg) * sinpa

end function radec2yz

!===============================================================================

subroutine radec2xy(ra, dec, x, y, x0, y0, cdeltinv0, ra0, dec0, pa0)
 
 use precision, only : p

 implicit none
 real(kind=p), intent(in)  :: y(*), z(*)     ! in degrees
 real(kind=p), intent(in)  :: ra0, dec0, pa0 ! PA in degrees
 real(kind=p), intent(out) :: ra(*), dec(*)
 real(kind=p)              :: cospa, sinpa
 real(kind=p), parameter   :: pi = 4.0_p * atan(1.0_p)
 real(kind=p), parameter   :: radeg = pi / 180_p

 cospa = cos(pa0*radeg)
 sinpa = sin(pa0*radeg)

 !y * sinpa + z * cospa = dec - dec0
 !y * cospa - z * sinpa = (ra - ra0) * cos(dec * radeg)

 x = x0 + cdeltinv0 * (dec - dec0) * sinpa + (ra - ra0) * cos(dec * radeg) * cospa
 y = y0 - cdeltinv0 * (dec - dec0) * cospa - (ra - ra0) * cos(dec * radeg) * sinpa

end function radec2xy


end module projection
