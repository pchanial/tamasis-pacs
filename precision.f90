module precision

   integer, parameter      :: p = selected_real_kind(12)
   real(kind=p), parameter :: zero = 0.0_p
   real(kind=p), parameter :: half = 0.5_p
   real(kind=p), parameter :: one  = 1.0_p

contains
   
   elemental function test_real_eq(a, b, n)
      implicit none
      real(kind=p), intent(in) :: a, b
      integer, intent(in)      :: n
      logical                  :: test_real_eq
      real(kind=p)             :: epsilon

      ! check for NaN values
      if (a /= a) then
          test_real_eq = b /= b
          return
      end if
      if (b /= b) then
          test_real_eq = .false.
          return
      end if

      epsilon = 10**(-real(n, kind=p))
      test_real_eq = abs(a-b) <= epsilon * abs(a)
   end function test_real_eq

end module precision
