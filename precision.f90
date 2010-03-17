module precision

   integer, parameter      :: sp = selected_real_kind(6)
   integer, parameter      :: dp = selected_real_kind(12)
   integer, parameter      :: p  = dp
   real(kind=p), parameter :: zero = 0.0_p
   real(kind=p), parameter :: half = 0.5_p
   real(kind=p), parameter :: one  = 1.0_p
   
end module precision
