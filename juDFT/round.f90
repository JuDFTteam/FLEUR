module m_juDFT_round
contains
   function round_to_deci(x, n) result(rounded_x)
      !rounds to the n-th decimal point
      implicit none
      real, intent(in)    :: x
      integer, intent(in) :: n
      real                :: deci_shift, rounded_x

      deci_shift = 10**n
      rounded_x = real(NINT(x * deci_shift ) / deci_shift)
   end function round_to_deci
end module m_juDFT_round
