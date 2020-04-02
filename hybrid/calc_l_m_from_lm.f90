module m_calc_l_m_from_lm

contains
   subroutine calc_l_m_from_lm(lm, l, m)
      use m_juDFT
      implicit none
      integer, intent(in)   :: lm
      integer, intent(out)  :: l, m
      if (lm <= 0) call judft_error("We define lm such that goes from 1..lmax**2")
      l = floor(sqrt(lm - 1.0))
      m = lm - (l**2 + l + 1)
   end subroutine calc_l_m_from_lm
end module m_calc_l_m_from_lm
