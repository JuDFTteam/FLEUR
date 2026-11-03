!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_calc_l_m_from_lm

   use m_juDFT
   implicit none
   private
   public :: calc_l_m_from_lm
contains
   pure subroutine calc_l_m_from_lm(lm, l, m)
      !$acc routine seq
      implicit none
      integer, intent(in)   :: lm
      integer, intent(out)  :: l, m
      !We define lm such that goes from 1..lmax**2
      l = floor(sqrt(lm - 1.0))
      m = lm - (l**2 + l + 1)
   end subroutine calc_l_m_from_lm
end module m_calc_l_m_from_lm
