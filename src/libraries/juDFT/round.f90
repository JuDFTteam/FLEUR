!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_juDFT_round
   implicit none
   private
   public :: round_to_deci
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
