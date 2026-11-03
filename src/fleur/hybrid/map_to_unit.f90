!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_map_to_unit
   implicit none
   private
   public :: map_to_unit
contains    
   function map_to_unit(x) result(y)
    implicit none 
    real, intent(in) :: x(:)
    real             :: y(size(x))

    y = x

    ! everything close to 0 or 1 get's mapped to 0 and 1
    where (abs(y - anint(y)) < 1e-6) y = anint(y)
    ! map to 0 -> 1 interval
    y = y - floor(y)
 end function map_to_unit
end module m_map_to_unit
