!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
! hidden first argument instead of in registers. gfortran's native convention
! disagrees, so such a BLAS silently returns garbage for these functions while
! everything else (zgemm, zheev, ...) works fine.
!
! This program calls zdotc with known inputs and compares against a reference
! computed in-line. It exits 0 if the convention matches, 1 otherwise.
program test_zdotc
   implicit none
   integer, parameter :: dp = kind(0.0d0)
   complex(dp) :: x(3), y(3), res, ref
   complex(dp), external :: zdotc
   integer :: i

   x(1) = (1.0_dp,  2.0_dp); x(2) = (3.0_dp, -1.0_dp); x(3) = (0.0_dp,  1.0_dp)
   y(1) = (2.0_dp,  1.0_dp); y(2) = (1.0_dp,  1.0_dp); y(3) = (1.0_dp, -1.0_dp)

   ! zdotc(n,x,1,y,1) = sum_i conjg(x_i) * y_i
   ref = (0.0_dp, 0.0_dp)
   do i = 1, 3
      ref = ref + conjg(x(i)) * y(i)
   end do

   res = zdotc(3, x, 1, y, 1)

   if (abs(res - ref) > 1.0e-10_dp) stop 1
end program test_zdotc
