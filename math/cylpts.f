!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_cylpts
      CONTAINS
      SUBROUTINE cylpts(x,n,rad)
c     *********************************************************
c     YM: generate random points on the cylindrical vacuum boundary
c     *********************************************************
      USE m_constants, ONLY : tpi_const
      USE m_qranf
      IMPLICIT NONE

      INTEGER n
      REAL    rad
      REAL    x(3,n)
      REAL    phi,t,tc,y1,x1,xr,yr
      INTEGER r,j,i

      INTRINSIC sqrt

     
      yr = sqrt(7.)

      j = 0

      DO 10 i = 1,n
         x(3,i) = qranf(yr,j) - 0.5
   10 CONTINUE

      j = 0
      xr = sqrt(13.e0)
      yr = sqrt(7.e0)
      DO 20 i = 1,n
         tc = 2.e0*qranf(xr,j) - 1.e0
         phi = tpi_const*qranf(yr,j)
         t = sqrt(1.e0-tc*tc)
         x1 = cos(phi)
         y1 = sin(phi)
         x(1,i) = rad*x1
         x(2,i) = rad*y1
 20   CONTINUE

      END SUBROUTINE cylpts
      END MODULE m_cylpts
