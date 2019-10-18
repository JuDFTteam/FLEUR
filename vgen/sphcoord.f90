!--------------------------------------------------------------------------------
! Copyright (c) 2019 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_sphcoord
   
   !--------------------------------------------------------------------------------
   ! Calculates the polar angles theta and phi of a vector with components
   ! x, y and z.
   !
   ! Based on pol_angle by P. Kurz.  
   !--------------------------------------------------------------------------------
   
   CONTAINS
      SUBROUTINE sphcoord(x,y,z,theta,phi)

      USE m_constants

      IMPLICIT NONE

      REAL, INTENT(IN)  :: x, y, z
      REAL, INTENT(OUT) :: theta, phi
      
      REAL, PARAMETER   :: eps = 1.0e-15

      REAL r, rho

      rho  = SQRT(x**2+y**2)
      r = SQRT(x**2+y**2+z**2)

      IF (r.LT.eps) THEN
         theta = 0.0
         phi   = 0.0
      ELSE
         theta = ACOS(z/r)
         IF (rho.LT.eps) THEN
            phi = 0.0
         ELSE IF ((x/rho).LT.eps) THEN
            phi = SIGN(1.0, y)*pi_const/2.0
         ELSE IF (x.GT.0.0) THEN
            phi = ATAN(y/x)!pi/2 - SIGN(1.0, vx)*pi/2
         ELSE IF (y.GE.0.0) THEN
            phi = ATAN(y/x)+pi_const
         ELSE 
            phi = ATAN(y/x)-pi_const
         
         ENDIF
      ENDIF

   END SUBROUTINE sphcoord
END MODULE m_sphcoord
