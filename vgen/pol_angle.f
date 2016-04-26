      MODULE m_polangle
c***********************************************************************
c calculates the polar angle theta and phi of a vector with components
c vx, vy and vz.
c Philipp Kurz 2000-02-08
c***********************************************************************
      CONTAINS
      SUBROUTINE pol_angle(
     >                     vx,vy,vz,
     <                     theta,phi)

      USE m_constants, ONLY : pimach
      IMPLICIT NONE

C     .. Scalar Arguments ..
      REAL, INTENT    (IN) :: vx,vy,vz
      REAL, INTENT   (OUT) :: theta,phi
C     ..
C     .. Local Scalars ..
      REAL eps,vxyz,vxy,pi
C     ..

      pi = pimach()
      eps = 1.0e-8

      vxy  = sqrt(vx**2 + vy**2)
      vxyz = sqrt(vx**2 + vy**2 + vz**2)

      IF ( (vxyz.LT.eps) .OR. (vxy.LT.eps) ) THEN
         theta = 0.0
         phi   = 0.0
      ELSE
c--->    due to rounding errors vxy/vxyz can become >1, if vz is very
c--->    small. therefore, make sure that vz is not to small.
         IF (abs(vz).LT.eps) THEN
            theta = pi/2
         ELSE
            theta = asin(vxy/vxyz)
         ENDIF
         IF ( vz.LT.0 ) theta = pi - theta

c--->    due to rounding errors vy/vxy can become >1, if vx is very
c--->    small. therefore, make sure that vx is not to small.
         IF (abs(vx).LT.eps) THEN
            phi = pi/2
         ELSE
            phi   = asin(abs(vy)/vxy)
         ENDIF
         IF ( vx.LT.0 ) phi = pi - phi
         IF ( vy.LT.0 ) phi = -phi
      ENDIF

      END SUBROUTINE pol_angle
      END MODULE m_polangle
