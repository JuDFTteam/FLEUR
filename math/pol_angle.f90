!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_polangle
   !-----------------------------------------------------------------------------
   ! Calculates the polar angles theta and phi of a vector with components
   ! vx, vy and vz.
   ! Philipp Kurz 2000-02-08
   ! Modernised A.N. 2020
   !-----------------------------------------------------------------------------
CONTAINS

    SUBROUTINE pol_angle(vx, vy, vz, theta, phi,l_minimize)
        USE m_constants, ONLY : pimach
        IMPLICIT NONE

        REAL, INTENT(IN)  :: vx, vy, vz
        REAL, INTENT(OUT) :: theta, phi
        LOGICAL,INTENT(IN),OPTIONAL :: l_minimize

        REAL :: eps, r, rho, pi

        pi = pimach()
        eps = 1.0e-9

        rho = SQRT(vx**2 + vy**2)
        r   = SQRT(vx**2 + vy**2 + vz**2)

        ! Zero vector:
        IF ( (r.LT.eps) .AND. (rho.LT.eps) ) THEN
           theta = 0.0
           phi   = 0.0
        ! v on positive z-axis:
        ELSE IF ( (rho.LT.eps) .AND. vz.GT.0) THEN
           theta = 0.0
           phi   = 0.0
        ! v on negative z-axis:
        ELSE IF ( (rho.LT.eps) .AND. vz.LT.0) THEN
           theta = pi
           phi   = 0.0
        ELSE
           ! v in xy-plane:
           IF (ABS(vz).LT.eps) THEN
              theta = pi/2
           ELSE
              theta = ASIN(rho/r)
           END IF

           IF ( vz.LT.0 ) THEN
              theta = pi - theta
           END IF

           ! v in yz-plane
           IF (ABS(vx).LT.eps) THEN
              phi = pi/2

           ELSE
              phi = ASIN(ABS(vy)/rho)
           END IF

           IF ( vx.LT.0 ) phi = pi - phi
           IF ( vy.LT.0 ) phi = -phi

           if (present(l_minimize))Then
             if (l_minimize) THEN
               !Make phi and theta minimal
               if (abs(phi)>pi/2) THEN
                 phi=phi-sign(pi,phi)
                 theta=-theta
               endif
             ENDIF
           endif
        END IF



    END SUBROUTINE pol_angle

    SUBROUTINE sphericaltocart(r, theta, phi, x, y, z)
       IMPLICIT NONE

       REAL, INTENT(IN)  :: r, theta, phi
       REAL, INTENT(OUT) :: x, y, z

       x=r*SIN(theta)*COS(phi)
       y=r*SIN(theta)*SIN(phi)
       z=r*COS(theta)
    END SUBROUTINE sphericaltocart

END MODULE m_polangle
