!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_visp5_0
      CONTAINS
      SUBROUTINE visp5_0(
     >             nmzxyd,nmzxy,delz,m,ivfft1,ivfft2,l,
     >             rxy,ani1,ani2,z1,amat,
     >             pvac,pint,tpi,jspd,val,
     <             vis_help)
      
      USE m_qsf
      USE m_angle
      IMPLICIT NONE

c---- this subroutine generates the gz=0 components of the 
c---- \tilde\tilde{V} potential on the 
c --- uniform grid in the unit cell \tilde{\Omega} for further 2-
c---- dimensional fast fourier transform 
c---- for the further understanding look at the vvacxy_4.F      
c----                           Y.Mokrousov

      INTEGER, INTENT (IN) :: nmzxyd,nmzxy,m,ivfft1,ivfft2,jspd,l
      REAL,    INTENT (IN) :: ani1,ani2,delz,z1,tpi
      COMPLEX, INTENT (IN) :: val
      REAL,    INTENT (IN) :: amat(3,3)
      COMPLEX, INTENT (IN) :: rxy(nmzxyd)
      COMPLEX, INTENT (IN) :: pvac(nmzxyd)
      COMPLEX, INTENT (IN) :: pint(nmzxyd)
      COMPLEX, INTENT (OUT)::
     &            vis_help(0:ivfft1-1,0:ivfft2-1)

      INTEGER imz,ix,iy,im,mult
      REAL    x,y,r,phi,zf,q

      mult = (-1)**(l+1)

      DO ix = 0,ivfft1 - 1
         DO iy = 0,ivfft2 - 1
             x = ix*ani1
             IF (x.GT.0.5) x = x - 1.
             y = iy*ani2
             IF (y.GT.0.5) y = y - 1.
             r = sqrt((x*amat(1,1) + y*amat(1,2))**2 +
     +                (x*amat(2,1) + y*amat(2,2))**2)
             phi = angle(x*amat(1,1) + y*amat(1,2),
     >                   x*amat(2,1) + y*amat(2,2))
            IF (r.GT.z1) THEN
               
               zf = (r-z1)/delz + 1.0
               im = zf
               q = zf - im
               vis_help(ix,iy) = (0.5* (q-1.)*
     *              (q-2.)*(pvac(im) + pint(im)) -
     +              q* (q-2.)*(pvac(im+1) + 
     +              pint(im+1)) +
     +              0.5*q* (q-1.)*(pvac(im+2) + 
     +              pint(im+2)))*
     *              exp(cmplx(0.,mult*m*phi))
               
            ELSE
               
               vis_help(ix,iy) = val*exp(cmplx(0.,mult*m*phi))*
     *                ((r**m)/(z1**m))

c------ here the contribution from the zero-fouirer component of
c------ the charge density is coming

               
            END IF
         END DO                 ! ix
      END DO                    ! iy
      
      END SUBROUTINE visp5_0
      END MODULE m_visp5_0
