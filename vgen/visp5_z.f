      MODULE m_visp5_z
      CONTAINS
      SUBROUTINE visp5_z(
     >             nmzxyd,nmzxy,delz,m,ivfft1,ivfft2,IIIR,
     >             rxy,ani1,ani2,z1,amat,pvac,pint,tpi,l,
     >             fpi,val_help,III,m1,gz,rmap,rr,
     <             vis_help)

      USE m_qsf
      USE m_angle
      USE m_modcyli
      IMPLICIT NONE

c---- this subroutine generates the gz.ne.0 components of the 
c---- \tilde\tilde{V} potential on the 
c --- uniform grid in the unit cell \tilde{\Omega} for further 2-
c---- dimensional fast fourier transform 
c---- for the further understanding look at the vvacxy_4.F      
c----                        Y.Mokrousov

      INTEGER, INTENT (IN) :: nmzxyd,nmzxy,m,ivfft1,ivfft2,l,m1
      REAL,    INTENT (IN) :: ani1,ani2,delz,z1,IIIR,tpi,fpi
      REAL,    INTENT (IN) :: III(ivfft1*ivfft2)
      REAL,    INTENT (IN) :: amat(3,3),gz
      COMPLEX, INTENT (IN) :: rxy(nmzxyd)
      COMPLEX, INTENT (IN) :: pvac(nmzxyd)
      COMPLEX, INTENT (IN) :: pint(nmzxyd)
      COMPLEX, INTENT (IN) :: val_help
      INTEGER, INTENT (IN) :: rmap(0:ivfft1-1,0:ivfft2-1)
      REAL,    INTENT (IN) :: rr(ivfft1*ivfft2)
      COMPLEX, INTENT (OUT)::
     &           vis_help(0:ivfft1-1,0:ivfft2-1)

      INTEGER imz,ix,iy,im,mult,irc
      REAL    x,y,r,phi,zf,q

      INTRINSIC real,aimag

      mult = (-1)**(l+1)

      DO ix = 0,ivfft1 - 1
         DO iy = 0,ivfft2 - 1
            x = ix*ani1
            IF (x.GT.0.5) x = x - 1.
            y = iy*ani2
            IF (y.GT.0.5) y = y - 1.
               r = sqrt((x*amat(1,1) + y*amat(1,2))**2 +
     +                  (x*amat(2,1) + y*amat(2,2))**2)
               phi = angle(x*amat(1,1) + y*amat(1,2),
     >                     x*amat(2,1) + y*amat(2,2))
            irc = rmap(ix,iy)
            r = rr(irc)
            IF (r.GT.z1) THEN
               zf = (r-z1)/delz + 1.0
               im = zf
               q = zf - im
               vis_help(ix,iy) = 
     +              (0.5* (q-1.)*
     *              (q-2.)*(pvac(im)+pint(im)) -
     +              q* (q-2.)*(pvac(im+1)
     +              +pint(im+1))+
     +              0.5*q* (q-1.)*(pvac(im+2)
     +              +pint(im+2)))*
     +              exp(cmplx(0.,mult*m*phi))
            ELSE  
               vis_help(ix,iy) = 
     *              III(irc)*(val_help)*
     *              exp(cmplx(0.,m1*phi))/IIIR

            END IF   
            
         END DO                 ! ix
      END DO                    ! iy
      
      END SUBROUTINE visp5_z
      END MODULE m_visp5_z
