      MODULE m_vacp5_0
      CONTAINS
      SUBROUTINE vacp5_0(
     >     nmzxyd,nmzxy,z1,tpi,rxy,m,delz,
     <     pvac,fact)     

c     calculates a part of a vacuum potential caused 
c     by the vacuum charge density for the g_z=0 component
c     based on the using of the Green function
c                      Y. Mokrousov
      
      USE m_qsf
      IMPLICIT NONE

      INTEGER, INTENT (IN) :: nmzxy,nmzxyd,m
      REAL,    INTENT (IN) :: z1,tpi,delz
      COMPLEX, INTENT (IN) :: rxy(nmzxyd)
      COMPLEX, INTENT (OUT):: pvac(nmzxyd)
      REAL,    INTENT (OUT):: fact(nmzxyd)

      INTEGER imz,imz1,l
      REAL    z,zp,vr(nmzxyd),vi(nmzxyd),fpi
      REAL    grfr(nmzxyd),grfi(nmzxyd)

      INTRINSIC real,aimag
      

      fpi = 2.*tpi

      DO 150 imz = 1,nmzxy
               
         z = z1 + delz*(imz-1)
         
c----------> goes to the interstitial contribution afterwards

         fact(imz) = tpi*z1/(m*((z/z1)**m))
         
c---------------------------------------------------------
         
         DO 250 imz1 = 1,nmzxy
            
c-----------------------------------> g(z,z') started
            
            zp = z1 + delz*(imz1-1)
            
            IF (imz1.LE.imz) THEN
               
c---------------------------------------> z'< z
          
               grfr(imz1) = tpi*zp*((zp/z)**m)*
     *              real(rxy(imz1))/m
               grfi(imz1) = tpi*zp*((zp/z)**m)*
     *              aimag(rxy(imz1))/m


c---------------------------------------->
            ELSE
c---------------------------------------> z'> z
            
               grfr(imz1) = tpi*zp*((z/zp)**m)*
     *              real(rxy(imz1))/m
               grfi(imz1) = tpi*zp*((z/zp)**m)*
     *              aimag(rxy(imz1))/m

c----------------------------------------> 
            END IF
            
 250     CONTINUE               ! imz1
               
c-----------------------------------> g(z,z') finished

c------>  obtaining the vacuum contribution for the potential
         
         CALL qsf(delz,grfr(1),vr,nmzxy,0)
         CALL qsf(delz,grfi(1),vi,nmzxy,0)
         
         pvac(imz) = cmplx(vr(1),vi(1))
               
 150  CONTINUE                  !  imz 
      
      END SUBROUTINE vacp5_0 
      END MODULE m_vacp5_0
