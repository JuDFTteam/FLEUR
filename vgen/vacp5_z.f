!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_vacp5_z
      CONTAINS
      SUBROUTINE vacp5_z(
     >     nmzxyd,nmzxy,z1,delz,fpi,II,KK,rxy,m,
     <     pvac)

c     calculates a part of a vacuum potential caused 
c     by the vacuum charge density for the g_z.ne.0 component
c     based on the using of the Green function
c                      Y. Mokrousov
      
      USE m_modcyli
      USE m_modcylk
      USE m_qsf
      IMPLICIT NONE

      INTEGER, INTENT (IN) :: nmzxy,nmzxyd,m
      REAL,    INTENT (IN) :: z1,delz,fpi
      REAL,    INTENT (IN) :: II(nmzxyd),KK(nmzxyd)
      COMPLEX, INTENT (IN) :: rxy(nmzxyd)
      COMPLEX, INTENT (OUT):: pvac(nmzxyd)

      INTEGER imz,imz1
      REAL    z,zp,vr(nmzxyd),vi(nmzxyd)
      REAL    grfr(nmzxyd),grfi(nmzxyd)

      DO 200 imz = 1,nmzxy
         
         z = z1 + delz*(imz-1)
         
         DO 300 imz1 = 1,nmzxy
            
            zp = z1 + delz*(imz1-1)


            IF (imz1.LE.imz) THEN

               grfr(imz1) = fpi*II(imz1)*KK(imz)*zp*
     *              real(rxy(imz1))
               grfi(imz1) = fpi*II(imz1)*KK(imz)*zp*
     *              aimag(rxy(imz1))
            
            ELSE
               
               grfr(imz1) = fpi*II(imz)*KK(imz1)*zp*
     *              real(rxy(imz1))
               grfi(imz1) = fpi*II(imz)*KK(imz1)*zp*
     *              aimag(rxy(imz1))

            END IF

 300     CONTINUE               ! imz1
         
         CALL qsf(delz,grfr(1),vr,nmzxy,0)
         CALL qsf(delz,grfi(1),vi,nmzxy,0)

         pvac(imz) = cmplx(vr(1),vi(1))
      
 200  CONTINUE                  ! imz

      END SUBROUTINE vacp5_z
      END MODULE m_vacp5_z
