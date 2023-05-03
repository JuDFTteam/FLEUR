!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_doswt
   !
   !     calculates the weights for each k-point for integrating functions
   !     of k.  the array w has beeen cleared before entering.
   !
   USE m_trisrt
   USE m_types

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE doswt(ei,nemax,jspins,kpts,eig,w)

      INTEGER,       INTENT(IN) :: jspins
      TYPE(t_kpts),  INTENT(IN) :: kpts
      REAL,          INTENT(IN) :: ei
      INTEGER,       INTENT(IN) :: nemax(:)
      REAL,          INTENT(IN) :: eig(:,:,:)      !(neig,nkpt,jspins)
      REAL,          INTENT(OUT):: w(:,:,:)        !(neig,nkpt,jspins)

      INTEGER :: jspin,iBand,itria
      INTEGER :: k1,k2,k3
      INTEGER :: neig
      REAL    :: e1,e2,e3
      REAl    :: ee,e32,e31,e21,s
      w=0.0 !init was missing
      DO jspin = 1,jspins
         neig = nemax(jspin)
         DO iBand = 1,neig
            DO itria = 1,kpts%ntet
               !Get the k-points and eigenvalues
               !of the current triangle
               k1 = kpts%ntetra(1,itria)
               k2 = kpts%ntetra(2,itria)
               k3 = kpts%ntetra(3,itria)
               e1 = eig(iBand,k1,jspin)
               e2 = eig(iBand,k2,jspin)
               e3 = eig(iBand,k3,jspin)
               !Sort by ascending eigenvalues
               CALL trisrt(e1,e2,e3,k1,k2,k3)
               IF (e1.LE.-9999.0) CYCLE !Not all eigenvalues available
               IF (ei.LE.e1) CYCLE !triangle not occupied
               IF (ei.GE.e3) THEN
                  !---> e3<e
                  s = kpts%voltet(itria)/kpts%ntet/3.0
                  w(iBand,k1,jspin) = w(iBand,k1,jspin) + s
                  w(iBand,k2,jspin) = w(iBand,k2,jspin) + s
                  w(iBand,k3,jspin) = w(iBand,k3,jspin) + s
               ELSEIF (ei.GT.e2) THEN
                  !---> e2<ei<e3
                  ee = e3 - ei
                  e31 = ee/ (e3-e1)
                  e32 = ee/ (e3-e2)
                  s = kpts%voltet(itria)/kpts%ntet/3.0
                  w(iBand,k1,jspin) = w(iBand,k1,jspin) + s* (1.-e31*e31*e32)
                  w(iBand,k2,jspin) = w(iBand,k2,jspin) + s* (1.-e31*e32*e32)
                  w(iBand,k3,jspin) = w(iBand,k3,jspin) + s* (1.-e31*e32*(3.-e31-e32))
               ELSE
                  !---> e1<ei<e2
                  ee = ei - e1
                  e31 = ee/ (e3-e1)
                  e21 = ee/ (e2-e1)
                  s = kpts%voltet(itria)/kpts%ntet*e31*e21/3.0
                  w(iBand,k1,jspin) = w(iBand,k1,jspin) + s* (3.0-e21-e31)
                  w(iBand,k2,jspin) = w(iBand,k2,jspin) + s*e21
                  w(iBand,k3,jspin) = w(iBand,k3,jspin) + s*e31
               ENDIF

            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE doswt
END MODULE m_doswt
