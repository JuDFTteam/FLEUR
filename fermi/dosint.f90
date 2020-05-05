!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dosint
   !
   !     integrated dos to ei
   !
   USE m_trisrt
   USE m_types

   IMPLICIT NONE

   CONTAINS
   SUBROUTINE dosint(ei,nemax,jspins,kpts,sfac,eig,ct)

      INTEGER,       INTENT(IN)  :: jspins
      REAL,          INTENT(IN)  :: ei,sfac
      TYPE(t_kpts),  INTENT(IN)  :: kpts
      INTEGER,       INTENT(IN)  :: nemax(:)
      REAL,          INTENT(IN)  :: eig(:,:,:)    !(neig,nkpt,jspins)
      REAL,          INTENT(OUT) :: ct

      INTEGER :: jspin,iBand,itria
      INTEGER :: k1,k2,k3
      INTEGER :: neig
      REAL    :: e1,e2,e3
      REAl    :: ee,e32,e31,e21,s

      s = 0.0
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
                  s = s + kpts%voltet(itria)/kpts%ntet
               ELSEIF (ei.GT.e2) THEN
                  e31 = e3 - e1
                  e32 = e3 - e2
                  ee = e3 - ei
                  s = s + kpts%voltet(itria)/kpts%ntet &
                         * (1.0-ee*ee/ (e31*e32))
               ELSE
                  e21 = e2 - e1
                  e31 = e3 - e1
                  ee = ei - e1
                  s = s + kpts%voltet(itria)/kpts%ntet &
                         * ee*ee/ (e21*e31)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      !Take into account spin-degeneracy
!jr   ct=2.*s
!gb   ct = (2./jspins)*s
      ct = sfac * s

   END SUBROUTINE dosint
END MODULE m_dosint
