!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dosef
   !
   !---  >    obtain dos at ei (here: ef)
   !
   USE m_constants
   USE m_types
   USE m_trisrt

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE dosef(ei,nemax,jspins,kpts,sfac,eig,l_output)

      INTEGER,       INTENT(IN) :: jspins
      TYPE(t_kpts),  INTENT(IN) :: kpts
      REAL,          INTENT(IN) :: ei,sfac
      INTEGER,       INTENT(IN) :: nemax(:)
      REAL,          INTENT(IN) :: eig(:,:,:) !(neig,nkpt,jspins)
      LOGICAL,INTENT(IN) :: l_output

      REAL     :: e1,e2,e21,e3,e31,e32,s
      INTEGER  :: iBand,jspin,k1,k2,k3,itria,neig

      DO jspin = 1,jspins
         neig = nemax(jspin)
         s = 0.0
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
               IF (e1.LE.-9999.0) cycle !Not all eigenvalues available
               IF ((ei.LT.e1) .OR. (ei.GE.e3)) cycle !triangle not contributing
               IF (ei.GT.e2) THEN
                  !--->  e2<ei<e3
                  e31 = e3 - e1
                  e32 = e3 - e2
                  s = s + 2.0*kpts%voltet(itria)/kpts%ntet &
                         * (e3-ei)/ (e31*e32)
               ELSE
                  !--->  e1<ei<e2
                  e31 = e3 - e1
                  e21 = e2 - e1
                  s = s + 2.0*kpts%voltet(itria)/kpts%ntet &
                         * (ei-e1)/ (e31*e21)
               ENDIF
            ENDDO
         ENDDO
         !gb  s = (2./jspins)*s
         s = sfac * s
         if (l_output) then
            WRITE (oUnit,FMT=8000) ei,jspin,s
8000        FORMAT (/,10x,'density of states at',f12.6,&
                    ' har for spin',i2,'=',e20.8,' states/har')
         end if
      ENDDO

   END SUBROUTINE dosef
END MODULE m_dosef
