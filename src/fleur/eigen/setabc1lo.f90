!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_setabc1lo
   !*********************************************************************
   ! Calculates the (lower case) a, b and c coefficients of the local
   ! orbitals: the LO basis function a*u + b*udot + c*ulo has zero value
   ! and zero derivative at the muffin-tin boundary and is normalized.
   ! An APW LO (l_dulo) lives on the udot slot, a*u + c*udot, and only
   ! its value vanishes at the boundary.
   !*********************************************************************
   USE m_judft
   USE m_types_atoms
   USE m_types_radfun
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: setabc1lo
CONTAINS
   SUBROUTINE setabc1lo(atoms,ntyp,rf,usp,alo1,blo1,clo1)
      TYPE(t_atoms), INTENT(IN)  :: atoms
      TYPE(t_radfun),INTENT(IN)  :: rf
      INTEGER,       INTENT(IN)  :: ntyp,usp
      REAL,          INTENT(OUT) :: alo1(:,:),blo1(:,:),clo1(:,:)

      REAL    :: ws,v(SIZE(rf%integral,1))
      INTEGER :: l,lo,i

      DO lo = 1,atoms%nlo(ntyp)
         l = atoms%llo(lo,ntyp)
         i = atoms%slot_of_lo(lo,ntyp)
         ASSOCIATE(b => rf%bnd(:,:,l,usp))
            v    = 0.0
            IF (atoms%l_dulo(lo,ntyp)) THEN
               ! a*b(1,1) + b(1,2) = 0
               v(1) = -b(1,2)/b(1,1)
               v(2) = 1.0
            ELSE
               ! a*b(:,1) + b*b(:,2) = -b(:,i)
               ws   = b(1,2)*b(2,1) - b(1,1)*b(2,2)
               v(1) = (b(2,2)*b(1,i) - b(1,2)*b(2,i))/ws
               v(2) = (b(1,1)*b(2,i) - b(2,1)*b(1,i))/ws
               v(i) = 1.0
            END IF
         END ASSOCIATE
         clo1(lo,usp) = 1.0/SQRT(DOT_PRODUCT(v,MATMUL(rf%integral(:,:,l,usp,usp),v)))
         alo1(lo,usp) = v(1)*clo1(lo,usp)
         blo1(lo,usp) = MERGE(0.0,v(2)*clo1(lo,usp),atoms%l_dulo(lo,ntyp))
      END DO
   END SUBROUTINE setabc1lo
END MODULE m_setabc1lo
