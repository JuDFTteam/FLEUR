MODULE m_lattHarmsSphHarmsConv

   USE m_types

CONTAINS

   SUBROUTINE lattHarmsRepToSphHarms(sym, atoms, lattHarms, iType, funcLattHarms, funcSphHarms)

      IMPLICIT NONE

      TYPE(t_sym),    INTENT(IN)    :: sym
      TYPE(t_atoms),  INTENT(IN)    :: atoms
      TYPE(t_sphhar), INTENT(IN)    :: lattHarms
      INTEGER,        INTENT(IN)    :: iType
      REAL,           INTENT(IN)    :: funcLattHarms(:,0:) ! (iR,iLH)
      COMPLEX,        INTENT(INOUT) :: funcSphHarms(:,:) ! (iR,lm)

      INTEGER :: iAtom, iLH, ns, l, iM, m, lm, iR

      iAtom = atoms%firstAtom(iType)
      ns = sym%ntypsy(iAtom)

      funcSphHarms = CMPLX(0.0,0.0)

      DO iLH = 0, lattHarms%nlh(ns)
         l = lattHarms%llh(iLH,ns)
         DO iM = 1, lattHarms%nmem(iLH,ns)
            m = lattHarms%mlh(iM,iLH,ns)
            lm = l*(l+1) + m + 1
            DO iR = 1, atoms%jri(iType)
               funcSphHarms(iR,lm) = funcSphHarms(iR,lm) + funcLattHarms(iR,iLH) * lattHarms%clnu(iM,iLH,ns)
            END DO
         END DO
      END DO

   END SUBROUTINE lattHarmsRepToSphHarms

SUBROUTINE sphHarmsRepToLattHarms(sym, atoms, lattHarms, iType, funcSphHarms, funcLattHarms)

   IMPLICIT NONE

      TYPE(t_sym),    INTENT(IN)    :: sym
      TYPE(t_atoms),  INTENT(IN)    :: atoms
      TYPE(t_sphhar), INTENT(IN)    :: lattHarms
      INTEGER,        INTENT(IN)    :: iType
      COMPLEX,        INTENT(IN)    :: funcSphHarms(:,:) ! (iR,lm)
      REAL,           INTENT(INOUT) :: funcLattHarms(:,0:) ! (iR,iLH)

      INTEGER :: iAtom, iLH, ns, l, iM, m, lm, iR

      iAtom = atoms%firstAtom(iType)
      ns = sym%ntypsy(iAtom)

      funcLattHarms = 0.0

      DO iLH = 0, lattHarms%nlh(ns)
         l = lattHarms%llh(iLH,ns)
         DO iM = 1, lattHarms%nmem(iLH,ns)
            m = lattHarms%mlh(iM,iLH,ns)
            lm = l*(l+1) + m + 1
            DO iR = 1, atoms%jri(iType)
               funcLattHarms(iR,iLH) = funcLattHarms(iR,iLH) + REAL(funcSphHarms(iR,lm) * CONJG(lattHarms%clnu(iM,iLH,ns)))
            END DO
         END DO
      END DO

   END SUBROUTINE sphHarmsRepToLattHarms

END MODULE m_lattHarmsSphHarmsConv
