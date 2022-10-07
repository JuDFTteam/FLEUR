!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_invert_HepsS

   IMPLICIT NONE

CONTAINS
   SUBROUTINE invert_HepsS(fmpi, atoms, noco, juPhon, lapwkpr, zMatkpr, eignukpr, eignuk, nekpr, nocck, l_real, invE, matE)
      !! Subroutine to calculate \((H-\epsilon_{\nu k} S)^{-1}\) as
      !! \(z_{k'}(\epsilon_{k'}-\epsilon_{\nu k})^(-1)z_{k'}^H\), i.e.
      !! in the spectral representation.

      USE m_types
      USE m_types_mpimat

      TYPE(t_mpi),INTENT(IN)       :: fmpi
      TYPE(t_atoms),   INTENT(IN) :: atoms
      TYPE(t_noco),    INTENT(IN) :: noco
      TYPE(t_juPhon),  INTENT(IN) :: juPhon
      TYPE(t_lapw),    INTENT(IN) :: lapwkpr
      CLASS(t_mat),    INTENT(IN) :: zMatkpr

      INTEGER,         INTENT(IN) :: nekpr, nocck
      REAL,            INTENT(IN) :: eignukpr(:), eignuk(:)
      LOGICAL,         INTENT(IN) :: l_real

      CLASS(t_mat), ALLOCATABLE, INTENT(OUT) :: invE(:)
      CLASS(t_mat), ALLOCATABLE, INTENT(OUT) :: matE(:)

      INTEGER :: nbasfcn, nu, iGpr, iG, nupr
      REAL    :: deps, invdeps

      IF (fmpi%n_size == 1) THEN
         ALLOCATE (t_mat::invE(nocck))
         ALLOCATE (t_mat::matE(nocck))
      ELSE
         ALLOCATE (t_mpimat::invE(nocck))
         ALLOCATE (t_mpimat::matE(nocck))
      END IF

      nbasfcn = MERGE(lapwkpr%nv(1)+lapwkpr%nv(2)+2*atoms%nlotot,lapwkpr%nv(1)+atoms%nlotot,noco%l_noco)
      DO nu = 1, nocck
         CALL invE(nu)%init(.TRUE., nekpr, nekpr)
         CALL matE(nu)%init(.TRUE., nekpr, nekpr)
         !DO nupr = 1, nocck
         DO nupr = 1, nekpr
         !DO nupr = nocck + 1, nekpr
            deps = eignukpr(nupr)-eignuk(nu)
            matE(nu)%data_r(nupr,nupr) = deps
            IF (ABS(deps)<juPhon%eDiffcut) CYCLE
            !IF (ABS(deps)<1e-12) deps = 1e-12
            invdeps = 1.0 / deps
            invE(nu)%data_r(nupr,nupr) = invdeps
         END DO
      END DO
   END SUBROUTINE invert_HepsS
 END MODULE m_invert_HepsS
