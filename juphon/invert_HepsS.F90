!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_invert_HepsS

   IMPLICIT NONE

CONTAINS
   SUBROUTINE invert_HepsS(fmpi, atoms, noco, juPhon, lapwkpr, zMatkpr, eignukpr, eignuk, nekpr, nocck, l_real, invE, nocckq, wk, &
                           wkq, matOcc, ikpt)
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

      INTEGER,         INTENT(IN) :: nekpr, nocck, nocckq, ikpt
      REAL,            INTENT(IN) :: eignukpr(:), eignuk(:), wk(:), wkq(:)
      LOGICAL,         INTENT(IN) :: l_real

      CLASS(t_mat), ALLOCATABLE, INTENT(OUT) :: invE(:)
      CLASS(t_mat), ALLOCATABLE, INTENT(OUT) :: matOcc(:)

      INTEGER :: nbasfcn, nu, iGpr, iG, nupr
      REAL    :: deps, invdeps, dwks

      IF (fmpi%n_size == 1) THEN
         ALLOCATE (t_mat::invE(nocck))
         ALLOCATE (t_mat::matOcc(nocck))
      ELSE
         ALLOCATE (t_mpimat::invE(nocck))
         ALLOCATE (t_mpimat::matOcc(nocck))
      END IF

      nbasfcn = MERGE(lapwkpr%nv(1)+lapwkpr%nv(2)+2*atoms%nlotot,lapwkpr%nv(1)+atoms%nlotot,noco%l_noco)
      DO nu = 1, nocck
         CALL invE(nu)%init(.TRUE., nekpr, nekpr)
         CALL matOcc(nu)%init(.TRUE., nocckq, nocckq)
         DO nupr = 1, nekpr
            deps = eignukpr(nupr)-eignuk(nu)
            write(6666,*) eignukpr(nupr)
            write(6666,*) eignuk(nu)
            write(6666,*) deps
            write(6666,*) "-------"
            IF (ABS(deps)<juPhon%eDiffcut) CYCLE
            !IF (ikpt==1) CYCLE
            !IF (ABS(deps)<1e-7) CYCLE
            invdeps = 1.0 / deps
            invE(nu)%data_r(nupr,nupr) = invdeps
            IF (nupr<=nocck) THEN
               dwks = wkq(nupr) - wk(nu)
               !IF (ABS(dwks)>=juPhon%fDiffcut) matOcc(nu)%data_r(nupr,nupr) = 1.0
            END IF
         END DO
      END DO
   END SUBROUTINE invert_HepsS
 END MODULE m_invert_HepsS
