!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_invert_HepsS

   IMPLICIT NONE

CONTAINS
   SUBROUTINE invert_HepsS(fmpi, atoms, noco, juPhon, lapwkpr, zMatkpr, eignukpr, eignuk, nekpr, nocck, l_real, invHepsS, invE)
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

      CLASS(t_mat), ALLOCATABLE, INTENT(OUT) :: invHepsS(:)
      CLASS(t_mat), ALLOCATABLE, INTENT(OUT) :: invE(:)

      INTEGER :: nbasfcn, nu, iGpr, iG, nupr
      REAL    :: deps, invdeps

      IF (fmpi%n_size == 1) THEN
         ALLOCATE (t_mat::invHepsS(nocck))
         ALLOCATE (t_mat::invE(nocck))
      ELSE
         ALLOCATE (t_mpimat::invHepsS(nocck))
         ALLOCATE (t_mpimat::invE(nocck))
      END IF

      nbasfcn = MERGE(lapwkpr%nv(1)+lapwkpr%nv(2)+2*atoms%nlotot,lapwkpr%nv(1)+atoms%nlotot,noco%l_noco)
      DO nu = 1, nocck
         CALL invHepsS(nu)%init(l_real, nbasfcn, nbasfcn)
         CALL invE(nu)%init(.TRUE., nekpr, nekpr)
         DO nupr = 1, nekpr
            deps = eignukpr(nupr)-eignuk(nu)
            IF (ABS(deps)<juPhon%eDiffcut) CYCLE
            invdeps = 1.0 / deps
            invE(nu)%data_r(nupr,nupr) = invdeps
            DO iG = 1, nbasfcn
               DO iGpr = 1, nbasfcn
                  IF (l_real) THEN
                     invHepsS(nu)%data_r(iGpr,iG) = invHepsS(nu)%data_r(iGpr,iG) &
                                                  + zMatkpr%data_r(iGpr, nupr) * &
                                                  & zMatkpr%data_r(iG, nupr) * invdeps
                  ELSE
                     invHepsS(nu)%data_c(iGpr,iG) = invHepsS(nu)%data_c(iGpr,iG) &
                                                  + zMatkpr%data_c(iGpr, nupr) * &
                                                  & CONJG(zMatkpr%data_c(iG, nupr)) * invdeps
                  END IF
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE invert_HepsS
 END MODULE m_invert_HepsS
