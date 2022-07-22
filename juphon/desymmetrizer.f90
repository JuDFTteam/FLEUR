!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_desymmetrizer
   USE m_types

   IMPLICIT NONE

CONTAINS
   SUBROUTINE desymmetrize_pw(sym, stars, stars_nosym, rhopw, rhopw_nosym)
      USE m_spgrot

      TYPE(t_sym),   INTENT(IN) :: sym
      TYPE(t_stars), INTENT(IN) :: stars, stars_nosym

      COMPLEX, INTENT(IN)     :: rhopw(:,:)
      COMPLEX, INTENT(INOUT)  :: rhopw_nosym(:,:)

      INTEGER :: iStar, iStar_nosym, iSym
      INTEGER :: kr(3,sym%nop)

      DO iStar = 1, stars%ng3
         CALL spgrot(sym%nop, sym%symor, sym%mrot, sym%tau, sym%invtab, stars%kv3(:, iStar), kr)
         DO iSym = 1, sym%nop
            iStar_nosym = stars_nosym%ig(kr(1,iSym),kr(2,iSym),kr(3,iSym))
            rhopw_nosym(iStar_nosym,:) = rhopw(iStar,:) * stars%rgphs(kr(1,iSym),kr(2,iSym),kr(3,iSym))
         END DO
      END DO

   END SUBROUTINE

   SUBROUTINE desymmetrize_mt(sym, sym_nosym, cell, atoms, atoms_nosym, sphhar, sphhar_nosym, rhomt, rhomt_nosym)
      USE m_dwigner

      TYPE(t_sym),    INTENT(IN) :: sym, sym_nosym
      TYPE(t_cell),   INTENT(IN) :: cell
      TYPE(t_atoms),  INTENT(IN) :: atoms, atoms_nosym
      TYPE(t_sphhar), INTENT(IN) :: sphhar, sphhar_nosym

      REAL, INTENT(IN)    :: rhomt(:,0:,:,:)
      REAL, INTENT(INOUT) :: rhomt_nosym(:,0:,:,:)

      INTEGER :: iAtom_new, iAtom_old, i_test, iType_old, nd_old, nd_new, iOp, m_wigner
      INTEGER :: iLH_new, llh_new, iMem_new, mlh_new, iLH_old, llh_old, iMem_old, mlh_old
      REAL    :: tau_new(3), tau_old(3)
      COMPLEX :: clnu_new, clnu_old, d_wigner_elem

      COMPLEX :: d_wigner_full(-atoms%lmaxd:atoms%lmaxd, -atoms%lmaxd:atoms%lmaxd, 0:atoms%lmaxd, sym%nop)

      CALL d_wigner(sym%nop, sym%mrot, cell%bmat, atoms%lmaxd, d_wigner_full(:, :, 1:, :sym%nop))
      d_wigner_full(:, :, 0, :) = 1

      DO iAtom_new = 1, atoms_nosym%ntype ! Same as atoms_nosym%nat
         tau_new = atoms_nosym%pos(:, iAtom_new) ! Position of this atom in the unsymmetrized system

         DO iAtom_old = 1, atoms%nat
            tau_old = atoms%pos(:, iAtom_old)
            i_test = iAtom_old
            IF (norm2(tau_new-tau_old)<1e-5) EXIT
         END DO

         iType_old = atoms%itype(iAtom_old)

         nd_old = sym%ntypsy(iAtom_old)
         nd_new = sym_nosym%ntypsy(iAtom_new)
         iOp    = sym%ngopr(iAtom_old)

         DO iLH_new = 0, sphhar_nosym%nlh(nd_new)
            llh_new = sphhar_nosym%llh(iLH_new,nd_new)
            DO iMem_new = 1, sphhar_nosym%nmem(iLH_new,nd_new)
               mlh_new = sphhar_nosym%mlh(iMem_new,iLH_new,nd_new)
               clnu_new = sphhar_nosym%clnu(iMem_new,iLH_new,nd_new)
               DO iLH_old = 0, sphhar%nlh(nd_old)
                  llh_old = sphhar%llh(iLH_old,nd_old)
                  DO iMem_old = 1, sphhar%nmem(iLH_old,nd_old)
                     mlh_old = sphhar%mlh(iMem_old,iLH_old,nd_old)
                     clnu_old = sphhar%clnu(iMem_old,iLH_old,nd_old)
                     DO m_wigner = -llh_old, llh_old
                        IF (llh_old==llh_new.AND.m_wigner==mlh_new) THEN
                           d_wigner_elem = d_wigner_full(mlh_old, m_wigner, llh_old, iOp)
                           rhomt_nosym(:atoms%jri(iType_old),iLH_new,iAtom_new,:) = &
                           rhomt_nosym(:atoms%jri(iType_old),iLH_new,iAtom_new,:) + &
                           CONJG(clnu_new) * clnu_old * CONJG(d_wigner_elem) * &
                           rhomt(:atoms%jri(iType_old),iLH_old,iType_old,:)
                        END IF ! L'=L, m''=m'(L'M')
                     END DO ! m_wigner
                  END DO ! iMem_old
               END DO ! iLH_old
            END DO ! iMem_new
         END DO ! iLH_new
      END DO ! iAtom_new

   END SUBROUTINE
END MODULE
