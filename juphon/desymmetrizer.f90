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

      INTEGER :: iAtom_new, iAtom_old, iType_old, nd_old, nd_new, iOp, m_wigner
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

   SUBROUTINE desymmetrize_types(input, input_nosym, atoms, atoms_nosym, noco, nococonv, nococonv_nosym, enpara, enpara_nosym, results, results_nosym)
      USE m_types_lapw

      TYPE(t_input),    INTENT(IN) :: input, input_nosym
      TYPE(t_atoms),    INTENT(IN) :: atoms, atoms_nosym
      TYPE(t_noco),     INTENT(IN) :: noco
      TYPE(t_nococonv), INTENT(IN) :: nococonv
      TYPE(t_enpara),   INTENT(IN) :: enpara
      TYPE(t_results),  INTENT(IN) :: results
      TYPE(t_nococonv), INTENT(INOUT) :: nococonv_nosym
      TYPE(t_enpara),   INTENT(INOUT) :: enpara_nosym
      TYPE(t_results),  INTENT(INOUT) :: results_nosym

      INTEGER :: neigd2, neigd2_nosym, iAtom_new, iAtom_old, iType_old
      REAL    :: tau_new(3), tau_old(3)

      ! TODO: Thes two should be identical!
      neigd2       = MIN(input%neig,lapw_dim_nbasfcn)
      neigd2_nosym = MIN(input_nosym%neig,lapw_dim_nbasfcn)
      IF (neigd2/=neigd2_nosym) WRITE(*,*) "neigd2 /= itself!!"

      IF (noco%l_soc.AND.(.NOT.noco%l_noco)) neigd2 = 2*neigd2

      ! Scalar/presized array quantities:
      nococonv_nosym%theta = nococonv%theta
      nococonv_nosym%phi   = nococonv%phi
      nococonv_nosym%qss   = nococonv%qss

      enpara_nosym%evac      = enpara%evac
      enpara_nosym%evac1     = enpara%evac1
      enpara_nosym%enmix     = enpara%enmix
      enpara_nosym%lchg_v    = enpara%lchg_v
      enpara_nosym%epara_min = enpara%epara_min
      enpara_nosym%ready     = enpara%ready
      enpara_nosym%floating  = enpara%floating

      results_nosym%ef       = results%ef
      results_nosym%seigc    = results%seigc
      results_nosym%seigv    = results%seigv
      results_nosym%ts       = results%ts
      results_nosym%te_vcoul = results%te_vcoul
      results_nosym%te_veff  = results%te_veff
      results_nosym%te_exc   = results%te_exc
      results_nosym%e_ldau   = results%e_ldau
      results_nosym%e_ldaopc = results%e_ldaopc
      results_nosym%e_vdw    = results%e_vdw
      results_nosym%tote     = results%tote
      results_nosym%bandgap  = results%bandgap
      results_nosym%te_hfex  = results%te_hfex

      results_nosym%te_hfex_loc         = results%te_hfex_loc
      results_nosym%last_distance       = results%last_distance
      results_nosym%last_mmpMatdistance = results%last_mmpMatdistance
      results_nosym%last_occdistance    = results%last_occdistance

      ! Allocated arrays:
      results_nosym%unfolding_weights = results%unfolding_weights
      results_nosym%w_iks             = results%w_iks
      results_nosym%eig               = results%eig
      results_nosym%neig              = results%neig
      IF(input%l_rdmft) THEN
         results_nosym%w_iksRDMFT = results_nosym%w_iksRDMFT
      END IF

      ! Atom loop:
      DO iAtom_new = 1, atoms_nosym%ntype ! Same as atoms_nosym%nat
         tau_new = atoms_nosym%pos(:, iAtom_new) ! Position of this atom in the unsymmetrized system

         DO iAtom_old = 1, atoms%nat
            tau_old = atoms%pos(:, iAtom_old)
            IF (norm2(tau_new-tau_old)<1e-5) EXIT
         END DO

         iType_old = atoms%itype(iAtom_old)

         enpara_nosym%el0(:,iAtom_new,:)   = enpara%el0(:,iType_old,:)
         enpara_nosym%el1(:,iAtom_new,:)   = enpara%el1(:,iType_old,:)
         enpara_nosym%ello0(:,iAtom_new,:) = enpara%ello0(:,iType_old,:)
         enpara_nosym%ello1(:,iAtom_new,:) = enpara%ello1(:,iType_old,:)

         enpara_nosym%skiplo(iAtom_new,:)    = enpara%skiplo(iType_old,:)
         enpara_nosym%lchange(:,iAtom_new,:) = enpara%lchange(:,iType_old,:)
         enpara_nosym%llochg(:,iAtom_new,:)  = enpara%llochg(:,iType_old,:)

         ! TODO: This is most DEFINITELY faulty, but we shouldn't fix it until
         !       the noco rotation logic itself is 100% cleaned up.
         IF (noco%l_noco) THEN
            nococonv_nosym%alph(iAtom_new)     = nococonv%alph(iType_old)
            nococonv_nosym%alphRlx(iAtom_new)  = nococonv%alphRlx(iType_old)
            nococonv_nosym%alphPrev(iAtom_new) = nococonv%alphPrev(iType_old)
            nococonv_nosym%beta(iAtom_new)     = nococonv%beta(iType_old)
            nococonv_nosym%betaRlx(iAtom_new)  = nococonv%betaRlx(iType_old)
            nococonv_nosym%betaPrev(iAtom_new) = nococonv%betaPrev(iType_old)

            nococonv_nosym%b_con(2,iAtom_new) = nococonv%b_con(2,iType_old)
         END IF
      END DO

      ! Omitted:
      ! results%force already exists as a desymmetrization function in force_w(?)
      ! results%force_old/_vdw as above

   END SUBROUTINE
END MODULE
