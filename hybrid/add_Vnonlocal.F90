!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_add_vnonlocal
   USE m_judft
! c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c
!     This module is the driver routine for the calculation of the Hartree    c
!     Fock exchange term by using the mixed basis set.                        c
!                                                                             c
!     hsfock                                                                  c
!         |                                                                   c
!         |- symm.F:                                                          c
!         |  calculates the irreducible representation                        c
!         |                                                                   c
!         |- wavefproducts.F:                 s      s*                       c
!         |  computes the repsentation of phi    phi       in the mixed basis c
!         |                                  n,k    n',k+q                    c
!         |                                                                   c
!         |- exchange.F:                                                      c
!         |  calculates valence-valence part of the exchange matrix (mat_ex), c
!         |                                                                   c
!         |- exchange_core.F                                                  c
!         |  calculate valence-core contribution                              c
!                                                                             c
!     variables:                                                              c
!         kpts%nkptf   :=   number of kpoints                                      c
!         kpts%nkpt   :=   number of irreducible kpoints                          c
!         nbands  :=   number of bands for which the exchange matrix (mat_ex) c
!                      in the space of the wavefunctions is calculated        c
!         te_hfex :=   hf exchange contribution to the total energy           c
!         mnobd   :=   maximum number of occupied bands                       c
!         parent  :=   parent(ikpt) points to the symmetry equivalent point   c
!                      under the little group of kpoint nk                    c
!         symop   :=   symop(ikpt) points to the symmetry operation, which    c
!                      maps parent(ikpt) on ikpt                              c
!                                                                             c
!                                                                             c
!                                               M.Betzinger (09/07)           c
! c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c
CONTAINS

   SUBROUTINE add_vnonlocal(nk, lapw, atoms, hybrid, dimension, kpts, jsp, results, xcpot, noco, hmat)

      USE m_symm_hf, ONLY: symm_hf
      USE m_util, ONLY: intgrf, intgrf_init
      USE m_exchange_valence_hf
      USE m_exchange_core
      USE m_symmetrizeh
      USE m_wrapper
      USE m_hsefunctional, ONLY: exchange_vccvHSE, exchange_ccccHSE
      USE m_types
      USE m_io_hybrid

      IMPLICIT NONE

      TYPE(t_results), INTENT(INOUT) :: results
      CLASS(t_xcpot), INTENT(IN)    :: xcpot
      TYPE(t_dimension), INTENT(IN)    :: dimension
      TYPE(t_hybrid), INTENT(INOUT) :: hybrid
      TYPE(t_kpts), INTENT(IN)    :: kpts
      TYPE(t_lapw), INTENT(IN)    :: lapw
      TYPE(t_atoms), INTENT(IN)    :: atoms
      TYPE(t_noco), INTENT(IN)    :: noco
      TYPE(t_mat), INTENT(INOUT) :: hmat

      INTEGER, INTENT(IN)    :: jsp
      INTEGER, INTENT(IN)    :: nk

      ! local scalars
      INTEGER                 :: n, nn, iband, nbasfcn
      REAL                    :: a_ex
      TYPE(t_mat)             :: olap, tmp, v_x, z
      COMPLEX                 :: exch(dimension%neigd, dimension%neigd)

      ! initialize weighting factor for HF exchange part
      a_ex = xcpot%get_exchange_weight()

      nbasfcn = MERGE(lapw%nv(1) + lapw%nv(2) + 2*atoms%nlotot, lapw%nv(1) + atoms%nlotot, noco%l_noco)
      CALL v_x%init(hmat%l_real, nbasfcn, nbasfcn)

      CALL read_v_x(v_x, kpts%nkpt*(jsp - 1) + nk)
      ! add non-local x-potential to the hamiltonian hmat
      DO n = 1, v_x%matsize1
         DO nn = 1, n
            IF (hmat%l_real) THEN
               hmat%data_r(nn, n) = hmat%data_r(nn, n) - a_ex*v_x%data_r(nn, n)
            ELSE
               hmat%data_c(nn, n) = hmat%data_c(nn, n) - a_ex*v_x%data_c(nn, n)
            ENDIF
         END DO
      END DO
      ! calculate HF energy
      IF (hybrid%l_calhf) THEN
         WRITE (6, '(A)') new_line('n')//new_line('n')//' ###     '//'        diagonal HF exchange elements (eV)              ###'

         WRITE (6, '(A)') new_line('n')//'         k-point      '//'band          tail           pole       total(valence+core)'
      END IF

      ! read in lower triangle part of overlap matrix from direct acces file olap
      CALL olap%init(hmat%l_real, nbasfcn, nbasfcn)
      CALL read_olap(olap, kpts%nkpt*(jsp - 1) + nk)
      IF (.NOT. olap%l_real) olap%data_c = conjg(olap%data_c)

      CALL z%init(olap%l_real, nbasfcn, dimension%neigd)

      CALL read_z(z, kpts%nkpt*(jsp - 1) + nk)

      ! calculate exchange contribution of current k-point nk to total energy (te_hfex)
      ! in the case of a spin-unpolarized calculation the factor 2 is added in eigen.F90
      IF (.NOT. v_x%l_real) v_x%data_c = conjg(v_x%data_c)
      exch = 0
      z%matsize1 = MIN(z%matsize1, v_x%matsize2)

      CALL v_x%multiply(z, tmp)

      DO iband = 1, hybrid%nbands(nk)
         IF (z%l_real) THEN
            exch(iband, iband) = dot_product(z%data_r(:z%matsize1, iband), tmp%data_r(:, iband))
         ELSE
            exch(iband, iband) = dot_product(z%data_c(:z%matsize1, iband), tmp%data_c(:, iband))
         END IF
         IF (iband <= hybrid%nobd(nk)) THEN
            results%te_hfex%valence = results%te_hfex%valence - a_ex*results%w_iks(iband, nk, jsp)*exch(iband, iband)
         END IF
         IF (hybrid%l_calhf) THEN
            WRITE (6, '(      ''  ('',F5.3,'','',F5.3,'','',F5.3,'')'',I4,4X,3F15.5)') &
               kpts%bkf(:, nk), iband, (REAL(exch(iband, iband)) - hybrid%div_vv(iband, nk, jsp))*(-27.211608), &
               hybrid%div_vv(iband, nk, jsp)*(-27.211608), REAL(exch(iband, iband))*(-27.211608)
         END IF
      END DO
   END SUBROUTINE add_vnonlocal

END MODULE m_add_vnonlocal
