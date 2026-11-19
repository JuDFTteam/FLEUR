!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hs_int_direct
#ifdef _OPENACC
   USE openacc
#endif
   USE m_types_mat
   USE m_types_mpi
   USE m_types_stars
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: hs_int_direct
CONTAINS
   SUBROUTINE hs_int_direct(fmpi, stars, bbmat, gvecPr, gvec, kvecPr, kvec, nvPr, nv, &
                          & iTkin, fact, l_smat, l_fullj, vpw, hmat, smat, theta_alt, vtau_pw)
      ! Calculates matrix elements of the form
      ! <\phi_{k'G'}|M|\phi_{kG}>
      ! for different use cases in the DFT/DFPT scf loop and operators M.
      !
      ! The spin loop and distinction is found on a higher level and translates
      ! into new switches.
      !
      ! DFT:
      ! M = \Theta_{IR} * (T + V) goes into hmat, M = \Theta_{IR} into smat.
      ! [vpw = \Theta_{IR} * V]
      ! [stars%ustep = \Theta_{IR}]
      ! [l_smat = F for offdiags, l_fullj = F]
      ! [iTkin = 0 for offdiags, 1 for l_useapw, 2 else]
      !
      ! DFPT:
      ! M = \Theta_{IR}^{(1)} * (T + V) + \Theta_{IR} * V^{(1)} goes into hmat,
      ! M = \Theta_{IR}^{(1)} into smat.
      ! [vpw = \Theta_{IR}^{(1)} * V + \Theta_{IR} * V^{(1)}]
      ! [stars%ustep = \Theta_{IR}^{(1)}]
      ! [l_smat = F for offdiags, l_fullj = T]
      ! [iTkin = 0 for offdiags, 1 else]


      ! If hmat and smat are mapped to the device, the elements are computed there.

      IMPLICIT NONE

      TYPE(t_mpi),   INTENT(IN)    :: fmpi
      TYPE(t_stars), INTENT(IN)    :: stars
      REAL,          INTENT(IN)    :: bbmat(3, 3)
      INTEGER,       INTENT(IN)    :: gvecPr(:, :), gvec(:, :)
      REAL,          INTENT(IN)    :: kvecPr(3), kvec(3)
      INTEGER,       INTENT(IN)    :: nvPr, nv, iTkin, fact
      LOGICAL,       INTENT(IN)    :: l_smat, l_fullj
      COMPLEX,       INTENT(IN)    :: vpw(:)

      CLASS(t_mat),  INTENT(INOUT) :: hmat, smat

      COMPLEX, OPTIONAL, INTENT(IN) :: theta_alt(:)
      ! Optional V_tau star coefficients for MetaGGA interstitial contribution
      COMPLEX, OPTIONAL, INTENT(IN) :: vtau_pw(:)

      INTEGER :: ikGPr, ikG, ikG0, gShift(3), mx(3), gp(3), gk(3), ns(3)
      COMPLEX :: th, ts
      LOGICAL :: l_vtau, l_hit, l_gpu
      COMPLEX, ALLOCATABLE :: theta(:), vtau(:)

      gShift = nint(kvecPr - kvec - stars%center)
      mx = [stars%mx1, stars%mx2, stars%mx3]
      IF (PRESENT(theta_alt)) THEN
         theta = theta_alt
      ELSE
         theta = stars%ustep
      END IF
      l_vtau = PRESENT(vtau_pw)
      IF (l_vtau) THEN
         vtau = vtau_pw
      ELSE
         ALLOCATE (vtau(1))
         vtau = 0.0
      END IF
      ns = [SIZE(vpw), SIZE(theta), SIZE(vtau)]

      l_gpu = .FALSE.
#ifdef _OPENACC
      IF (hmat%l_real) THEN
         l_gpu = acc_is_present(hmat%data_r) .AND. acc_is_present(smat%data_r)
      ELSE
         l_gpu = acc_is_present(hmat%data_c) .AND. acc_is_present(smat%data_c)
      END IF
      IF (l_gpu) THEN
         IF (hmat%l_real) THEN
            CALL priv_gpu_r(hmat%data_r, smat%data_r, fmpi%n_rank, fmpi%n_size, nvPr, nv, l_fullj, gvecPr, gvec, &
                            kvecPr, kvec, gShift, fact, mx, stars%ig, stars%rgphs, ns, vpw, theta, vtau, bbmat, iTkin, l_smat, l_vtau)
         ELSE
            CALL priv_gpu_c(hmat%data_c, smat%data_c, fmpi%n_rank, fmpi%n_size, nvPr, nv, l_fullj, gvecPr, gvec, &
                            kvecPr, kvec, gShift, fact, mx, stars%ig, stars%rgphs, ns, vpw, theta, vtau, bbmat, iTkin, l_smat, l_vtau)
         END IF
         RETURN
      END IF
#endif

      !$OMP PARALLEL DO SCHEDULE(dynamic) DEFAULT(none) &
      !$OMP SHARED(fmpi, stars, bbmat, gvecPr, gvec, kvecPr, kvec, gShift, mx, ns) &
      !$OMP SHARED(nvPr, nv, iTkin, fact, l_smat, l_fullj, l_vtau, vpw, theta, vtau, hmat, smat) &
      !$OMP PRIVATE(ikGPr, ikG, ikG0, gp, gk, th, ts, l_hit)
      DO ikG = fmpi%n_rank + 1, nv, fmpi%n_size
         ikG0 = (ikG-1) / fmpi%n_size + 1
         DO  ikGPr = 1, MERGE(nvPr, MIN(ikG, nvPr), l_fullj)
            gp = gvecPr(:, ikGPr)
            gk = gvec(:, ikG)
            CALL priv_element(gp, gk, kvecPr, kvec, gShift, fact, mx, stars%ig, stars%rgphs, &
                              ns, vpw, theta, vtau, bbmat, iTkin, l_smat, l_vtau, th, ts, l_hit)
            IF (.NOT.l_hit) CYCLE
            IF (hmat%l_real) THEN
               hmat%data_r(ikGPr, ikG0) = REAL(th)
               smat%data_r(ikGPr, ikG0) = REAL(ts)
            ELSE
               hmat%data_c(ikGPr, ikG0) = th
               smat%data_c(ikGPr, ikG0) = ts
            END IF
         END DO
      END DO
      !$OMP END PARALLEL DO
   END SUBROUTINE hs_int_direct

   PURE SUBROUTINE priv_element(gPr, g, kvecPr, kvec, gShift, fact, mx, ig, rgphs, ns, vpw, theta, vtau, bbmat, &
                                iTkin, l_smat, l_vtau, th, ts, l_hit)
      ! One matrix element <k'+G'|M|k+G>; l_hit is false if G'-G is no star member.
      !$acc routine seq
      INTEGER, INTENT(IN)  :: gPr(3), g(3), gShift(3), fact, mx(3), iTkin, ns(3)
      INTEGER, INTENT(IN)  :: ig(-mx(1):mx(1), -mx(2):mx(2), -mx(3):mx(3))
      COMPLEX, INTENT(IN)  :: rgphs(-mx(1):mx(1), -mx(2):mx(2), -mx(3):mx(3))
      COMPLEX, INTENT(IN)  :: vpw(ns(1)), theta(ns(2)), vtau(ns(3))
      REAL,    INTENT(IN)  :: kvecPr(3), kvec(3), bbmat(3, 3)
      LOGICAL, INTENT(IN)  :: l_smat, l_vtau
      COMPLEX, INTENT(OUT) :: th, ts
      LOGICAL, INTENT(OUT) :: l_hit

      INTEGER :: gPrG(3), gInd
      COMPLEX :: phase
      REAL    :: bvecPr(3), bvec(3), r2

      th = 0.0
      ts = 0.0
      gPrG = fact * (gPr - g + gShift)
      l_hit = .NOT.(ABS(gPrG(1)) > mx(1) .OR. ABS(gPrG(2)) > mx(2) .OR. ABS(gPrG(3)) > mx(3))
      IF (.NOT.l_hit) RETURN
      gInd = ig(gPrG(1), gPrG(2), gPrG(3))
      l_hit = gInd /= 0
      IF (.NOT.l_hit) RETURN

      phase = rgphs(gPrG(1), gPrG(2), gPrG(3))
      th = phase * vpw(gInd)

      IF (iTkin > 0) THEN
         bvecPr = kvecPr + gPr
         bvec = kvec + g
         IF (iTkin == 1) THEN ! Symmetric Dirac form
            r2 = 0.5 * priv_bdot(bvecPr, bvec, bbmat)
         ELSE IF (iTkin == 2) THEN ! Symmetrized Laplace form
            r2 =      0.25 * priv_bdot(bvecPr, bvecPr, bbmat)
            r2 = r2 + 0.25 * priv_bdot(bvec, bvec, bbmat)
         ELSE ! Pure Laplace form
            r2 = 0.5 * priv_bdot(bvec, bvec, bbmat)
         END IF
         th = th + phase * r2 * theta(gInd)

         ! MetaGGA V_tau interstitial contribution:
         ! H_tau(G',G) = (1/2) V_tau(G'-G) * (k+G') . bbmat . (k+G)
         ! hs_int only passes vtau_pw for iSpinPr == iSpin, which implies iTkin > 0.
         IF (l_vtau) th = th + phase * 0.5 * priv_bdot(bvecPr, bvec, bbmat) * vtau(gInd)
      END IF

      IF (l_smat) ts = phase * theta(gInd)
   END SUBROUTINE priv_element

   PURE REAL FUNCTION priv_bdot(a, b, bbmat)
      ! a . bbmat . b
      !$acc routine seq
      REAL, INTENT(IN) :: a(3), b(3), bbmat(3, 3)
      INTEGER :: i, j
      priv_bdot = 0.0
      DO j = 1, 3
         DO i = 1, 3
            priv_bdot = priv_bdot + a(i)*bbmat(i, j)*b(j)
         END DO
      END DO
   END FUNCTION priv_bdot

#ifdef _OPENACC
   SUBROUTINE priv_gpu_r(h, s, n_rank, n_size, nvPr, nv, l_fullj, gvecPr, gvec, kvecPr, kvec, gShift, fact, mx, &
                         ig, rgphs, ns, vpw, theta, vtau, bbmat, iTkin, l_smat, l_vtau)
      REAL,    INTENT(INOUT) :: h(:, :), s(:, :)
      INTEGER, INTENT(IN)    :: n_rank, n_size, nvPr, nv, gShift(3), fact, mx(3), iTkin, ns(3)
      LOGICAL, INTENT(IN)    :: l_fullj, l_smat, l_vtau
      INTEGER, INTENT(IN)    :: gvecPr(:, :), gvec(:, :)
      INTEGER, INTENT(IN)    :: ig(-mx(1):mx(1), -mx(2):mx(2), -mx(3):mx(3))
      COMPLEX, INTENT(IN)    :: rgphs(-mx(1):mx(1), -mx(2):mx(2), -mx(3):mx(3))
      COMPLEX, INTENT(IN)    :: vpw(ns(1)), theta(ns(2)), vtau(ns(3))
      REAL,    INTENT(IN)    :: kvecPr(3), kvec(3), bbmat(3, 3)

      INTEGER :: ikG, ikGPr, ikG0, gp(3), gk(3)
      COMPLEX :: th, ts
      LOGICAL :: l_hit

      !$acc parallel loop gang present(h, s) copyin(gvecPr, gvec, ig, rgphs, vpw, theta, vtau, kvecPr, kvec, bbmat, gShift, mx, ns) &
      !$acc private(ikG0)
      DO ikG = n_rank + 1, nv, n_size
         ikG0 = (ikG-1) / n_size + 1
         !$acc loop vector private(gp, gk, th, ts, l_hit)
         DO ikGPr = 1, MERGE(nvPr, MIN(ikG, nvPr), l_fullj)
            gp = gvecPr(:, ikGPr)
            gk = gvec(:, ikG)
            CALL priv_element(gp, gk, kvecPr, kvec, gShift, fact, mx, ig, rgphs, &
                              ns, vpw, theta, vtau, bbmat, iTkin, l_smat, l_vtau, th, ts, l_hit)
            IF (l_hit) THEN
               h(ikGPr, ikG0) = REAL(th)
               s(ikGPr, ikG0) = REAL(ts)
            END IF
         END DO
      END DO
   END SUBROUTINE priv_gpu_r

   SUBROUTINE priv_gpu_c(h, s, n_rank, n_size, nvPr, nv, l_fullj, gvecPr, gvec, kvecPr, kvec, gShift, fact, mx, &
                         ig, rgphs, ns, vpw, theta, vtau, bbmat, iTkin, l_smat, l_vtau)
      COMPLEX, INTENT(INOUT) :: h(:, :), s(:, :)
      INTEGER, INTENT(IN)    :: n_rank, n_size, nvPr, nv, gShift(3), fact, mx(3), iTkin, ns(3)
      LOGICAL, INTENT(IN)    :: l_fullj, l_smat, l_vtau
      INTEGER, INTENT(IN)    :: gvecPr(:, :), gvec(:, :)
      INTEGER, INTENT(IN)    :: ig(-mx(1):mx(1), -mx(2):mx(2), -mx(3):mx(3))
      COMPLEX, INTENT(IN)    :: rgphs(-mx(1):mx(1), -mx(2):mx(2), -mx(3):mx(3))
      COMPLEX, INTENT(IN)    :: vpw(ns(1)), theta(ns(2)), vtau(ns(3))
      REAL,    INTENT(IN)    :: kvecPr(3), kvec(3), bbmat(3, 3)

      INTEGER :: ikG, ikGPr, ikG0, gp(3), gk(3)
      COMPLEX :: th, ts
      LOGICAL :: l_hit

      !$acc parallel loop gang present(h, s) copyin(gvecPr, gvec, ig, rgphs, vpw, theta, vtau, kvecPr, kvec, bbmat, gShift, mx, ns) &
      !$acc private(ikG0)
      DO ikG = n_rank + 1, nv, n_size
         ikG0 = (ikG-1) / n_size + 1
         !$acc loop vector private(gp, gk, th, ts, l_hit)
         DO ikGPr = 1, MERGE(nvPr, MIN(ikG, nvPr), l_fullj)
            gp = gvecPr(:, ikGPr)
            gk = gvec(:, ikG)
            CALL priv_element(gp, gk, kvecPr, kvec, gShift, fact, mx, ig, rgphs, &
                              ns, vpw, theta, vtau, bbmat, iTkin, l_smat, l_vtau, th, ts, l_hit)
            IF (l_hit) THEN
               h(ikGPr, ikG0) = th
               s(ikGPr, ikG0) = ts
            END IF
         END DO
      END DO
   END SUBROUTINE priv_gpu_c
#endif
END MODULE m_hs_int_direct
