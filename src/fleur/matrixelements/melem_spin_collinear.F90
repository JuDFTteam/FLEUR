!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The combined 2N spin operator of a collinear jspins=2 calculation.
!>
!>  It is an operator, not orchestration, but it cannot run in the coarse pass: the two spin
!>  channels are separate eigenproblems, and the cross-spin overlap <up|dn> only becomes a
!>  spin matrix once BOTH channels have been wannierised, since it has to be rotated with
!>  both gauges. That is why it is called at the end of the wannierization instead.
MODULE m_melem_spin_collinear
   USE m_juDFT
   USE m_constants, ONLY: ImagUnit, oUnit
   USE m_types_atoms
   USE m_types_cell
   USE m_types_input
   USE m_types_kpts
   USE m_types_lapw
   USE m_types_noco
   USE m_types_nococonv
   USE m_types_sym
   USE m_types_mpi
   USE m_types_stars
   USE m_types_usdus
   USE m_types_mat
   USE m_types_radfun
   USE m_types_abc
   USE m_melem_spin, ONLY: melem_spin_mt_block
   USE m_melem_ft, ONLY: melem_ft_to_real_reduce
   USE m_melem_overlap, ONLY: melem_overlap_interstitial
   USE m_melem_get_z, ONLY: melem_get_z
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: melem_rspauli_collinear

CONTAINS

   SUBROUTINE melem_rspauli_collinear(num_bands, num_wann, min_band, max_band, atoms, input, sym, cell, noco, nococonv, kpts, &
                                      stars, usdus, radfun, eig_id, l_real_wann, distk, fmpi, v_ch)
      INTEGER, INTENT(IN) :: num_bands, num_wann     !> sizes of the Bloch and Wannier manifolds
      INTEGER, INTENT(IN) :: min_band, max_band      !> band window, forwarded to get_z
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_sym), INTENT(IN) :: sym
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_noco), INTENT(IN) :: noco
      TYPE(t_nococonv), INTENT(IN) :: nococonv
      TYPE(t_kpts), INTENT(IN) :: kpts
      TYPE(t_stars), INTENT(IN) :: stars
      TYPE(t_usdus), INTENT(IN) :: usdus
      TYPE(t_radfun), INTENT(IN) :: radfun(:)
      INTEGER, INTENT(IN) :: eig_id
      LOGICAL, INTENT(IN) :: l_real_wann
      INTEGER, INTENT(IN) :: distk(:)
      TYPE(t_mpi), INTENT(IN) :: fmpi
      !> The Wannier gauge V = u_opt.u_matrix of each spin channel, over the whole mesh:
      !> (num_bands, num_wann, nkptf, channel). Both channels have to be wannierised before
      !> the cross-spin overlap becomes a spin matrix, which is why they arrive together and
      !> already assembled rather than being built here.
      COMPLEX, INTENT(IN) :: v_ch(:, :, :, :)

      TYPE(t_lapw) :: lapw_u, lapw_d
      TYPE(t_mat)  :: zMat_u, zMat_d
      TYPE(t_abc), ALLOCATABLE :: abc_both(:, :)   ! (ntype,2): 1=up, 2=dn
      INTEGER :: nb, nw, n2, nkl, kl, gk, itype, ikpt, iu, irpt, i, j, kk, nrpts, gb(3)
      INTEGER, ALLOCATABLE :: gk_loc(:), irvec(:, :), ndegen(:)
      COMPLEX, ALLOCATABLE :: o_uu(:, :), o_dd(:, :), o_ud(:, :), o_du(:, :), Xk(:, :), tmp(:, :)
      COMPLEX, ALLOCATABLE :: sig_loc(:, :, :, :), s1(:, :, :), sr(:, :, :, :)

      nb = num_bands; nw = num_wann; n2 = 2*nw; gb = 0

      !> The gauges index the same manifolds the overlap does, and a mismatch would be a
      !> silently wrong rotation rather than a failure.
      IF (SIZE(v_ch, 1) /= nb .OR. SIZE(v_ch, 2) /= nw .OR. SIZE(v_ch, 4) /= 2) &
         CALL juDFT_error("melem_rspauli_collinear: the gauges do not match the manifold", &
                          calledby="melem_rspauli_collinear")
      nkl = COUNT(distk == fmpi%irank); ALLOCATE (gk_loc(nkl)); j = 0
      DO ikpt = 1, SIZE(distk)
         IF (distk(ikpt) == fmpi%irank) THEN; j = j + 1; gk_loc(j) = ikpt; END IF
      END DO
      ALLOCATE (abc_both(atoms%ntype, 2))
      ALLOCATE (o_uu(nb, nb), o_dd(nb, nb), o_ud(nb, nb), o_du(nb, nb), Xk(nw, nw), tmp(nb, nw))
      ALLOCATE (sig_loc(n2, n2, 3, MAX(1, nkl)), source=CMPLX(0.0, 0.0))

      DO kl = 1, nkl
         gk = gk_loc(kl)
         ! up + down eigenvectors + local coefficients at this k (same basis, different eigenvectors)
         CALL melem_get_z(min_band, max_band, eig_id, input, atoms, noco, nococonv, kpts, sym, cell, gk, 1, l_real_wann, lapw_u, zMat_u)
         DO itype = 1, atoms%ntype
            CALL abc_both(itype, 1)%init(input, atoms, nb, itype)
            CALL abc_both(itype, 1)%calc_abc(input, atoms, sym, cell, lapw_u, nb, usdus, noco, nococonv, 1, itype, zMat_u)
         END DO
         CALL melem_get_z(min_band, max_band, eig_id, input, atoms, noco, nococonv, kpts, sym, cell, gk, 2, l_real_wann, lapw_d, zMat_d)
         DO itype = 1, atoms%ntype
            CALL abc_both(itype, 2)%init(input, atoms, nb, itype)
            CALL abc_both(itype, 2)%calc_abc(input, atoms, sym, cell, lapw_d, nb, usdus, noco, nococonv, MERGE(1, 2, input%jspins == 1), itype, zMat_d)
         END DO
         ! cross-spin overlap o_ud = <up|dn>: interstitial (b=0) + muffin-tin (both channels' abc)
         o_uu = CMPLX(0.0, 0.0); o_dd = CMPLX(0.0, 0.0); o_ud = CMPLX(0.0, 0.0); o_du = CMPLX(0.0, 0.0)
         CALL melem_overlap_interstitial(stars, lapw_u, lapw_d, zMat_u, zMat_d, 0, 0, o_ud)
         CALL melem_spin_mt_block(atoms, abc_both, radfun, o_uu, o_dd, o_ud, o_du)
         ! rotate to the WF gauge: X = V_up^dagger o_ud V_dn
         tmp = MATMUL(o_ud, v_ch(:, :, gk, 2))
         Xk = MATMUL(CONJG(TRANSPOSE(v_ch(:, :, gk, 1))), tmp)
         ! assemble the 2N Pauli in the WF gauge (sigma_z block-diagonal +/- I by orthonormality)
         DO i = 1, nw
            sig_loc(i, i, 3, kl) = CMPLX(1.0, 0.0)
            sig_loc(nw + i, nw + i, 3, kl) = CMPLX(-1.0, 0.0)
         END DO
         sig_loc(1:nw, nw + 1:n2, 1, kl) = Xk
         sig_loc(nw + 1:n2, 1:nw, 1, kl) = CONJG(TRANSPOSE(Xk))
         sig_loc(1:nw, nw + 1:n2, 2, kl) = -ImagUnit*Xk
         sig_loc(nw + 1:n2, 1:nw, 2, kl) = ImagUnit*CONJG(TRANSPOSE(Xk))
      END DO
      DEALLOCATE (o_uu, o_dd, o_ud, o_du, Xk, tmp, abc_both)

      ! FT-reduce each of the 3 components (collective), rank 0 writes rspauli.1
      DO kk = 1, 3
         CALL melem_ft_to_real_reduce(cell, kpts, sig_loc(:, :, kk, :), gk_loc, fmpi%mpi_comm, s1, irvec, ndegen, nrpts)
         IF (kk == 1) ALLOCATE (sr(n2, n2, nrpts, 3))
         sr(:, :, :, kk) = s1; DEALLOCATE (s1)
      END DO
      IF (fmpi%irank == 0) THEN
         OPEN (newunit=iu, file='rspauli.1', status='replace')
         DO irpt = 1, nrpts
            DO j = 1, n2
               DO i = 1, n2
                  DO kk = 1, 3
                     WRITE (iu, '(i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,f20.8,1x,f20.8)') &
                        irvec(1, irpt), irvec(2, irpt), irvec(3, irpt), i, j, kk, REAL(sr(i, j, irpt, kk)), AIMAG(sr(i, j, irpt, kk))
                  END DO
               END DO
            END DO
         END DO
         CLOSE (iu)
         WRITE (oUnit, '(a,i0,a)') 'wannierlib: wrote rspauli.1 (combined 2N collinear spin, ', nrpts, ' R-vectors, distributed FT)'
      END IF
      DEALLOCATE (sig_loc, gk_loc)
      IF (ALLOCATED(sr)) DEALLOCATE (sr)
      IF (ALLOCATED(irvec)) DEALLOCATE (irvec, ndegen)
   END SUBROUTINE melem_rspauli_collinear

END MODULE m_melem_spin_collinear
