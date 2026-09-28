!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The 2N operators of a collinear jspins=2 calculation, in real space: the ones that
!>  cannot be assembled until BOTH spin channels have been wannierised, because they need
!>  both gauges at once. Everything separable per channel is written by wgauge_operators_r.
!>
!>  For the spin, the ingredient is the cross-spin overlap <up|dn> on the coarse mesh, a
!>  Bloch quantity produced with the other coarse matrices; what is left here is rotating
!>  it with both gauges and assembling the 2N Pauli matrices. sigma_z is +/-1 on the
!>  diagonal by the orthonormality of each channel, so it is written down rather than
!>  computed, and the transverse components follow from the single rotated block.
!>
!>  For the orbital moment there is nothing to rotate across channels: L acts on the
!>  spatial part alone, so <up|L|dn> carries the spin overlap <up|dn> = 0 as a factor and
!>  the cross block vanishes identically. The 2N matrix is block-diagonal, and the two
!>  blocks are each channel's own L in its own gauge.
MODULE m_wgauge_spin_collinear
   USE m_juDFT
   USE m_constants, ONLY: ImagUnit, oUnit
   USE m_types_cell
   USE m_types_kpts
   USE m_types_mpi
   USE m_wgauge_ft, ONLY: wgauge_ft_to_real_reduce
   USE m_wgauge_io, ONLY: wgauge_write_realspace
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: wgauge_rspauli_collinear, wgauge_anglmom_collinear, wgauge_soc_collinear

CONTAINS

   SUBROUTINE wgauge_rspauli_collinear(num_wann, x0, v_ch, cell, kpts, distk, fmpi)
      INTEGER, INTENT(IN) :: num_wann
      !> The cross-spin overlap <up|dn> on this rank's k-slice, in the Bloch basis:
      !> (num_bands, num_bands, nk_loc), in ascending global-k order.
      COMPLEX, INTENT(IN) :: x0(:, :, :)
      !> The Wannier gauge V = u_opt.u_matrix of each spin channel over the whole mesh:
      !> (num_bands, num_wann, nkptf, channel). Both channels are needed at once, which is
      !> why this step waits until neither wannierization is still running.
      COMPLEX, INTENT(IN) :: v_ch(:, :, :, :)
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_kpts), INTENT(IN) :: kpts
      INTEGER, INTENT(IN) :: distk(:)
      TYPE(t_mpi), INTENT(IN) :: fmpi

      INTEGER :: nb, nw, n2, nkl, kl, gk, iu, irpt, i, j, kk, nrpts
      INTEGER, ALLOCATABLE :: gk_loc(:), irvec(:, :), ndegen(:)
      COMPLEX, ALLOCATABLE :: Xk(:, :), tmp(:, :)
      COMPLEX, ALLOCATABLE :: sig_loc(:, :, :, :), s1(:, :, :), sr(:, :, :, :)

      nb = SIZE(x0, 1); nw = num_wann; n2 = 2*nw

      !> The gauges index the same manifolds the overlap does, and a mismatch would be a
      !> silently wrong rotation rather than a failure.
      IF (SIZE(v_ch, 1) /= nb .OR. SIZE(v_ch, 2) /= nw .OR. SIZE(v_ch, 4) /= 2) &
         CALL juDFT_error("wgauge_rspauli_collinear: the gauges do not match the manifold", &
                          calledby="wgauge_rspauli_collinear")

      nkl = COUNT(distk == fmpi%irank); ALLOCATE (gk_loc(nkl)); j = 0
      DO i = 1, SIZE(distk)
         IF (distk(i) == fmpi%irank) THEN; j = j + 1; gk_loc(j) = i; END IF
      END DO
      IF (SIZE(x0, 3) < nkl) CALL juDFT_error( &
         "wgauge_rspauli_collinear: fewer overlap slices than k-points on this rank", &
         calledby="wgauge_rspauli_collinear")

      ALLOCATE (Xk(nw, nw), tmp(nb, nw))
      ALLOCATE (sig_loc(n2, n2, 3, MAX(1, nkl)), source=CMPLX(0.0, 0.0))

      DO kl = 1, nkl
         gk = gk_loc(kl)
         ! rotate to the WF gauge: X = V_up^dagger o_ud V_dn
         tmp = MATMUL(x0(:, :, kl), v_ch(:, :, gk, 2))
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
      DEALLOCATE (Xk, tmp)

      ! FT-reduce each of the 3 components (collective), rank 0 writes rspauli.1
      DO kk = 1, 3
         CALL wgauge_ft_to_real_reduce(cell, kpts, sig_loc(:, :, kk, :), gk_loc, fmpi%mpi_comm, s1, irvec, ndegen, nrpts)
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
   END SUBROUTINE wgauge_rspauli_collinear

   !> The site-summed orbital moment of both channels as one 2N block-diagonal matrix,
   !> written as anglmomrs.1 -- the same single-file contract the 2N spin uses, so that a
   !> reader of a collinear run does not have to know how many channels produced it.
   SUBROUTINE wgauge_anglmom_collinear(num_wann, l0, v_ch, cell, kpts, distk, fmpi)
      INTEGER, INTENT(IN) :: num_wann
      !> Site-resolved L on this rank's k-slice, in the Bloch basis, for both channels:
      !> (num_bands, num_bands, 3, natoms, 2, nk_loc). The sum over sites is taken here.
      COMPLEX, INTENT(IN) :: l0(:, :, :, :, :, :)
      !> The Wannier gauge V = u_opt.u_matrix of each channel: (num_bands, num_wann, nkptf, 2).
      COMPLEX, INTENT(IN) :: v_ch(:, :, :, :)
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_kpts), INTENT(IN) :: kpts
      INTEGER, INTENT(IN) :: distk(:)
      TYPE(t_mpi), INTENT(IN) :: fmpi

      INTEGER :: nb, nw, n2, nkl, kl, gk, iu, irpt, i, j, kk, ic, nrpts
      INTEGER, ALLOCATABLE :: gk_loc(:), irvec(:, :), ndegen(:)
      COMPLEX, ALLOCATABLE :: lb(:, :), tmp(:, :)
      COMPLEX, ALLOCATABLE :: lop_loc(:, :, :, :), l1(:, :, :), lr(:, :, :, :)

      nb = SIZE(l0, 1); nw = num_wann; n2 = 2*nw

      IF (SIZE(v_ch, 1) /= nb .OR. SIZE(v_ch, 2) /= nw .OR. SIZE(v_ch, 4) /= 2) &
         CALL juDFT_error("wgauge_anglmom_collinear: the gauges do not match the manifold", &
                          calledby="wgauge_anglmom_collinear")
      IF (SIZE(l0, 5) /= 2) CALL juDFT_error( &
         "wgauge_anglmom_collinear: L was not stored for both spin channels", &
         calledby="wgauge_anglmom_collinear")

      nkl = COUNT(distk == fmpi%irank); ALLOCATE (gk_loc(nkl)); j = 0
      DO i = 1, SIZE(distk)
         IF (distk(i) == fmpi%irank) THEN; j = j + 1; gk_loc(j) = i; END IF
      END DO
      IF (SIZE(l0, 6) < nkl) CALL juDFT_error( &
         "wgauge_anglmom_collinear: fewer L slices than k-points on this rank", &
         calledby="wgauge_anglmom_collinear")

      ALLOCATE (lb(nb, nb), tmp(nb, nw))
      ALLOCATE (lop_loc(n2, n2, 3, MAX(1, nkl)), source=CMPLX(0.0, 0.0))

      DO kl = 1, nkl
         gk = gk_loc(kl)
         DO kk = 1, 3
            DO ic = 1, 2
               lb = SUM(l0(:, :, kk, :, ic, kl), DIM=3)          ! sum over sites
               tmp = MATMUL(lb, v_ch(:, :, gk, ic))
               i = (ic - 1)*nw
               lop_loc(i + 1:i + nw, i + 1:i + nw, kk, kl) = &
                  MATMUL(CONJG(TRANSPOSE(v_ch(:, :, gk, ic))), tmp)
            END DO
         END DO
      END DO
      DEALLOCATE (lb, tmp)

      DO kk = 1, 3
         CALL wgauge_ft_to_real_reduce(cell, kpts, lop_loc(:, :, kk, :), gk_loc, fmpi%mpi_comm, l1, irvec, ndegen, nrpts)
         IF (kk == 1) ALLOCATE (lr(n2, n2, nrpts, 3))
         lr(:, :, :, kk) = l1; DEALLOCATE (l1)
      END DO
      IF (fmpi%irank == 0) THEN
         OPEN (newunit=iu, file='anglmomrs.1', status='replace')
         DO irpt = 1, nrpts
            DO j = 1, n2
               DO i = 1, n2
                  DO kk = 1, 3
                     WRITE (iu, '(i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,f20.8,1x,f20.8)') &
                        irvec(1, irpt), irvec(2, irpt), irvec(3, irpt), i, j, kk, REAL(lr(i, j, irpt, kk)), AIMAG(lr(i, j, irpt, kk))
                  END DO
               END DO
            END DO
         END DO
         CLOSE (iu)
         WRITE (oUnit, '(a,i0,a)') 'wannierlib: wrote anglmomrs.1 (combined 2N collinear orbital, ', nrpts, ' R-vectors, distributed FT)'
      END IF
      DEALLOCATE (lop_loc, gk_loc)
      IF (ALLOCATED(lr)) DEALLOCATE (lr)
      IF (ALLOCATED(irvec)) DEALLOCATE (irvec, ndegen)
   END SUBROUTINE wgauge_anglmom_collinear

   !> Spin-orbit coupling as an OPERATOR in the 2N collinear Wannier basis, written as
   !> rssocmat.1. This is the piece that makes a collinear wannierisation usable for SOC
   !> physics: the coupling never enters the eigenproblem, so the Wannier functions stay
   !> spin-pure and every operator stays localised, and whoever consumes the files adds
   !> H_SOC to H(R) themselves.
   !>
   !> Layout. Rows and columns already carry the spin (1..N is up, N+1..2N is down), so the
   !> 2x2 index the file also has is redundant: an entry is non-zero only in the component
   !> matching its own quadrant, which is how the reference files are shaped too. The
   !> consumer sums the four components per (m,n) without looking at the indices, so what
   !> the sum has to come out to is the single block value -- and it does, the other three
   !> being zero.
   SUBROUTINE wgauge_soc_collinear(num_wann, soc4, v_ch, cell, kpts, distk, fmpi)
      INTEGER, INTENT(IN) :: num_wann
      !> The four spin blocks of H_SOC on this rank's k-slice, in the Bloch basis:
      !> (num_bands, num_bands, 4, nk_loc), ordered (1,1) (1,2) (2,1) (2,2).
      COMPLEX, INTENT(IN) :: soc4(:, :, :, :)
      COMPLEX, INTENT(IN) :: v_ch(:, :, :, :)
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_kpts), INTENT(IN) :: kpts
      INTEGER, INTENT(IN) :: distk(:)
      TYPE(t_mpi), INTENT(IN) :: fmpi

      INTEGER :: nb, nw, n2, nkl, kl, gk, iu, irpt, i, j, kk, is, js, ib, jb, ic, nrpts
      INTEGER, ALLOCATABLE :: gk_loc(:), irvec(:, :), ndegen(:)
      COMPLEX, ALLOCATABLE :: tmp(:, :), blk(:, :)
      COMPLEX, ALLOCATABLE :: soc_loc(:, :, :, :), s1(:, :, :), sr(:, :, :, :)

      nb = SIZE(soc4, 1); nw = num_wann; n2 = 2*nw

      IF (SIZE(v_ch, 1) /= nb .OR. SIZE(v_ch, 2) /= nw .OR. SIZE(v_ch, 4) /= 2) &
         CALL juDFT_error("wgauge_soc_collinear: the gauges do not match the manifold", &
                          calledby="wgauge_soc_collinear")
      IF (SIZE(soc4, 3) /= 4) CALL juDFT_error( &
         "wgauge_soc_collinear: the SOC operator did not keep its four spin blocks", &
         calledby="wgauge_soc_collinear")

      nkl = COUNT(distk == fmpi%irank); ALLOCATE (gk_loc(nkl)); j = 0
      DO i = 1, SIZE(distk)
         IF (distk(i) == fmpi%irank) THEN; j = j + 1; gk_loc(j) = i; END IF
      END DO
      IF (SIZE(soc4, 4) < nkl) CALL juDFT_error( &
         "wgauge_soc_collinear: fewer SOC slices than k-points on this rank", &
         calledby="wgauge_soc_collinear")

      ALLOCATE (tmp(nb, nw), blk(nw, nw))
      ALLOCATE (soc_loc(n2, n2, 4, MAX(1, nkl)), source=CMPLX(0.0, 0.0))

      DO kl = 1, nkl
         gk = gk_loc(kl)
         DO is = 1, 2                    ! spin of the bra, i.e. of the row block
            DO js = 1, 2                 ! spin of the ket
               tmp = MATMUL(soc4(:, :, (is - 1)*2 + js, kl), v_ch(:, :, gk, js))
               blk = MATMUL(CONJG(TRANSPOSE(v_ch(:, :, gk, is))), tmp)
               ib = (is - 1)*nw; jb = (js - 1)*nw
               !> The writer builds its component index as (ii-1)*2 + jj and prints the pair
               !> as (jj, ii), so this is the component that comes out labelled (is, js).
               ic = (js - 1)*2 + is
               soc_loc(ib + 1:ib + nw, jb + 1:jb + nw, ic, kl) = blk
            END DO
         END DO
      END DO
      DEALLOCATE (tmp, blk)

      DO kk = 1, 4
         CALL wgauge_ft_to_real_reduce(cell, kpts, soc_loc(:, :, kk, :), gk_loc, fmpi%mpi_comm, s1, irvec, ndegen, nrpts)
         IF (kk == 1) ALLOCATE (sr(n2, n2, nrpts, 4))
         sr(:, :, :, kk) = s1; DEALLOCATE (s1)
      END DO
      IF (fmpi%irank == 0) THEN
         CALL wgauge_write_realspace(sr, irvec, ndegen, nrpts, n2, 4, 'spinor2x2', 'rssocmat.1', fmpi%irank)
         WRITE (oUnit, '(a,i0,a)') 'wannierlib: wrote rssocmat.1 (2N collinear spin-orbit, ', nrpts, ' R-vectors, distributed FT)'
      END IF
      DEALLOCATE (soc_loc, gk_loc)
      IF (ALLOCATED(sr)) DEALLOCATE (sr)
      IF (ALLOCATED(irvec)) DEALLOCATE (irvec, ndegen)
   END SUBROUTINE wgauge_soc_collinear

END MODULE m_wgauge_spin_collinear
