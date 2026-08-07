!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Real-space export of the Wannier matrix elements O(R) -- "Fourier step 3".
!>
!>  Given a Bloch-basis coarse operator slice O^0(k) on this rank's k-points and the
!>  Wannier gauge V(k) = u_opt(k).u_matrix(k), each operator is rotated into the Wannier
!>  gauge and Fourier transformed to real space with the distributed reduce in m_melem_ft;
!>  rank 0 writes the plain-text file. No band interpolation happens here.
!>
!>  Files written (legacy formats, consumed by external transport post-processing):
!>    WF<n>_hr.dat   H(R) in eV, Wannier90 hr format
!>    WF<n>_r.dat    A(R) = <0n|r|Rm> in Angstrom, Wannier90 r format
!>    rspauli.1      spin (3 components)
!>    anglmomrs.<n>  orbital moment (3 components)
!>    rssocmat.1     spin-orbit (2x2 spinor blocks)
!>    wig_vectors    the Wigner-Seitz R-mesh
!>
!>  Extracted from m_wannierlib_w90_adapter. The only Wannier90 coupling that came with
!>  it -- the b-shell weights and the reference Wannier centres -- now arrives as a plain
!>  t_melem_bmesh, so this module is free of the Wannier90 library.
MODULE m_melem_operators_r
  USE m_juDFT
  USE m_constants, ONLY : oUnit, hartree_to_ev_const
  USE m_types_cell
  USE m_types_kpts
  USE m_types_melem_manifold, ONLY : t_melem_manifold
  USE m_types_melem_request, ONLY : t_melem_request
  USE m_types_melem_bmesh
  USE m_melem_ft, ONLY : melem_ft_to_real, melem_ws_vectors, melem_ft_to_real_reduce
  USE m_melem_io, ONLY : melem_write_realspace
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: melem_write_operators_r, melem_op_rs_distributed, melem_build_berry_aw_r
  PUBLIC :: melem_check_berry_centres, melem_build_hamk, melem_write_hr, melem_write_ar
  PUBLIC :: melem_write_wig_once
CONTAINS

  ! Build the Wannier-gauge Berry connection in real space A^(W)_alpha(R), distributed:
  ! each rank forms A^(W)_alpha(k) = i sum_b w_b b_alpha (M^(W,b)(k) - delta) on its OWN k-slice
  ! (M^(W,b)(k) = V(k)^dagger M^(0,b)(k) V(k_b), V = u_opt.u_matrix, k_b = nnlist(b,k)) from the
  ! local overlaps mmn_loc (global k = gk_loc), then reduces coarse -> R with the distributed
  ! FT-reduce (collective). The full-mesh overlaps are never gathered. A^(W)(R) is exactly the
  ! position operator A(R) = <0n|r|Rm> and the R->k' interpolant used by the velocity/Berry curvature.
  ! kmesh (wb/bk/nnlist) arrives as a plain t_melem_bmesh -> no Wannier90 coupling here.
  ! Collective over mpi_comm; result valid on all ranks.
  SUBROUTINE melem_build_berry_aw_r(this, cell, kpts, mmn_loc, gk_loc, u_opt, u_matrix, bmesh, mpi_comm, &
                                         aw_r, irvec, ndegen, nrpts)
    TYPE(t_melem_manifold), INTENT(IN) :: this
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    COMPLEX, INTENT(IN) :: mmn_loc(:, :, :, :)    ! (nb,nb,nntot,nk_loc) this rank's overlap slice
    INTEGER, INTENT(IN) :: gk_loc(:)              ! (nk_loc) global k index of each slice entry
    COMPLEX, INTENT(IN) :: u_opt(:, :, :)         ! (nb,nw,nk) full mesh
    COMPLEX, INTENT(IN) :: u_matrix(:, :, :)      ! (nw,nw,nk) full mesh
    TYPE(t_melem_bmesh), INTENT(IN) :: bmesh      ! b-shell weights / neighbour list of the coarse mesh
    INTEGER, INTENT(IN) :: mpi_comm
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: aw_r(:, :, :, :)   ! (nw,nw,nrpts,3) reduced Berry connection / A(R)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: irvec(:, :), ndegen(:)
    INTEGER, INTENT(OUT) :: nrpts
    INTEGER :: nb, nw, nk, nnt, k, kb, nn, a, i, kl, nkl
    REAL :: wb, b(3)
    COMPLEX, ALLOCATABLE :: Vk(:, :), Vkb(:, :), Mw(:, :), tmp(:, :), aw_loc(:, :, :, :), a1(:, :, :)

    IF (bmesh%nntot < 1 .OR. .NOT. ALLOCATED(bmesh%nnlist)) THEN
      ! no b-mesh available (e.g. built without the Wannier90 library) -> nothing to build
      nrpts = 0
      ALLOCATE(aw_r(1, 1, 1, 3), irvec(3, 1), ndegen(1))
      RETURN
    END IF
    nb = this%num_bands; nw = this%num_wann; nk = kpts%nkptf
    nnt = bmesh%nntot
    nkl = SIZE(gk_loc)
    ALLOCATE(aw_loc(nw, nw, 3, MAX(1, nkl)), source=CMPLX(0.0, 0.0))
    ALLOCATE(Vk(nb, nw), Vkb(nb, nw), Mw(nw, nw), tmp(nb, nw))
    DO kl = 1, nkl
      k = gk_loc(kl)
      Vk = MATMUL(u_opt(:, :, k), u_matrix(:, :, k))
      DO nn = 1, nnt
        kb = bmesh%nnlist(k, nn)                       ! shape: nnlist(num_kpts, nntot)
        IF (kb < 1 .OR. kb > nk) CALL juDFT_error('wannierlib: nnlist neighbour index out of range', &
                                                  calledby='melem_build_berry_aw_r')
        wb = bmesh%wb(nn)
        b  = bmesh%bk(:, nn, k)
        Vkb = MATMUL(u_opt(:, :, kb), u_matrix(:, :, kb))
        tmp = MATMUL(mmn_loc(:, :, nn, kl), Vkb)            ! (nb x nw)
        Mw  = MATMUL(CONJG(TRANSPOSE(Vk)), tmp)             ! (nw x nw) = M^(W,b)(k)
        DO i = 1, nw
          Mw(i, i) = Mw(i, i) - CMPLX(1.0, 0.0)            ! subtract delta on the diagonal
        END DO
        DO a = 1, 3
          aw_loc(:, :, a, kl) = aw_loc(:, :, a, kl) + CMPLX(0.0, wb*b(a)) * Mw
        END DO
      END DO
    END DO
    DEALLOCATE(Vk, Vkb, Mw, tmp)
    DO a = 1, 3
      CALL melem_ft_to_real_reduce(cell, kpts, aw_loc(:, :, a, :), gk_loc, mpi_comm, a1, irvec, ndegen, nrpts)
      IF (a == 1) ALLOCATE(aw_r(nw, nw, nrpts, 3))
      aw_r(:, :, :, a) = a1; DEALLOCATE(a1)
    END DO
    DEALLOCATE(aw_loc)
  END SUBROUTINE melem_build_berry_aw_r

  ! Validation: the diagonal of A^(W)_alpha at R=0 is the Wannier centre <r_alpha>_n
  ! (Marzari-Vanderbilt). Since A^(W)_alpha(R=0) = (1/Nk) sum_k A^(W)_alpha(k), it is exactly
  ! the R=0 entry of the reduced aw_r. Compare to the reference centres carried by the
  ! b-mesh bundle (filled from w90_get_centres on the wannierization side) to calibrate the
  ! conj/sign convention of the overlaps. Writes berry_centre_check.dat (rank 0).
  ! No reference centres available -> the check is silently skipped.
  SUBROUTINE melem_check_berry_centres(this, aw_r, irvec, nrpts, bmesh)
    TYPE(t_melem_manifold), INTENT(IN) :: this
    COMPLEX, INTENT(IN) :: aw_r(:, :, :, :)     ! (nw,nw,nrpts,3)
    INTEGER, INTENT(IN) :: irvec(:, :), nrpts
    TYPE(t_melem_bmesh), INTENT(IN) :: bmesh
    INTEGER :: nw, a, n, iu, irpt, irpt0
    COMPLEX :: aR0

    IF (.NOT. ALLOCATED(bmesh%centres)) RETURN
    nw = this%num_wann
    irpt0 = 0
    DO irpt = 1, nrpts
      IF (ALL(irvec(:, irpt) == 0)) THEN; irpt0 = irpt; EXIT; END IF
    END DO
    IF (irpt0 == 0) RETURN   ! no R=0 vector (should not happen) -> skip the check
    OPEN(newunit=iu, file='berry_centre_check.dat', status='replace')
    WRITE(iu,'(a)') '# n  alpha    Re[A_nn(R=0)]        -Re[A_nn(R=0)]       w90_centre(Bohr)'
    DO n = 1, nw
      DO a = 1, 3
        aR0 = aw_r(n, n, irpt0, a)
        WRITE(iu,'(2i4,3(2x,f18.10))') n, a, REAL(aR0), -REAL(aR0), bmesh%centres(a, n)
      END DO
    END DO
    CLOSE(iu)
    WRITE(oUnit,'(a)') 'wannierlib: wrote berry_centre_check.dat (A^(W) R=0 diag vs w90 centres)'
  END SUBROUTINE melem_check_berry_centres

  ! Build the Wannier-gauge Hamiltonian ham_k = U^dag diag(eigval2) U (same as m_melem_interpolate_ham).
  SUBROUTINE melem_build_hamk(this, eig, u_matrix, u_opt, ham_k)
    TYPE(t_melem_manifold), INTENT(IN) :: this
    REAL,    INTENT(IN) :: eig(:, :)
    COMPLEX, INTENT(IN) :: u_matrix(:, :, :), u_opt(:, :, :)
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: ham_k(:, :, :)
    INTEGER :: nw, nb, nk, k, i, j, m, counter
    LOGICAL :: have_dis
    REAL, ALLOCATABLE :: eigval2(:, :), eigval_opt(:)
    nw = this%num_wann; nb = this%num_bands; nk = SIZE(u_matrix, 3); have_dis = (nb > nw)
    ALLOCATE(eigval2(nw, nk), source=0.0)
    IF (have_dis) THEN
      ALLOCATE(eigval_opt(nb))
      DO k = 1, nk
        counter = 0; eigval_opt = 0.0
        DO j = 1, nb
          IF (eig(j, k) >= this%dis_win_min .AND. eig(j, k) <= this%dis_win_max) THEN
            counter = counter + 1; eigval_opt(counter) = eig(j, k)
          END IF
        END DO
        DO m = 1, nw
          DO i = 1, counter
            eigval2(m, k) = eigval2(m, k) + eigval_opt(i) * ABS(u_opt(i, m, k))**2
          END DO
        END DO
      END DO
      DEALLOCATE(eigval_opt)
    ELSE
      eigval2(1:nw, :) = eig(1:nw, :)
    END IF
    ALLOCATE(ham_k(nw, nw, nk), source=CMPLX(0.0, 0.0))
    DO k = 1, nk
      DO j = 1, nw
        DO i = 1, nw
          DO m = 1, nw
            ham_k(i, j, k) = ham_k(i, j, k) + eigval2(m, k) * CONJG(u_matrix(m, i, k)) * u_matrix(m, j, k)
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(eigval2)
  END SUBROUTINE melem_build_hamk

  ! Write H(R) in Wannier90 seedname_hr.dat format (energies in eV). Rank-0 only.
  SUBROUTINE melem_write_hr(this, cell, kpts, ham_k, wfpref)
    TYPE(t_melem_manifold), INTENT(IN) :: this
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    COMPLEX, INTENT(IN) :: ham_k(:, :, :)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: wfpref   ! seedname prefix 'WF1'/'WF2' (collinear jspins=2 channel); default 'WF1'
    COMPLEX, ALLOCATABLE :: hr(:, :, :)
    COMPLEX, ALLOCATABLE :: hr4(:, :, :, :)   ! one-component view for the shared writer
    INTEGER, ALLOCATABLE :: irvec(:, :), ndegen(:)
    INTEGER :: nrpts, nw, irpt, i, j, iu, c
    CHARACTER(LEN=64) :: fn
    CALL melem_ft_to_real(cell, ham_k, kpts, hr, irvec, ndegen, nrpts)
    nw = this%num_wann
    fn = 'WF1_hr'
    IF (PRESENT(wfpref)) fn = TRIM(wfpref)//'_hr'
    ALLOCATE(hr4(nw, nw, nrpts, 1)); hr4(:, :, :, 1) = hr(1:nw, 1:nw, 1:nrpts)
    CALL melem_write_realspace(hr4, irvec, ndegen, nrpts, nw, 1, 'hr', TRIM(fn)//'.dat', 0)
    DEALLOCATE(hr4)
    WRITE(oUnit,'(a)') 'wannierlib: wrote '//TRIM(fn)//'.dat (H(R), eV)'
    DEALLOCATE(hr, irvec, ndegen)
  END SUBROUTINE melem_write_hr

  ! Write A(R) = <0n|r_alpha|Rm> in Wannier90 seedname_r.dat format (positions in Angstrom). Rank-0.
  ! aw_r is the already-reduced Berry connection A^(W)(R) (= A(R)), so no Fourier transform here.
  SUBROUTINE melem_write_ar(this, aw_r, irvec, nrpts, wfpref)
    TYPE(t_melem_manifold), INTENT(IN) :: this
    COMPLEX, INTENT(IN) :: aw_r(:, :, :, :)          ! (nw,nw,nrpts,3) reduced A(R)
    INTEGER, INTENT(IN) :: irvec(:, :), nrpts
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: wfpref   ! seedname prefix 'WF1'/'WF2' (collinear jspins=2 channel); default 'WF1'
    INTEGER :: nw, irpt, i, j, a, iu
    CHARACTER(LEN=64) :: fn
    REAL, PARAMETER :: bohr2ang = 0.5291772109
    nw = this%num_wann
    fn = 'WF1_r'
    IF (PRESENT(wfpref)) fn = TRIM(wfpref)//'_r'
    CALL melem_write_realspace(aw_r, irvec, [(0, i = 1, nrpts)], nrpts, nw, 3, 'r', TRIM(fn)//'.dat', 0)
    WRITE(oUnit,'(a)') 'wannierlib: wrote '//TRIM(fn)//'.dat (A(R), Ang)'
  END SUBROUTINE melem_write_ar

  ! ---------------------------------------------------------------------------
  ! <operators_r>: real-space Wannier operator matrices O(R) (Fourier step 3),
  ! in Wannier90/FLEUR standalone format for external transport post-processing.
  ! No band interpolation. Rank 0 writes; spin/orbital/spin_orbit use a distributed FT-reduce
  ! over MPI ranks. Files: WF1_hr.dat (H, eV), WF1_r.dat (position, Ang), rspauli.1 (spin),
  ! anglmomrs.1 (orbital), rssocmat.1 (SOC), wig_vectors.
  ! ---------------------------------------------------------------------------
  SUBROUTINE melem_write_operators_r(this, request, cell, kpts, eig, u_matrix, u_opt, &
                                          s0_loc, l0_loc, soc4_loc, bmesh, distk, mpi_comm, mmn_loc, irank, wf_channel, l_collinear, l0col_loc)
    TYPE(t_melem_manifold), INTENT(IN) :: this
    TYPE(t_melem_request), INTENT(IN) :: request
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    REAL,    INTENT(IN) :: eig(:, :)
    COMPLEX, INTENT(IN) :: u_matrix(:, :, :)         ! (nw,nw,nk) MLWF gauge
    COMPLEX, INTENT(IN) :: u_opt(:, :, :)            ! (nb,nw,nk) disentangled (amn_local)
    COMPLEX, INTENT(IN) :: s0_loc(:, :, :, :), l0_loc(:, :, :, :, :), soc4_loc(:, :, :, :) ! per-rank coarse slices
    TYPE(t_melem_bmesh), INTENT(IN) :: bmesh         ! b-shell weights (position/Berry operator)
    COMPLEX, INTENT(IN) :: mmn_loc(:, :, :, :)       ! (nb,nb,nntot,nk_loc) this rank's overlap slice (position/Berry)
    INTEGER, INTENT(IN) :: distk(:), mpi_comm, irank
    INTEGER, INTENT(IN) :: wf_channel               ! collinear jspins=2 spin channel (1/2); 1 in the spinor case
    LOGICAL, INTENT(IN) :: l_collinear              ! collinear jspins=2 (spin channels separable, per-channel WFn files)
    COMPLEX, INTENT(IN), OPTIONAL :: l0col_loc(:, :, :, :)  ! (nb,nb,3,nk_loc) collinear per-channel Bloch orbital L
    INTEGER :: iop, nb, ik, j, nkl, aw_nrpts
    LOGICAL :: l_wig_done
    CHARACTER(LEN=8) :: wfpref                       ! 'WF1'/'WF2' seedname prefix for H(R)/position
    INTEGER, ALLOCATABLE :: gk_loc(:), aw_irvec(:, :), aw_ndegen(:)
    COMPLEX, ALLOCATABLE :: ham_k(:, :, :), aw_r(:, :, :, :), o0l(:, :, :, :), vloc(:, :, :)
    IF (.NOT. request%l_operators_r .OR. request%n_op_r < 1) RETURN   ! all ranks (reduce is collective)
    CALL timestart('melem_write_operators_r')
    nb = this%num_bands
    nkl = COUNT(distk == irank); ALLOCATE(gk_loc(nkl)); j = 0
    DO ik = 1, SIZE(distk)
      IF (distk(ik) == irank) THEN; j = j + 1; gk_loc(j) = ik; END IF
    END DO
    ! precompute the gauge V(gk)=u_opt.u_matrix for this rank's k-slice ONCE (reused by spin/orbital/soc)
    ALLOCATE(vloc(this%num_bands, this%num_wann, MAX(1, nkl)), source=CMPLX(0.0,0.0))
    DO j = 1, nkl
      vloc(:, :, j) = MATMUL(u_opt(:, :, gk_loc(j)), u_matrix(:, :, gk_loc(j)))
    END DO
    ! collinear jspins=2: per-channel seedname WF1/WF2 for H(R)/position; spinor case keeps WF1.
    WRITE(wfpref, '(a,i0)') 'WF', wf_channel
    IF (irank == 0) THEN; l_wig_done = .FALSE.; CALL melem_write_wig_once(cell, kpts, l_wig_done); END IF
    DO iop = 1, request%n_op_r
      SELECT CASE (TRIM(request%op_r_name(iop)))
      CASE ('hamiltonian')   ! cheap, no coarse array -> rank 0 serial
        IF (irank == 0) THEN
          CALL melem_build_hamk(this, eig, u_matrix, u_opt, ham_k)
          CALL melem_write_hr(this, cell, kpts, ham_k, TRIM(wfpref))
          DEALLOCATE(ham_k)
        END IF
      CASE ('position')      ! Berry connection A^(W)(R)=A(R): distributed reduce over the local overlaps
        CALL melem_build_berry_aw_r(this, cell, kpts, mmn_loc, gk_loc, u_opt, u_matrix, bmesh, mpi_comm, &
                                         aw_r, aw_irvec, aw_ndegen, aw_nrpts)
        IF (irank == 0) CALL melem_write_ar(this, aw_r, aw_irvec, aw_nrpts, TRIM(wfpref))
        IF (ALLOCATED(aw_r)) DEALLOCATE(aw_r, aw_irvec, aw_ndegen)
      CASE ('spin')          ! distributed reduce over the coarse spin slice
        ! collinear: the spin operator is combined (2N) across both channels, so it is written once
        ! after both wannierisations by melem_rspauli_collinear (main), not per channel here.
        IF (l_collinear) THEN
          CONTINUE
        ELSE
          CALL melem_op_rs_distributed(this, cell, kpts, vloc, s0_loc, gk_loc, 3, mpi_comm, irank, .FALSE., 'rspauli.1')
        END IF
      CASE ('orbital')
        IF (l_collinear) THEN
          ! per-channel Bloch L built in main's mmn loop (single spin) -> reduce -> anglmomrs.{1,2}
          IF (PRESENT(l0col_loc)) THEN
            CALL melem_op_rs_distributed(this, cell, kpts, vloc, l0col_loc, gk_loc, 3, mpi_comm, irank, .FALSE., &
                                              'anglmomrs.'//ACHAR(48+wf_channel))
          ELSE IF (irank == 0) THEN
            WRITE(oUnit,'(a)') 'wannierlib operators_r: collinear orbital slice missing -> skipped'
          END IF
        ELSE
          ALLOCATE(o0l(nb, nb, 3, SIZE(l0_loc, 5))); o0l = SUM(l0_loc, DIM=4)
          CALL melem_op_rs_distributed(this, cell, kpts, vloc, o0l, gk_loc, 3, mpi_comm, irank, .FALSE., 'anglmomrs.1')
          DEALLOCATE(o0l)
        END IF
      CASE ('spin_orbit')
        IF (l_collinear) THEN
          IF (irank == 0) WRITE(oUnit,'(a)') 'wannierlib operators_r: spin_orbit has no collinear (no-SOC) meaning -> skipped'
        ELSE
          CALL melem_op_rs_distributed(this, cell, kpts, vloc, soc4_loc, gk_loc, 4, mpi_comm, irank, .TRUE., 'rssocmat.1')
        END IF
      CASE DEFAULT
        !> Accepted by the request, so the gap is this dispatch and not the input.
        CALL judft_bug('melem_operators_r: "'//TRIM(request%op_r_name(iop))// &
                       '" is an accepted operator with no branch in this pass')
      END SELECT
    END DO
    DEALLOCATE(gk_loc, vloc)
    CALL timestop('melem_write_operators_r')
  END SUBROUTINE melem_write_operators_r

  ! Distributed real-space export of a coarse operator (spin/orbital/soc): each rank builds
  ! ow_loc = V(gk)^dagger o0_loc V(gk) for its k-slice, FT-reduces to O(R), rank 0 writes.
  ! is_soc=.TRUE. -> rssocmat.1 format (R i j jj ii); else rspauli/anglmomrs (R i j comp).
  SUBROUTINE melem_op_rs_distributed(this, cell, kpts, vloc, o0_loc, gk_loc, ncomp, mpi_comm, irank, is_soc, fname)
    TYPE(t_melem_manifold), INTENT(IN) :: this
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    COMPLEX, INTENT(IN) :: vloc(:, :, :)                      ! (nb,nw,>=nk_loc) precomputed gauge V(gk)
    COMPLEX, INTENT(IN) :: o0_loc(:, :, :, :)                 ! (nb,nb,ncomp,>=nk_loc) local Bloch op
    INTEGER, INTENT(IN) :: gk_loc(:), ncomp, mpi_comm, irank
    LOGICAL, INTENT(IN) :: is_soc
    CHARACTER(LEN=*), INTENT(IN) :: fname
    COMPLEX, ALLOCATABLE :: ow_loc(:, :, :, :), or_(:, :, :, :), o1(:, :, :), tmp(:, :)
    INTEGER, ALLOCATABLE :: irvec(:, :), ndegen(:)
    INTEGER :: nb, nw, nkl, kl, a, nrpts, irpt, i, j, kk, ii, jj, c, iu
    nb = SIZE(o0_loc, 1); nw = this%num_wann; nkl = SIZE(gk_loc)
    ALLOCATE(ow_loc(nw, nw, ncomp, MAX(1, nkl)), source=CMPLX(0.0, 0.0))
    ALLOCATE(tmp(nb, nw))
    DO kl = 1, nkl
      DO a = 1, ncomp
        tmp = MATMUL(o0_loc(:, :, a, kl), vloc(:, :, kl))
        ow_loc(:, :, a, kl) = MATMUL(CONJG(TRANSPOSE(vloc(:, :, kl))), tmp)
      END DO
    END DO
    DEALLOCATE(tmp)
    DO a = 1, ncomp
      CALL melem_ft_to_real_reduce(cell, kpts, ow_loc(:, :, a, :), gk_loc, mpi_comm, o1, irvec, ndegen, nrpts)
      IF (a == 1) ALLOCATE(or_(nw, nw, nrpts, ncomp))
      or_(:, :, :, a) = o1; DEALLOCATE(o1)
    END DO
    DEALLOCATE(ow_loc)
    IF (irank == 0) THEN
      CALL melem_write_realspace(or_, irvec, ndegen, nrpts, nw, ncomp, &
                                 MERGE('soc    ', 'generic', is_soc), TRIM(fname), irank)
      WRITE(oUnit, '(a,i0,a)') 'wannierlib: wrote '//TRIM(fname)//' (', nrpts, ' R-vectors, distributed FT)'
    END IF
    IF (ALLOCATED(or_)) DEALLOCATE(or_)
    IF (ALLOCATED(irvec)) DEALLOCATE(irvec, ndegen)
  END SUBROUTINE melem_op_rs_distributed

  ! wig_vectors (once): idx R1 R2 R3 ndegen (list-directed; = orbitrans rpts).
  SUBROUTINE melem_write_wig_once(cell, kpts, done)
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    LOGICAL, INTENT(INOUT) :: done
    INTEGER, ALLOCATABLE :: irvec(:, :), ndegen(:)
    INTEGER :: nrpts, ii, iu
    IF (done) RETURN
    CALL melem_ws_vectors(cell, kpts%nkpt3, irvec, ndegen, nrpts)
    OPEN(newunit=iu, file='wig_vectors', status='replace')
    DO ii = 1, nrpts
      WRITE(iu, *) ii, irvec(1, ii), irvec(2, ii), irvec(3, ii), ndegen(ii)
    END DO
    CLOSE(iu)
    WRITE(oUnit, '(a,i0,a)') 'wannierlib: wrote wig_vectors (', nrpts, ' R-vectors)'
    DEALLOCATE(irvec, ndegen)
    done = .TRUE.
  END SUBROUTINE melem_write_wig_once

END MODULE m_melem_operators_r
