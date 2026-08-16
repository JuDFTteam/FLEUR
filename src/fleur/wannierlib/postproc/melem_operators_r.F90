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
!>  The b-shell weights and the reference Wannier centres arrive as a plain t_melem_bmesh,
!>  so this module is free of the Wannier90 library.
MODULE m_melem_operators_r
  USE m_juDFT
  USE m_constants, ONLY : oUnit, hartree_to_ev_const, tpi_const
  USE m_types_cell
  USE m_types_kpts
  USE m_types_melem_manifold, ONLY : t_melem_manifold
  USE m_types_melem_request, ONLY : t_melem_request
  USE m_types_melem_bmesh
  USE m_melem_ft, ONLY : melem_ft_to_real, melem_ws_vectors, melem_ft_to_real_reduce
  USE m_melem_io, ONLY : melem_write_realspace
  USE m_melem_hamk, ONLY : melem_build_hamk
  USE m_melem_coeff_a, ONLY : melem_build_berry_aw_r, melem_check_berry_centres, melem_write_ar
  USE m_melem_coeff_b, ONLY : melem_write_bmn
  USE m_melem_coeff_f, ONLY : melem_write_fmn
  USE m_melem_coeff_c, ONLY : melem_write_cmn
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: melem_write_operators_r, melem_op_rs_distributed
  PUBLIC :: melem_write_hr, melem_write_wig_once
CONTAINS

  ! Write H(R) in Wannier90 seedname_hr.dat format (energies in eV). Rank-0 only.
  SUBROUTINE melem_write_hr(this, cell, kpts, hr, irvec, ndegen, nrpts, wfpref)
    TYPE(t_melem_manifold), INTENT(IN) :: this
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    COMPLEX, INTENT(IN) :: hr(:, :, :)                 ! H(R), already transformed
    INTEGER, INTENT(IN) :: irvec(:, :), ndegen(:), nrpts
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: wfpref   ! seedname prefix 'WF1'/'WF2' (collinear jspins=2 channel); default 'WF1'
    COMPLEX, ALLOCATABLE :: hr4(:, :, :, :)   ! one-component view for the shared writer
    INTEGER :: nw, irpt, i, j, iu, c
    CHARACTER(LEN=64) :: fn
    nw = this%num_wann
    fn = 'WF1_hr'
    IF (PRESENT(wfpref)) fn = TRIM(wfpref)//'_hr'
    ALLOCATE(hr4(nw, nw, nrpts, 1)); hr4(:, :, :, 1) = hr(1:nw, 1:nw, 1:nrpts)
    CALL melem_write_realspace(hr4, irvec, ndegen, nrpts, nw, 1, 'hr', TRIM(fn)//'.dat', 0)
    DEALLOCATE(hr4)
    WRITE(oUnit,'(a)') 'wannierlib: wrote '//TRIM(fn)//'.dat (H(R), eV)'
  END SUBROUTINE melem_write_hr

  ! ---------------------------------------------------------------------------
  ! <operators_r>: real-space Wannier operator matrices O(R) (Fourier step 3),
  ! in Wannier90/FLEUR standalone format for external transport post-processing.
  ! No band interpolation. Rank 0 writes; spin/orbital/spin_orbit use a distributed FT-reduce
  ! over MPI ranks. Files: WF1_hr.dat (H, eV), WF1_r.dat (position, Ang), rspauli.1 (spin),
  ! anglmomrs.1 (orbital), rssocmat.1 (SOC), wig_vectors.
  ! ---------------------------------------------------------------------------
  SUBROUTINE melem_write_operators_r(this, request, cell, kpts, eig, u_matrix, u_opt, &
                                          s0_loc, l0_loc, soc4_loc, f0_loc, c0_loc, bmesh, distk, mpi_comm, mmn_loc, irank, wf_channel, l_collinear)
    TYPE(t_melem_manifold), INTENT(IN) :: this
    TYPE(t_melem_request), INTENT(IN) :: request
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    REAL,    INTENT(IN) :: eig(:, :)
    COMPLEX, INTENT(IN) :: u_matrix(:, :, :)         ! (nw,nw,nk) MLWF gauge
    COMPLEX, INTENT(IN) :: u_opt(:, :, :)            ! (nb,nw,nk) disentangled (amn_local)
    !> Per-rank coarse slices. ALLOCATABLE because a slice nobody asked for is not
    !> built at all: every use below sits behind the request that would have built it.
    COMPLEX, ALLOCATABLE, INTENT(IN) :: s0_loc(:, :, :, :), l0_loc(:, :, :, :, :, :), soc4_loc(:, :, :, :)
    !> (nw,nw,3,3,nk_loc) the geometric tensor, already contracted over the pairs of
    !> neighbours and already gauged, because that has to happen where the wavefunctions are.
    COMPLEX, ALLOCATABLE, INTENT(IN) :: f0_loc(:, :, :, :, :)
    !> (nw,nw,3,3,nk_loc) el mismo tensor con el hamiltoniano dentro. Sin asignar
    !> salvo que se pida, por lo mismo que f0_loc.
    COMPLEX, ALLOCATABLE, INTENT(IN) :: c0_loc(:, :, :, :, :)
    TYPE(t_melem_bmesh), INTENT(IN) :: bmesh         ! b-shell weights (position/Berry operator)
    COMPLEX, INTENT(IN) :: mmn_loc(:, :, :, :)       ! (nb,nb,nntot,nk_loc) this rank's overlap slice (position/Berry)
    INTEGER, INTENT(IN) :: distk(:), mpi_comm, irank
    INTEGER, INTENT(IN) :: wf_channel               ! collinear jspins=2 spin channel (1/2); 1 in the spinor case
    LOGICAL, INTENT(IN) :: l_collinear              ! collinear jspins=2 (spin channels separable, per-channel WFn files)
    INTEGER :: iop, nb, ik, j, nkl, aw_nrpts
    LOGICAL :: l_wig_done
    CHARACTER(LEN=8) :: wfpref                       ! 'WF1'/'WF2' seedname prefix for H(R)/position
    INTEGER, ALLOCATABLE :: gk_loc(:), aw_irvec(:, :), aw_ndegen(:)
    INTEGER, ALLOCATABLE :: hr_irvec(:, :), hr_ndegen(:)
    COMPLEX, ALLOCATABLE :: ham_r(:, :, :)
    INTEGER :: hr_nrpts
    LOGICAL :: l_need_hr
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
    !> H(R) is wanted by the Hamiltonian export, and its transform is what gives B(R) the
    !> Wigner-Seitz R set it sums over. Built once, and on EVERY rank
    !> rather than only where it is written: melem_write_bmn is collective, so all of them
    !> reach it, and handing an unallocated actual to a non-allocatable dummy is undefined
    !> behaviour that a serial run can never show. It needs only eig, u_matrix and u_opt,
    !> which every rank already holds, so the duplicated work is the cheap way out.
    l_need_hr = request%has_op_r('hamiltonian') .OR. request%has_op_r('bmn') &
                .OR. request%has_op_r('fmn') .OR. request%has_op_r('cmn')
    IF (l_need_hr) THEN
      CALL melem_build_hamk(this, eig, u_matrix, u_opt, ham_k)
      CALL melem_ft_to_real(cell, ham_k, kpts, ham_r, hr_irvec, hr_ndegen, hr_nrpts)
      DEALLOCATE(ham_k)
    END IF
    DO iop = 1, request%n_op_r
      SELECT CASE (TRIM(request%op_r_name(iop)))
      CASE ('hamiltonian')   ! H(R) already built above
        IF (irank == 0) CALL melem_write_hr(this, cell, kpts, ham_r, hr_irvec, hr_ndegen, &
                                            hr_nrpts, TRIM(wfpref))
      CASE ('bmn')           ! collective: the k-sum is distributed like the others
        CALL melem_write_bmn(this, kpts, eig, u_matrix, u_opt, mmn_loc, gk_loc, &
                             bmesh, hr_irvec, hr_nrpts, mpi_comm, irank, TRIM(wfpref))
      CASE ('fmn')           ! collective: only the sum over R is left, the gauge is already in
        CALL melem_write_fmn(this, kpts, f0_loc, gk_loc, hr_irvec, hr_nrpts, mpi_comm, &
                             irank, TRIM(wfpref))
      CASE ('cmn')           ! collective: como fmn, el gauge ya viene aplicado
        CALL melem_write_cmn(this, kpts, c0_loc, gk_loc, hr_irvec, hr_nrpts, mpi_comm, &
                             irank, TRIM(wfpref))
      CASE ('position')      ! Berry connection A^(W)(R)=A(R): distributed reduce over the local overlaps
        CALL melem_build_berry_aw_r(this, cell, kpts, mmn_loc, gk_loc, u_opt, u_matrix, bmesh, mpi_comm, &
                                         aw_r, aw_irvec, aw_ndegen, aw_nrpts)
        !> The R=0 diagonal is the Wannier centre, so it can be checked against the
        !> centres the b-mesh carries whenever A is built -- not only when a band
        !> interpolation happens to be asked for as well.
        IF (irank == 0) CALL melem_check_berry_centres(this, aw_r, aw_irvec, aw_nrpts, bmesh)
        IF (irank == 0) CALL melem_write_ar(this, aw_r, aw_irvec, aw_nrpts, TRIM(wfpref))
        IF (ALLOCATED(aw_r)) DEALLOCATE(aw_r, aw_irvec, aw_ndegen)
      CASE ('position_pw90') ! A(R) in the postw90 convention: Eq. (44) on the diagonal too,
                             ! plus A <- (A + A^dagger)/2 at each k. `position` is NOT touched:
                             ! orbitrans and the Wannier-centre check need the _r.dat form, and
                             ! berry.F90 refuses that form for the orbital magnetisation
                             ! ("transl_inv=T disabled for morb"). Written as WF*_rpw.dat.
        CALL melem_build_berry_aw_r(this, cell, kpts, mmn_loc, gk_loc, u_opt, u_matrix, bmesh, mpi_comm, &
                                         aw_r, aw_irvec, aw_ndegen, aw_nrpts, l_pw90=.TRUE.)
        !> No centre check here: the R=0 diagonal of this variant is not the Wannier centre.
        IF (irank == 0) CALL melem_write_ar(this, aw_r, aw_irvec, aw_nrpts, TRIM(wfpref), '_rpw')
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
        ! Site-summed L of this channel -> reduce -> anglmomrs.<channel>. The spinor case is
        ! the same thing with one channel, so there is nothing to tell apart here.
        ALLOCATE(o0l(nb, nb, 3, SIZE(l0_loc, 6)))
        o0l = SUM(l0_loc(:, :, :, :, wf_channel, :), DIM=4)
        CALL melem_op_rs_distributed(this, cell, kpts, vloc, o0l, gk_loc, 3, mpi_comm, irank, .FALSE., &
                                          'anglmomrs.'//ACHAR(48+wf_channel))
        DEALLOCATE(o0l)
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
    IF (ALLOCATED(ham_r)) DEALLOCATE(ham_r, hr_irvec, hr_ndegen)
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
