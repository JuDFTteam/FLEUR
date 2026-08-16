!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  A(R) = <0n|r|Rm>, the Wannier-gauge Berry connection in real space.
!>
!>  Off the diagonal this is Eq. (44) of Wang, Yates, Souza and Vanderbilt, PRB 74,
!>  195118 (2006); on it, Eqs. (27,29,32) of Marzari and Vanderbilt, PRB 56, 12847
!>  (1997). The R=0 diagonal is the Wannier centre, which is what the check asserts.
MODULE m_melem_coeff_a
  USE m_juDFT
  USE m_constants, ONLY : oUnit
  USE m_types_cell
  USE m_types_kpts
  USE m_types_melem_manifold, ONLY : t_melem_manifold
  USE m_types_melem_bmesh
  USE m_melem_ft, ONLY : melem_ft_to_real_reduce
  USE m_melem_io, ONLY : melem_write_realspace
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: melem_build_berry_aw_r, melem_check_berry_centres, melem_write_ar
CONTAINS

  ! Build the Wannier-gauge Berry connection in real space A^(W)_alpha(R), distributed:
  ! each rank forms A^(W)_alpha(k) on its OWN k-slice, off the diagonal as
  ! i sum_b w_b b_alpha M^(W,b)_nm(k) and on it as -sum_b w_b b_alpha Im log M^(W,b)_nn(k)
  ! (M^(W,b)(k) = V(k)^dagger M^(0,b)(k) V(k_b), V = u_opt.u_matrix, k_b = nnlist(b,k)) from the
  ! local overlaps mmn_loc (global k = gk_loc), then reduces coarse -> R with the distributed
  ! FT-reduce (collective). The full-mesh overlaps are never gathered. A^(W)(R) is exactly the
  ! position operator A(R) = <0n|r|Rm> and the R->k' interpolant used by the velocity/Berry curvature.
  ! kmesh (wb/bk/nnlist) arrives as a plain t_melem_bmesh -> no Wannier90 coupling here.
  ! Collective over mpi_comm; result valid on all ranks.
  SUBROUTINE melem_build_berry_aw_r(this, cell, kpts, mmn_loc, gk_loc, u_opt, u_matrix, bmesh, mpi_comm, &
                                         aw_r, irvec, ndegen, nrpts, l_pw90)
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
    !> Diagonal convention. Absent or .FALSE. gives the `_r.dat` of Wannier90
    !> (hamiltonian_write_rmn): Eq. (44) WYSV06 off the diagonal, log branch on it, no
    !> hermitisation. That is what the Wannier centres and orbitrans expect, so it must
    !> stay the default. .TRUE. gives the postw90 object (get_AA_R): Eq. (44) everywhere
    !> plus A <- (A + A^dagger)/2 at each k, which is what berry.F90 uses for the orbital
    !> magnetisation -- and it REFUSES the log branch there ("transl_inv=T disabled for
    !> morb"). The two agree to first order only while |M_nn| stays near 1, i.e. for
    !> well localised functions; they part company exactly where ours are worst.
    LOGICAL, INTENT(IN), OPTIONAL :: l_pw90
    LOGICAL :: pw90
    INTEGER :: nb, nw, nk, nnt, k, kb, nn, a, i, kl, nkl
    REAL :: wb, b(3)
    COMPLEX, ALLOCATABLE :: Vk(:, :), Vkb(:, :), Mw(:, :), tmp(:, :), aw_loc(:, :, :, :), a1(:, :, :)
    COMPLEX, ALLOCATABLE :: Mb(:, :)

    IF (bmesh%nntot < 1 .OR. .NOT. ALLOCATED(bmesh%nnlist)) THEN
      ! no b-mesh available (e.g. built without the Wannier90 library) -> nothing to build
      nrpts = 0
      ALLOCATE(aw_r(1, 1, 1, 3), irvec(3, 1), ndegen(1))
      RETURN
    END IF
    nb = this%num_bands; nw = this%num_wann; nk = kpts%nkptf
    pw90 = .FALSE.
    IF (PRESENT(l_pw90)) pw90 = l_pw90
    nnt = bmesh%nntot
    nkl = SIZE(gk_loc)
    ALLOCATE(aw_loc(nw, nw, 3, MAX(1, nkl)), source=CMPLX(0.0, 0.0))
    ALLOCATE(Vk(nb, nw), Vkb(nb, nw), Mw(nw, nw), Mb(nw, nw), tmp(nb, nw))
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
        ! The diagonal takes the phase of the overlap and not its linear part, so it
        ! stays real and reproduces the Wannier centre to the accuracy of the b-shell.
        ! It needs M^(W,b)_nn /= 0, which a converged wannierisation gives.
        Mb = CMPLX(0.0, 1.0) * Mw
        IF (.NOT. pw90) THEN
          DO i = 1, nw
            Mb(i, i) = CMPLX(-AIMAG(LOG(Mw(i, i))), 0.0)
          END DO
        END IF
        DO a = 1, 3
          aw_loc(:, :, a, kl) = aw_loc(:, :, a, kl) + wb*b(a) * Mb
        END DO
      END DO
      !> Eq. (44) does not preserve hermiticity, so postw90 restores it by hand at each k
      !> (get_oper.F90). Done here, before the Fourier transform, so that A(k) is hermitian
      !> for every k and not merely on average.
      IF (pw90) THEN
        DO a = 1, 3
          aw_loc(:, :, a, kl) = 0.5 * (aw_loc(:, :, a, kl) &
                                       + CONJG(TRANSPOSE(aw_loc(:, :, a, kl))))
        END DO
      END IF
    END DO
    DEALLOCATE(Vk, Vkb, Mw, Mb, tmp)
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

  ! Write A(R) = <0n|r_alpha|Rm> in Wannier90 seedname_r.dat format (positions in Angstrom). Rank-0.
  ! aw_r is the already-reduced Berry connection A^(W)(R) (= A(R)), so no Fourier transform here.
  SUBROUTINE melem_write_ar(this, aw_r, irvec, nrpts, wfpref, suffix)
    TYPE(t_melem_manifold), INTENT(IN) :: this
    COMPLEX, INTENT(IN) :: aw_r(:, :, :, :)          ! (nw,nw,nrpts,3) reduced A(R)
    INTEGER, INTENT(IN) :: irvec(:, :), nrpts
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: wfpref   ! seedname prefix 'WF1'/'WF2' (collinear jspins=2 channel); default 'WF1'
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: suffix   ! '_r' (default) or '_rpw' for the postw90 variant
    INTEGER :: nw, irpt, i, j, a, iu
    CHARACTER(LEN=64) :: fn
    REAL, PARAMETER :: bohr2ang = 0.5291772109
    nw = this%num_wann
    fn = 'WF1_r'
    IF (PRESENT(wfpref)) fn = TRIM(wfpref)//'_r'
    IF (PRESENT(suffix)) THEN
      fn = 'WF1'//TRIM(suffix)
      IF (PRESENT(wfpref)) fn = TRIM(wfpref)//TRIM(suffix)
    END IF
    CALL melem_write_realspace(aw_r, irvec, [(0, i = 1, nrpts)], nrpts, nw, 3, 'r', TRIM(fn)//'.dat', 0)
    WRITE(oUnit,'(a)') 'wannierlib: wrote '//TRIM(fn)//'.dat (A(R), Ang)'
  END SUBROUTINE melem_write_ar

END MODULE m_melem_coeff_a
