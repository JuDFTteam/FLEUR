!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Operator-agnostic Wannier-gauge interpolation driver (framework Class A/B/C).
!>
!>  Given the Bloch-basis operator matrices O0_alpha(k) (nb x nb, alpha=1..ncomp),
!>  it does the shared pipeline for ANY operator:
!>    O_W,alpha(k) = V^dagger O0_alpha(k) V ,   V = U_dis U
!>    O_alpha(k')  = FT[ O_W,alpha ]            (shared core m_melem_ft)
!>    H(k')        = FT[ H_W ] -> diag -> E_n(k'), C(k')
!>    <O_alpha>_n(k') = [ C^dagger O_alpha(k') C ]_nn
!>  and writes <outfile>.dat: kdist, [ E_n(eV), <O_1>_n, ..., <O_ncomp>_n ] per band.
!>
!>  A new operator only supplies its O0(k) (a provider) and calls this with the
!>  right ncomp/outfile -- steps above are never rewritten. Master rank only.
MODULE m_melem_interpolate_op
  USE m_juDFT
  USE m_constants, ONLY : oUnit, hartree_to_ev_const
  USE m_types_cell
  USE m_types_kpts
  USE m_types_melem_manifold, ONLY: t_melem_manifold
  USE m_melem_hamk, ONLY : melem_build_hamk
  USE m_melem_ft, ONLY : melem_ft_interpolate, melem_ft_to_real_reduce, melem_ft_rtok
  USE m_melem_interp_util, ONLY : melem_kpath, melem_zheev_workspace
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: melem_interpolate_operator
CONTAINS

  SUBROUTINE melem_interpolate_operator(this, cell, kpts, eig, u_matrix, u_opt, o0_loc, gk_loc, ncomp, kfrac, outfile, irank, mpicm)
    TYPE(t_melem_manifold), INTENT(IN) :: this
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    REAL,    INTENT(IN) :: eig(:, :)              ! (num_bands, nk)
    COMPLEX, INTENT(IN) :: u_matrix(:, :, :)      ! (num_wann, num_wann, nk)  MLWF gauge (full mesh)
    COMPLEX, INTENT(IN) :: u_opt(:, :, :)         ! (num_bands, num_wann, nk) disentangled (full mesh)
    COMPLEX, INTENT(IN) :: o0_loc(:, :, :, :)     ! (num_bands, num_bands, ncomp, nk_loc) this rank's Bloch slice
    INTEGER, INTENT(IN) :: gk_loc(:)              ! (nk_loc) global k index of each slice entry
    INTEGER, INTENT(IN) :: ncomp
    CHARACTER(LEN=*), INTENT(IN) :: outfile
    !> The domain's k-set and the names its files take, both decided by the caller: the
    !> k-points come from a named kPointList and the names from the exposure table plus the
    !> domain suffix. Unallocated off rank 0, which never reaches them.
    REAL, ALLOCATABLE, INTENT(IN) :: kfrac(:, :)          !> (3, np) fractional mesh
    INTEGER, INTENT(IN) :: irank, mpicm

    INTEGER :: num_wann, num_bands, m, ip, np, iu, info, lwork, a
    INTEGER :: nkl, kl, nrpts
    REAL,    ALLOCATABLE :: kdist(:), evals(:), rwork(:), oexp(:)
    COMPLEX, ALLOCATABLE :: ham_k(:, :, :), H_interp(:, :, :), o_interp(:, :, :, :)
    COMPLEX, ALLOCATABLE :: hk(:, :), work(:), vloc(:, :, :), tmp(:, :), cvec(:, :), oc(:, :, :)
    COMPLEX, ALLOCATABLE :: ow_loc(:, :, :, :), o_r(:, :, :, :), o1(:, :, :)
    INTEGER, ALLOCATABLE :: irvec(:, :), ndegen(:)

    num_wann  = this%num_wann
    num_bands = this%num_bands
    CALL timestart('melem_interpolate_operator')

    ! ---- PHASE A (ALL ranks): O_W,alpha(k) = V(gk)^dagger O0_alpha V(gk) on this rank's k-slice,
    !      then coarse -> real space O_alpha(R) via the distributed FT-reduce (collective). ----
    nkl = SIZE(gk_loc)
    ALLOCATE(vloc(num_bands, num_wann, MAX(1, nkl)), source=CMPLX(0.0, 0.0))
    DO kl = 1, nkl
      vloc(:, :, kl) = MATMUL(u_opt(:, :, gk_loc(kl)), u_matrix(:, :, gk_loc(kl)))
    END DO
    ALLOCATE(ow_loc(num_wann, num_wann, ncomp, MAX(1, nkl)), source=CMPLX(0.0, 0.0))
    ALLOCATE(tmp(num_bands, num_wann))
    DO kl = 1, nkl
      DO a = 1, ncomp
        tmp = MATMUL(o0_loc(:, :, a, kl), vloc(:, :, kl))
        ow_loc(:, :, a, kl) = MATMUL(CONJG(TRANSPOSE(vloc(:, :, kl))), tmp)
      END DO
    END DO
    DEALLOCATE(tmp, vloc)
    DO a = 1, ncomp
      CALL melem_ft_to_real_reduce(cell, kpts, ow_loc(:, :, a, :), gk_loc, mpicm, o1, irvec, ndegen, nrpts)
      IF (a == 1) ALLOCATE(o_r(num_wann, num_wann, nrpts, ncomp))
      o_r(:, :, :, a) = o1; DEALLOCATE(o1)
    END DO
    DEALLOCATE(ow_loc)

    ! only rank 0 does the R -> fine-path interpolation, diagonalization and write
    IF (irank /= 0) THEN
      IF (ALLOCATED(o_r)) DEALLOCATE(o_r)
      IF (ALLOCATED(irvec)) DEALLOCATE(irvec, ndegen)
      CALL timestop('melem_interpolate_operator'); RETURN
    END IF

    np = SIZE(kfrac, 2)   ! the caller resolved the domain; there is nothing to skip

    CALL melem_kpath(cell, kfrac, kdist)   ! abscissa of the output, from the mesh just read

    ! ---- H_W(k) via eigval2 (same construction as the validated band driver), full mesh on rank 0 ----
    CALL melem_build_hamk(this, eig, u_matrix, u_opt, ham_k)

    ! ---- interpolate H (full mesh -> R -> k', shared core) and the operator (R -> k' only:
    !      O_alpha(R) is already assembled by the distributed reduce above) ----
    CALL melem_ft_interpolate(cell, ham_k, kpts, kfrac, H_interp)
    ALLOCATE(o_interp(num_wann, num_wann, ncomp, np))
    BLOCK
      COMPLEX, ALLOCATABLE :: o_one(:, :, :)
      DO a = 1, ncomp
        CALL melem_ft_rtok(o_r(:, :, :, a), irvec, ndegen, nrpts, kfrac, o_one)
        o_interp(:, :, a, :) = o_one
      END DO
    END BLOCK
    DEALLOCATE(o_r, irvec, ndegen)

    ! ---- diagonalize H(k') with eigenvectors, project operator, write ----
    ALLOCATE(evals(num_wann), hk(num_wann, num_wann), cvec(num_wann, num_wann), &
             oc(num_wann, num_wann, ncomp), oexp(ncomp))
    CALL melem_zheev_workspace('V', num_wann, work, rwork, lwork)

    OPEN(newunit=iu, file=TRIM(outfile)//'.dat', status='replace')
    WRITE(iu,'(a,i0,a)') '# kdist   [ E_n(eV)  <O_1>_n .. <O_', ncomp, '>_n ] for n=1..num_wann'
    DO ip = 1, np
      hk = H_interp(:, :, ip)
      CALL zheev('V', 'U', num_wann, hk, num_wann, evals, work, lwork, rwork, info)
      IF (info /= 0) CALL juDFT_error('zheev failed', calledby='melem_interpolate_operator')
      cvec = hk
      DO a = 1, ncomp
        oc(:, :, a) = MATMUL(o_interp(:, :, a, ip), cvec)
      END DO
      WRITE(iu,'(f12.6)', advance='no') kdist(ip)
      DO m = 1, num_wann
        DO a = 1, ncomp
          oexp(a) = REAL(DOT_PRODUCT(cvec(:, m), oc(:, m, a)))
        END DO
        WRITE(iu,'(2x,f14.8)', advance='no') hartree_to_ev_const*evals(m)
        DO a = 1, ncomp
          WRITE(iu,'(2x,f12.6)', advance='no') oexp(a)
        END DO
      END DO
      WRITE(iu,'(a)') ''
    END DO
    CLOSE(iu)
    WRITE(oUnit,'(a,i0,a)') 'wannierlib operator interpolation: wrote '//TRIM(outfile)//'.dat (', np, ' k-points)'
    CALL timestop('melem_interpolate_operator')
  END SUBROUTINE melem_interpolate_operator

END MODULE m_melem_interpolate_op
