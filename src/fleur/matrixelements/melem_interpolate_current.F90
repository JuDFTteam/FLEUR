!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Spin / orbital current operator (framework Class C: algebra of already
!>  interpolated quantities). For a Cartesian transport direction alpha and a
!>  spin/orbital component beta:
!>
!>      j_{alpha,beta}(k') = 1/2 { v_alpha(k') , O_beta(k') }
!>
!>  with v = dH/dk (velocity variant of the FT core) and O = sigma (spin current)
!>  or L (orbital current), both built in the Wannier gauge and interpolated. All
!>  9 components (alpha,beta = x,y,z) are projected on the band eigenvectors and
!>  written: <j_{alpha,beta}>_n = [ C^dagger j_{alpha,beta}(k') C ]_nn.
!>
!>  NOTE (gauge): v and O are taken in the Wannier gauge; the diagonal band
!>  velocity is exact but the off-diagonal velocity omits the Berry-connection
!>  (position-operator) term. The projected current is therefore the leading
!>  (intraband) contribution -- exact for the diagonal v, approximate where the
!>  anticommutator mixes interband elements. A rigorous version needs the position
!>  operator A(R) (Class B), added in a later stage.
!>  Output: kdist, [ E_n(eV), j_xx, j_xy, j_xz, j_yx, ... j_zz (eV*bohr) ] per band.
!>
!>  Parallelism: the operator part O_beta(k) is distributed exactly like the generic
!>  operator driver (m_melem_interpolate_op) -- each rank builds O_W,beta on its
!>  own k-slice and reduces to O_beta(R) (collective); rank 0 then does R -> fine path,
!>  the velocity build, diagonalization and write. The full-mesh coarse operator is
!>  never materialized.
MODULE m_melem_interpolate_current
  USE m_juDFT
  USE m_constants, ONLY : oUnit, hartree_to_ev_const
  USE m_types_cell
  USE m_types_kpts
  USE m_types_melem_manifold, ONLY: t_melem_manifold
  USE m_melem_hamk, ONLY : melem_build_hamk
  USE m_melem_domains, ONLY : melem_read_kset
  USE m_melem_ft, ONLY : melem_ft_to_real, melem_ft_rtok_velocity, &
                              melem_ft_to_real_reduce, melem_ft_rtok
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: melem_interpolate_current
CONTAINS

  SUBROUTINE melem_interpolate_current(this, cell, kpts, eig, u_matrix, u_opt, o0_loc, gk_loc, outfile, irank, mpicm)
    TYPE(t_melem_manifold), INTENT(IN) :: this
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    REAL,    INTENT(IN) :: eig(:, :)              ! (num_bands, nk)
    COMPLEX, INTENT(IN) :: u_matrix(:, :, :)      ! (num_wann, num_wann, nk)  MLWF gauge (full mesh)
    COMPLEX, INTENT(IN) :: u_opt(:, :, :)         ! (num_bands, num_wann, nk) disentangled (full mesh)
    COMPLEX, INTENT(IN) :: o0_loc(:, :, :, :)     ! (num_bands, num_bands, 3, nk_loc) this rank's Bloch spin OR orbital slice
    INTEGER, INTENT(IN) :: gk_loc(:)              ! (nk_loc) global k index of each slice entry
    CHARACTER(LEN=*), INTENT(IN) :: outfile
    INTEGER, INTENT(IN) :: irank, mpicm

    INTEGER :: num_wann, num_bands, nk, k, i, j, m, counter, ip, np, iu, info, lwork, al, be, c
    INTEGER :: nkl, kl, nrpts
    REAL    :: kpath, dk(3), dkc(3)
    REAL,    ALLOCATABLE :: kfrac(:, :), evals(:), rwork(:), jexp(:)
    COMPLEX, ALLOCATABLE :: ham_k(:, :, :), H_interp(:, :, :), o_interp(:, :, :, :)
    COMPLEX, ALLOCATABLE :: v_interp(:, :, :, :), hk(:, :), work(:), vloc(:, :, :), tmp(:, :), cvec(:, :)
    COMPLEX, ALLOCATABLE :: jmat(:, :, :), jc(:, :), ow_loc(:, :, :, :), o_r(:, :, :, :), o1(:, :, :)
    INTEGER, ALLOCATABLE :: irvec(:, :), ndegen(:)
    COMPLEX, ALLOCATABLE :: ham_r(:, :, :)
    INTEGER, ALLOCATABLE :: h_irvec(:, :), h_ndegen(:)
    INTEGER :: h_nrpts
    COMPLEX :: wq(1)

    num_wann  = this%num_wann
    num_bands = this%num_bands
    nk        = kpts%nkptf
    CALL timestart('melem_interpolate_current')

    ! ---- PHASE A (ALL ranks): O_W,beta(k) = V(gk)^dagger O0_beta V(gk) on this rank's k-slice,
    !      then coarse -> real space O_beta(R) via the distributed FT-reduce (collective). ----
    nkl = SIZE(gk_loc)
    ALLOCATE(vloc(num_bands, num_wann, MAX(1, nkl)), source=CMPLX(0.0, 0.0))
    DO kl = 1, nkl
      vloc(:, :, kl) = MATMUL(u_opt(:, :, gk_loc(kl)), u_matrix(:, :, gk_loc(kl)))
    END DO
    ALLOCATE(ow_loc(num_wann, num_wann, 3, MAX(1, nkl)), source=CMPLX(0.0, 0.0))
    ALLOCATE(tmp(num_bands, num_wann))
    DO kl = 1, nkl
      DO be = 1, 3
        tmp = MATMUL(o0_loc(:, :, be, kl), vloc(:, :, kl))
        ow_loc(:, :, be, kl) = MATMUL(CONJG(TRANSPOSE(vloc(:, :, kl))), tmp)
      END DO
    END DO
    DEALLOCATE(tmp, vloc)
    DO be = 1, 3
      CALL melem_ft_to_real_reduce(cell, kpts, ow_loc(:, :, be, :), gk_loc, mpicm, o1, irvec, ndegen, nrpts)
      IF (be == 1) ALLOCATE(o_r(num_wann, num_wann, nrpts, 3))
      o_r(:, :, :, be) = o1; DEALLOCATE(o1)
    END DO
    DEALLOCATE(ow_loc)

    ! only rank 0 does the R -> fine-path interpolation, velocity build, diagonalization and write
    IF (irank /= 0) THEN
      IF (ALLOCATED(o_r)) DEALLOCATE(o_r)
      IF (ALLOCATED(irvec)) DEALLOCATE(irvec, ndegen)
      CALL timestop('melem_interpolate_current'); RETURN
    END IF

    IF (.NOT. melem_read_kset(kfrac, np, 'melem_interpolate_current')) THEN
      WRITE(oUnit,'(a)') 'wannierlib current interpolation: no kpts_interpol file -> skipped'
      IF (ALLOCATED(o_r)) DEALLOCATE(o_r)
      IF (ALLOCATED(irvec)) DEALLOCATE(irvec, ndegen)
      CALL timestop('melem_interpolate_current'); RETURN
    END IF

    ! ---- H_W(k) via eigval2 (same construction as the validated band driver), full mesh on rank 0 ----
    CALL melem_build_hamk(this, eig, u_matrix, u_opt, ham_k)

    ! ---- interpolate H (eigenvectors), v = dH/dk (3), and O_beta (3, R -> k' from the reduced O(R)) ----
    !> One transform of H_W to real space, then both the interpolant and its derivative
    !> off the same H(R); they used to be two calls that each rebuilt it.
    CALL melem_ft_to_real(cell, ham_k, kpts, ham_r, h_irvec, h_ndegen, h_nrpts)
    CALL melem_ft_rtok(ham_r, h_irvec, h_ndegen, h_nrpts, kfrac, H_interp)
    CALL melem_ft_rtok_velocity(cell, ham_r, h_irvec, h_ndegen, h_nrpts, kfrac, v_interp)
    DEALLOCATE(ham_r, h_irvec, h_ndegen)
    ALLOCATE(o_interp(num_wann, num_wann, 3, np))
    BLOCK
      COMPLEX, ALLOCATABLE :: o_one(:, :, :)
      DO be = 1, 3
        CALL melem_ft_rtok(o_r(:, :, :, be), irvec, ndegen, nrpts, kfrac, o_one)
        o_interp(:, :, be, :) = o_one
      END DO
    END BLOCK
    DEALLOCATE(o_r, irvec, ndegen)

    ! ---- diagonalize H(k'); j_{al,be} = 1/2 { v_al, O_be }; project diagonal; write 9 comps ----
    ALLOCATE(evals(num_wann), hk(num_wann, num_wann), cvec(num_wann, num_wann), &
             jmat(num_wann, num_wann, 9), jc(num_wann, num_wann), jexp(9), rwork(MAX(1, 3*num_wann-2)))
    CALL zheev('V', 'U', num_wann, hk, num_wann, evals, wq, -1, rwork, info)
    lwork = MAX(1, NINT(REAL(wq(1)))); ALLOCATE(work(lwork))

    OPEN(newunit=iu, file=outfile, status='replace')
    WRITE(iu,'(a)') '# kdist   [ E_n(eV)  j_ab (a,b=x,y,z: xx xy xz yx yy yz zx zy zz, eV*bohr) ] for n=1..num_wann'
    kpath = 0.0
    DO ip = 1, np
      hk = H_interp(:, :, ip)
      CALL zheev('V', 'U', num_wann, hk, num_wann, evals, work, lwork, rwork, info)
      IF (info /= 0) CALL juDFT_error('zheev failed', calledby='melem_interpolate_current')
      cvec = hk
      IF (ip > 1) THEN
        dk = kfrac(:, ip) - kfrac(:, ip-1); dkc = MATMUL(cell%bmat, dk)
        kpath = kpath + SQRT(DOT_PRODUCT(dkc, dkc))
      END IF
      WRITE(iu,'(f12.6)', advance='no') kpath
      !> The nine anticommutators belong to this k-point, not to a band. Formed here, once:
      !> inside the band loop below each of them was built num_wann times out of the same two
      !> matrices. Component order stays xx xy xz yx yy yz zx zy zz.
      c = 0
      DO al = 1, 3
        DO be = 1, 3
          c = c + 1
          ! j = 1/2 (v_al O_be + O_be v_al)
          jmat(:, :, c) = 0.5 * ( MATMUL(v_interp(:, :, al, ip), o_interp(:, :, be, ip)) &
                                + MATMUL(o_interp(:, :, be, ip), v_interp(:, :, al, ip)) )
        END DO
      END DO
      DO m = 1, num_wann
        DO c = 1, 9
          jc(:, m) = MATMUL(jmat(:, :, c), cvec(:, m))
          jexp(c) = hartree_to_ev_const * REAL(DOT_PRODUCT(cvec(:, m), jc(:, m)))
        END DO
        WRITE(iu,'(2x,f14.8)', advance='no') hartree_to_ev_const*evals(m)
        DO c = 1, 9
          WRITE(iu,'(2x,f12.6)', advance='no') jexp(c)
        END DO
      END DO
      WRITE(iu,'(a)') ''
    END DO
    CLOSE(iu)
    WRITE(oUnit,'(a,i0,a)') 'wannierlib current interpolation: wrote '//TRIM(outfile)//' (', np, ' k-points)'
    CALL timestop('melem_interpolate_current')
  END SUBROUTINE melem_interpolate_current

END MODULE m_melem_interpolate_current
