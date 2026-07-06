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
!>  operator A(R) (Class B), added in a later stage. Master rank only.
!>  Output: kdist, [ E_n(eV), j_xx, j_xy, j_xz, j_yx, ... j_zz (eV*bohr) ] per band.
MODULE m_wannierlib_interpolate_current
  USE m_juDFT
  USE m_constants, ONLY : oUnit, hartree_to_ev_const
  USE m_types_cell
  USE m_types_kpts
  USE m_types_wannierlib
  USE m_wannierlib_ft, ONLY : wannierlib_ft_interpolate, wannierlib_ft_velocity
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: wannierlib_interpolate_current
CONTAINS

  SUBROUTINE wannierlib_interpolate_current(this, cell, kpts, eig, u_matrix, u_opt, o0, outfile, irank)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    REAL,    INTENT(IN) :: eig(:, :)              ! (num_bands, nk)
    COMPLEX, INTENT(IN) :: u_matrix(:, :, :)      ! (num_wann, num_wann, nk)  MLWF gauge
    COMPLEX, INTENT(IN) :: u_opt(:, :, :)         ! (num_bands, num_wann, nk) disentangled
    COMPLEX, INTENT(IN) :: o0(:, :, :, :)         ! (num_bands, num_bands, 3, nk) Bloch spin OR orbital
    CHARACTER(LEN=*), INTENT(IN) :: outfile
    INTEGER, INTENT(IN) :: irank

    INTEGER :: num_wann, num_bands, nk, k, i, j, m, counter, ip, np, ios, iu, info, lwork, al, be, c
    LOGICAL :: have_dis, lexist
    REAL    :: kpath, dk(3), dkc(3)
    REAL,    ALLOCATABLE :: eigval2(:, :), eigval_opt(:), kfrac(:, :), evals(:), rwork(:), jexp(:)
    COMPLEX, ALLOCATABLE :: ham_k(:, :, :), ow_k(:, :, :, :), H_interp(:, :, :), o_interp(:, :, :, :)
    COMPLEX, ALLOCATABLE :: v_interp(:, :, :, :), hk(:, :), work(:), vmat(:, :), tmp(:, :), cvec(:, :)
    COMPLEX, ALLOCATABLE :: jmat(:, :), jc(:, :)
    COMPLEX :: wq(1)

    IF (irank /= 0) RETURN
    num_wann  = this%num_wann
    num_bands = this%num_bands
    nk        = kpts%nkptf
    have_dis  = (num_bands > num_wann)
    CALL timestart('wannierlib_interpolate_current')

    INQUIRE(file='kpts_interpol', exist=lexist)
    IF (.NOT. lexist) THEN
      WRITE(oUnit,'(a)') 'wannierlib current interpolation: no kpts_interpol file -> skipped'
      CALL timestop('wannierlib_interpolate_current'); RETURN
    END IF
    OPEN(newunit=iu, file='kpts_interpol', status='old')
    READ(iu, *, iostat=ios) np
    IF (ios /= 0 .OR. np <= 0) CALL juDFT_error('bad kpts_interpol header', calledby='wannierlib_interpolate_current')
    ALLOCATE(kfrac(3, np))
    DO ip = 1, np
      READ(iu, *, iostat=ios) kfrac(:, ip)
      IF (ios /= 0) CALL juDFT_error('bad kpts_interpol line', calledby='wannierlib_interpolate_current')
    END DO
    CLOSE(iu)

    ! ---- H_W(k) via eigval2 (same construction as the validated band driver) ----
    ALLOCATE(eigval2(num_wann, nk), source=0.0)
    IF (have_dis) THEN
      ALLOCATE(eigval_opt(num_bands))
      DO k = 1, nk
        counter = 0; eigval_opt = 0.0
        DO j = 1, num_bands
          IF (eig(j, k) >= this%dis_win_min .AND. eig(j, k) <= this%dis_win_max) THEN
            counter = counter + 1; eigval_opt(counter) = eig(j, k)
          END IF
        END DO
        DO m = 1, num_wann
          DO i = 1, counter
            eigval2(m, k) = eigval2(m, k) + eigval_opt(i) * ABS(u_opt(i, m, k))**2
          END DO
        END DO
      END DO
      DEALLOCATE(eigval_opt)
    ELSE
      eigval2(1:num_wann, :) = eig(1:num_wann, :)
    END IF

    ALLOCATE(ham_k(num_wann, num_wann, nk), source=CMPLX(0.0, 0.0))
    DO k = 1, nk
      DO j = 1, num_wann
        DO i = 1, num_wann
          DO m = 1, num_wann
            ham_k(i, j, k) = ham_k(i, j, k) + eigval2(m, k) * CONJG(u_matrix(m, i, k)) * u_matrix(m, j, k)
          END DO
        END DO
      END DO
    END DO

    ! ---- O_W,beta(k) = V^dagger O0_beta V,   V = u_opt . u_matrix ----
    ALLOCATE(vmat(num_bands, num_wann), tmp(num_bands, num_wann))
    ALLOCATE(ow_k(num_wann, num_wann, 3, nk), source=CMPLX(0.0, 0.0))
    DO k = 1, nk
      vmat = MATMUL(u_opt(:, :, k), u_matrix(:, :, k))
      DO be = 1, 3
        tmp = MATMUL(o0(:, :, be, k), vmat)
        ow_k(:, :, be, k) = MATMUL(CONJG(TRANSPOSE(vmat)), tmp)
      END DO
    END DO

    ! ---- interpolate H (eigenvectors), v = dH/dk (3), and O_beta (3) ----
    CALL wannierlib_ft_interpolate(cell, ham_k, kpts, kfrac, H_interp)
    CALL wannierlib_ft_velocity(cell, ham_k, kpts, kfrac, v_interp)     ! (nw,nw,3,np)
    ALLOCATE(o_interp(num_wann, num_wann, 3, np))
    BLOCK
      COMPLEX, ALLOCATABLE :: o_one(:, :, :)
      DO be = 1, 3
        CALL wannierlib_ft_interpolate(cell, ow_k(:, :, be, :), kpts, kfrac, o_one)
        o_interp(:, :, be, :) = o_one
      END DO
    END BLOCK

    ! ---- diagonalize H(k'); j_{al,be} = 1/2 { v_al, O_be }; project diagonal; write 9 comps ----
    ALLOCATE(evals(num_wann), hk(num_wann, num_wann), cvec(num_wann, num_wann), &
             jmat(num_wann, num_wann), jc(num_wann, num_wann), jexp(9), rwork(MAX(1, 3*num_wann-2)))
    CALL zheev('V', 'U', num_wann, hk, num_wann, evals, wq, -1, rwork, info)
    lwork = MAX(1, NINT(REAL(wq(1)))); ALLOCATE(work(lwork))

    OPEN(newunit=iu, file=outfile, status='replace')
    WRITE(iu,'(a)') '# kdist   [ E_n(eV)  j_ab (a,b=x,y,z: xx xy xz yx yy yz zx zy zz, eV*bohr) ] for n=1..num_wann'
    kpath = 0.0
    DO ip = 1, np
      hk = H_interp(:, :, ip)
      CALL zheev('V', 'U', num_wann, hk, num_wann, evals, work, lwork, rwork, info)
      IF (info /= 0) CALL juDFT_error('zheev failed', calledby='wannierlib_interpolate_current')
      cvec = hk
      IF (ip > 1) THEN
        dk = kfrac(:, ip) - kfrac(:, ip-1); dkc = MATMUL(cell%bmat, dk)
        kpath = kpath + SQRT(DOT_PRODUCT(dkc, dkc))
      END IF
      WRITE(iu,'(f12.6)', advance='no') kpath
      DO m = 1, num_wann
        c = 0
        DO al = 1, 3
          DO be = 1, 3
            c = c + 1
            ! j = 1/2 (v_al O_be + O_be v_al)
            jmat = 0.5 * ( MATMUL(v_interp(:, :, al, ip), o_interp(:, :, be, ip)) &
                         + MATMUL(o_interp(:, :, be, ip), v_interp(:, :, al, ip)) )
            jc(:, m) = MATMUL(jmat, cvec(:, m))
            jexp(c) = hartree_to_ev_const * REAL(DOT_PRODUCT(cvec(:, m), jc(:, m)))
          END DO
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
    CALL timestop('wannierlib_interpolate_current')
  END SUBROUTINE wannierlib_interpolate_current

END MODULE m_wannierlib_interpolate_current
