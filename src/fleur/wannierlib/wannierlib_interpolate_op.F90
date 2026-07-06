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
!>    O_alpha(k')  = FT[ O_W,alpha ]            (shared core m_wannierlib_ft)
!>    H(k')        = FT[ H_W ] -> diag -> E_n(k'), C(k')
!>    <O_alpha>_n(k') = [ C^dagger O_alpha(k') C ]_nn
!>  and writes 'outfile': kdist, [ E_n(eV), <O_1>_n, ..., <O_ncomp>_n ] per band.
!>
!>  A new operator only supplies its O0(k) (a provider) and calls this with the
!>  right ncomp/outfile -- steps above are never rewritten. Master rank only.
MODULE m_wannierlib_interpolate_op
  USE m_juDFT
  USE m_constants, ONLY : oUnit, hartree_to_ev_const
  USE m_types_cell
  USE m_types_kpts
  USE m_types_wannierlib
  USE m_wannierlib_ft, ONLY : wannierlib_ft_interpolate
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: wannierlib_interpolate_operator
CONTAINS

  SUBROUTINE wannierlib_interpolate_operator(this, cell, kpts, eig, u_matrix, u_opt, o0, ncomp, outfile, irank)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    REAL,    INTENT(IN) :: eig(:, :)              ! (num_bands, nk)
    COMPLEX, INTENT(IN) :: u_matrix(:, :, :)      ! (num_wann, num_wann, nk)  MLWF gauge
    COMPLEX, INTENT(IN) :: u_opt(:, :, :)         ! (num_bands, num_wann, nk) disentangled
    COMPLEX, INTENT(IN) :: o0(:, :, :, :)         ! (num_bands, num_bands, ncomp, nk) Bloch operator
    INTEGER, INTENT(IN) :: ncomp
    CHARACTER(LEN=*), INTENT(IN) :: outfile
    INTEGER, INTENT(IN) :: irank

    INTEGER :: num_wann, num_bands, nk, k, i, j, m, counter, ip, np, ios, iu, info, lwork, a
    LOGICAL :: have_dis, lexist
    REAL    :: kpath, dk(3), dkc(3)
    REAL,    ALLOCATABLE :: eigval2(:, :), eigval_opt(:), kfrac(:, :), evals(:), rwork(:), oexp(:)
    COMPLEX, ALLOCATABLE :: ham_k(:, :, :), ow_k(:, :, :, :), H_interp(:, :, :), o_interp(:, :, :, :)
    COMPLEX, ALLOCATABLE :: hk(:, :), work(:), vmat(:, :), tmp(:, :), cvec(:, :), oc(:, :, :)
    COMPLEX :: wq(1)

    IF (irank /= 0) RETURN

    num_wann  = this%num_wann
    num_bands = this%num_bands
    nk        = kpts%nkptf
    have_dis  = (num_bands > num_wann)
    CALL timestart('wannierlib_interpolate_operator')

    INQUIRE(file='kpts_interpol', exist=lexist)
    IF (.NOT. lexist) THEN
      WRITE(oUnit,'(a)') 'wannierlib operator interpolation: no kpts_interpol file -> skipped'
      CALL timestop('wannierlib_interpolate_operator'); RETURN
    END IF
    OPEN(newunit=iu, file='kpts_interpol', status='old')
    READ(iu, *, iostat=ios) np
    IF (ios /= 0 .OR. np <= 0) CALL juDFT_error('bad kpts_interpol header', calledby='wannierlib_interpolate_operator')
    ALLOCATE(kfrac(3, np))
    DO ip = 1, np
      READ(iu, *, iostat=ios) kfrac(:, ip)
      IF (ios /= 0) CALL juDFT_error('bad kpts_interpol line', calledby='wannierlib_interpolate_operator')
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

    ! ---- O_W,alpha(k) = V^dagger O0_alpha V,   V = u_opt . u_matrix ----
    ALLOCATE(vmat(num_bands, num_wann), tmp(num_bands, num_wann))
    ALLOCATE(ow_k(num_wann, num_wann, ncomp, nk), source=CMPLX(0.0, 0.0))
    DO k = 1, nk
      vmat = MATMUL(u_opt(:, :, k), u_matrix(:, :, k))
      DO a = 1, ncomp
        tmp = MATMUL(o0(:, :, a, k), vmat)
        ow_k(:, :, a, k) = MATMUL(CONJG(TRANSPOSE(vmat)), tmp)
      END DO
    END DO

    ! ---- interpolate H and the ncomp operator components (shared core) ----
    CALL wannierlib_ft_interpolate(cell, ham_k, kpts, kfrac, H_interp)
    ALLOCATE(o_interp(num_wann, num_wann, ncomp, np))
    BLOCK
      COMPLEX, ALLOCATABLE :: o_one(:, :, :)
      DO a = 1, ncomp
        CALL wannierlib_ft_interpolate(cell, ow_k(:, :, a, :), kpts, kfrac, o_one)
        o_interp(:, :, a, :) = o_one
      END DO
    END BLOCK

    ! ---- diagonalize H(k') with eigenvectors, project operator, write ----
    ALLOCATE(evals(num_wann), hk(num_wann, num_wann), cvec(num_wann, num_wann), &
             oc(num_wann, num_wann, ncomp), oexp(ncomp), rwork(MAX(1, 3*num_wann-2)))
    CALL zheev('V', 'U', num_wann, hk, num_wann, evals, wq, -1, rwork, info)
    lwork = MAX(1, NINT(REAL(wq(1)))); ALLOCATE(work(lwork))

    OPEN(newunit=iu, file=outfile, status='replace')
    WRITE(iu,'(a,i0,a)') '# kdist   [ E_n(eV)  <O_1>_n .. <O_', ncomp, '>_n ] for n=1..num_wann'
    kpath = 0.0
    DO ip = 1, np
      hk = H_interp(:, :, ip)
      CALL zheev('V', 'U', num_wann, hk, num_wann, evals, work, lwork, rwork, info)
      IF (info /= 0) CALL juDFT_error('zheev failed', calledby='wannierlib_interpolate_operator')
      cvec = hk
      IF (ip > 1) THEN
        dk = kfrac(:, ip) - kfrac(:, ip-1); dkc = MATMUL(cell%bmat, dk)
        kpath = kpath + SQRT(DOT_PRODUCT(dkc, dkc))
      END IF
      DO a = 1, ncomp
        oc(:, :, a) = MATMUL(o_interp(:, :, a, ip), cvec)
      END DO
      WRITE(iu,'(f12.6)', advance='no') kpath
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
    WRITE(oUnit,'(a,i0,a)') 'wannierlib operator interpolation: wrote '//TRIM(outfile)//' (', np, ' k-points)'
    CALL timestop('wannierlib_interpolate_operator')
  END SUBROUTINE wannierlib_interpolate_operator

END MODULE m_wannierlib_interpolate_op
