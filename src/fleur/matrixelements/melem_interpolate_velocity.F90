!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Band velocity v_alpha = dE_n/dk_alpha (framework Class C: no new Bloch matrix
!>  element -- built entirely from the interpolated Wannier-gauge Hamiltonian).
!>
!>  Pipeline:
!>    H_W(k)      : same Wannier-gauge Hamiltonian as the band/operator drivers
!>    v_W,alpha(k') = FT[ i R_cart(alpha) H_W ]     (m_melem_ft: velocity variant)
!>    H(k')       = FT[ H_W ] -> diag -> E_n(k'), C(k')
!>    <v_alpha>_n = [ C^dagger v_W,alpha(k') C ]_nn   (diagonal band velocity, exact)
!>
!>  The diagonal <n|v|n> = dE_n/dk needs no gauge (Berry-connection) correction, so
!>  it is exact here. Output bands_wann_velocity.dat: kdist, [ E_n(eV), vx, vy, vz ]
!>  per band, with v in eV*bohr (dE/dk). Master rank only.
MODULE m_melem_interpolate_velocity
  USE m_juDFT
  USE m_constants, ONLY : oUnit, hartree_to_ev_const
  USE m_types_cell
  USE m_types_kpts
  USE m_types_wannierlib
  USE m_melem_ft, ONLY : melem_ft_interpolate, melem_ft_velocity, melem_ft_rtok
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: melem_interpolate_velocity
CONTAINS

  SUBROUTINE melem_interpolate_velocity(this, cell, kpts, eig, u_matrix, u_opt, aw_r, irvec, ndegen, nrpts, irank)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    REAL,    INTENT(IN) :: eig(:, :)              ! (num_bands, nk)
    COMPLEX, INTENT(IN) :: u_matrix(:, :, :)      ! (num_wann, num_wann, nk)  MLWF gauge
    COMPLEX, INTENT(IN) :: u_opt(:, :, :)         ! (num_bands, num_wann, nk) disentangled
    COMPLEX, INTENT(IN) :: aw_r(:, :, :, :)       ! (num_wann,num_wann,nrpts,3) Berry connection A^(W)_a(R), reduced
    INTEGER, INTENT(IN) :: irvec(:, :), ndegen(:), nrpts   ! Wigner-Seitz R-mesh from the reduce (rank 0 only used)
    INTEGER, INTENT(IN) :: irank

    INTEGER :: num_wann, num_bands, nk, k, i, j, m, n, counter, ip, np, ios, iu, iuc, info, lwork, a
    INTEGER :: ax(3), ay(3)
    LOGICAL :: have_dis, lexist, l_berry
    REAL    :: kpath, dk(3), dkc(3), de
    REAL,    ALLOCATABLE :: eigval2(:, :), eigval_opt(:), kfrac(:, :), evals(:), rwork(:), vexp(:), omega(:, :)
    COMPLEX, ALLOCATABLE :: ham_k(:, :, :), H_interp(:, :, :), v_interp(:, :, :, :), A_interp(:, :, :, :)
    COMPLEX, ALLOCATABLE :: hk(:, :), work(:), cvec(:, :), vc(:, :, :), Hbar(:, :, :), Abar(:, :, :), vfull(:, :, :)
    COMPLEX :: wq(1), acc

    IF (irank /= 0) RETURN
    num_wann  = this%num_wann
    num_bands = this%num_bands
    nk        = kpts%nkptf
    have_dis  = (num_bands > num_wann)
    CALL timestart('melem_interpolate_velocity')

    INQUIRE(file='kpts_interpol', exist=lexist)
    IF (.NOT. lexist) THEN
      WRITE(oUnit,'(a)') 'wannierlib velocity interpolation: no kpts_interpol file -> skipped'
      CALL timestop('melem_interpolate_velocity'); RETURN
    END IF
    OPEN(newunit=iu, file='kpts_interpol', status='old')
    READ(iu, *, iostat=ios) np
    IF (ios /= 0 .OR. np <= 0) CALL juDFT_error('bad kpts_interpol header', calledby='melem_interpolate_velocity')
    ALLOCATE(kfrac(3, np))
    DO ip = 1, np
      READ(iu, *, iostat=ios) kfrac(:, ip)
      IF (ios /= 0) CALL juDFT_error('bad kpts_interpol line', calledby='melem_interpolate_velocity')
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

    ! ---- interpolate H (for eigenvectors) and v = dH/dk (velocity variant of the core) ----
    CALL melem_ft_interpolate(cell, ham_k, kpts, kfrac, H_interp)
    CALL melem_ft_velocity(cell, ham_k, kpts, kfrac, v_interp)   ! (nw,nw,3,np)  = dH/dk

    ! ---- interband part: R -> k' of the reduced Wannier Berry connection A^(W)_a(R) -> A^(W)_a(k') ----
    l_berry = (nrpts > 0 .AND. SIZE(aw_r, 1) == num_wann .AND. SIZE(aw_r, 4) == 3)
    IF (l_berry) THEN
      ALLOCATE(A_interp(num_wann, num_wann, 3, np))
      BLOCK
        COMPLEX, ALLOCATABLE :: a_one(:, :, :)
        DO a = 1, 3
          CALL melem_ft_rtok(aw_r(:, :, :, a), irvec, ndegen, nrpts, kfrac, a_one)
          A_interp(:, :, a, :) = a_one
        END DO
      END BLOCK
    END IF

    ! ---- diagonalize H(k'), project the diagonal band velocity, write ----
    ALLOCATE(evals(num_wann), hk(num_wann, num_wann), cvec(num_wann, num_wann), &
             vc(num_wann, num_wann, 3), vexp(3), rwork(MAX(1, 3*num_wann-2)))
    IF (l_berry) ALLOCATE(Hbar(num_wann, num_wann, 3), Abar(num_wann, num_wann, 3), &
                          vfull(num_wann, num_wann, 3), omega(3, num_wann))
    ax = (/ 2, 3, 1 /)   ! Omega_gamma = eps_{gamma,alpha,beta}: (Ox<-yz, Oy<-zx, Oz<-xy)
    ay = (/ 3, 1, 2 /)
    CALL zheev('V', 'U', num_wann, hk, num_wann, evals, wq, -1, rwork, info)
    lwork = MAX(1, NINT(REAL(wq(1)))); ALLOCATE(work(lwork))

    OPEN(newunit=iu, file='bands_wann_velocity.dat', status='replace')
    WRITE(iu,'(a)') '# kdist   [ E_n(eV)  vx vy vz (eV*bohr, dE/dk) ] for n=1..num_wann'
    IF (l_berry) THEN
      OPEN(newunit=iuc, file='bands_wann_berrycurv.dat', status='replace')
      WRITE(iuc,'(a)') '# kdist   [ E_n(eV)  Omega_x Omega_y Omega_z (bohr^2) ] for n=1..num_wann'
    END IF
    kpath = 0.0
    DO ip = 1, np
      hk = H_interp(:, :, ip)
      CALL zheev('V', 'U', num_wann, hk, num_wann, evals, work, lwork, rwork, info)
      IF (info /= 0) CALL juDFT_error('zheev failed', calledby='melem_interpolate_velocity')
      cvec = hk
      IF (ip > 1) THEN
        dk = kfrac(:, ip) - kfrac(:, ip-1); dkc = MATMUL(cell%bmat, dk)
        kpath = kpath + SQRT(DOT_PRODUCT(dkc, dkc))
      END IF
      DO a = 1, 3
        vc(:, :, a) = MATMUL(v_interp(:, :, a, ip), cvec)         ! v_interp_a . C
      END DO
      ! ---- diagonal band velocity <n|v|n> = dE_n/dk (unchanged, byte-identical) ----
      WRITE(iu,'(f12.6)', advance='no') kpath
      DO m = 1, num_wann
        DO a = 1, 3
          vexp(a) = hartree_to_ev_const * REAL(DOT_PRODUCT(cvec(:, m), vc(:, m, a)))
        END DO
        WRITE(iu,'(2x,f14.8)', advance='no') hartree_to_ev_const*evals(m)
        DO a = 1, 3
          WRITE(iu,'(2x,f12.6)', advance='no') vexp(a)
        END DO
      END DO
      WRITE(iu,'(a)') ''

      ! ---- interband velocity matrix + Berry curvature (a.u.: v in Ha*bohr, Omega in bohr^2) ----
      IF (l_berry) THEN
        DO a = 1, 3
          Hbar(:, :, a) = MATMUL(CONJG(TRANSPOSE(cvec)), vc(:, :, a))                       ! C^dag dH/dk_a C
          Abar(:, :, a) = MATMUL(CONJG(TRANSPOSE(cvec)), MATMUL(A_interp(:, :, a, ip), cvec)) ! C^dag A^(W)_a C
        END DO
        DO a = 1, 3
          DO n = 1, num_wann
            DO m = 1, num_wann
              IF (n == m) THEN
                vfull(n, m, a) = Hbar(n, m, a)
              ELSE
                vfull(n, m, a) = Hbar(n, m, a) + CMPLX(0.0, evals(n) - evals(m)) * Abar(n, m, a)
              END IF
            END DO
          END DO
        END DO
        DO n = 1, num_wann
          DO a = 1, 3               ! a = gamma (x,y,z); curvature from the (ax,ay) plane
            acc = CMPLX(0.0, 0.0)
            DO m = 1, num_wann
              IF (m == n) CYCLE
              de = evals(n) - evals(m)
              IF (ABS(de) < 1.0e-8) CYCLE
              acc = acc + vfull(n, m, ax(a)) * vfull(m, n, ay(a)) / CMPLX(de*de, 0.0)
            END DO
            omega(a, n) = -2.0 * AIMAG(acc)
          END DO
        END DO
        WRITE(iuc,'(f12.6)', advance='no') kpath
        DO m = 1, num_wann
          WRITE(iuc,'(2x,f14.8)', advance='no') hartree_to_ev_const*evals(m)
          DO a = 1, 3
            WRITE(iuc,'(2x,es14.6)', advance='no') omega(a, m)
          END DO
        END DO
        WRITE(iuc,'(a)') ''
      END IF
    END DO
    CLOSE(iu)
    IF (l_berry) CLOSE(iuc)
    WRITE(oUnit,'(a,i0,a)') 'wannierlib velocity interpolation: wrote bands_wann_velocity.dat (', np, ' k-points)'
    CALL timestop('melem_interpolate_velocity')
  END SUBROUTINE melem_interpolate_velocity

END MODULE m_melem_interpolate_velocity
