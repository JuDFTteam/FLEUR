!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Export of the Wannier-Hamiltonian eigenstates C(k') along the interpolation
!>  mesh. C(k') is the unitary that diagonalizes the interpolated Hamiltonian,
!>      H(k') C(k') = C(k') E(k') ,
!>  so its columns are the band eigenvectors expressed in the Wannier basis -- the
!>  Wannier-gauge -> Hamiltonian-gauge rotation U^(H) of Wang-Yates-Souza-Vanderbilt.
!>
!>  Unlike the band operators this is NOT a projected expectation value <O>_n but a
!>  full num_wann x num_wann matrix per k', so it is written as a matrix (like the
!>  tight-binding exports). It lets any operator be reconstructed in the band basis
!>  in post-processing (e.g. <O>_n = [C^dagger O(k') C]_nn) without re-diagonalizing.
!>
!>  Builds H_W(k) exactly as the band driver (m_melem_interpolate_ham),
!>  Fourier-interpolates to H(k') via the shared core, diagonalizes WITH eigenvectors,
!>  and writes C(k'). Master rank only.
MODULE m_melem_interpolate_eigenstates
  USE m_juDFT
  USE m_constants, ONLY : oUnit
  USE m_types_cell
  USE m_types_kpts
  USE m_types_melem_manifold, ONLY: t_melem_manifold
  USE m_melem_ft, ONLY : melem_ft_interpolate
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: melem_interpolate_eigenstates
CONTAINS

  SUBROUTINE melem_interpolate_eigenstates(this, cell, kpts, eig, u_matrix, u_opt, irank)
    TYPE(t_melem_manifold), INTENT(IN) :: this
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    REAL,    INTENT(IN) :: eig(:, :)          ! (num_bands, nk)
    COMPLEX, INTENT(IN) :: u_matrix(:, :, :)  ! (num_wann, num_wann, nk)   MLWF gauge
    COMPLEX, INTENT(IN) :: u_opt(:, :, :)     ! (num_bands, num_wann, nk)  disentangled
    INTEGER, INTENT(IN) :: irank

    INTEGER :: num_wann, num_bands, nk, k, i, j, m, counter, ip, np, ios, iu, info, lwork
    LOGICAL :: have_dis, lexist
    REAL    :: kpath, dk(3), dkc(3)
    REAL,    ALLOCATABLE :: eigval2(:, :), eigval_opt(:), kfrac(:, :), evals(:), rwork(:)
    COMPLEX, ALLOCATABLE :: ham_k(:, :, :), H_interp(:, :, :), cvec(:, :), work(:)
    COMPLEX :: wq(1)

    IF (irank /= 0) RETURN                      ! only the master holds the full U(k)

    num_wann  = this%num_wann
    num_bands = this%num_bands
    nk        = kpts%nkptf
    have_dis  = (num_bands > num_wann)
    CALL timestart('melem_interpolate_eigenstates')

    ! ---- k-mesh from kpts_interpol ----
    INQUIRE(file='kpts_interpol', exist=lexist)
    IF (.NOT. lexist) THEN
      WRITE(oUnit,'(a)') 'wannierlib eigenstates: no kpts_interpol file -> skipped'
      CALL timestop('melem_interpolate_eigenstates'); RETURN
    END IF
    OPEN(newunit=iu, file='kpts_interpol', status='old')
    READ(iu, *, iostat=ios) np
    IF (ios /= 0 .OR. np <= 0) CALL juDFT_error('bad kpts_interpol header', calledby='melem_interpolate_eigenstates')
    ALLOCATE(kfrac(3, np))
    DO ip = 1, np
      READ(iu, *, iostat=ios) kfrac(:, ip)
      IF (ios /= 0) CALL juDFT_error('bad kpts_interpol line', calledby='melem_interpolate_eigenstates')
    END DO
    CLOSE(iu)

    ! ---- H_W(k) via eigval2 (identical construction to the band driver) ----
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

    ! ---- Fourier-interpolate H_W(k) to the fine mesh (shared core) ----
    CALL melem_ft_interpolate(cell, ham_k, kpts, kfrac, H_interp)

    ! ---- diagonalize H(k') WITH eigenvectors; write C(k') ----
    ALLOCATE(evals(num_wann), cvec(num_wann, num_wann), rwork(MAX(1, 3*num_wann-2)))
    CALL zheev('V', 'U', num_wann, cvec, num_wann, evals, wq, -1, rwork, info)
    lwork = MAX(1, NINT(REAL(wq(1)))); ALLOCATE(work(lwork))

    OPEN(newunit=iu, file='bands_wann_eigenstates.dat', status='replace')
    WRITE(iu,'(a)') '# Wannier-Hamiltonian eigenstates C(k): H(k) C = C E, columns of C = band'
    WRITE(iu,'(a)') '# eigenvectors in the Wannier (tight-binding) basis (Hamiltonian-gauge U^(H)).'
    WRITE(iu,'(a,i0,a,i0)') '# num_wann = ', num_wann, ' ,  n_kpts = ', np
    WRITE(iu,'(a)') '# per k: "k <ip> <kx> <ky> <kz> <kdist>", then num_wann^2 lines "i j Re(C_ij) Im(C_ij)"'
    kpath = 0.0
    DO ip = 1, np
      cvec = H_interp(:, :, ip)
      CALL zheev('V', 'U', num_wann, cvec, num_wann, evals, work, lwork, rwork, info)
      IF (info /= 0) CALL juDFT_error('zheev failed', calledby='melem_interpolate_eigenstates')
      IF (ip > 1) THEN
        dk = kfrac(:, ip) - kfrac(:, ip-1); dkc = MATMUL(cell%bmat, dk)
        kpath = kpath + SQRT(DOT_PRODUCT(dkc, dkc))
      END IF
      WRITE(iu,'(a,i0,4(1x,f12.8))') 'k ', ip, kfrac(:, ip), kpath
      DO j = 1, num_wann          ! column j = eigenstate j
        DO i = 1, num_wann        ! row i = Wannier index
          WRITE(iu,'(2i5,2(2x,es18.10))') i, j, REAL(cvec(i, j)), AIMAG(cvec(i, j))
        END DO
      END DO
    END DO
    CLOSE(iu)
    WRITE(oUnit,'(a,i0,a)') 'wannierlib eigenstates: wrote bands_wann_eigenstates.dat (', np, ' k-points)'
    CALL timestop('melem_interpolate_eigenstates')
  END SUBROUTINE melem_interpolate_eigenstates

END MODULE m_melem_interpolate_eigenstates
