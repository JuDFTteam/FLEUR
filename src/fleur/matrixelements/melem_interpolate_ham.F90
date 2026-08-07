!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Wannier-gauge band-structure interpolation for the library-mode wannierization.
!>
!>  This is the Hamiltonian-specific driver on top of the operator-agnostic
!>  Fourier core m_melem_ft. It reproduces Wannier90's get_hr to build the
!>  Wannier-gauge Hamiltonian H_W(k) (num_wann x num_wann):
!>    1) slim eig to the outer window (lwindow from dis_win_min/max),
!>    2) eigval2(m,k) = sum_win eigval_opt * |u_opt(:,m)|^2,
!>    3) H_W(i,j,k) = sum_m eigval2(m,k) conjg(u(m,i)) u(m,j),
!>  delegates the Fourier interpolation H_W(k) -> H(R) -> H(k') to
!>  melem_ft_interpolate, diagonalizes H(k') and writes the bands.
!>
!>  To interpolate a different operator, write an analogous driver that builds
!>  its own operator matrix O_W(k) and reuses melem_ft_interpolate unchanged.
!>
!>  Runs only on the master rank (irank==0), where the W90 U matrices are complete.
MODULE m_melem_interpolate_ham
  USE m_juDFT
  USE m_constants, ONLY : oUnit, hartree_to_ev_const
  USE m_types_cell
  USE m_types_kpts
  USE m_types_melem_manifold, ONLY: t_melem_manifold
  USE m_melem_hamk, ONLY : melem_build_hamk
  USE m_melem_ft, ONLY : melem_ft_interpolate
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: melem_interpolate_ham
CONTAINS

  SUBROUTINE melem_interpolate_ham(this, cell, kpts, eig, u_matrix, u_opt, irank)
    TYPE(t_melem_manifold), INTENT(IN) :: this
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    REAL,    INTENT(IN) :: eig(:, :)          ! (num_bands, nk)
    COMPLEX, INTENT(IN) :: u_matrix(:, :, :)  ! (num_wann, num_wann, nk)   MLWF gauge
    COMPLEX, INTENT(IN) :: u_opt(:, :, :)     ! (num_bands, num_wann, nk)  disentangled
    INTEGER, INTENT(IN) :: irank              ! MPI rank (interpolate only on master)

    INTEGER :: num_wann, num_bands, nk, k, i, j, m, counter, ip, np, ios, iu, iuev, info, lwork
    LOGICAL :: lexist
    REAL    :: kpath, dk(3), dkc(3)
    REAL,    ALLOCATABLE :: kfrac(:, :), evals(:), rwork(:)
    COMPLEX, ALLOCATABLE :: ham_k(:, :, :), H_interp(:, :, :), hk(:, :), work(:)
    COMPLEX :: wq(1)

    IF (irank /= 0) RETURN                      ! only the master holds the full U(k)

    num_wann  = this%num_wann
    num_bands = this%num_bands
    nk        = kpts%nkptf
    CALL timestart('melem_interpolate_ham')

    ! ---- k-path from kpts_interpol ----
    INQUIRE(file='kpts_interpol', exist=lexist)
    IF (.NOT. lexist) THEN
      WRITE(oUnit,'(a)') 'wannierlib interpolation: no kpts_interpol file -> skipped'
      CALL timestop('melem_interpolate_ham'); RETURN
    END IF
    OPEN(newunit=iu, file='kpts_interpol', status='old')
    READ(iu, *, iostat=ios) np
    IF (ios /= 0 .OR. np <= 0) CALL juDFT_error('bad kpts_interpol header', calledby='melem_interpolate_ham')
    ALLOCATE(kfrac(3, np))
    DO ip = 1, np
      READ(iu, *, iostat=ios) kfrac(:, ip)
      IF (ios /= 0) CALL juDFT_error('bad kpts_interpol line', calledby='melem_interpolate_ham')
    END DO
    CLOSE(iu)

    ! ---- steps 1+2: eigval2 in the optimal (num_wann) subspace ----
    CALL melem_build_hamk(this, eig, u_matrix, u_opt, ham_k)

    ! ---- generic core: Fourier-interpolate H_W(k) to the fine path ----
    CALL melem_ft_interpolate(cell, ham_k, kpts, kfrac, H_interp)

    ! ---- diagonalize interpolated H(k') (complex Hermitian) ----
    ALLOCATE(evals(num_wann), hk(num_wann, num_wann), rwork(MAX(1, 3*num_wann-2)))
    CALL zheev('N', 'U', num_wann, hk, num_wann, evals, wq, -1, rwork, info)
    lwork = MAX(1, NINT(REAL(wq(1)))); ALLOCATE(work(lwork))

    OPEN(newunit=iu,   file='bands_wann_interpol.dat',    status='replace')
    OPEN(newunit=iuev, file='bands_wann_interpol_ev.dat', status='replace')
    WRITE(iu,  '(a)') '# kdist   E_1..E_numwann  (Htr)'
    WRITE(iuev,'(a)') '# kdist   E_1..E_numwann  (eV, absolute)'
    kpath = 0.0
    DO ip = 1, np
      hk = H_interp(:, :, ip)
      CALL zheev('N', 'U', num_wann, hk, num_wann, evals, work, lwork, rwork, info)
      IF (info /= 0) CALL juDFT_error('zheev failed', calledby='melem_interpolate_ham')
      IF (ip > 1) THEN
        dk = kfrac(:, ip) - kfrac(:, ip-1); dkc = MATMUL(cell%bmat, dk)
        kpath = kpath + SQRT(DOT_PRODUCT(dkc, dkc))
      END IF
      WRITE(iu,  '(f12.6,*(2x,f14.8))') kpath, evals(:)
      WRITE(iuev,'(f12.6,*(2x,f14.8))') kpath, hartree_to_ev_const * evals(:)
    END DO
    CLOSE(iu); CLOSE(iuev)
    WRITE(oUnit,'(a,i0,a)') 'wannierlib interpolation: wrote bands_wann_interpol.dat (', np, ' k-points)'
    CALL timestop('melem_interpolate_ham')
  END SUBROUTINE melem_interpolate_ham

END MODULE m_melem_interpolate_ham
