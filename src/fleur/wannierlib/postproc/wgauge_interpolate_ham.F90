!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Wannier-gauge band-structure interpolation for the library-mode wannierization.
!>
!>  This is the Hamiltonian-specific driver on top of the operator-agnostic
!>  Fourier core m_wgauge_ft. It reproduces Wannier90's get_hr to build the
!>  Wannier-gauge Hamiltonian H_W(k) (num_wann x num_wann):
!>    1) slim eig to the outer window (lwindow from dis_win_min/max),
!>    2) eigval2(m,k) = sum_win eigval_opt * |u_opt(:,m)|^2,
!>    3) H_W(i,j,k) = sum_m eigval2(m,k) conjg(u(m,i)) u(m,j),
!>  delegates the Fourier interpolation H_W(k) -> H(R) -> H(k') to
!>  wgauge_ft_interpolate, diagonalizes H(k') and writes the bands.
!>
!>  To interpolate a different operator, write an analogous driver that builds
!>  its own operator matrix O_W(k) and reuses wgauge_ft_interpolate unchanged.
!>
!>  Runs only on the master rank (irank==0), where the W90 U matrices are complete.
MODULE m_wgauge_interpolate_ham
  USE m_juDFT
  USE m_constants, ONLY : oUnit, hartree_to_ev_const
  USE m_types_cell
  USE m_types_kpts
  USE m_types_wgauge_manifold, ONLY: t_wgauge_manifold
  USE m_wgauge_hamk, ONLY : wgauge_build_hamk
  USE m_wgauge_ft, ONLY : wgauge_ft_interpolate
  USE m_wgauge_interp_util, ONLY : wgauge_kpath, wgauge_zheev_workspace
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: wgauge_interpolate_ham
CONTAINS

  SUBROUTINE wgauge_interpolate_ham(this, cell, kpts, eig, u_matrix, u_opt, kfrac, out1, out2, irank)
    TYPE(t_wgauge_manifold), INTENT(IN) :: this
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    REAL,    INTENT(IN) :: eig(:, :)          ! (num_bands, nk)
    COMPLEX, INTENT(IN) :: u_matrix(:, :, :)  ! (num_wann, num_wann, nk)   MLWF gauge
    COMPLEX, INTENT(IN) :: u_opt(:, :, :)     ! (num_bands, num_wann, nk)  disentangled
    !> The domain's k-set and the names its files take, both decided by the caller: the
    !> k-points come from a named kPointList and the names from the exposure table plus the
    !> domain suffix. Unallocated off rank 0, which never reaches them.
    REAL, ALLOCATABLE, INTENT(IN) :: kfrac(:, :)          !> (3, np) fractional mesh
    CHARACTER(LEN=*), INTENT(IN) :: out1
    CHARACTER(LEN=*), INTENT(IN) :: out2
    INTEGER, INTENT(IN) :: irank              ! MPI rank (interpolate only on master)

    INTEGER :: num_wann, ip, np, iu, iuev, info, lwork
    REAL,    ALLOCATABLE :: kdist(:), evals(:), rwork(:)
    COMPLEX, ALLOCATABLE :: ham_k(:, :, :), H_interp(:, :, :), hk(:, :), work(:)

    IF (irank /= 0) RETURN                      ! only the master holds the full U(k)

    num_wann  = this%num_wann
    CALL timestart('wgauge_interpolate_ham')

    ! ---- the domain's k-set, already resolved by the caller ----
    np = SIZE(kfrac, 2)   ! the caller resolved the domain; there is nothing to skip

    CALL wgauge_kpath(cell, kfrac, kdist)   ! abscissa of the output, from the mesh just read

    ! ---- steps 1+2: eigval2 in the optimal (num_wann) subspace ----
    CALL wgauge_build_hamk(this, eig, u_matrix, u_opt, ham_k)

    ! ---- generic core: Fourier-interpolate H_W(k) to the fine path ----
    CALL wgauge_ft_interpolate(cell, ham_k, kpts, kfrac, H_interp)

    ! ---- diagonalize interpolated H(k') (complex Hermitian) ----
    ALLOCATE(evals(num_wann), hk(num_wann, num_wann))
    CALL wgauge_zheev_workspace('N', num_wann, work, rwork, lwork)

    OPEN(newunit=iu,   file=TRIM(out1)//'.dat', status='replace')
    OPEN(newunit=iuev, file=TRIM(out2)//'.dat', status='replace')
    WRITE(iu,  '(a)') '# kdist   E_1..E_numwann  (Htr)'
    WRITE(iuev,'(a)') '# kdist   E_1..E_numwann  (eV, absolute)'
    DO ip = 1, np
      hk = H_interp(:, :, ip)
      CALL zheev('N', 'U', num_wann, hk, num_wann, evals, work, lwork, rwork, info)
      IF (info /= 0) CALL juDFT_error('zheev failed', calledby='wgauge_interpolate_ham')
      WRITE(iu,  '(f12.6,*(2x,f14.8))') kdist(ip), evals(:)
      WRITE(iuev,'(f12.6,*(2x,f14.8))') kdist(ip), hartree_to_ev_const * evals(:)
    END DO
    CLOSE(iu); CLOSE(iuev)
    WRITE(oUnit,'(a,i0,a)') 'wannierlib interpolation: wrote '//TRIM(out1)//'.dat (', np, ' k-points)'
    CALL timestop('wgauge_interpolate_ham')
  END SUBROUTINE wgauge_interpolate_ham

END MODULE m_wgauge_interpolate_ham
