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
  USE m_melem_hamk, ONLY : melem_build_hamk
  USE m_melem_domains, ONLY : melem_read_kset
  USE m_melem_ft, ONLY : melem_ft_interpolate
  USE m_melem_interp_util, ONLY : melem_kpath, melem_zheev_workspace
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

    INTEGER :: num_wann, i, j, ip, np, iu, info, lwork
    REAL,    ALLOCATABLE :: kfrac(:, :), kdist(:), evals(:), rwork(:)
    COMPLEX, ALLOCATABLE :: ham_k(:, :, :), H_interp(:, :, :), cvec(:, :), work(:)

    IF (irank /= 0) RETURN                      ! only the master holds the full U(k)

    num_wann  = this%num_wann
    CALL timestart('melem_interpolate_eigenstates')

    ! ---- k-mesh from kpts_interpol ----
    IF (.NOT. melem_read_kset(kfrac, np, 'melem_interpolate_eigenstates')) THEN
      WRITE(oUnit,'(a)') 'wannierlib eigenstates: no kpts_interpol file -> skipped'
      CALL timestop('melem_interpolate_eigenstates'); RETURN
    END IF

    CALL melem_kpath(cell, kfrac, kdist)   ! abscissa of the output, from the mesh just read

    ! ---- H_W(k) via eigval2 (identical construction to the band driver) ----
    CALL melem_build_hamk(this, eig, u_matrix, u_opt, ham_k)

    ! ---- Fourier-interpolate H_W(k) to the fine mesh (shared core) ----
    CALL melem_ft_interpolate(cell, ham_k, kpts, kfrac, H_interp)

    ! ---- diagonalize H(k') WITH eigenvectors; write C(k') ----
    ALLOCATE(evals(num_wann), cvec(num_wann, num_wann))
    CALL melem_zheev_workspace('V', num_wann, work, rwork, lwork)

    OPEN(newunit=iu, file='bands_wann_eigenstates.dat', status='replace')
    WRITE(iu,'(a)') '# Wannier-Hamiltonian eigenstates C(k): H(k) C = C E, columns of C = band'
    WRITE(iu,'(a)') '# eigenvectors in the Wannier (tight-binding) basis (Hamiltonian-gauge U^(H)).'
    WRITE(iu,'(a,i0,a,i0)') '# num_wann = ', num_wann, ' ,  n_kpts = ', np
    WRITE(iu,'(a)') '# per k: "k <ip> <kx> <ky> <kz> <kdist>", then num_wann^2 lines "i j Re(C_ij) Im(C_ij)"'
    DO ip = 1, np
      cvec = H_interp(:, :, ip)
      CALL zheev('V', 'U', num_wann, cvec, num_wann, evals, work, lwork, rwork, info)
      IF (info /= 0) CALL juDFT_error('zheev failed', calledby='melem_interpolate_eigenstates')
      WRITE(iu,'(a,i0,4(1x,f12.8))') 'k ', ip, kfrac(:, ip), kdist(ip)
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
