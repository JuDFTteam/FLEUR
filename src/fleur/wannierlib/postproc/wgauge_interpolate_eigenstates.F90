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
!>  Takes the eigenvectors the shared H-gauge step (m_wgauge_hamk) already produced for
!>  this domain and writes them. Master rank only.
MODULE m_wgauge_interpolate_eigenstates
  USE m_juDFT
  USE m_constants, ONLY : oUnit
  USE m_wgauge_hamk, ONLY : t_wgauge_hgauge
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: wgauge_interpolate_eigenstates
CONTAINS

  SUBROUTINE wgauge_interpolate_eigenstates(hg, kfrac, out1, irank)
    TYPE(t_wgauge_hgauge), INTENT(IN) :: hg   !> the domain's C(k'), already diagonalized
    !> The domain's k-set and the name its file takes, both decided by the caller: the
    !> k-points come from a named kPointList and the name from the exposure table plus the
    !> domain suffix. Unallocated off rank 0, which never reaches them.
    REAL, ALLOCATABLE, INTENT(IN) :: kfrac(:, :)          !> (3, np) fractional mesh
    CHARACTER(LEN=*), INTENT(IN) :: out1
    INTEGER, INTENT(IN) :: irank

    INTEGER :: num_wann, i, j, ip, np, iu

    IF (irank /= 0) RETURN                      ! only the master holds the full U(k)
    CALL timestart('wgauge_interpolate_eigenstates')
    num_wann = SIZE(hg%evals, 1)
    np = SIZE(hg%kdist)

    OPEN(newunit=iu, file=TRIM(out1)//'.dat', status='replace')
    WRITE(iu,'(a)') '# Wannier-Hamiltonian eigenstates C(k): H(k) C = C E, columns of C = band'
    WRITE(iu,'(a)') '# eigenvectors in the Wannier (tight-binding) basis (Hamiltonian-gauge U^(H)).'
    WRITE(iu,'(a,i0,a,i0)') '# num_wann = ', num_wann, ' ,  n_kpts = ', np
    WRITE(iu,'(a)') '# per k: "k <ip> <kx> <ky> <kz> <kdist>", then num_wann^2 lines "i j Re(C_ij) Im(C_ij)"'
    DO ip = 1, np
      WRITE(iu,'(a,i0,4(1x,f12.8))') 'k ', ip, kfrac(:, ip), hg%kdist(ip)
      DO j = 1, num_wann          ! column j = eigenstate j
        DO i = 1, num_wann        ! row i = Wannier index
          WRITE(iu,'(2i5,2(2x,es18.10))') i, j, REAL(hg%cvec(i, j, ip)), AIMAG(hg%cvec(i, j, ip))
        END DO
      END DO
    END DO
    CLOSE(iu)
    WRITE(oUnit,'(a,i0,a)') 'wannierlib eigenstates: wrote '//TRIM(out1)//'.dat (', np, ' k-points)'
    CALL timestop('wgauge_interpolate_eigenstates')
  END SUBROUTINE wgauge_interpolate_eigenstates

END MODULE m_wgauge_interpolate_eigenstates
