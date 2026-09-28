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
!>    v_W,alpha(k') = FT[ i R_cart(alpha) H_W ]     (m_wgauge_ft: velocity variant)
!>    H(k')       = FT[ H_W ] -> diag -> E_n(k'), C(k')
!>    <v_alpha>_n = [ C^dagger v_W,alpha(k') C ]_nn   (diagonal band velocity, exact)
!>
!>  The diagonal <n|v|n> = dE_n/dk needs no gauge (Berry-connection) correction, so
!>  it is exact here. Output bands_wann_velocity.dat: kdist, [ E_n(eV), vx, vy, vz ]
!>  per band, with v in eV*bohr (dE/dk). Master rank only.
MODULE m_wgauge_interpolate_velocity
  USE m_juDFT
  USE m_wgauge_bands_io, ONLY: wgauge_bands_open, wgauge_bands_row
  USE m_constants, ONLY : oUnit, hartree_to_ev_const
  USE m_types_cell
  USE m_wgauge_hamk, ONLY : t_wgauge_hgauge
  USE m_wgauge_ft, ONLY : wgauge_ft_rtok_velocity, wgauge_ft_rtok
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: wgauge_interpolate_velocity
CONTAINS

  SUBROUTINE wgauge_interpolate_velocity(cell, ham_r, h_irvec, h_ndegen, h_nrpts, &
                                        aw_r, irvec, ndegen, nrpts, kfrac, hg, out1, out2, irank)
    TYPE(t_cell), INTENT(IN) :: cell
    !> H(R) and its Wigner-Seitz mesh, transformed once by the caller: the derivative is
    !> taken off the same real-space Hamiltonian the bands were interpolated from, so the
    !> velocity and the energies it is written next to cannot drift apart.
    COMPLEX, ALLOCATABLE, INTENT(IN) :: ham_r(:, :, :)     ! (num_wann,num_wann,h_nrpts)
    INTEGER, ALLOCATABLE, INTENT(IN) :: h_irvec(:, :), h_ndegen(:)
    INTEGER, INTENT(IN) :: h_nrpts
    COMPLEX, INTENT(IN) :: aw_r(:, :, :, :)       ! (num_wann,num_wann,nrpts,3) Berry connection A^(W)_a(R), reduced
    INTEGER, INTENT(IN) :: irvec(:, :), ndegen(:), nrpts   ! Wigner-Seitz R-mesh from the reduce (rank 0 only used)
    !> The domain's k-set and the names its files take, both decided by the caller: the
    !> k-points come from a named kPointList and the names from the exposure table plus the
    !> domain suffix. Unallocated off rank 0, which never reaches them.
    REAL, ALLOCATABLE, INTENT(IN) :: kfrac(:, :)          !> (3, np) fractional mesh
    TYPE(t_wgauge_hgauge), INTENT(IN) :: hg   !> the domain's bands, already diagonalized
    CHARACTER(LEN=*), INTENT(IN) :: out1
    CHARACTER(LEN=*), INTENT(IN) :: out2
    INTEGER, INTENT(IN) :: irank

    INTEGER :: num_wann, m, n, ip, np, iu, iuc, a
    INTEGER :: ax(3), ay(3)
    LOGICAL :: l_berry
    REAL    :: de
    REAL,    ALLOCATABLE :: vrow(:, :), omega(:, :)
    COMPLEX, ALLOCATABLE :: v_interp(:, :, :, :), A_interp(:, :, :, :)
    COMPLEX, ALLOCATABLE :: cvec(:, :), vc(:, :, :), Hbar(:, :, :), Abar(:, :, :), vfull(:, :, :)
    COMPLEX :: acc

    IF (irank /= 0) RETURN
    num_wann = SIZE(hg%evals, 1)
    CALL timestart('wgauge_interpolate_velocity')

    np = SIZE(kfrac, 2)   ! the caller resolved the domain; there is nothing to skip

    ! ---- v = dH/dk off the H(R) the caller already transformed ----
    CALL wgauge_ft_rtok_velocity(cell, ham_r, h_irvec, h_ndegen, h_nrpts, kfrac, v_interp)

    ! ---- interband part: R -> k' of the reduced Wannier Berry connection A^(W)_a(R) -> A^(W)_a(k') ----
    l_berry = (nrpts > 0 .AND. SIZE(aw_r, 1) == num_wann .AND. SIZE(aw_r, 4) == 3)
    IF (l_berry) THEN
      ALLOCATE(A_interp(num_wann, num_wann, 3, np))
      BLOCK
        COMPLEX, ALLOCATABLE :: a_one(:, :, :)
        DO a = 1, 3
          CALL wgauge_ft_rtok(aw_r(:, :, :, a), irvec, ndegen, nrpts, kfrac, a_one)
          A_interp(:, :, a, :) = a_one
        END DO
      END BLOCK
    END IF

    ! ---- project the diagonal band velocity on the bands of this domain, write ----
    ALLOCATE(cvec(num_wann, num_wann), vc(num_wann, num_wann, 3), vrow(3, num_wann))
    IF (l_berry) ALLOCATE(Hbar(num_wann, num_wann, 3), Abar(num_wann, num_wann, 3), &
                          vfull(num_wann, num_wann, 3), omega(3, num_wann))
    ax = (/ 2, 3, 1 /)   ! Omega_gamma = eps_{gamma,alpha,beta}: (Ox<-yz, Oy<-zx, Oz<-xy)
    ay = (/ 3, 1, 2 /)
    CALL wgauge_bands_open(iu, out1, &
      '# kdist   [ E_n(eV)  vx vy vz (eV*bohr, dE/dk) ] for n=1..num_wann')
    IF (l_berry) CALL wgauge_bands_open(iuc, out2, &
      '# kdist   [ E_n(eV)  Omega_x Omega_y Omega_z (bohr^2) ] for n=1..num_wann')
    DO ip = 1, np
      cvec = hg%cvec(:, :, ip)
      DO a = 1, 3
        vc(:, :, a) = MATMUL(v_interp(:, :, a, ip), cvec)         ! v_interp_a . C
      END DO
      ! ---- diagonal band velocity <n|v|n> = dE_n/dk ----
      DO m = 1, num_wann
        DO a = 1, 3
          vrow(a, m) = hartree_to_ev_const * REAL(DOT_PRODUCT(cvec(:, m), vc(:, m, a)))
        END DO
      END DO
      CALL wgauge_bands_row(iu, hg%kdist(ip), hartree_to_ev_const*hg%evals(:, ip), vrow, '2x,f12.6')

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
                vfull(n, m, a) = Hbar(n, m, a) + CMPLX(0.0, hg%evals(n, ip) - hg%evals(m, ip)) * Abar(n, m, a)
              END IF
            END DO
          END DO
        END DO
        DO n = 1, num_wann
          DO a = 1, 3               ! a = gamma (x,y,z); curvature from the (ax,ay) plane
            acc = CMPLX(0.0, 0.0)
            DO m = 1, num_wann
              IF (m == n) CYCLE
              de = hg%evals(n, ip) - hg%evals(m, ip)
              IF (ABS(de) < 1.0e-8) CYCLE
              acc = acc + vfull(n, m, ax(a)) * vfull(m, n, ay(a)) / CMPLX(de*de, 0.0)
            END DO
            omega(a, n) = -2.0 * AIMAG(acc)
          END DO
        END DO
        CALL wgauge_bands_row(iuc, hg%kdist(ip), hartree_to_ev_const*hg%evals(:, ip), omega, '2x,es14.6')
      END IF
    END DO
    CLOSE(iu)
    IF (l_berry) CLOSE(iuc)
    WRITE(oUnit,'(a,i0,a)') 'wannierlib velocity interpolation: wrote '//TRIM(out1)//'.dat (', np, ' k-points)'
    CALL timestop('wgauge_interpolate_velocity')
  END SUBROUTINE wgauge_interpolate_velocity

END MODULE m_wgauge_interpolate_velocity
