!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Operator-agnostic Wannier (slow) Fourier interpolation core.
!>
!>  Given any k-space matrix mat_k(nw,nw,nk) sampled on the coarse uniform
!>  Wannier mesh (in the Wannier gauge), interpolates it onto an arbitrary
!>  fractional k-path via the standard Wannier90 real-space transform:
!>
!>    mat_r(R) = (1/Nk) sum_k e^{-i 2pi k.R} mat_k(k)           (coarse mesh -> R)
!>    mat(k')  = sum_R  e^{+i 2pi k'.R} / ndegen(R) * mat_r(R)  (R -> fine path)
!>
!>  over the Wigner-Seitz supercell of R vectors with their degeneracies.
!>  This routine is agnostic to the physical meaning of mat: feeding it the
!>  Wannier-gauge Hamiltonian H_W(k) interpolates the band structure; feeding
!>  it a spin / orbital-moment / position operator interpolates that operator.
!>  A new operator interpolator only has to build its own mat_k(k) and reuse
!>  this core unchanged (see m_wannierlib_interpolate for the H example).
MODULE m_wannierlib_ft
  USE m_constants, ONLY : tpi_const
  USE m_types_cell
  USE m_types_kpts
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: wannierlib_ft_interpolate, wannierlib_ft_velocity, wannierlib_ft_to_real, wannierlib_ws_vectors, wannierlib_ft_to_real_reduce, wannierlib_ft_rtok
  ! Wigner-Seitz R-vectors depend only on the mesh (operator-independent) -> compute once, cache, reuse.
  INTEGER, ALLOCATABLE :: ws_irvec_c(:, :), ws_ndegen_c(:)
  INTEGER :: ws_nrpts_c = 0, ws_mp_c(3) = 0
CONTAINS

  !> Coarse-mesh Wannier-gauge matrix -> real space: mat_r(R) = (1/Nk) sum_k e^{-i2pi k.R} mat_k(k),
  !> over the Wigner-Seitz R-mesh (same irvec/ndegen as the interpolation). Used to export the
  !> tight-binding operators H(R), A(R) in Wannier90 _hr.dat / _r.dat format for external post-proc.
  SUBROUTINE wannierlib_ft_to_real(cell, mat_k, kpts, mat_r, irvec, ndegen, nrpts)
    TYPE(t_cell), INTENT(IN) :: cell
    COMPLEX, INTENT(IN) :: mat_k(:, :, :)     ! (nw, nw, nk)  Wannier-gauge matrix on coarse mesh
    TYPE(t_kpts), INTENT(IN) :: kpts
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: mat_r(:, :, :)      ! (nw, nw, nrpts)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: irvec(:, :), ndegen(:)
    INTEGER, INTENT(OUT) :: nrpts
    INTEGER :: nw, nk, k, irpt, mp_grid(3)
    REAL :: rdotk
    nw = SIZE(mat_k, 1); nk = SIZE(mat_k, 3); mp_grid = kpts%nkpt3
    CALL ws_get(cell, mp_grid, irvec, ndegen, nrpts)
    ALLOCATE(mat_r(nw, nw, nrpts), source=CMPLX(0.0, 0.0))
    DO irpt = 1, nrpts
      DO k = 1, nk
        rdotk = tpi_const * DOT_PRODUCT(kpts%bkf(:, k), REAL(irvec(:, irpt)))
        mat_r(:, :, irpt) = mat_r(:, :, irpt) + (EXP(CMPLX(0.0, -rdotk)) / REAL(nk)) * mat_k(:, :, k)
      END DO
    END DO
  END SUBROUTINE wannierlib_ft_to_real

  !> Velocity interpolation v_alpha(k') = dH/dk_alpha in the Wannier gauge:
  !>    v_alpha(k') = sum_R  i R_cart(alpha) e^{+i2pi k'.R} / ndegen(R) * mat_r(R)
  !> with mat_r the Wannier-gauge Hamiltonian in real space (same mat_r as the core)
  !> and R_cart = amat . irvec the Cartesian lattice vector (Bohr). Returns the three
  !> Cartesian components. Projected on the band eigenvectors, the DIAGONAL <n|v|n> =
  !> dE_n/dk is exact; off-diagonal elements omit the Berry-connection (gauge) term,
  !> so use the diagonal (band velocity) only until the position operator is available.
  SUBROUTINE wannierlib_ft_velocity(cell, mat_k, kpts, kfrac, vel_interp)
    TYPE(t_cell), INTENT(IN) :: cell
    COMPLEX, INTENT(IN) :: mat_k(:, :, :)     ! (nw, nw, nk)  Wannier-gauge H on coarse mesh
    TYPE(t_kpts), INTENT(IN) :: kpts
    REAL,    INTENT(IN) :: kfrac(:, :)        ! (3, nfine)
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: vel_interp(:, :, :, :)  ! (nw, nw, 3, nfine)

    INTEGER :: nw, nk, nfine, k, irpt, ip, nrpts, mp_grid(3), alpha
    INTEGER, ALLOCATABLE :: irvec(:, :), ndegen(:)
    COMPLEX, ALLOCATABLE :: mat_r(:, :, :)
    REAL :: rdotk, rcart(3)
    COMPLEX :: base

    nw = SIZE(mat_k, 1); nk = SIZE(mat_k, 3); nfine = SIZE(kfrac, 2)
    mp_grid = kpts%nkpt3
    CALL ws_get(cell, mp_grid, irvec, ndegen, nrpts)

    ! coarse mesh -> R  (identical to the core)
    ALLOCATE(mat_r(nw, nw, nrpts), source=CMPLX(0.0, 0.0))
    DO irpt = 1, nrpts
      DO k = 1, nk
        rdotk = tpi_const * DOT_PRODUCT(kpts%bkf(:, k), REAL(irvec(:, irpt)))
        mat_r(:, :, irpt) = mat_r(:, :, irpt) + (EXP(CMPLX(0.0, -rdotk)) / REAL(nk)) * mat_k(:, :, k)
      END DO
    END DO

    ! R -> fine path, weighted by i R_cart(alpha)
    ALLOCATE(vel_interp(nw, nw, 3, nfine), source=CMPLX(0.0, 0.0))
    DO ip = 1, nfine
      DO irpt = 1, nrpts
        rcart = MATMUL(cell%amat, REAL(irvec(:, irpt)))
        rdotk = tpi_const * DOT_PRODUCT(kfrac(:, ip), REAL(irvec(:, irpt)))
        base = EXP(CMPLX(0.0, rdotk)) / REAL(ndegen(irpt))
        DO alpha = 1, 3
          vel_interp(:, :, alpha, ip) = vel_interp(:, :, alpha, ip) + &
              CMPLX(0.0, rcart(alpha)) * base * mat_r(:, :, irpt)
        END DO
      END DO
    END DO
    DEALLOCATE(mat_r, irvec, ndegen)
  END SUBROUTINE wannierlib_ft_velocity

  !> Fourier-interpolate a coarse-mesh k-space matrix onto a fine k-path.
  SUBROUTINE wannierlib_ft_interpolate(cell, mat_k, kpts, kfrac, mat_interp)
    TYPE(t_cell), INTENT(IN) :: cell
    COMPLEX, INTENT(IN) :: mat_k(:, :, :)     ! (nw, nw, nk)  Wannier-gauge matrix on coarse mesh
    TYPE(t_kpts), INTENT(IN) :: kpts
    REAL,    INTENT(IN) :: kfrac(:, :)        ! (3, nfine)    fractional coords of the fine path
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: mat_interp(:, :, :)  ! (nw, nw, nfine)

    INTEGER :: nw, nk, nfine, k, irpt, ip, nrpts, mp_grid(3)
    INTEGER, ALLOCATABLE :: irvec(:, :), ndegen(:)
    COMPLEX, ALLOCATABLE :: mat_r(:, :, :)
    REAL :: rdotk
    COMPLEX :: fac

    nw = SIZE(mat_k, 1); nk = SIZE(mat_k, 3); nfine = SIZE(kfrac, 2)
    mp_grid = kpts%nkpt3

    CALL ws_get(cell, mp_grid, irvec, ndegen, nrpts)

    ! coarse mesh -> R
    ALLOCATE(mat_r(nw, nw, nrpts), source=CMPLX(0.0, 0.0))
    DO irpt = 1, nrpts
      DO k = 1, nk
        rdotk = tpi_const * DOT_PRODUCT(kpts%bkf(:, k), REAL(irvec(:, irpt)))
        fac = EXP(CMPLX(0.0, -rdotk)) / REAL(nk)
        mat_r(:, :, irpt) = mat_r(:, :, irpt) + fac * mat_k(:, :, k)
      END DO
    END DO

    ! R -> fine path
    ALLOCATE(mat_interp(nw, nw, nfine), source=CMPLX(0.0, 0.0))
    DO ip = 1, nfine
      DO irpt = 1, nrpts
        rdotk = tpi_const * DOT_PRODUCT(kfrac(:, ip), REAL(irvec(:, irpt)))
        fac = EXP(CMPLX(0.0, rdotk)) / REAL(ndegen(irpt))
        mat_interp(:, :, ip) = mat_interp(:, :, ip) + fac * mat_r(:, :, irpt)
      END DO
    END DO
    DEALLOCATE(mat_r, irvec, ndegen)
  END SUBROUTINE wannierlib_ft_interpolate

  !> Wigner-Seitz supercell R vectors + degeneracies (replicates W90
  !> hamiltonian_wigner_seitz). Kept local to preserve W90's exact ndegen
  !> convention, on which the interpolation weights depend; FLEUR's generic
  !> cell%calculate_WSweight uses a different convention and is not a drop-in.
  SUBROUTINE wigner_seitz(cell, mp_grid, irvec, ndegen, nrpts)
    TYPE(t_cell), INTENT(IN) :: cell
    INTEGER, INTENT(IN)  :: mp_grid(3)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: irvec(:, :), ndegen(:)
    INTEGER, INTENT(OUT) :: nrpts
    INTEGER, PARAMETER :: ws = 2
    REAL,    PARAMETER :: tol = 1.0e-5
    INTEGER :: n1, n2, n3, i1, i2, i3, icnt, i, jj, ndeg, pass, ndiff(3), dist_dim, imid
    REAL    :: metric(3, 3), dist_min
    REAL, ALLOCATABLE :: dist(:)

    metric = MATMUL(TRANSPOSE(cell%amat), cell%amat)
    dist_dim = (2*ws + 3)**3; imid = (dist_dim + 1) / 2
    ALLOCATE(dist(dist_dim))
    DO pass = 1, 2
      nrpts = 0
      DO n1 = -ws*mp_grid(1), ws*mp_grid(1)
        DO n2 = -ws*mp_grid(2), ws*mp_grid(2)
          DO n3 = -ws*mp_grid(3), ws*mp_grid(3)
            icnt = 0
            DO i1 = -ws-1, ws+1
              DO i2 = -ws-1, ws+1
                DO i3 = -ws-1, ws+1
                  icnt = icnt + 1
                  ndiff = (/ n1 - i1*mp_grid(1), n2 - i2*mp_grid(2), n3 - i3*mp_grid(3) /)
                  dist(icnt) = 0.0
                  DO i = 1, 3
                    DO jj = 1, 3
                      dist(icnt) = dist(icnt) + REAL(ndiff(i))*metric(i, jj)*REAL(ndiff(jj))
                    END DO
                  END DO
                END DO
              END DO
            END DO
            dist_min = MINVAL(dist)
            IF (ABS(dist(imid) - dist_min) < tol**2) THEN
              nrpts = nrpts + 1
              IF (pass == 2) THEN
                ndeg = 0
                DO i = 1, dist_dim
                  IF (ABS(dist(i) - dist_min) < tol**2) ndeg = ndeg + 1
                END DO
                ndegen(nrpts) = ndeg
                irvec(:, nrpts) = (/ n1, n2, n3 /)
              END IF
            END IF
          END DO
        END DO
      END DO
      IF (pass == 1) ALLOCATE(irvec(3, nrpts), ndegen(nrpts))
    END DO
    DEALLOCATE(dist)
  END SUBROUTINE wigner_seitz

  ! Cached Wigner-Seitz accessor: computes irvec/ndegen once per (cell,mp_grid), then returns
  ! copies. Eliminates the ~15x recomputation across operators/components (costly at 16^3+).
  SUBROUTINE ws_get(cell, mp_grid, irvec, ndegen, nrpts)
    TYPE(t_cell), INTENT(IN) :: cell
    INTEGER, INTENT(IN)  :: mp_grid(3)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: irvec(:, :), ndegen(:)
    INTEGER, INTENT(OUT) :: nrpts
    IF (ws_nrpts_c == 0 .OR. ANY(ws_mp_c /= mp_grid)) THEN
      IF (ALLOCATED(ws_irvec_c)) DEALLOCATE(ws_irvec_c, ws_ndegen_c)
      CALL wigner_seitz(cell, mp_grid, ws_irvec_c, ws_ndegen_c, ws_nrpts_c)
      ws_mp_c = mp_grid
    END IF
    nrpts = ws_nrpts_c
    ALLOCATE(irvec(3, nrpts), ndegen(nrpts))
    irvec = ws_irvec_c; ndegen = ws_ndegen_c
  END SUBROUTINE ws_get

  !> Distributed coarse-mesh -> real space: each rank sums its OWN k-slice (mat_loc, global
  !> indices gk_loc), then MPI_ALLREDUCE the small mat_r(nw,nw,nrpts). Same result as the
  !> serial wannierlib_ft_to_real but never materializes the full-mesh coarse matrix.
  SUBROUTINE wannierlib_ft_to_real_reduce(cell, kpts, mat_loc, gk_loc, commw, mat_r, irvec, ndegen, nrpts)
#ifdef CPP_MPI
    use mpi
#endif
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_kpts), INTENT(IN) :: kpts
    COMPLEX, INTENT(IN) :: mat_loc(:, :, :)   ! (nw,nw,>=nk_loc) Wannier-gauge, this rank's slice
    INTEGER, INTENT(IN) :: gk_loc(:)          ! (nk_loc) global k index of each slice entry
    INTEGER, INTENT(IN) :: commw
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: mat_r(:, :, :)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: irvec(:, :), ndegen(:)
    INTEGER, INTENT(OUT) :: nrpts
    INTEGER :: nw, kl, irpt, gk, mp_grid(3)
    REAL :: rdotk
#ifdef CPP_MPI
    INTEGER :: ierr
#endif
    nw = SIZE(mat_loc, 1); mp_grid = kpts%nkpt3
    CALL ws_get(cell, mp_grid, irvec, ndegen, nrpts)
    ALLOCATE(mat_r(nw, nw, nrpts), source=CMPLX(0.0, 0.0))
    DO irpt = 1, nrpts
      DO kl = 1, SIZE(gk_loc)
        gk = gk_loc(kl)
        rdotk = tpi_const * DOT_PRODUCT(kpts%bkf(:, gk), REAL(irvec(:, irpt)))
        mat_r(:, :, irpt) = mat_r(:, :, irpt) + (EXP(CMPLX(0.0, -rdotk)) / REAL(kpts%nkptf)) * mat_loc(:, :, kl)
      END DO
    END DO
#ifdef CPP_MPI
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, mat_r, SIZE(mat_r), MPI_DOUBLE_COMPLEX, MPI_SUM, commw, ierr)
#endif
  END SUBROUTINE wannierlib_ft_to_real_reduce

  !> R -> fine path only: mat(k') = sum_R e^{+i2pi k'.R} / ndegen(R) * mat_r(R).
  !> Second half of wannierlib_ft_interpolate, split out so a coarse->R matrix already
  !> assembled by the distributed reduce (wannierlib_ft_to_real_reduce) can be interpolated
  !> onto the fine path on rank 0 without rebuilding it from the full coarse mesh.
  SUBROUTINE wannierlib_ft_rtok(mat_r, irvec, ndegen, nrpts, kfrac, mat_interp)
    COMPLEX, INTENT(IN) :: mat_r(:, :, :)     ! (nw, nw, nrpts)  real-space matrix
    INTEGER, INTENT(IN) :: irvec(:, :), ndegen(:), nrpts
    REAL,    INTENT(IN) :: kfrac(:, :)        ! (3, nfine)    fractional coords of the fine path
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: mat_interp(:, :, :)  ! (nw, nw, nfine)
    INTEGER :: nw, nfine, ip, irpt
    REAL :: rdotk
    COMPLEX :: fac
    nw = SIZE(mat_r, 1); nfine = SIZE(kfrac, 2)
    ALLOCATE(mat_interp(nw, nw, nfine), source=CMPLX(0.0, 0.0))
    DO ip = 1, nfine
      DO irpt = 1, nrpts
        rdotk = tpi_const * DOT_PRODUCT(kfrac(:, ip), REAL(irvec(:, irpt)))
        fac = EXP(CMPLX(0.0, rdotk)) / REAL(ndegen(irpt))
        mat_interp(:, :, ip) = mat_interp(:, :, ip) + fac * mat_r(:, :, irpt)
      END DO
    END DO
  END SUBROUTINE wannierlib_ft_rtok

  !> Public wrapper exposing the Wigner-Seitz R-vectors + degeneracies (W90 convention)
  !> so the real-space operator export (operators_r) can write wig_vectors.
  SUBROUTINE wannierlib_ws_vectors(cell, mp_grid, irvec, ndegen, nrpts)
    TYPE(t_cell), INTENT(IN) :: cell
    INTEGER, INTENT(IN)  :: mp_grid(3)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: irvec(:, :), ndegen(:)
    INTEGER, INTENT(OUT) :: nrpts
    CALL ws_get(cell, mp_grid, irvec, ndegen, nrpts)
  END SUBROUTINE wannierlib_ws_vectors

END MODULE m_wannierlib_ft
