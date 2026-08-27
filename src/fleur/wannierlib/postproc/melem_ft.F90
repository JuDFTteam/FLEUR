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
!>  this core unchanged.
MODULE m_melem_ft
  USE m_juDFT
  USE m_constants, ONLY : tpi_const, oUnit, bohr_to_angstrom_const
  USE m_types_cell
  USE m_types_kpts
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: melem_ft_interpolate, melem_ft_rtok_velocity, melem_ft_to_real, melem_ws_vectors, melem_ft_to_real_reduce, melem_ft_rtok, melem_ws_distance
  PUBLIC :: melem_mdrs_set, melem_mdrs_clear
  ! Wigner-Seitz R-vectors depend only on the mesh (operator-independent) -> compute once, cache, reuse.
  INTEGER, ALLOCATABLE :: ws_irvec_c(:, :), ws_ndegen_c(:)
  INTEGER :: ws_nrpts_c = 0, ws_mp_c(3) = 0
  !> MDRS replica table for the wannierization currently being post-processed. It is a
  !> property of the pair (mesh, Wannier centres) exactly as irvec is a property of the mesh,
  !> so it is cached here rather than threaded through five call chains -- and, more to the
  !> point, every interpolation of one run has to share it: MDRS on the Hamiltonian and off
  !> on an operator would evaluate the two in different gauges, and their product would be
  !> meaningless. Set by melem_mdrs_set, which the interpolation driver calls once.
  LOGICAL :: mdrs_on = .FALSE.
  INTEGER, ALLOCATABLE :: mdrs_ndeg(:, :, :), mdrs_irdist(:, :, :, :, :)
  INTEGER :: mdrs_nw = 0, mdrs_nrpts = 0
CONTAINS

  !> Coarse-mesh Wannier-gauge matrix -> real space: mat_r(R) = (1/Nk) sum_k e^{-i2pi k.R} mat_k(k),
  !> over the Wigner-Seitz R-mesh (same irvec/ndegen as the interpolation). Used to export the
  !> tight-binding operators H(R), A(R) in Wannier90 _hr.dat / _r.dat format for external post-proc.
  SUBROUTINE melem_ft_to_real(cell, mat_k, kpts, mat_r, irvec, ndegen, nrpts)
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
  END SUBROUTINE melem_ft_to_real

  !> Velocity from an ALREADY transformed Wannier-gauge Hamiltonian:
  !>    v_alpha(k') = sum_R  i R_cart(alpha) e^{+i2pi k'.R} / ndegen(R) * mat_r(R)
  !> with R_cart = amat . irvec (Bohr). Takes mat_r rather than the coarse matrix so that a
  !> caller wanting both H(k') and dH/dk transforms the coarse mesh once.
  !>
  !> Projected on the band eigenvectors the DIAGONAL <n|v|n> = dE_n/dk is exact; off-diagonal
  !> elements omit the Berry-connection term.
  SUBROUTINE melem_ft_rtok_velocity(cell, mat_r, irvec, ndegen, nrpts, kfrac, vel_interp)
    TYPE(t_cell), INTENT(IN) :: cell
    COMPLEX, INTENT(IN) :: mat_r(:, :, :)     ! (nw, nw, nrpts) Wannier-gauge H in real space
    INTEGER, INTENT(IN) :: irvec(:, :), ndegen(:), nrpts
    REAL,    INTENT(IN) :: kfrac(:, :)        ! (3, nfine)
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: vel_interp(:, :, :, :)  ! (nw, nw, 3, nfine)

    INTEGER :: nw, nfine, irpt, ip, alpha, m, n, j
    REAL :: rdotk, rcart(3)
    COMPLEX :: base, facv(3)

    nw = SIZE(mat_r, 1); nfine = SIZE(kfrac, 2)
    CALL mdrs_assert(nw, nrpts, 'melem_ft_rtok_velocity')
    ALLOCATE(vel_interp(nw, nw, 3, nfine), source=CMPLX(0.0, 0.0))
    DO ip = 1, nfine
      DO irpt = 1, nrpts
        IF (mdrs_on) THEN
          !> The same replica average melem_ft_rtok does, with the derivative taken INSIDE
          !> it: the factor i R_cart belongs to the replica that carries the phase and not
          !> to the nominal R. Pulling it out would interpolate v in a different gauge from
          !> H, and the identity <n|v|n> = dE_n/dk would stop holding.
          DO n = 1, nw
            DO m = 1, nw
              facv = CMPLX(0.0, 0.0)
              DO j = 1, mdrs_ndeg(m, n, irpt)
                rcart = MATMUL(cell%amat, REAL(mdrs_irdist(:, j, m, n, irpt)))
                rdotk = tpi_const * DOT_PRODUCT(kfrac(:, ip), REAL(mdrs_irdist(:, j, m, n, irpt)))
                DO alpha = 1, 3
                  facv(alpha) = facv(alpha) + CMPLX(0.0, rcart(alpha)) * EXP(CMPLX(0.0, rdotk))
                END DO
              END DO
              facv = facv / (REAL(ndegen(irpt)) * REAL(mdrs_ndeg(m, n, irpt)))
              DO alpha = 1, 3
                vel_interp(m, n, alpha, ip) = vel_interp(m, n, alpha, ip) + facv(alpha) * mat_r(m, n, irpt)
              END DO
            END DO
          END DO
        ELSE
          rcart = MATMUL(cell%amat, REAL(irvec(:, irpt)))
          rdotk = tpi_const * DOT_PRODUCT(kfrac(:, ip), REAL(irvec(:, irpt)))
          base = EXP(CMPLX(0.0, rdotk)) / REAL(ndegen(irpt))
          DO alpha = 1, 3
            vel_interp(:, :, alpha, ip) = vel_interp(:, :, alpha, ip) + &
                CMPLX(0.0, rcart(alpha)) * base * mat_r(:, :, irpt)
          END DO
        END IF
      END DO
    END DO
  END SUBROUTINE melem_ft_rtok_velocity

  !> Fourier-interpolate a coarse-mesh k-space matrix onto a fine k-path.
  SUBROUTINE melem_ft_interpolate(cell, mat_k, kpts, kfrac, mat_interp)
    TYPE(t_cell), INTENT(IN) :: cell
    COMPLEX, INTENT(IN) :: mat_k(:, :, :)     ! (nw, nw, nk)  Wannier-gauge matrix on coarse mesh
    TYPE(t_kpts), INTENT(IN) :: kpts
    REAL,    INTENT(IN) :: kfrac(:, :)        ! (3, nfine)    fractional coords of the fine path
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: mat_interp(:, :, :)  ! (nw, nw, nfine)

    INTEGER :: nrpts
    INTEGER, ALLOCATABLE :: irvec(:, :), ndegen(:)
    COMPLEX, ALLOCATABLE :: mat_r(:, :, :)

    !> The two halves, in the order they have to happen. Kept as one entry point because most
    !> callers want exactly this; a caller that also needs the derivative uses the halves.
    CALL melem_ft_to_real(cell, mat_k, kpts, mat_r, irvec, ndegen, nrpts)
    CALL melem_ft_rtok(mat_r, irvec, ndegen, nrpts, kfrac, mat_interp)
    DEALLOCATE(mat_r, irvec, ndegen)
  END SUBROUTINE melem_ft_interpolate

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

  !> MDRS -- minimum distance replica selection. For each Wannier pair (m,n) and each R of the
  !> Wigner-Seitz mesh, find the supercell translation T = (j1,j2,j3)*mp_grid that minimises the
  !> REAL distance |R + T + tau_n - tau_m|, and count the replicas that tie.
  !>
  !> Why it exists: the plain transform phases the hopping <0m|H|Rn> as e^{ikR}, i.e. as if it
  !> connected cell ORIGINS. It does not -- m sits at tau_m and n at tau_n, so the physical
  !> separation is R + tau_n - tau_m. The plain assignment is discontinuous at the Wigner-Seitz
  !> boundary (the same hopping lands on different R for different pairs), H(R) then decays
  !> slowly, and a slowly-decaying Fourier series rings. Averaging over the tied minimum-distance
  !> replicas restores the real distance and the ringing softens.
  !>
  !> It matters when the centres sit far from the cell origin -- bond-centred, or pushed into the
  !> vacuum of a film -- and when the cell is anisotropic. Both hold here.
  !>
  !> Modelled on Wannier90's ws_distance.F90 (Paulatto, Gibertini, Gresch, Pizzi), reimplemented
  !> rather than called: that module is W90-internal and our interpolation never enters W90.
  SUBROUTINE melem_ws_distance(cell, mp_grid, irvec, nrpts, cfrac, ndeg, irdist)
    TYPE(t_cell), INTENT(IN) :: cell
    INTEGER, INTENT(IN) :: mp_grid(3), irvec(:, :), nrpts
    REAL,    INTENT(IN) :: cfrac(:, :)        !> (3, nw) Wannier centres, FRACTIONAL
    INTEGER, ALLOCATABLE, INTENT(OUT) :: ndeg(:, :, :)       !> (nw, nw, nrpts)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: irdist(:, :, :, :, :) !> (3, ndegmax, nw, nw, nrpts)
    !> One MORE than the half-width wigner_seitz searches, and Wannier90 says why: a centre
    !> that wandered out of its reference cell needs the extra ring to find its own minimum.
    !> A film is where that happens -- the centres are pushed into the vacuum.
    INTEGER, PARAMETER :: ws = 3
    INTEGER, PARAMETER :: ndegmax = 8            !> a point can touch at most 8 cells at a vertex
    REAL,    PARAMETER :: tol = 1.0e-5
    INTEGER :: m, n, irpt, j1, j2, j3, nd, i, jj
    REAL    :: metric(3, 3), d, dmin, v(3), shift(3)
    INTEGER :: cand(3, (2*ws+1)**3), nc
    REAL    :: dcand((2*ws+1)**3)

    metric = MATMUL(TRANSPOSE(cell%amat), cell%amat)
    ALLOCATE(ndeg(SIZE(cfrac, 2), SIZE(cfrac, 2), nrpts), source=1)
    ALLOCATE(irdist(3, ndegmax, SIZE(cfrac, 2), SIZE(cfrac, 2), nrpts), source=0)

    DO irpt = 1, nrpts
      DO n = 1, SIZE(cfrac, 2)
        DO m = 1, SIZE(cfrac, 2)
          !> The separation carried by this hopping, in fractional coordinates.
          shift = REAL(irvec(:, irpt)) + cfrac(:, n) - cfrac(:, m)
          nc = 0; dmin = HUGE(1.0)
          DO j1 = -ws, ws
            DO j2 = -ws, ws
              DO j3 = -ws, ws
                v = shift + (/ REAL(j1*mp_grid(1)), REAL(j2*mp_grid(2)), REAL(j3*mp_grid(3)) /)
                d = 0.0
                DO i = 1, 3
                  DO jj = 1, 3
                    d = d + v(i)*metric(i, jj)*v(jj)
                  END DO
                END DO
                nc = nc + 1
                cand(:, nc) = (/ j1*mp_grid(1), j2*mp_grid(2), j3*mp_grid(3) /)
                dcand(nc) = d
                dmin = MIN(dmin, d)
              END DO
            END DO
          END DO
          !> Keep every replica within tol of the minimum -- ties happen when the partner sits on
          !> the Wigner-Seitz boundary, and dropping them would break the point-group symmetry
          !> of the interpolation.
          !>
          !> The comparison is on the DISTANCE and not on its square, which is what Wannier90
          !> compares and is the whole point of the tolerance: 1e-5 on d is ~1e-11 on d^2 for
          !> the separations here, and at that width rounding alone splits a genuine tie.
          nd = 0
          DO i = 1, nc
            IF (ABS(SQRT(dcand(i)) - SQRT(dmin)) < tol) THEN
              nd = nd + 1
              IF (nd > ndegmax) CALL juDFT_error( &
                'melem_ws_distance: more than 8 minimum-distance replicas, which no point can have', &
                hint='the Wannier centres or the Bravais matrix are not what they should be', &
                calledby='melem_ws_distance')
              irdist(:, nd, m, n, irpt) = irvec(:, irpt) + cand(:, i)
            END IF
          END DO
          ndeg(m, n, irpt) = MAX(nd, 1)
          IF (nd == 0) irdist(:, 1, m, n, irpt) = irvec(:, irpt)
        END DO
      END DO
    END DO
  END SUBROUTINE melem_ws_distance

  !> Switch MDRS on for every interpolation that follows, from the Wannier centres Wannier90
  !> reported. Idempotent: it rebuilds the table, so a second wannierization (the second spin
  !> channel of a collinear run) replaces the first channel's centres instead of reusing them.
  !>
  !> The centres arrive CARTESIAN and in ANGSTROM -- that is what w90_get_centres returns,
  !> because the adapter hands Wannier90 the lattice as bohr_to_angstrom_const*cell%amat --
  !> while the replica search works in fractional coordinates on a lattice in Bohr. That
  !> conversion is the one step here where a unit or a transpose can go wrong without any
  !> symptom, so it solves amat.x = r outright instead of trusting a stored inverse, checks
  !> itself by mapping back, and prints what it got: the centres it echoes are in Angstrom
  !> and must reproduce the ones Wannier90 printed.
  SUBROUTINE melem_mdrs_set(cell, mp_grid, centres_ang, irank)
    TYPE(t_cell), INTENT(IN) :: cell
    INTEGER, INTENT(IN) :: mp_grid(3)
    REAL,    INTENT(IN) :: centres_ang(:, :)   !> (3, nw) cartesian, Angstrom
    INTEGER, INTENT(IN) :: irank
    INTEGER :: nw, n, nrpts
    INTEGER, ALLOCATABLE :: irvec(:, :), ndegen(:)
    REAL,    ALLOCATABLE :: cfrac(:, :)
    REAL :: rc(3), back(3), err, scal

    CALL melem_mdrs_clear()
    nw = SIZE(centres_ang, 2)
    scal = MAXVAL(ABS(cell%amat))
    IF (ABS(det3(cell%amat)) < 1.0e-10 * scal**3) CALL juDFT_error( &
      'melem_mdrs_set: the Bravais matrix is singular', calledby='melem_mdrs_set')

    ALLOCATE(cfrac(3, nw))
    err = 0.0
    DO n = 1, nw
      rc = centres_ang(:, n) / bohr_to_angstrom_const
      cfrac(:, n) = solve3(cell%amat, rc)
      back = MATMUL(cell%amat, cfrac(:, n))
      err = MAX(err, MAXVAL(ABS(back - rc)))
    END DO
    IF (err > 1.0e-8 * scal) CALL juDFT_error( &
      'melem_mdrs_set: cartesian -> fractional does not invert', &
      hint='the Bravais matrix is near-singular or amat is not the column convention', &
      calledby='melem_mdrs_set')

    CALL ws_get(cell, mp_grid, irvec, ndegen, nrpts)
    CALL melem_ws_distance(cell, mp_grid, irvec, nrpts, cfrac, mdrs_ndeg, mdrs_irdist)
    mdrs_nw = nw; mdrs_nrpts = nrpts; mdrs_on = .TRUE.

    IF (irank == 0) THEN
      WRITE(oUnit, '(a,i0,a,f6.3)') 'wannierlib: MDRS on -- ', nrpts, &
        ' R vectors, mean replicas per pair ', REAL(SUM(mdrs_ndeg)) / REAL(nw*nw*nrpts)
      WRITE(oUnit, '(a)') '  n    centre (fractional)              centre (Angstrom, echoed back)'
      DO n = 1, nw
        WRITE(oUnit, '(i4,3f10.5,4x,3f10.5)') n, cfrac(:, n), &
          MATMUL(cell%amat, cfrac(:, n)) * bohr_to_angstrom_const
      END DO
    END IF
    DEALLOCATE(cfrac, irvec, ndegen)
  END SUBROUTINE melem_mdrs_set

  !> Back to the plain transform. Called at the end of the interpolation pass so the state
  !> cannot outlive the wannierization it describes.
  SUBROUTINE melem_mdrs_clear()
    IF (ALLOCATED(mdrs_ndeg)) DEALLOCATE(mdrs_ndeg)
    IF (ALLOCATED(mdrs_irdist)) DEALLOCATE(mdrs_irdist)
    mdrs_on = .FALSE.; mdrs_nw = 0; mdrs_nrpts = 0
  END SUBROUTINE melem_mdrs_clear

  !> What every consumer of the MDRS state assumes about it. Cached state that outlives one
  !> call is worth exactly the check that it still describes the caller.
  SUBROUTINE mdrs_assert(nw, nrpts, who)
    INTEGER, INTENT(IN) :: nw, nrpts
    CHARACTER(LEN=*), INTENT(IN) :: who
    IF (.NOT. mdrs_on) RETURN
    IF (nw /= mdrs_nw .OR. nrpts /= mdrs_nrpts) CALL judft_bug( &
      TRIM(who)//': the MDRS replica table was built for a different manifold or R-mesh')
  END SUBROUTINE mdrs_assert

  !> x of a.x = b for a 3x3, by Cramer. Local rather than cell%bmat so that no convention of
  !> a stored reciprocal matrix -- factor of 2pi, transpose -- can enter unnoticed: the only
  !> thing assumed is r_cart = amat . r_frac, which is how the cell is built.
  PURE FUNCTION solve3(a, b) RESULT(x)
    REAL, INTENT(IN) :: a(3, 3), b(3)
    REAL :: x(3), m(3, 3), d
    INTEGER :: i
    d = det3(a)
    DO i = 1, 3
      m = a; m(:, i) = b
      x(i) = det3(m) / d
    END DO
  END FUNCTION solve3

  PURE FUNCTION det3(a) RESULT(d)
    REAL, INTENT(IN) :: a(3, 3)
    REAL :: d
    d = a(1, 1)*(a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)) &
      - a(1, 2)*(a(2, 1)*a(3, 3) - a(2, 3)*a(3, 1)) &
      + a(1, 3)*(a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1))
  END FUNCTION det3

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
  !> serial melem_ft_to_real but never materializes the full-mesh coarse matrix.
  SUBROUTINE melem_ft_to_real_reduce(cell, kpts, mat_loc, gk_loc, commw, mat_r, irvec, ndegen, nrpts)
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
  END SUBROUTINE melem_ft_to_real_reduce

  !> R -> fine path only: mat(k') = sum_R e^{+i2pi k'.R} / ndegen(R) * mat_r(R).
  !> Second half of melem_ft_interpolate, split out so a coarse->R matrix already
  !> assembled by the distributed reduce (melem_ft_to_real_reduce) can be interpolated
  !> onto the fine path on rank 0 without rebuilding it from the full coarse mesh.
  !> MDRS is off unless melem_mdrs_set switched it on, in which case the phase stops being
  !> one scalar per R and becomes a MATRIX per R, because the minimum-distance replica depends
  !> on WHICH pair of Wannier centres the hopping connects. See melem_ws_distance for why.
  !> With it off this is byte for byte the transform this routine has always done.
  SUBROUTINE melem_ft_rtok(mat_r, irvec, ndegen, nrpts, kfrac, mat_interp)
    COMPLEX, INTENT(IN) :: mat_r(:, :, :)     ! (nw, nw, nrpts)  real-space matrix
    INTEGER, INTENT(IN) :: irvec(:, :), ndegen(:), nrpts
    REAL,    INTENT(IN) :: kfrac(:, :)        ! (3, nfine)    fractional coords of the fine path
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: mat_interp(:, :, :)  ! (nw, nw, nfine)
    INTEGER :: nw, nfine, ip, irpt, m, n, j
    REAL :: rdotk
    COMPLEX :: fac
    nw = SIZE(mat_r, 1); nfine = SIZE(kfrac, 2)
    CALL mdrs_assert(nw, nrpts, 'melem_ft_rtok')
    ALLOCATE(mat_interp(nw, nw, nfine), source=CMPLX(0.0, 0.0))
    DO ip = 1, nfine
      DO irpt = 1, nrpts
        IF (mdrs_on) THEN
          DO n = 1, nw
            DO m = 1, nw
              fac = CMPLX(0.0, 0.0)
              DO j = 1, mdrs_ndeg(m, n, irpt)
                rdotk = tpi_const * DOT_PRODUCT(kfrac(:, ip), REAL(mdrs_irdist(:, j, m, n, irpt)))
                fac = fac + EXP(CMPLX(0.0, rdotk))
              END DO
              fac = fac / (REAL(ndegen(irpt)) * REAL(mdrs_ndeg(m, n, irpt)))
              mat_interp(m, n, ip) = mat_interp(m, n, ip) + fac * mat_r(m, n, irpt)
            END DO
          END DO
        ELSE
          rdotk = tpi_const * DOT_PRODUCT(kfrac(:, ip), REAL(irvec(:, irpt)))
          fac = EXP(CMPLX(0.0, rdotk)) / REAL(ndegen(irpt))
          mat_interp(:, :, ip) = mat_interp(:, :, ip) + fac * mat_r(:, :, irpt)
        END IF
      END DO
    END DO
  END SUBROUTINE melem_ft_rtok

  !> Public wrapper exposing the Wigner-Seitz R-vectors + degeneracies (W90 convention)
  !> so the real-space operator export (operators_r) can write wig_vectors.
  SUBROUTINE melem_ws_vectors(cell, mp_grid, irvec, ndegen, nrpts)
    TYPE(t_cell), INTENT(IN) :: cell
    INTEGER, INTENT(IN)  :: mp_grid(3)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: irvec(:, :), ndegen(:)
    INTEGER, INTENT(OUT) :: nrpts
    CALL ws_get(cell, mp_grid, irvec, ndegen, nrpts)
  END SUBROUTINE melem_ws_vectors

END MODULE m_melem_ft
