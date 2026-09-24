!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Scaffolding shared by the Wannier-gauge interpolation drivers
!>  (m_melem_interpolate_ham / _eigenstates / _op / _velocity).
!>
!>  Nothing here is operator-specific: it is the bookkeeping every driver repeats
!>  around its own loop over the interpolation mesh. The Fourier transforms
!>  themselves live in m_melem_ft; the k-set arrives as an argument from m_melem_run.
MODULE m_melem_interp_util
  USE m_juDFT
  USE m_types_cell
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: melem_kpath, melem_zheev_workspace
CONTAINS

  !>  Cumulative distance along the interpolation mesh: the abscissa written as
  !>  the first column of every bands_wann_*.dat, in the units of cell%bmat.
  !>
  !>  Only a path gives this a physical meaning. On a domain that is not one -- a
  !>  plane, a grid -- the number still grows monotonically but is not a coordinate,
  !>  and of the drivers only the eigenstates one writes kx,ky,kz beside it. The
  !>  others write the abscissa alone, so such an output is matched to its input
  !>  kPointList by row order and by nothing else.
  SUBROUTINE melem_kpath(cell, kfrac, kdist)
    TYPE(t_cell), INTENT(IN) :: cell
    REAL, INTENT(IN) :: kfrac(:, :)                 ! (3, np) fractional mesh, in path order
    REAL, ALLOCATABLE, INTENT(OUT) :: kdist(:)      ! (np)

    INTEGER :: ip, np
    REAL    :: dkc(3)

    np = SIZE(kfrac, 2)
    ALLOCATE(kdist(MAX(1, np)))
    kdist = 0.0
    DO ip = 2, np
      dkc = MATMUL(cell%bmat, kfrac(:, ip) - kfrac(:, ip-1))
      kdist(ip) = kdist(ip-1) + SQRT(DOT_PRODUCT(dkc, dkc))
    END DO
  END SUBROUTINE melem_kpath

  !>  Allocate the zheev workspace for an n x n Hermitian problem, sized by
  !>  LAPACK own query instead of by a hardcoded rule.
  !>
  !>  jobz has to match the jobz of the solve that follows: the optimal lwork is
  !>  not the same for eigenvalues only and for eigenvectors. Meant to be called
  !>  once outside the k loop, since the size depends on n alone.
  SUBROUTINE melem_zheev_workspace(jobz, n, work, rwork, lwork)
    CHARACTER(LEN=1), INTENT(IN) :: jobz
    INTEGER, INTENT(IN) :: n
    COMPLEX, ALLOCATABLE, INTENT(OUT) :: work(:)
    REAL,    ALLOCATABLE, INTENT(OUT) :: rwork(:)
    INTEGER, INTENT(OUT) :: lwork

    INTEGER :: info
    COMPLEX :: wq(1)
    REAL,    ALLOCATABLE :: evals(:)
    COMPLEX, ALLOCATABLE :: a(:, :)

    !> The query does not read a or evals, but LAPACK still expects arrays of the
    !> declared shape, so they are allocated here rather than faked.
    ALLOCATE(rwork(MAX(1, 3*n - 2)), evals(MAX(1, n)), a(MAX(1, n), MAX(1, n)))
    CALL zheev(jobz, "U", n, a, MAX(1, n), evals, wq, -1, rwork, info)
    IF (info /= 0) CALL juDFT_error("zheev workspace query failed", calledby="melem_zheev_workspace")
    lwork = MAX(1, NINT(REAL(wq(1))))
    ALLOCATE(work(lwork))
  END SUBROUTINE melem_zheev_workspace

END MODULE m_melem_interp_util
