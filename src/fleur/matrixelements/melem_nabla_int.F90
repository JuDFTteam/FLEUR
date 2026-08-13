!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The interstitial half of <u_a| e^{-i b.r} p |u_b>, the momentum rather than the
!>  gradient: acting on a plane wave the momentum returns the plane wave times its own
!>  wavevector, so what the overlap needs is the same sum with one real factor more, and
!>  that factor depends on the ket alone.
!>
!>  Which is why the ket coefficients are the ones scaled: the step function between the
!>  two basis sets is untouched and the two matrix products that follow are the ones the
!>  plain overlap already does.
!>
!>  Returns p and not nabla, so the caller has no factor to add here; the muffin-tin half
!>  returns nabla and is the one that carries the -i. The two halves are only comparable
!>  once that is applied.
MODULE m_melem_nabla_int
  USE m_juDFT
  USE m_types
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: melem_nabla_int
CONTAINS

  !> bkpt_b is the ket k-point, and the wavevector that multiplies each of its plane waves
  !> is k_b + G, taken to cartesian coordinates because the component asked for is.
  !>
  !> Accumulates, like the overlap it mirrors.
  SUBROUTINE melem_nabla_int(stars, cell, lapw, lapw_b, jspin, jspin_b, zMat, zMat_b, &
                             bkpt_b, gb, pnk, ioff, ioff_b)
    TYPE(t_stars), INTENT(IN) :: stars
    TYPE(t_cell), INTENT(IN) :: cell
    TYPE(t_lapw), INTENT(IN) :: lapw
    TYPE(t_lapw), INTENT(IN) :: lapw_b
    INTEGER, INTENT(IN) :: jspin
    INTEGER, INTENT(IN) :: jspin_b
    TYPE(t_mat), INTENT(IN) :: zMat
    TYPE(t_mat), INTENT(IN) :: zMat_b
    REAL, INTENT(IN) :: bkpt_b(3)
    INTEGER, INTENT(IN) :: gb(3)
    COMPLEX, INTENT(INOUT) :: pnk(:, :, :)   !> (nbnd, nbnd, 3), cartesian
    INTEGER, INTENT(IN), OPTIONAL :: ioff, ioff_b

    INTEGER :: nbnd, nv, nv_b, i, j, j1, j2, j3, i1, i2, i3, in, io, io_b, a
    COMPLEX :: phasefac
    COMPLEX, ALLOCATABLE :: stepf(:, :), zb(:, :), zbs(:, :), za(:, :), tmp(:, :), res(:, :)
    REAL, ALLOCATABLE :: qg(:, :)
    REAL :: gvec(3)

    CALL timestart("melem_nabla_int")

    nbnd = SIZE(pnk, 1)
    nv = lapw%nv(jspin)
    nv_b = lapw_b%nv(jspin_b)
    io = 0;   IF (PRESENT(ioff))   io = ioff
    io_b = 0; IF (PRESENT(ioff_b)) io_b = ioff_b

    ALLOCATE (stepf(nv_b, nv), source=CMPLX(0.0, 0.0))
    DO i = 1, nv
      j1 = -lapw%k1(i, jspin) - gb(1)
      j2 = -lapw%k2(i, jspin) - gb(2)
      j3 = -lapw%k3(i, jspin) - gb(3)
      DO j = 1, nv_b
        i1 = j1 + lapw_b%k1(j, jspin_b)
        i2 = j2 + lapw_b%k2(j, jspin_b)
        i3 = j3 + lapw_b%k3(j, jspin_b)
        in = stars%ig(i1, i2, i3)
        IF (in == 0) CYCLE
        phasefac = stars%rgphs(i1, i2, i3)*stars%ustep(in)
        stepf(j, i) = CONJG(phasefac)
      END DO
    END DO

    !> The wavevector of each of the ket's plane waves, in cartesian coordinates.
    ALLOCATE (qg(nv_b, 3))
    DO j = 1, nv_b
      gvec = [REAL(lapw_b%k1(j, jspin_b)), REAL(lapw_b%k2(j, jspin_b)), REAL(lapw_b%k3(j, jspin_b))]
      qg(j, :) = MATMUL(bkpt_b + gvec, cell%bmat)
    END DO

    !> Real eigenvectors are promoted rather than given their own path: the momentum has
    !> no reason to stay real, and one arithmetic here is easier to trust than two.
    ALLOCATE (za(nv, nbnd), zb(nv_b, nbnd), zbs(nv_b, nbnd), tmp(nv, nbnd), res(nbnd, nbnd))
    IF (zMat%l_real) THEN
      za = CMPLX(zMat%data_r(1 + io:nv + io, 1:nbnd), 0.0)
      zb = CMPLX(zMat_b%data_r(1 + io_b:nv_b + io_b, 1:nbnd), 0.0)
    ELSE
      za = zMat%data_c(1 + io:nv + io, 1:nbnd)
      zb = zMat_b%data_c(1 + io_b:nv_b + io_b, 1:nbnd)
    END IF

    DO a = 1, 3
      DO j = 1, nv_b
        zbs(j, :) = qg(j, a)*zb(j, :)
      END DO
      CALL zgemm('T', 'N', nv, nbnd, nv_b, CMPLX(1.0, 0.0), &
                 stepf, nv_b, zbs, nv_b, CMPLX(0.0, 0.0), tmp, nv)
      tmp = CONJG(tmp)
      CALL zgemm('T', 'N', nbnd, nbnd, nv, CMPLX(1.0, 0.0), &
                 za, nv, tmp, nv, CMPLX(0.0, 0.0), res, nbnd)
      pnk(:, :, a) = pnk(:, :, a) + res
    END DO

    DEALLOCATE (stepf, qg, za, zb, zbs, tmp, res)
    CALL timestop("melem_nabla_int")
  END SUBROUTINE melem_nabla_int

END MODULE m_melem_nabla_int
