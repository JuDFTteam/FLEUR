!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The interstitial half of <u_a| e^{-i b.r} |u_b>: the step function between the two
!>  basis sets, contracted with the coefficients of bra and ket by the two matrix products
!>  that follow it.
!>
!>  The two sides are independent -- their own lapw, their own zMat, their own spin -- and
!>  gb is the reciprocal lattice vector that brings the ket's k-point back into the zone.
!>  That the caller happens to use them with k and k+b is its business, not this routine's.
!>
!>  Accumulates into mmnk, so the muffin-tin half can be added on top.
MODULE m_melem_mmkb_int
  USE m_juDFT
  USE m_types
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: melem_mmkb_int
CONTAINS

  SUBROUTINE melem_mmkb_int(stars, lapw, lapw_b, jspin, jspin_b, zMat, zMat_b, gb, mmnk, ioff, ioff_b)
    TYPE(t_stars), INTENT(IN) :: stars
    TYPE(t_lapw), INTENT(IN) :: lapw
    TYPE(t_lapw), INTENT(IN) :: lapw_b
    INTEGER, INTENT(IN) :: jspin
    INTEGER, INTENT(IN) :: jspin_b
    TYPE(t_mat), INTENT(IN) :: zMat
    TYPE(t_mat), INTENT(IN) :: zMat_b
    INTEGER, INTENT(IN) :: gb(3)
    COMPLEX, INTENT(INOUT) :: mmnk(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: ioff, ioff_b   ! noco offset to the spin-down block (nv(1)+nlotot)

    INTEGER :: nbnd
    COMPLEX, ALLOCATABLE :: stepf_c(:, :)
    COMPLEX, ALLOCATABLE :: phasusbmat_c(:, :)
    REAL, ALLOCATABLE :: stepf_r(:, :)
    REAL, ALLOCATABLE :: phasusbmat_r(:, :)
    REAL, ALLOCATABLE :: mmnk_tmp(:, :)
    INTEGER :: i, j1, j2, j3, i1, i2, i3, j, in
    COMPLEX :: phasefac
    INTEGER :: io, io_b

    CALL timestart("melem_mmkb_int")

    nbnd=size(mmnk, 1)
    io = 0;   IF (PRESENT(ioff))   io   = ioff
    io_b = 0; IF (PRESENT(ioff_b)) io_b = ioff_b
    

    IF (zMat%l_real) THEN
      ALLOCATE(mmnk_tmp(nbnd, nbnd))
      ALLOCATE(phasusbmat_r(lapw%nv(jspin), nbnd))
      ALLOCATE(stepf_r(lapw_b%nv(jspin_b), lapw%nv(jspin)))
      stepf_r = 0.0
    ELSE
      ALLOCATE(phasusbmat_c(lapw%nv(jspin), nbnd))
      ALLOCATE(stepf_c(lapw_b%nv(jspin_b), lapw%nv(jspin)))
      stepf_c = CMPLX(0.0, 0.0)
    END IF

    DO i = 1, lapw%nv(jspin)
      j1 = -lapw%k1(i, jspin) - gb(1)
      j2 = -lapw%k2(i, jspin) - gb(2)
      j3 = -lapw%k3(i, jspin) - gb(3)
      DO j = 1, lapw_b%nv(jspin_b)
        i1 = j1 + lapw_b%k1(j, jspin_b)
        i2 = j2 + lapw_b%k2(j, jspin_b)
        i3 = j3 + lapw_b%k3(j, jspin_b)
        in = stars%ig(i1, i2, i3)
        IF (in == 0) CYCLE
        phasefac = stars%rgphs(i1, i2, i3) * stars%ustep(in)
        IF (zMat%l_real) THEN
          stepf_r(j, i) = REAL(phasefac)
        ELSE
          stepf_c(j, i) = CONJG(phasefac)
        END IF
      END DO
    END DO

    IF (zMat%l_real) THEN
      CALL dgemm('T', 'N', lapw%nv(jspin), nbnd, lapw_b%nv(jspin_b), REAL(1.0), &
                 stepf_r, lapw_b%nv(jspin_b), zMat_b%data_r(1+io_b, 1), zMat_b%matsize1, &
                 REAL(0.0), phasusbmat_r, lapw%nv(jspin))
      CALL dgemm('T', 'N', nbnd, nbnd, lapw%nv(jspin), REAL(1.0), &
                 zMat%data_r(1+io, 1), zMat%matsize1, phasusbmat_r, lapw%nv(jspin), &
                 REAL(0.0), mmnk_tmp, nbnd)
      mmnk(1:nbnd, 1:nbnd) = mmnk(1:nbnd, 1:nbnd) + mmnk_tmp(1:nbnd, 1:nbnd) 
    ELSE
      CALL zgemm('T', 'N', lapw%nv(jspin), nbnd, lapw_b%nv(jspin_b), CMPLX(1.0), &
                 stepf_c, lapw_b%nv(jspin_b), zMat_b%data_c(1+io_b, 1), zMat_b%matsize1, CMPLX(0.0), &
                 phasusbmat_c, lapw%nv(jspin))
      phasusbmat_c = CONJG(phasusbmat_c)
      CALL zgemm('T', 'N', nbnd, nbnd, lapw%nv(jspin), CMPLX(1.0, 0.0), &
                 zMat%data_c(1+io, 1), zMat%matsize1, phasusbmat_c, lapw%nv(jspin), &
                 CMPLX(1.0, 0.0), mmnk, nbnd)
    END IF


    CALL timestop("melem_mmkb_int")
  END SUBROUTINE melem_mmkb_int

END MODULE m_melem_mmkb_int
