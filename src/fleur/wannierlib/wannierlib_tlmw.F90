!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The angular half of a projection: the trial orbital the user asked for, written in the
!>  basis of complex spherical harmonics the wavefunction is expanded in.
!>
!>  tlm holds the real harmonics of Wannier90's convention -- p_x = (Y_1,-1 - Y_1,1)/sqrt(2),
!>  p_y = i(Y_1,-1 + Y_1,1)/sqrt(2), and so on -- and tlmwf then picks, per Wannier function,
!>  the combination its (l, mr) names.
!>
!>  A NEGATIVE l is a hybrid, and that is where the coefficients stop being one row of tlm:
!>  sp3 mixes l=0 and l=1 with 1/sqrt(3) and 1/sqrt(6), sp3d2 reaches l=2 as well. Those are
!>  the guesses that localise on a bond rather than on an atom.
!>
!>  m_wannierlib_rad_twd carries the same table for the spin-orbit case, where the orbital is
!>  named by (j, m_j) and the coefficients are Clebsch-Gordan instead.
MODULE m_wannierlib_tlmw
  USE m_juDFT
  USE m_types_wannierlib
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: wannierlib_tlmw
CONTAINS

  SUBROUTINE wannierlib_tlmw(wannierlib, nwfs, l_nocosoc, l_spinors, jspin, tlmwf)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: wannierlib
    INTEGER, INTENT(IN) :: nwfs
    LOGICAL, INTENT(IN) :: l_nocosoc
    !> True whenever the run carries spinors (noco OR soc). The column guard below
    !> needs THIS, not l_nocosoc: l_nocosoc is noco AND NOT soc, so under SOC it is
    !> false and the guard went inert -- both spinor components then filled all
    !> num_wann columns and the projection matrix came out rank num_wann/2.
    LOGICAL, INTENT(IN) :: l_spinors
    INTEGER, INTENT(IN) :: jspin
    COMPLEX, INTENT(OUT) :: tlmwf(0:3, -3:3, nwfs)

    INTEGER :: nwf, lr, mr
    COMPLEX :: tlm(0:3, -3:3, 1:7)
    COMPLEX :: ic

    IF (nwfs > SIZE(wannierlib%proj_l)) THEN
      CALL juDFT_error('wannierlib_tlmw: nwfs exceeds configured projections', calledby='wannierlib_tlmw')
    END IF

    CALL timestart('wannierlib_tlmw')
    ic = CMPLX(0.0, 1.0)
    tlm = CMPLX(0.0, 0.0)

    tlm(0, 0, 1) = CMPLX(1.0, 0.0)

    tlm(1, 0, 1) = CMPLX(1.0, 0.0)
    tlm(1, -1, 2) = CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(1, 1, 2) = -CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(1, -1, 3) = ic / SQRT(2.0)
    tlm(1, 1, 3) = ic / SQRT(2.0)

    tlm(2, 0, 1) = CMPLX(1.0, 0.0)
    tlm(2, 1, 2) = -CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(2, -1, 2) = CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(2, 1, 3) = ic / SQRT(2.0)
    tlm(2, -1, 3) = ic / SQRT(2.0)
    tlm(2, 2, 4) = CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(2, -2, 4) = CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(2, 2, 5) = -ic / SQRT(2.0)
    tlm(2, -2, 5) = ic / SQRT(2.0)

    tlm(3, 0, 1) = CMPLX(1.0, 0.0)
    tlm(3, 1, 2) = -CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(3, -1, 2) = CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(3, 1, 3) = ic / SQRT(2.0)
    tlm(3, -1, 3) = ic / SQRT(2.0)
    tlm(3, 2, 4) = CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(3, -2, 4) = CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(3, 2, 5) = -ic / SQRT(2.0)
    tlm(3, -2, 5) = ic / SQRT(2.0)
    tlm(3, 3, 6) = -CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(3, -3, 6) = CMPLX(1.0, 0.0) / SQRT(2.0)
    tlm(3, 3, 7) = ic / SQRT(2.0)
    tlm(3, -3, 7) = ic / SQRT(2.0)

    tlmwf = CMPLX(0.0, 0.0)

    DO nwf = 1, nwfs
      IF (l_spinors .AND. ((3 - 2 * jspin) /= wannierlib%proj_spin(nwf))) CYCLE

      lr = wannierlib%proj_l(nwf)
      mr = wannierlib%proj_m(nwf)

      IF (lr >= 0) THEN
        tlmwf(lr, :, nwf) = tlm(lr, :, mr)
      ELSE IF (lr == -1) THEN
        IF (mr == 1) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(2.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = (1.0 / SQRT(2.0)) * tlm(1, :, 2)
        ELSE IF (mr == 2) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(2.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = -(1.0 / SQRT(2.0)) * tlm(1, :, 2)
        END IF
      ELSE IF (lr == -2) THEN
        IF (mr == 1) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(3.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = -(1.0 / SQRT(6.0)) * tlm(1, :, 2) + (1.0 / SQRT(2.0)) * tlm(1, :, 3)
        ELSE IF (mr == 2) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(3.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = -(1.0 / SQRT(6.0)) * tlm(1, :, 2) - (1.0 / SQRT(2.0)) * tlm(1, :, 3)
        ELSE IF (mr == 3) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(3.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = (2.0 / SQRT(6.0)) * tlm(1, :, 2)
        END IF
      ELSE IF (lr == -3) THEN
        IF (mr == 1) THEN
          tlmwf(0, :, nwf) = 0.5 * tlm(0, :, 1)
          tlmwf(1, :, nwf) = 0.5 * (tlm(1, :, 2) + tlm(1, :, 3) + tlm(1, :, 1))
        ELSE IF (mr == 2) THEN
          tlmwf(0, :, nwf) = 0.5 * tlm(0, :, 1)
          tlmwf(1, :, nwf) = 0.5 * (tlm(1, :, 2) - tlm(1, :, 3) - tlm(1, :, 1))
        ELSE IF (mr == 3) THEN
          tlmwf(0, :, nwf) = 0.5 * tlm(0, :, 1)
          tlmwf(1, :, nwf) = 0.5 * (-tlm(1, :, 2) + tlm(1, :, 3) - tlm(1, :, 1))
        ELSE IF (mr == 4) THEN
          tlmwf(0, :, nwf) = 0.5 * tlm(0, :, 1)
          tlmwf(1, :, nwf) = 0.5 * (-tlm(1, :, 2) - tlm(1, :, 3) + tlm(1, :, 1))
        END IF
      ELSE IF (lr == -4) THEN
        IF (mr == 1) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(3.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = -(1.0 / SQRT(6.0)) * tlm(1, :, 2) + (1.0 / SQRT(2.0)) * tlm(1, :, 3)
        ELSE IF (mr == 2) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(3.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = -(1.0 / SQRT(6.0)) * tlm(1, :, 2) - (1.0 / SQRT(2.0)) * tlm(1, :, 3)
        ELSE IF (mr == 3) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(3.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = (2.0 / SQRT(6.0)) * tlm(1, :, 2)
        ELSE IF (mr == 4) THEN
          tlmwf(1, :, nwf) = (1.0 / SQRT(2.0)) * tlm(1, :, 1)
          tlmwf(2, :, nwf) = (1.0 / SQRT(2.0)) * tlm(2, :, 1)
        ELSE IF (mr == 5) THEN
          tlmwf(1, :, nwf) = -(1.0 / SQRT(2.0)) * tlm(1, :, 1)
          tlmwf(2, :, nwf) = (1.0 / SQRT(2.0)) * tlm(2, :, 1)
        END IF
      ELSE IF (lr == -5) THEN
        IF (mr == 1) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(6.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = -(1.0 / SQRT(2.0)) * tlm(1, :, 2)
          tlmwf(2, :, nwf) = -(1.0 / SQRT(12.0)) * tlm(2, :, 1) + 0.5 * tlm(2, :, 4)
        ELSE IF (mr == 2) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(6.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = (1.0 / SQRT(2.0)) * tlm(1, :, 2)
          tlmwf(2, :, nwf) = -(1.0 / SQRT(12.0)) * tlm(2, :, 1) + 0.5 * tlm(2, :, 4)
        ELSE IF (mr == 3) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(6.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = -(1.0 / SQRT(2.0)) * tlm(1, :, 3)
          tlmwf(2, :, nwf) = -(1.0 / SQRT(12.0)) * tlm(2, :, 1) - 0.5 * tlm(2, :, 4)
        ELSE IF (mr == 4) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(6.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = (1.0 / SQRT(2.0)) * tlm(1, :, 3)
          tlmwf(2, :, nwf) = -(1.0 / SQRT(12.0)) * tlm(2, :, 1) - 0.5 * tlm(2, :, 4)
        ELSE IF (mr == 5) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(6.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = -(1.0 / SQRT(2.0)) * tlm(1, :, 1)
          tlmwf(2, :, nwf) = (1.0 / SQRT(3.0)) * tlm(2, :, 1)
        ELSE IF (mr == 6) THEN
          tlmwf(0, :, nwf) = (1.0 / SQRT(6.0)) * tlm(0, :, 1)
          tlmwf(1, :, nwf) = (1.0 / SQRT(2.0)) * tlm(1, :, 1)
          tlmwf(2, :, nwf) = (1.0 / SQRT(3.0)) * tlm(2, :, 1)
        END IF
      ELSE
        CALL juDFT_error('no tlmw for this lr', calledby='wannierlib_tlmw')
      END IF
    END DO

    CALL timestop('wannierlib_tlmw')
  END SUBROUTINE wannierlib_tlmw

END MODULE m_wannierlib_tlmw
