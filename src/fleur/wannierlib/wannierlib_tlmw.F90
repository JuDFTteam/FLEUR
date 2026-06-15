!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
!--------------------------------------------------------------------------------
MODULE m_wannierlib_tlmw
  USE m_juDFT
  USE m_types_wannierlib
  IMPLICIT NONE
CONTAINS

  SUBROUTINE wannierlib_tlmw(wannierlib, nwfs, l_nocosoc, jspin, tlmwf)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: wannierlib
    INTEGER, INTENT(IN) :: nwfs
    LOGICAL, INTENT(IN) :: l_nocosoc
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
      IF (l_nocosoc .AND. ((3 - 2 * jspin) /= wannierlib%proj_spin(nwf))) CYCLE

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
! test commit 
END MODULE m_wannierlib_tlmw
