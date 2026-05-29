!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
!--------------------------------------------------------------------------------
MODULE m_wannierlib_amn_stub
  USE m_juDFT
  IMPLICIT NONE
CONTAINS

  SUBROUTINE wannierlib_build_amn_stub(num_bands, num_wann, num_kpts, amn, implemented)
    INTEGER, INTENT(IN) :: num_bands
    INTEGER, INTENT(IN) :: num_wann
    INTEGER, INTENT(IN) :: num_kpts
    COMPLEX, ALLOCATABLE, INTENT(INOUT) :: amn(:, :, :)
    LOGICAL, INTENT(OUT) :: implemented

    IF (.NOT.ALLOCATED(amn)) THEN
      IF ((num_bands > 0) .AND. (num_wann > 0) .AND. (num_kpts > 0)) THEN
        ALLOCATE(amn(num_bands, num_wann, num_kpts))
        amn = CMPLX(0.0, 0.0)
      END IF
    END IF

    implemented = .FALSE.
    CALL juDFT_error('wannierlib AMN stub reached: not implemented yet', calledby='wannierlib_build_amn_stub')
  END SUBROUTINE wannierlib_build_amn_stub

END MODULE m_wannierlib_amn_stub
