!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_kpts_kplib
  IMPLICIT NONE
CONTAINS
  SUBROUTINE kpts_kplib(cell,sym,kpts,minDistance)
    USE m_types_cell
    USE m_types_sym
    USE m_types_kpts
    USE m_judft
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_kpts),INTENT(INOUT):: kpts
    REAL,INTENT(in),OPTIONAL  :: minDistance

    CALL judft_error("You compiled without support for kplib")

  END SUBROUTINE kpts_kplib
END MODULE m_kpts_kplib
