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
