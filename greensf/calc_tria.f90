MODULE m_calc_tria

   CONTAINS

      SUBROUTINE calc_tria(kpts,cell,input,sym)

      USE m_types
      USE m_juDFT
      USE m_constants
      USE m_triang

      IMPLICIT NONE

      TYPE(t_kpts),  INTENT(INOUT) :: kpts 
      TYPE(t_input), INTENT(INOUT) :: input
      TYPE(t_sym),   INTENT(IN)    :: sym
      TYPE(t_cell),  INTENT(IN)    :: cell

      INTEGER :: itria(3,2*kpts%nkpt), ntria
      REAL    :: atr(2*kpts%nkpt), as

      LOGICAL l_tria, l_film


      l_film = .TRUE.
      CALL timestart("Calculation of Triangles")
      CALL triang(kpts%bk,kpts%nkpt,itria,ntria,atr,as,l_film)
      l_tria = .true.
      IF (sym%invs) THEN
         IF (abs(sym%nop2*as-0.5).GT.0.000001) l_tria=.false.
      ELSE
         IF (abs(sym%nop2*as-1.0).GT.0.000001) l_tria=.false.
      ENDIF

      !IF(.NOT.l_tria) CALL juDFT_warn("Triangles may not cover whol BZ",calledby="calc_tria")

      !Write to types_kpts
      kpts%ntria = ntria
      ALLOCATE(kpts%itria(3,kpts%ntria))
      kpts%itria = itria
      ALLOCATE(kpts%voltria(kpts%ntria))
      kpts%voltria = atr(:)/as

      CALL timestop("Calculation of Triangles")

   END SUBROUTINE calc_tria



END MODULE m_calc_tria