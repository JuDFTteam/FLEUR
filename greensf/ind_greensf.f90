MODULE m_ind_greensf

!--------------------------------------------------------------------
! The information on which index i_gf in the green's function arrays
! corresponds to which (l,lp,nType,nTypep) block of the system is stored
! in the gfelem array in atoms
! This fucntion returns the index i_gf by supplying the four indices
! If lp or nTypep are not given in the call they are assumed to be equal
! to l/nType
!--------------------------------------------------------------------

CONTAINS

   FUNCTION ind_greensf(atoms,l,nType,lp,nTypep) result(i_gf)

      USE m_types_setup
      USE m_juDFT

      !Maps between the four indices (l,lp,nType,nTypep) and the position in the
      !gf arrays

      IMPLICIT NONE

      TYPE(t_atoms),       INTENT(IN)  :: atoms
      INTEGER,             INTENT(IN)  :: l
      INTEGER,             INTENT(IN)  :: nType
      INTEGER, OPTIONAL,   INTENT(IN)  :: lp
      INTEGER, OPTIONAL,   INTENT(IN)  :: nTypep

      INTEGER i_gf

      LOGICAL search

      search = .TRUE.
      i_gf = 0

      DO WHILE(search)
         i_gf = i_gf + 1

         IF(i_gf.GT.atoms%n_gf) THEN
            !Something went wrong here
            CALL juDFT_error("Greens function element not found", calledby="ind_greensf")
         ENDIF
         !--------------------------------------------
         ! Check the current element
         !--------------------------------------------
         IF(atoms%gfelem(i_gf)%l.NE.l) CYCLE
         IF(atoms%gfelem(i_gf)%atomType.NE.nType) CYCLE
         IF(PRESENT(lp)) THEN
            IF(atoms%gfelem(i_gf)%lp.NE.lp) CYCLE
         ELSE
            IF(atoms%gfelem(i_gf)%lp.NE.l) CYCLE
         ENDIF
         IF(PRESENT(nTypep)) THEN
            IF(atoms%gfelem(i_gf)%atomTypep.NE.nTypep) CYCLE
         ELSE
            IF(atoms%gfelem(i_gf)%atomTypep.NE.nType) CYCLE
         ENDIF
         !If we are here we found the element
         search = .FALSE.
      ENDDO
   END FUNCTION ind_greensf

END MODULE m_ind_greensf