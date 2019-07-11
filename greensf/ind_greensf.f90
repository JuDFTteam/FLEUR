MODULE m_ind_greensf


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
         IF(atoms%gfelem(i_gf)%l.EQ.l.AND.atoms%gfelem(i_gf)%atomType.EQ.nType) THEN
            IF(PRESENT(lp).OR.PRESENT(nTypep)) THEN
               !Check second l-argument
               IF(PRESENT(lp)) THEN
                  IF(atoms%gfelem(i_gf)%lp.EQ.lp) THEN
                     search = .FALSE.
                  ELSE
                     search = .TRUE.
                  ENDIF
               ELSE IF(atoms%gfelem(i_gf)%lp.EQ.l) THEN
                  search = .FALSE.
               ELSE
                  search = .TRUE.
               ENDIF
               !check second type argument
               IF(PRESENT(nTypep)) THEN
                  IF(atoms%gfelem(i_gf)%atomTypep.EQ.nTypep) THEN
                     search = .FALSE.
                  ELSE
                     search = .TRUE.
                  ENDIF
               ELSE IF(atoms%gfelem(i_gf)%atomTypep.EQ.nType) THEN
                  search = .FALSE.
               ELSE
                  search = .TRUE.
               ENDIF
            ELSE IF(atoms%gfelem(i_gf)%lp.EQ.l.AND.atoms%gfelem(i_gf)%atomTypep.EQ.nType) THEN
               search = .FALSE.
            ENDIF
         ENDIF
         IF(search.AND.i_gf.EQ.atoms%n_gf) THEN
            CALL juDFT_error("Greens function element not found", calledby="ind_greensf")
         ENDIF
      ENDDO 
   END FUNCTION ind_greensf

END MODULE m_ind_greensf