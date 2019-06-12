MODULE m_denmat_dist
   
   CONTAINS
   SUBROUTINE n_mmp_dist(n_mmp_in,n_mmp_out,hub1,results,jspins)

      USE m_types
      USE m_constants

      TYPE(t_hub1ham),     INTENT(IN)     :: hub1
      COMPLEX,             INTENT(IN)     :: n_mmp_in(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,hub1%n_hia,jspins)
      COMPLEX,             INTENT(IN)     :: n_mmp_out(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,hub1%n_hia,jspins)
      TYPE(t_results),     INTENT(INOUT)  :: results
      INTEGER,             INTENT(IN)     :: jspins

      INTEGER ispin,i_hia,j,k
      REAL n_in,n_out
      
      !Calculates the distance for two density matrices (maximum distance between two elements)
      n_out = 0.0
      n_in = 0.0
      results%last_mmpMatdistance = 0.0
      DO ispin = 1, jspins
         DO i_hia = 1, hub1%n_hia
            DO j = -3,3
               DO k = -3,3
                  IF((ABS(n_mmp_out(k,j,i_hia,ispin) - n_mmp_in(k,j,i_hia,ispin))).GT.results%last_mmpMatdistance) THEN
                     results%last_mmpMatdistance = ABS(n_mmp_out(k,j,i_hia,ispin) - n_mmp_in(k,j,i_hia,ispin))
                  ENDIF
                  IF(j.EQ.k) THEN
                     n_out = n_out + REAL(n_mmp_out(k,j,i_hia,ispin))
                     n_in = n_in + REAL(n_mmp_in(k,j,i_hia,ispin))
                  ENDIF
               END DO
            END DO
         END DO
      ENDDO
      results%last_occdistance = ABS(n_out-n_in)
      WRITE(6,*) "Occupation distance: ", results%last_occdistance
      WRITE(6,*) "Density matrix distance: ", results%last_mmpMatdistance
   
   END SUBROUTINE n_mmp_dist
END MODULE m_denmat_dist