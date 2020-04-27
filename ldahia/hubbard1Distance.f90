MODULE m_hubbard1Distance

   USE m_types
   USE m_constants

   IMPLICIT NONE

   CONTAINS
   SUBROUTINE hubbard1Distance(n_mmp_in,n_mmp_out,input,gfinp,results)

      TYPE(t_input),       INTENT(IN)     :: input
      TYPE(t_gfinp),       INTENT(IN)     :: gfinp
      COMPLEX,             INTENT(IN)     :: n_mmp_in(-lmaxU_const:,-lmaxU_const:,:)
      COMPLEX,             INTENT(IN)     :: n_mmp_out(-lmaxU_const:,-lmaxU_const:,:)
      TYPE(t_results),     INTENT(INOUT)  :: results

      REAL    :: elementDistance(MERGE(3,input%jspins,gfinp%l_mperp))
      INTEGER :: ispin,j,k
      REAL    :: n_in,n_out

      !Calculates the distance for two density matrices (maximum distance between two elements)
      n_out = 0.0
      n_in = 0.0
      elementDistance = cmplx_0
      DO ispin = 1, MERGE(3,input%jspins,gfinp%l_mperp)
         DO j = -lmaxU_const,lmaxU_const
            DO k = -lmaxU_const,lmaxU_const
               !Distance between individual elements
               IF((ABS(n_mmp_out(k,j,ispin) - n_mmp_in(k,j,ispin))).GT.elementDistance(ispin)) THEN
                  elementDistance(ispin) = ABS(n_mmp_out(k,j,ispin) - n_mmp_in(k,j,ispin))
               ENDIF
               !Distance of the spin diagonal trace
               IF(j.EQ.k.AND.ispin<3) THEN
                  n_out = n_out + REAL(n_mmp_out(k,k,ispin))
                  n_in = n_in + REAL(n_mmp_in(k,k,ispin))
               ENDIF
            END DO
         END DO
      ENDDO
      results%last_occdistance = results%last_occdistance + ABS(n_out-n_in)
      results%last_mmpMatdistance = MAXVAL(elementDistance)

      !IO to out file
      WRITE(oUnit,'(A)') "Hubbard 1 Distances:"
      WRITE(oUnit,'(A)') "-------------------------------"
      WRITE(oUnit,*) "  Occupation distance: ", results%last_occdistance
      WRITE(oUnit,*) "  Density matrix distance: ", results%last_mmpMatdistance
      DO ispin = 1, MERGE(3,input%jspins,gfinp%l_mperp)
         WRITE(oUnit,*)
         WRITE(oUnit,'(A11,I1,A2,f14.8)') "Delta Spin ", ispin, ": ", elementDistance(ispin)
         WRITE(oUnit,'(14f12.7)') n_mmp_out(-lmaxU_const:,-lmaxU_const:,ispin) - n_mmp_in(-lmaxU_const:,-lmaxU_const:,ispin)
      ENDDO

   END SUBROUTINE hubbard1Distance
END MODULE m_hubbard1Distance
