!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hubbard1Distance

   USE m_constants
   USE m_types_misc

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: hubbard1distance

   CONTAINS
   SUBROUTINE hubbard1Distance(n_mmp_in,n_mmp_out,results)

      COMPLEX,             INTENT(IN)     :: n_mmp_in (-lmaxU_const:,-lmaxU_const:,:)
      COMPLEX,             INTENT(IN)     :: n_mmp_out(-lmaxU_const:,-lmaxU_const:,:)
      TYPE(t_results),     INTENT(INOUT)  :: results

      REAL    :: elementDistance(SIZE(n_mmp_in,3))
      INTEGER :: ispin,m,mp
      REAL    :: n_in,n_out

      !Calculates the distance for two density matrices (maximum distance between two elements)
      n_out = 0.0
      n_in = 0.0
      elementDistance = 0.0
      DO ispin = 1, SIZE(n_mmp_in,3)
         DO m = -lmaxU_const,lmaxU_const
            DO mp = -lmaxU_const,lmaxU_const
               !Distance between individual elements
               IF((ABS(n_mmp_out(m,mp,ispin) - n_mmp_in(m,mp,ispin))).GT.elementDistance(ispin)) THEN
                  elementDistance(ispin) = ABS(n_mmp_out(m,mp,ispin) - n_mmp_in(m,mp,ispin))
               ENDIF
               !Distance of the spin diagonal trace
               IF(m.EQ.mp.AND.ispin<3) THEN
                  n_out = n_out + REAL(n_mmp_out(m,m,ispin))
                  n_in  = n_in  + REAL(n_mmp_in (m,m,ispin))
               ENDIF
            END DO
         END DO
      ENDDO
      results%last_occdistance    = MAX(results%last_occdistance   ,ABS(n_out-n_in))
      results%last_mmpMatdistance = MAX(results%last_mmpMatdistance,MAXVAL(elementDistance))

      !IO to out file
      WRITE(oUnit,'(A)') "Hubbard 1 Distances:"
      WRITE(oUnit,'(A)') "-------------------------------"
      WRITE(oUnit,*) "  Occupation distance:     ", ABS(n_out-n_in)
      WRITE(oUnit,*) "  Density matrix distance: ", MAXVAL(elementDistance)
      DO ispin = 1, SIZE(n_mmp_in,3)
         WRITE(oUnit,'(/,A,I1,A,f14.8)') "Delta Spin ", ispin, ": ", elementDistance(ispin)
         WRITE(oUnit,'(7f12.7)') REAL(  n_mmp_out(-lmaxU_const:,-lmaxU_const:,ispin) &
                                      - n_mmp_in (-lmaxU_const:,-lmaxU_const:,ispin)  )
      ENDDO

   END SUBROUTINE hubbard1Distance
END MODULE m_hubbard1Distance
