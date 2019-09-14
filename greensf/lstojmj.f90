MODULE m_lsTOjmj

   USE m_clebsch
   USE m_types
   USE m_juDFT


   IMPLICIT NONE

   CONTAINS

   SUBROUTINE lsTOjmj(cmat,l)

      TYPE(t_mat),      INTENT(INOUT)  :: cmat
      INTEGER,          INTENT(IN)     :: l

      INTEGER jj,j,mj,s,ml,i,k

      IF(.NOT.ALLOCATED(cmat%data_r)) THEN
         CALL cmat%init(.TRUE.,2*(2*l+1),2*(2*l+1))
      ENDIF

      !Calculate the matrix of CG-coefficients to transform from |l,ml,ms> to |j=l\pm1/2,mjz>
      !compare utils/occup.f90 from libedsolver
      cmat%data_r = 0.0
      DO jj = -1, 1, 2
         j = 2*l+jj
         DO mj = 1, j+1
            k = mj+(jj+1)*l 
            DO s = -1, 1, 2
               DO ml = -l, l
                  IF(ml-s*0.5.NE.mj-1-j*0.5) CYCLE
                  !In libedsolver spin up and down are flipped
                  i = ml+l+1+(1+s)/2.0*(2*l+1)
                  !The minus sign in contrast to occup.f90 stems from the fact 
                  !that the spins are stored in reversed order in the solver
                  cmat%data_r(i,k) = clebsch(1.0*l,0.5,1.0*ml,-s*0.5,j*0.5,mj-1-j*0.5)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE lsTOjmj


END MODULE m_lsTOjmj