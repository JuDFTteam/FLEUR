MODULE m_gfDOS

   USE m_types
   USE m_lsTOjmj
   USE m_constants

   IMPLICIT NONE

   CONTAINS


   SUBROUTINE gfDOS(g,l,nType,jobID,atoms,input)

      !Provides the spin up/down DOS as well as the high/low J DOS 
      !calculated from the greens function in gfDOS.jobID

      TYPE(t_greensf),     INTENT(IN)  :: g
      TYPE(t_atoms),       INTENT(IN)  :: atoms
      TYPE(t_input),       INTENT(IN)  :: input
      INTEGER,             INTENT(IN)  :: l
      INTEGER,             INTENT(IN)  :: nType
      INTEGER,             INTENT(IN)  :: jobID

      INTEGER iz,ipm,i,ns
      INTEGER io_error
      COMPLEX dos(4) !up,down,low,high (only at the current energy point)
      TYPE(t_mat) :: gmat,cmat,jmat
      CHARACTER(len=10) :: filename


      IF(g%mode.NE.3) WRITE(*,*) "You are using an energy contour where the gfDOS might not be very informative"

9000  FORMAT("gfDOS."I1)
      WRITE(filename,9000) jobID

      OPEN(unit=3456,file=filename,status="replace",action="write",iostat=io_error)
      IF(io_error) CALL juDFT_error("IO-error",calledby="gfDOS")

      !Calculate the transformation between |L,S> and |J,mj> basis
      ns = 2*l+1
      CALL cmat%init(.TRUE.,2*ns,2*ns)
      CALL jmat%init(.FALSE.,2*ns,2*ns)
      CALL lsTOjmj(cmat,l)

      DO iz = 1, g%nz
         dos = 0.0
         DO ipm = 1, 2  !Sum over G^+ and G^-
            !Get the full gf matrix at the energy point
            CALL g%get_gf(gmat,atoms,input,iz,l,nType,ipm.EQ.2)
            !Calculate up/down dos
            DO i = 1, ns
               dos(1) = dos(1) - 1/(2*pi_const) * (-1)**(ipm-1) * gmat%data_c(i,i)
            ENDDO
            DO i = ns+1, 2*ns
               dos(2) = dos(2) - 1/(2*pi_const) * (-1)**(ipm-1) * gmat%data_c(i,i)
            ENDDO
            !Transform to |J,mj> basis
            jmat%data_c = matmul(gmat%data_c,cmat%data_r)
            jmat%data_c = matmul(transpose(cmat%data_r),jmat%data_c)
            !Calculate low/high dos
            DO i = 1, ns-1
               dos(3) = dos(3) - 1/(2*pi_const) * (-1)**(ipm-1) * jmat%data_c(i,i)
            ENDDO
            DO i = ns, 2*ns
               dos(4) = dos(4) - 1/(2*pi_const) * (-1)**(ipm-1) * jmat%data_c(i,i)
            ENDDO
            CALL gmat%free()
         ENDDO 
         WRITE(3456,"(5f14.8)") REAL(g%e(iz)), (AIMAG(dos(i)),i=1,4)
      ENDDO
      CLOSE(3456,iostat=io_error)
      IF(io_error) CALL juDFT_error("IO-error",calledby="gfDOS")

      CALL cmat%free()
      CALL jmat%free()

   END SUBROUTINE gfDOS

END MODULE m_gfDOS