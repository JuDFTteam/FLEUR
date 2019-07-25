MODULE m_hybridization


   USE m_types
   USE m_constants
   USE m_juDFT

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE hybridization(l,nType,gf,atoms,input)

      !------------------------------------------------------
      ! Evaluates the hybridization function 
      ! acc. to. sci. rep. 5, 15429 (2015)
      !------------------------------------------------------
      ! \Delta(E) = -1/pi Im TR[G^-1_{DFT}(E+i\delta)]
      !------------------------------------------------------

      INTEGER,         INTENT(IN) :: l
      INTEGER,         INTENT(IN) :: nType
      TYPE(t_greensf), INTENT(IN) :: gf 
      TYPE(t_atoms),   INTENT(IN) :: atoms
      TYPE(t_input),   INTENT(IN) :: input

      TYPE(t_mat) :: gmat
      INTEGER io_error,iz,i
      COMPLEX tr

      !Open file
      OPEN(unit=1337,file="hybridization.dat",status="replace",action="write",iostat=io_error)
      IF(io_error.NE.0) CALL juDFT_error("IO error",calledby="hybridization")

      DO iz = 1, gf%nz  
         !Get the full Greens function matrix for the current energy point
         CALL gf%get_gf(gmat,atoms,input,iz,l,nType,.FALSE.)
         !--------------------------------------------------
         !Invert the matrix using the routines in types_mat
         !--------------------------------------------------
         CALL gmat%inverse()
         !Compute the trace
         tr = 0.0
         DO i = 1, gmat%matsize1
            tr = tr + gmat%data_c(i,i)
         ENDDO
         WRITE(1337,"(2f14.8)") REAL(gf%e(iz)), -1/pi_const * AIMAG(tr)
         !Free up the gmat matrix (it is initialized in gf%get_gf)
         CALL gmat%free()
      ENDDO

      CLOSE(unit=1337,iostat=io_error)
      IF(io_error.NE.0) CALL juDFT_error("IO error",calledby="hybridization")

   END SUBROUTINE hybridization

END MODULE m_hybridization