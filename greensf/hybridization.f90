MODULE m_hybridization


   USE m_types
   USE m_constants
   USE m_juDFT

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE hybridization(gf,i_gf,gfinp,input,ef)

      !------------------------------------------------------
      ! Evaluates the hybridization function
      ! acc. to. sci. rep. 5, 15429 (2015)
      !------------------------------------------------------
      ! \Delta(E) = -1/(pi*N_l) Im TR[G^-1_{DFT}(E+i\delta)]
      !------------------------------------------------------

      INTEGER,         INTENT(IN) :: i_gf
      TYPE(t_greensf), INTENT(IN) :: gf
      TYPE(t_gfinp),   INTENT(IN) :: gfinp
      TYPE(t_input),   INTENT(IN) :: input
      REAL,            INTENT(IN) :: ef

      TYPE(t_mat) :: gmat
      INTEGER io_error,iz,i,ipm
      REAL v_low,v_high,nf,avg_delta,ellow,elup
      COMPLEX tr
      REAL Delta(gf%contour%nz) !Hybridization function
      INTEGER contourShape

      !Open file
      OPEN(unit=1337,file="hybridization.dat",status="replace",action="write",iostat=io_error)
      IF(io_error.NE.0) CALL juDFT_error("IO error",calledby="hybridization")

      contourShape = gfinp%contour(gfinp%elem(i_gf)%iContour)%shape

      Delta = 0.0
      DO iz = 1, gf%contour%nz
         tr = 0.0
         DO ipm = 1, 2
            !--------------------------------------------------
            ! Get the full Greens function matrix for the current energy point
            !--------------------------------------------------
            CALL gf%get(i_gf,gmat,gfinp,input,iz,ipm.EQ.2)
            !--------------------------------------------------
            ! Invert the matrix using the routines in types_mat
            !--------------------------------------------------
            CALL gmat%inverse()
            !Compute the trace
            DO i = 1, gmat%matsize1
               tr = tr + (-1)**(ipm-1) * gmat%data_c(i,i)
            ENDDO
         ENDDO
         Delta(iz) = -1/(tpi_const*gmat%matsize1) * AIMAG(tr)
         WRITE(1337,"(2f14.8)") REAL(gf%contour%e(iz)-ef)*hartree_to_ev_const, Delta(iz)
         !Free up the gmat matrix (it is initialized in gf%get)
         CALL gmat%free()
      ENDDO

      CLOSE(unit=1337,iostat=io_error)
      IF(io_error.NE.0) CALL juDFT_error("IO error",calledby="hybridization")

      IF(contourShape.EQ.CONTOUR_DOS_CONST) THEN
         !-----------------------------------------------------
         ! Try to fit a bath state to this configuration
         !-----------------------------------------------------
         !Average over this interval (relative to e_fermi)
         ellow = -0.5/hartree_to_ev_const
         elup  = +0.5/hartree_to_ev_const
         !-----------------------------------------------------
         ! Average the hybridization function over this interval
         !-----------------------------------------------------
         avg_delta = 0.0
         DO iz = 1, gf%contour%nz
            IF(REAL(gf%contour%e(iz)).LT.ef+ellow.OR.REAL(gf%contour%e(iz)).GT.ef+elup) CYCLE
            avg_delta = avg_delta + Delta(iz)*hartree_to_ev_const*REAL(gf%contour%de(iz))
         ENDDO
         WRITE(*,*) avg_delta
         !low J
         nf = 2*gfinp%elem(i_gf)%l
         v_low = sqrt(-avg_delta/(nf))
         !high J
         nf = 2*gfinp%elem(i_gf)%l+2
         v_high = sqrt(-avg_delta/(nf))
         WRITE(*,*) v_low,v_high
      ENDIF


   END SUBROUTINE hybridization

END MODULE m_hybridization