MODULE m_types_greensfContourData
   USE m_juDFT
   USE m_types_gfinp
   USE m_constants
   IMPLICIT NONE
   PRIVATE

   TYPE t_greensfContourData
      !This type is used in types_greensf for the storage of the energy points and weights
      !It is defined here, because we want to use it in eContour_gfinp, without importing types_greensf
      INTEGER :: nz = 0      !Number of points
      COMPLEX, ALLOCATABLE :: e(:)  !energy points
      COMPLEX, ALLOCATABLE :: de(:) !integration weights

   CONTAINS
      PROCEDURE,PASS :: init       => init_greensfContourData
      PROCEDURE      :: eContour   => eContour_greensfContourData
      PROCEDURE      :: mpi_bc     => mpi_bc_greensfContourData
   END TYPE t_greensfContourData

   PUBLIC t_greensfContourData

   CONTAINS

   SUBROUTINE init_greensfContourData(this,contourInp,contour_in)

      CLASS(t_greensfContourData),           INTENT(INOUT)  :: this
      TYPE(t_contourInp),                    INTENT(IN)     :: contourInp
      TYPE(t_greensfContourData), OPTIONAL,  INTENT(IN)     :: contour_in
      !
      !Setting up parameters for the energy contour if one was passed
      !
      IF(PRESENT(contour_in)) THEN
         this%nz = contour_in%nz
         ALLOCATE(this%e(this%nz),source=cmplx_0)
         ALLOCATE(this%de(this%nz),source=cmplx_0)

         this%e  = contour_in%e
         this%de = contour_in%de
      ELSE
         SELECT CASE(contourInp%shape)
         CASE(CONTOUR_RECTANGLE_CONST)
            this%nz = contourInp%n1+contourInp%n2+contourInp%n3+contourInp%nmatsub
         CASE(CONTOUR_SEMICIRCLE_CONST)
            this%nz = contourInp%ncirc
         CASE(CONTOUR_DOS_CONST)
            this%nz = contourInp%nDOS
         CASE DEFAULT
            CALL juDFT_error("No valid energy contour mode",calledby="init_contourDataType")
         END SELECT

         ALLOCATE(this%e(this%nz),source=cmplx_0)
         ALLOCATE(this%de(this%nz),source=cmplx_0)
      ENDIF

   END SUBROUTINE init_greensfContourData

   SUBROUTINE mpi_bc_greensfContourData(this,mpi_comm,irank)
         USE m_mpi_bc_tool
         CLASS(t_greensfContourData), INTENT(INOUT)::this
         INTEGER, INTENT(IN):: mpi_comm
         INTEGER, INTENT(IN), OPTIONAL::irank
         INTEGER ::rank
         IF (PRESENT(irank)) THEN
            rank = irank
         ELSE
            rank = 0
         END IF

         CALL mpi_bc(this%nz,rank,mpi_comm)
         CALL mpi_bc(this%e,rank,mpi_comm)
         CALL mpi_bc(this%de,rank,mpi_comm)

   END SUBROUTINE mpi_bc_greensfContourData

   SUBROUTINE eContour_greensfContourData(this,contourInp,ef,irank)

      USE m_grule

      !Calculates the complex energy contour and
      !writes it into the corresponding arrays in gf

      CLASS(t_greensfContourData), INTENT(INOUT) :: this
      TYPE(t_contourInp),          INTENT(IN)    :: contourInp
      REAL,                        INTENT(IN)    :: ef
      INTEGER,                     INTENT(IN)    :: irank

      INTEGER iz,n
      REAL e1, e2, sigma
      COMPLEX del
      REAL r, xr, expo, ff
      REAL, ALLOCATABLE :: x(:), w(:)

      !Help arrays
      ALLOCATE(x(this%nz),source=0.0)
      ALLOCATE(w(this%nz),source=0.0)


      !Transform from relative to ef to absolute
      e1 = ef+contourInp%eb
      e2 = ef+contourInp%et

      SELECT CASE(contourInp%shape)

      CASE(CONTOUR_RECTANGLE_CONST)
         sigma = contourInp%sigma * pi_const
         IF(contourInp%nmatsub > 0) THEN
            n = 0

            !Left Vertical part (e1,0) -> (e1,sigma)
            del = contourInp%nmatsub * CMPLX(0.0,sigma)
            CALL grule(contourInp%n1,x(1:(contourInp%n1)/2),w(1:(contourInp%n1)/2))
            x = -x
            DO iz = 1, (contourInp%n1+3)/2-1
               x(contourInp%n1-iz+1) = -x(iz)
               w(contourInp%n1-iz+1) =  w(iz)
            ENDDO
            DO iz = 1, contourInp%n1
               n = n + 1
               IF(n.GT.this%nz) CALL juDFT_error("Dimension error in energy mesh",calledby="eContour_gfinp")
               this%e(n)  = e1 + del + del * x(iz)
               this%de(n) = w(iz)*del
            ENDDO

            !Horizontal Part (eb,sigma) -> (et,sigma)
            del = (ef-30*contourInp%sigma-e1)/2.0
            CALL grule(contourInp%n2,x(1:(contourInp%n2)/2),w(1:(contourInp%n2)/2))
            x = -x
            DO iz = 1, (contourInp%n2+3)/2-1
               x(contourInp%n2-iz+1) = -x(iz)
               w(contourInp%n2-iz+1) =  w(iz)
            ENDDO
            DO iz = 1, contourInp%n2
               n = n + 1
               IF(n.GT.this%nz) CALL juDFT_error("Dimension error in energy mesh",calledby="eContour_gfinp")
               this%e(n)  = del*x(iz) + del + e1 + 2 * contourInp%nmatsub * ImagUnit * sigma
               this%de(n) = del*w(iz)
            ENDDO

            !Right Vertical part (et,sigma) -> infty
            CALL grule(contourInp%n3,x(1:(contourInp%n3)/2),w(1:(contourInp%n3)/2))
            x = -x
            DO iz = 1, (contourInp%n3+3)/2-1
               x(contourInp%n3-iz+1) = -x(iz)
               w(contourInp%n3-iz+1) =  w(iz)
            ENDDO
            del = 30*contourInp%sigma
            DO iz = 1, contourInp%n3
               n = n + 1
               IF(n.GT.this%nz) CALL juDFT_error("Dimension error in energy mesh",calledby="eContour_gfinp")
               this%e(n)  = del*x(iz)+ef +  2 * contourInp%nmatsub * ImagUnit * sigma
               expo = -ABS(REAL(this%e(n))-ef)/contourInp%sigma
               expo = EXP(expo)
               IF(REAL(this%e(n))<ef) THEN
                  ff = 1.0/(expo+1.0)
               ELSE
                  ff = expo/(expo+1.0)
               ENDIF
               this%de(n) = w(iz)*del * ff
            ENDDO

            !Matsubara frequencies
            DO iz = contourInp%nmatsub , 1, -1
               n = n + 1
               IF(n.GT.this%nz) CALL juDFT_error("Dimension error in energy mesh",calledby="eContour_gfinp")
               this%e(n)  = ef + (2*iz-1) * ImagUnit *sigma
               this%de(n) = -2 * ImagUnit * sigma
            ENDDO
         ENDIF
      CASE(CONTOUR_SEMICIRCLE_CONST)

         !Radius
         r  = (e2-e1)*0.5
         !midpoint
         xr = (e2+e1)*0.5

         CALL grule(contourInp%ncirc,x(1:(contourInp%ncirc)/2),w(1:(contourInp%ncirc)/2))

         DO iz = 1, contourInp%ncirc/2
            x(contourInp%ncirc-iz+1) = -x(iz)
            w(contourInp%ncirc-iz+1) =  w(iz)
         ENDDO
         DO iz = 1, contourInp%ncirc
            this%e(iz)  = xr + ImagUnit * r * EXP(ImagUnit*pi_const/2.0 * x(iz))
            this%de(iz) = pi_const/2.0 * r * w(iz) * EXP(ImagUnit*pi_const/2.0 * x(iz))
            !Scale the imaginary part with the given factor alpha
            this%e(iz)  = REAL(this%e(iz))  + ImagUnit * contourInp%alpha * AIMAG(this%e(iz))
            this%de(iz) = REAL(this%de(iz)) + ImagUnit * contourInp%alpha * AIMAG(this%de(iz))
         ENDDO

      CASE(CONTOUR_DOS_CONST)

         !Equidistant contour (without vertical edges)
         del = (contourInp%et-contourInp%eb)/REAL(this%nz-1)
         DO iz = 1, this%nz
            this%e(iz) = (iz-1) * del + e1 + ImagUnit * contourInp%sigmaDOS
            IF(contourInp%l_dosfermi) THEN
               expo = -ABS(REAL(this%e(iz))-ef)/contourInp%sigmaDOS
               expo = EXP(expo)
               IF(REAL(this%e(iz))<ef) THEN
                  ff = 1.0/(expo+1.0)
               ELSE
                  ff = expo/(expo+1.0)
               ENDIF
               this%de(iz) = del * ff
            ELSE
               this%de(iz) = del
            ENDIF
         ENDDO

         !Not really important but for trapezian method
         !the weight is half at the edges
         this%de(1) = this%de(1)/2.0
         this%de(this%nz) = this%de(this%nz)/2.0

      CASE DEFAULT
         CALL juDFT_error("Invalid mode for energy contour in Green's function calculation", calledby="eContour_gfinp")
      END SELECT

      IF(irank.EQ.0) THEN
         !Write out the information about the energy contour
         WRITE(oUnit,"(A)") "---------------------------------------------"
         WRITE(oUnit,"(A)") " Green's function energy contour"
         WRITE(oUnit,"(A)") "---------------------------------------------"
         WRITE(oUnit,999)  TRIM(ADJUSTL(contourInp%label))
         WRITE(oUnit,1000) contourInp%shape

         SELECT CASE(contourInp%shape)

         CASE(CONTOUR_RECTANGLE_CONST)
            WRITE(oUnit,"(A)") "Rectangular Contour: "
            WRITE(oUnit,1010) this%nz, contourInp%nmatsub,contourInp%n1,contourInp%n2,contourInp%n3
            WRITE(oUnit,"(A)") "Energy limits (rel. to fermi energy): "
            WRITE(oUnit,1040) contourInp%eb,0.0
         CASE(CONTOUR_SEMICIRCLE_CONST)
            WRITE(oUnit,"(A)") "Semicircle Contour: "
            WRITE(oUnit,1020) this%nz, contourInp%alpha
            WRITE(oUnit,"(A)") "Energy limits (rel. to fermi energy): "
            WRITE(oUnit,1040) contourInp%eb,contourInp%et
         CASE(CONTOUR_DOS_CONST)
            WRITE(oUnit,"(A)") "Equidistant Contour for DOS calculations: "
            WRITE(oUnit,1030) this%nz, contourInp%sigmaDOS
            WRITE(oUnit,"(A)") "Energy limits (rel. to fermi energy): "
            WRITE(oUnit,1040) contourInp%eb,contourInp%et
         CASE DEFAULT

         END SELECT

         !Write out points and weights
         WRITE(oUnit,*)
         WRITE(oUnit,"(A)") " Energy points: "
         WRITE(oUnit,"(A)") "---------------------------------------------"
         DO iz = 1, this%nz
            WRITE(oUnit,1050) REAL(this%e(iz)), AIMAG(this%e(iz)), REAL(this%de(iz)), AIMAG(this%de(iz))
         ENDDO

999      FORMAT("Name of the contour:", A)
1000     FORMAT("Using energy contour mode: ", I1,/)
1010     FORMAT("nz: ", I5.1,"; nmatsub: ", I5.1,"; n1: ", I5.1,"; n2: ", I5.1,"; n3: ", I5.1)
1020     FORMAT("nz: ", I5.1," alpha: ", f8.4)
1030     FORMAT("n: ", I5.1,"; sigma: ", f8.4)
1040     FORMAT("eb: ", f8.4,"; et: ",f8.4)
1050     FORMAT(2f8.4,"      weight: ",2e15.4)
      ENDIF

   END SUBROUTINE eContour_greensfContourData


END MODULE m_types_greensfContourData
