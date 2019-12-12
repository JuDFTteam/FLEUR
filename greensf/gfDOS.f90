MODULE m_gfDOS

   USE m_types
   USE m_lsTOjmj
   USE m_constants

   IMPLICIT NONE

   CONTAINS


   SUBROUTINE gfDOS(g,l,nType,jobID,atoms,input,ef)

      !Provides the spin up/down DOS as well as the high/low J DOS
      !calculated from the greens function in gfDOS.jobID

      TYPE(t_greensf),     INTENT(IN)  :: g
      TYPE(t_atoms),       INTENT(IN)  :: atoms
      TYPE(t_input),       INTENT(IN)  :: input
      INTEGER,             INTENT(IN)  :: l
      INTEGER,             INTENT(IN)  :: nType
      INTEGER,             INTENT(IN)  :: jobID
      REAL, OPTIONAL,      INTENT(IN)  :: ef

      INTEGER iz,ipm,i,ns
      INTEGER io_error
      COMPLEX dos(4,2*l+2),re(2) !up,down,low,high (only at the current energy point)
      TYPE(t_mat) :: gmat,cmat,jmat
      CHARACTER(len=10) :: filename

9000  FORMAT("gfDOS.",I3.3)
9001  FORMAT("mgfDOS.",I3.3)
      WRITE(filename,9000) jobID

      OPEN(unit=3456,file=filename,status="replace",action="write",iostat=io_error)
      WRITE(filename,9001) jobID
      OPEN(unit=3457,file=filename,status="replace",action="write",iostat=io_error)
      IF(io_error.NE.0) CALL juDFT_error("IO-error",calledby="gfDOS")
      !Write out warnings
      IF(.NOT.PRESENT(ef)) WRITE(3456,"(A)") "This gfDOS is not corrected to have ef=0"
      IF(g%mode.NE.3) WRITE(3456,"(A)") "You are using an energy contour where the gfDOS might not be very informative"
      !Calculate the transformation between |L,S> and |J,mj> basis
      ns = 2*l+1
      CALL cmat%init(.TRUE.,2*ns,2*ns)
      CALL jmat%init(.FALSE.,2*ns,2*ns)
      CALL lsTOjmj(cmat,l)

      DO iz = 1, g%nz
         dos = 0.0
         re = 0.0
         DO ipm = 1, 2  !Sum over G^+ and G^-
            !Get the full gf matrix at the energy point
            CALL g%get(gmat,atoms,input,iz,l,nType,ipm.EQ.2)
            !Convert to eV^-1
            gmat%data_c = gmat%data_c/hartree_to_eV_const
            !Calculate up/down dos
            DO i = 1, ns
               dos(1,i) = dos(1,i) - 1.0/tpi_const * (-1)**(ipm-1) * gmat%data_c(i,i)
            ENDDO
            DO i = ns+1, 2*ns
               dos(2,i-ns) = dos(2,i-ns) - 1.0/tpi_const * (-1)**(ipm-1) * gmat%data_c(i,i)
            ENDDO
            !Transform to |J,mj> basis
            jmat%data_c = matmul(gmat%data_c,cmat%data_r)
            jmat%data_c = matmul(transpose(cmat%data_r),jmat%data_c)
            !Calculate low/high dos
            DO i = 1, ns-1
               dos(3,i) = dos(3,i) - 1.0/tpi_const * (-1)**(ipm-1) * jmat%data_c(i,i)
            ENDDO
            DO i = ns, 2*ns
               dos(4,i-ns+1) = dos(4,i-ns+1) - 1.0/tpi_const * (-1)**(ipm-1) * jmat%data_c(i,i)
            ENDDO
            !Real part
            DO i = 1, ns
               re(1) = re(1) - 1.0/tpi_const  * gmat%data_c(i,i)
            ENDDO
            DO i = ns+1, 2*ns
               re(2) = re(2) - 1.0/tpi_const* gmat%data_c(i,i)
            ENDDO
            CALL gmat%free()
         ENDDO
         WRITE(3456,"(7f14.8)") (REAL(g%e(iz))-MERGE(ef,0.0,PRESENT(ef)))*hartree_to_eV_const, &
                                SUM(AIMAG(dos(1,1:ns))),SUM(AIMAG(dos(2,1:ns))),&
                                SUM(AIMAG(dos(3,1:ns-1))),SUM(AIMAG(dos(4,1:ns+1))), REAL(re(1)), REAL(re(2))
         WRITE(3457,"(15f10.5)") (REAL(g%e(iz))-MERGE(ef,0.0,PRESENT(ef)))*hartree_to_eV_const, (AIMAG(dos(1,i)),i=1, ns),(AIMAG(dos(2,i)),i=1, ns)


      ENDDO
      CLOSE(3456,iostat=io_error)
      CLOSE(3457)
      IF(io_error.NE.0) CALL juDFT_error("IO-error",calledby="gfDOS")

      CALL cmat%free()
      CALL jmat%free()

   END SUBROUTINE gfDOS

END MODULE m_gfDOS