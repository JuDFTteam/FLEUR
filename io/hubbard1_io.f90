MODULE m_hubbard1_io
   
   !------------------------------------------------------------------------------
   !
   ! MODULE: m_hubbard1_io
   !
   !> @author
   !> Henning Janßen
   !
   ! DESCRIPTION: 
   !>  This module provides an interface with the Hubbard 1 Solver written by 
   !>  J. Kolorenč
   !
   ! REVISION HISTORY:
   ! 20 03 2019 - Initial Version
   ! TODO:  Replace Stop by calls to juDFT_error
   !------------------------------------------------------------------------------

   USE m_juDFT
   USE m_types

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE write_hubbard1_input(path,l,f0,f2,f4,f6,xi,bz,n_min,n_max,beta,mu)

      CHARACTER(len=*), INTENT(IN)  :: path
      INTEGER,          INTENT(IN)  :: l
      REAL,             INTENT(IN)  :: f0, f2, f4, f6
      REAL,             INTENT(IN)  :: xi
      REAL,             INTENT(IN)  :: bz 
      INTEGER,          INTENT(IN)  :: n_min,n_max
      REAL,             INTENT(IN)  :: beta
      REAL,             INTENT(IN)  :: mu 

      CHARACTER(len=300) :: fname
      INTEGER :: info, io_error,io_unit

      io_unit = 17

      fname = "hubbard1.cfg"

      OPEN(unit=io_unit, file=TRIM(ADJUSTL(path)) // TRIM(ADJUSTL(fname)), status="replace", action="write", iostat=io_error)

      IF(io_error.NE.0) CALL juDFT_error("IO-Error in Hubbard 1 IO", calledby="write_hubbard1_input")

      WRITE(io_unit,"(A)") "#**********************************************************",&
                           "#  Parameters for the atomic Hamiltonian"                   ,&                   
                           "#**********************************************************"
      
      WRITE(io_unit,"(A)") "#  Orbital quantum number"
      WRITE(io_unit,9010)  l
9010  FORMAT(TR3,'Lorb',I4.1)

      WRITE(io_unit,"(A)") "#  Slater Integrals"
      SELECT CASE(l)
      CASE(3) 
         WRITE(io_unit,9020) f0, f2, f4, f6
9020     FORMAT(TR3,'Fk',TR2,4f7.2)
      CASE(2)
         WRITE(io_unit,9030) f0, f2, f4, f6
9030     FORMAT(TR3,'Fk',TR2,3f7.2)
      END SELECT

      WRITE(io_unit,"(A)") "#  Spin-orbit-coupling parameter"
      WRITE(io_unit,9040) xi  
9040  FORMAT(TR3,'gfact',f7.3)

      WRITE(io_unit,"(A)") "#  External field"
      WRITE(io_unit,9050) bz 
9050  FORMAT(TR3,'Bz',f15.8)

      WRITE(io_unit,"(A)") "#  Chemical potential"
      WRITE(io_unit,9060) mu
9060  FORMAT(TR3,'mu',f15.8)

      WRITE(io_unit,"(A)") "#**********************************************************",&
                           "#  Parameters for the Solver"                               ,&                   
                           "#**********************************************************"
      WRITE(io_unit,"(A)") "#  Minimum and maximum occupation of the orbital"
      WRITE(io_unit,9070) n_min
9070  FORMAT(TR3,'Nap_min',I4.1)
      WRITE(io_unit,9080) n_max
9080  FORMAT(TR3,'Nap_max',I4.1)

      WRITE(io_unit,"(A)") "#  Setting the solver to use the power lanczos method"
      WRITE(io_unit,"(TR3,A12)") "method_lancz"

      WRITE(io_unit,"(A)") "#  Number of iterations"
      WRITE(io_unit,9090)  100
9090  FORMAT(TR3,'N_lancz_iter',I4.1)

      WRITE(io_unit,"(A)") "#  Number of eigenstates calculated"
      WRITE(io_unit,9100) 35
9100  FORMAT(TR3,'N_lancz_states',I3.1)

      WRITE(io_unit,"(A)") "#**********************************************************",&
                           "#  Parameters for the frequency/energy axis"                ,&                   
                           "#**********************************************************"
      WRITE(io_unit,*) "real_freq_axis {"
      WRITE(io_unit,*) "omegamin", -20.0
      WRITE(io_unit,*) "omegamax", 20.0
      WRITE(io_unit,*) "Nomega",    2000
      WRITE(io_unit,*) "eps",      1.0e-2
      WRITE(io_unit,*) "}"

      WRITE(io_unit,*) "matsub_freq_axis {"
      WRITE(io_unit,*) "Nmatsub",     200
      WRITE(io_unit,*) "}"

      CLOSE(unit=io_unit)
   END SUBROUTINE write_hubbard1_input

   SUBROUTINE write_onsite_gf(fname,gOnsite,i_gf)

      !writes out the onsite green's function calculated from the KS-eigenstates

      TYPE(t_greensf),  INTENT(IN)  :: gOnsite
      CHARACTER(len=*), INTENT(IN)  :: fname
      INTEGER,          INTENT(IN)  :: i_gf

      INTEGER  io_unit,io_error
      INTEGER  iz

      io_unit = 17

      OPEN(unit=io_unit, file=TRIM(ADJUSTL(fname)), status="replace", action="write", iostat=io_error)

      IF(io_error.NE.0) CALL juDFT_error("IO-Error in Hubbard 1 IO", calledby="write_onsite_gf")

      DO iz = 1, gOnsite%nz
         !Write out energy
         WRITE(io_unit,9020) gOnsite%e(iz)
         !WRITE(io_unit,"(A)") "Spin up"
         !WRITE(io_unit,"(A)") "   Real part"
         !WRITE(io_unit,9010)  REAL(gOnsite%gmmpMat(1,iz,i_gf,:,:,1))
         !WRITE(io_unit,"(A)") "   Imaginary part"
         !WRITE(io_unit,9010)  AIMAG(gOnsite%gmmpMat(1,iz,i_gf,:,:,1))
         !WRITE(io_unit,"(A)") "Spin down"
         !WRITE(io_unit,"(A)") "   Real part"
         !WRITE(io_unit,9010)  REAL(gOnsite%gmmpMat(1,iz,i_gf,:,:,2))
         !WRITE(io_unit,"(A)") "   Imaginary part"
         !WRITE(io_unit,9010)  AIMAG(gOnsite%gmmpMat(1,iz,i_gf,:,:,2))
      ENDDO
      CLOSE(unit=io_unit)

9010  FORMAT(7f14.8)
9020  FORMAT("Energy:"2f14.8)
   END SUBROUTINE write_onsite_gf

   SUBROUTINE read_selfen(path,selfen_re,selfen_im,ne,l,e)
      
      !This Subroutine reads in the self-energy
      !produced by the hubbard 1 solver

      REAL, ALLOCATABLE,   INTENT(OUT) :: selfen_re(:,:,:)
      REAL, ALLOCATABLE,   INTENT(OUT) :: selfen_im(:,:,:)
      CHARACTER(len=*),    INTENT(IN)  :: path
      INTEGER,             INTENT(IN)  :: ne
      INTEGER,             INTENT(IN)  :: l
      REAL, ALLOCATABLE,   INTENT(OUT) :: e(:)    
      
      INTEGER io_error,io_unit
      INTEGER n,m,i
      INTEGER matsize
      io_unit = 17
      !Open the selfenergy file
      OPEN(unit=io_unit, file=TRIM(ADJUSTL(path)) // "se.atom", status="replace", action="read", iostat=io_error)

      IF(io_error.NE.0) CALL juDFT_error("IO-Error in reading the self-energy", calledby="read_selfen")
          
      matsize = 2*(2*l+1)
      ALLOCATE(selfen_re(ne,matsize,matsize))
      ALLOCATE(selfen_im(ne,matsize,matsize))
      ALLOCATE(e(ne))
      
      DO i = 1, ne

         READ(io_unit,9010) e(i)
         READ(io_unit,9020) ((selfen_re(i,m,n), m= 1, matsize), n= 1, matsize)
         READ(io_unit,9020) ((selfen_im(i,m,n), m= 1, matsize), n= 1, matsize)

      ENDDO

9010  FORMAT(f10.5)
9020  FORMAT(7f10.5)
   END SUBROUTINE


END MODULE m_hubbard1_io

