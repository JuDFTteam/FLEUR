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
   USE m_constants

   IMPLICIT NONE
   !------------------------------------------------------------------
   !Here the keywords for the hubbard 1 solver input file are defined
   !------------------------------------------------------------------

   CHARACTER(len=300), PARAMETER :: input_filename ="hubbard1.cfg"
   INTEGER, PARAMETER            :: input_iounit = 17

   CONTAINS

   SUBROUTINE write_hubbard1_input(path,i_hia,l,f0,f2,f4,f6,hub1,mu,n,ne,nmatsub,e_min,e_max,sigma)

      USE m_generic_txtio

      IMPLICIT NONE

      CHARACTER(len=*), INTENT(IN)  :: path
      INTEGER,          INTENT(IN)  :: i_hia
      INTEGER,          INTENT(IN)  :: l
      REAL,             INTENT(IN)  :: f0,f2,f4,f6
      TYPE(t_hub1ham),  INTENT(IN)  :: hub1
      REAL,             INTENT(IN)  :: mu
      INTEGER,          INTENT(IN)  :: n
      INTEGER,          INTENT(IN)  :: ne
      INTEGER,          INTENT(IN)  :: nmatsub
      REAL,             INTENT(IN)  :: e_min
      REAL,             INTENT(IN)  :: e_max
      REAL,             INTENT(IN)  :: sigma 
      
      INTEGER :: info, io_error

      !Main input file
      OPEN(unit=input_iounit, file=TRIM(ADJUSTL(path)) // TRIM(ADJUSTL(input_filename)),&
          status="replace", action="write", iostat=io_error)
      !IF(io_error.NE.0) CALL juDFT_error("IO-Error in Hubbard 1 IO", calledby="write_hubbard1_input")

      CALL header(input_iounit,"Parameters for the atomic Hamiltonian in eV",1)

      CALL comment(input_iounit,"Orbital quantum number",1)
      CALL writeValue(input_iounit,"Lorb",l)

      CALL comment(input_iounit,"Slater Integrals",1)
      CALL writeValue(input_iounit,"Fk",(/f0,f2,f4,f6/))

      CALL comment(input_iounit,"Spin-orbit-coupling parameter",1)
      CALL writeValue(input_iounit,"gfact",hub1%xi(i_hia))

      CALL comment(input_iounit,"External field",1)
      CALL writeValue(input_iounit,"Bz",hub1%bz(i_hia))

      CALL comment(input_iounit,"Inverse temperature",1)
      CALL writeValue(input_iounit,"beta",hub1%beta)

      CALL comment(input_iounit,"Chemical potential",1)
      CALL writeValue(input_iounit,"mu",mu)

      CALL comment(input_iounit,"Crystal field factor",1)
      CALL writeValue(input_iounit,"ccf",1.0)


      CALL header(input_iounit,"Parameters for the Solver",1)

      CALL comment(input_iounit,"Minimum and maximum occupation of the orbital",1)
      CALL writeValue(input_iounit,"Nap_min",n-hub1%n_exc)
      CALL writeValue(input_iounit,"Nap_max",n+hub1%n_exc)

      CALL comment(input_iounit,"Setting the solver to use the power lanczos method",1)
      WRITE(input_iounit,"(TR3,A12)") "method_lancz"

      CALL comment(input_iounit,"Number of iterations",1)
      CALL writeValue(input_iounit,"N_lancz_iter",100)

      CALL comment(input_iounit,"Number of eigenstates calculated",1)
      CALL writeValue(input_iounit,"N_lancz_states",35)

      CALL comment(input_iounit,"Parameters for the frequency/energy axis",1)
      WRITE(input_iounit,*) "real_freq_axis {"
      WRITE(input_iounit,*) "omegamin", e_min*hartree_to_ev_const
      WRITE(input_iounit,*) "omegamax", e_max*hartree_to_ev_const
      WRITE(input_iounit,*) "Nomega",   ne
      WRITE(input_iounit,*) "eps",      sigma*hartree_to_ev_const
      WRITE(input_iounit,*) "}"

      WRITE(input_iounit,*) "matsub_freq_axis {"
      WRITE(input_iounit,*) "Nmatsub",     nmatsub
      WRITE(input_iounit,*) "}"

      CLOSE(unit=input_iounit)
   END SUBROUTINE write_hubbard1_input

   SUBROUTINE write_gf(app,gOnsite,i_gf)

      !writes out the onsite green's function calculated from the KS-eigenstates

      TYPE(t_greensf),  INTENT(IN)  :: gOnsite
      CHARACTER(len=*), INTENT(IN)  :: app
      INTEGER,          INTENT(IN)  :: i_gf

      INTEGER  io_unit,io_error
      INTEGER  iz,ipm
      CHARACTER(len=2) :: ret

      io_unit = 17
      DO ipm = 1, 2
         IF(ipm.EQ.1) THEN
            ret = "_R"
         ELSE
            ret = "_A"
         ENDIF
         OPEN(unit=io_unit, file=TRIM(ADJUSTL(app)) // TRIM(ADJUSTL(ret)), status="replace", action="write", iostat=io_error)

         IF(io_error.NE.0) CALL juDFT_error("IO-Error in Hubbard 1 IO", calledby="write_onsite_gf")

         DO iz = 1, gOnsite%nz-gOnsite%Nmatsub
            !Write out energy
            WRITE(io_unit,9020) gOnsite%e(iz)
            WRITE(io_unit,"(A)") "Spin up"
            WRITE(io_unit,"(A)") "   Real part"
            WRITE(io_unit,9010)  REAL(gOnsite%gmmpMat(iz,i_gf,:,:,1,ipm))
            WRITE(io_unit,"(A)") "   Imaginary part"
            WRITE(io_unit,9010)  AIMAG(gOnsite%gmmpMat(iz,i_gf,:,:,1,ipm))
            WRITE(io_unit,"(A)") "Spin down"
            WRITE(io_unit,"(A)") "   Real part"
            WRITE(io_unit,9010)  REAL(gOnsite%gmmpMat(iz,i_gf,:,:,2,ipm))
            WRITE(io_unit,"(A)") "   Imaginary part"
            WRITE(io_unit,9010)  AIMAG(gOnsite%gmmpMat(iz,i_gf,:,:,2,ipm))
         ENDDO
         CLOSE(unit=io_unit)
      ENDDO

9010  FORMAT(7f14.8)
9020  FORMAT("Energy:"2f14.8)
   END SUBROUTINE write_gf

   SUBROUTINE write_ccfmat(path,ccfmat,l)

      CHARACTER(len=*), INTENT(IN)  :: path
      REAL,             INTENT(IN)  :: ccfmat(-l:l,-l:l)
      INTEGER,          INTENT(IN)  :: l

      CHARACTER(len=300) :: fname
      INTEGER :: info, io_error,io_unit

      io_unit = 17

      fname = "ccf.dat"

      OPEN(unit=io_unit, file=TRIM(ADJUSTL(path)) // TRIM(ADJUSTL(fname)), status="replace", action="write", iostat=io_error)

      IF(l.EQ.2) THEN
         WRITE(io_unit,"(5f10.6)") ccfmat*hartree_to_ev_const
      ELSE IF(l.EQ.3) THEN
         WRITE(io_unit,"(7f10.6)") ccfmat*hartree_to_ev_const
      ENDIF

      CLOSE(io_unit)

   END SUBROUTINE write_ccfmat

   SUBROUTINE read_selfen(path,selfen,ne,matsize,l_matsub)
      
      USE m_constants
      !This Subroutine reads in the self-energy
      !produced by the hubbard 1 solver

      COMPLEX,          INTENT(OUT) :: selfen(:,:,:)
      CHARACTER(len=*), INTENT(IN)  :: path
      INTEGER,          INTENT(IN)  :: ne
      INTEGER,          INTENT(IN)  :: matsize
      LOGICAL,          INTENT(IN)  :: l_matsub  
      
      INTEGER io_error,io_unit
      INTEGER n,m,i
      REAL tmp(matsize,matsize)
      io_unit = 17
      !Open the selfenergy file
      IF(l_matsub) THEN
         OPEN(unit=io_unit, file=TRIM(ADJUSTL(path)) // "selfen_matsub_bundle.dat",status="old", action="read", iostat=io_error)

         IF(io_error.NE.0) CALL juDFT_error("IO-Error in reading the self-energy", calledby="read_selfen")
         READ(io_unit,*)
         DO i = 1, ne
            READ(io_unit,*) 
            DO m = 1, matsize
               READ(io_unit,*) selfen(i,1:matsize,m)
            ENDDO
         ENDDO
      ELSE
         OPEN(unit=io_unit, file=TRIM(ADJUSTL(path)) // "se.atom", status="old", action="read", iostat=io_error)

         IF(io_error.NE.0) CALL juDFT_error("IO-Error in reading the self-energy", calledby="read_selfen")
         
         DO i = 1, ne
            READ(io_unit,9010) 
            READ(io_unit,9020) ((tmp(m,n), m= 1, matsize), n= 1, matsize)
            selfen(i,1:matsize,1:matsize) = tmp(1:matsize,1:matsize)/hartree_to_ev_const 
            READ(io_unit,9020) ((tmp(m,n), m= 1, matsize), n= 1, matsize)
            selfen(i,1:matsize,1:matsize) = selfen(i,1:matsize,1:matsize) + ImagUnit * tmp(1:matsize,1:matsize)/hartree_to_ev_const 
         ENDDO
      ENDIF


      CLOSE(io_unit)

9010  FORMAT(f10.5)
9020  FORMAT(7f11.5)
   END SUBROUTINE read_selfen


END MODULE m_hubbard1_io

