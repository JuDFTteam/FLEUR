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

   !Filenames for input
   CHARACTER(*), PARAMETER :: cfg_file_main ="hubbard1.cfg"
   CHARACTER(*), PARAMETER :: cfg_file_bath ="bath.cfg"
   CHARACTER(*), PARAMETER :: cfg_file_hloc ="hloc.cfg"
   CHARACTER(*), PARAMETER :: cfg_file_ccf = "ccf.dat"
   CHARACTER(*), PARAMETER :: file_G0_fits  ="G0_fit_monitor.dat"
   CHARACTER(*), PARAMETER :: file_hybr_fits ="hyb_fit_monitor.dat"
   INTEGER, PARAMETER      :: input_iounit  = 17

   !Real freq axis parameters
   REAL, PARAMETER    :: emin = -13.0
   REAL, PARAMETER    :: emax =  13.0
   INTEGER, PARAMETER :: ne   =  2600
   REAL, PARAMETER    :: sigma = 0.0314
   INTEGER, PARAMETER :: nmats = 0


   CONTAINS

   SUBROUTINE hubbard1_input(path,i_hia,l,f0,f2,f4,f6,hub1,mu,n,l_bath,l_new)

      IMPLICIT NONE

      CHARACTER(len=*), INTENT(IN)  :: path
      INTEGER,          INTENT(IN)  :: i_hia
      INTEGER,          INTENT(IN)  :: l
      REAL,             INTENT(IN)  :: f0,f2,f4,f6
      TYPE(t_hub1ham),  INTENT(IN)  :: hub1
      REAL,             INTENT(IN)  :: mu
      INTEGER,          INTENT(IN)  :: n
      LOGICAL,          INTENT(IN)  :: l_bath
      LOGICAL,          INTENT(IN)  :: l_new

      !Old or new input format
      IF(l_new) THEN
         CALL write_hubbard1_input_new(path,i_hia,l,f0,f2,f4,f6,hub1,mu,n,l_bath)
      ELSE
         CALL write_hubbard1_input_old(path,i_hia,l,f0,f2,f4,f6,hub1,mu,n)
      ENDIF

   END SUBROUTINE hubbard1_input


   SUBROUTINE write_hubbard1_input_new(path,i_hia,l,f0,f2,f4,f6,hub1,mu,n,l_bath)

      USE m_generic_txtio

      IMPLICIT NONE

      CHARACTER(len=*), INTENT(IN)  :: path
      INTEGER,          INTENT(IN)  :: i_hia
      INTEGER,          INTENT(IN)  :: l
      REAL,             INTENT(IN)  :: f0,f2,f4,f6
      TYPE(t_hub1ham),  INTENT(IN)  :: hub1
      REAL,             INTENT(IN)  :: mu
      INTEGER,          INTENT(IN)  :: n
      LOGICAL,          INTENT(IN)  :: l_bath
      
      INTEGER :: info, io_error,i,j,k,ind1,ind2,i_exc,i_arg
      REAL exc
      TYPE(t_mat) :: cfmat

      !Main input file
      OPEN(unit=input_iounit, file=TRIM(ADJUSTL(path)) // TRIM(ADJUSTL(cfg_file_main)),&
          status="replace", action="write", iostat=io_error)
      IF(io_error.NE.0) CALL juDFT_error("IO-Error in Hubbard 1 IO", calledby="write_hubbard1_input_new")

      CALL startSection(input_iounit,"hamiltonian")
      CALL comment(input_iounit,"Slater Integrals",1)
      CALL writeValue(input_iounit,"Fk",(/f0,f2,f4,f6/))
      CALL writeValue(input_iounit, "include", cfg_file_hloc)
      IF(l_bath) CALL writeValue(input_iounit, "include", cfg_file_bath)
      CALL endSection(input_iounit)

      CALL startSection(input_iounit,"fock_space")
      CALL comment(input_iounit,"Min/Max Occupation",1)
      CALL writeValue(input_iounit,"Np_min",MAX(0,n-hub1%n_exc))
      CALL writeValue(input_iounit,"Np_max",MIN(2*(2*l+1),n+hub1%n_exc))
      CALL endSection(input_iounit)

      CALL startSection(input_iounit,"GC_ensemble")
      CALL comment(input_iounit,"Inverse temperature",1)
      CALL writeValue(input_iounit,"beta",hub1%beta)
      CALL comment(input_iounit,"States with smaller weight are dropped",1)
      CALL writeValue(input_iounit, "weight_limit",1.0e-4)
      CALL endSection(input_iounit)

      CALL startSection(input_iounit,"method")
      CALL writeValue(input_iounit, "lancz")
      CALL comment(input_iounit,"Number of iterations",1)
      CALL writeValue(input_iounit,"N_lancz_iter",100)
      CALL comment(input_iounit,"Number of eigenstates calculated",1)
      CALL writeValue(input_iounit,"N_lancz_states",80)
      CALL endSection(input_iounit)

      CALL startSection(input_iounit,"real_freq_axis")
      CALL writeValue(input_iounit, "omegamin", emin)
      CALL writeValue(input_iounit, "omegamax", emax)
      CALL writeValue(input_iounit, "Nomega", ne)
      CALL writeValue(input_iounit, "eps", sigma)
      CALL endSection(input_iounit)

      CLOSE(input_iounit)


      !local hamiltonian
      OPEN(unit=input_iounit, file=TRIM(ADJUSTL(path)) // TRIM(ADJUSTL(cfg_file_hloc)),&
      status="replace", action="write", iostat=io_error)
      IF(io_error.NE.0) CALL juDFT_error("IO-Error in Hubbard 1 IO", calledby="write_hubbard1_input_new")

      CALL comment(input_iounit,"Orbital quantum number",1)
      CALL writeValue(input_iounit,"Lorb",l)
      CALL comment(input_iounit,"Energy level of the atomic level",1)
      CALL writeValue(input_iounit,"ea",-mu)
      CALL comment(input_iounit,"Spin-orbit-coupling parameter",1)
      CALL writeValue(input_iounit,"xiSOC",hub1%xi(i_hia))
      !calculate the additional exchange splitting
      exc = 0.0
      DO i_exc = 1, hub1%n_exc_given(i_hia)
         exc = exc + hub1%exc(i_hia,i_exc)*hub1%mag_mom(i_hia,i_exc)
      ENDDO
      !Only write the exchange splitting her if its not zero to not conflict with possible additional args
      IF(exc.NE.0.0) THEN
         CALL comment(input_iounit,"Exchange splitting",1)
         !Sign??
         CALL writeValue(input_iounit,"Exc",-exc)
      ENDIF
      CALL comment(input_iounit,"Additional arguments",1)
      DO i_arg = 1, hub1%n_addArgs(i_hia)
         CALL writeValue(input_iounit, TRIM(ADJUSTL(hub1%arg_keys(i_hia,i_arg))),hub1%arg_vals(i_hia,i_arg))
      ENDDO
      CALL writeValue(input_iounit, "cf")

      CALL cfmat%init(.true.,2*(2*l+1),2*(2*l+1))
      cfmat%data_r= 0.0
      DO i = 1, 2
         DO j = 1, (2*l+1)
            DO k = 1, (2*l+1)
               ind1 = (i-1)*(2*l+1) + j
               ind2 = (i-1)*(2*l+1) + k
               cfmat%data_r(ind1,ind2) = hub1%ccfmat(i_hia,j-l-1,k-l-1)*hartree_to_ev_const*hub1%ccf(i_hia)
            ENDDO
         ENDDO
      ENDDO
      CALL writeValue(input_iounit, cfmat)

      CLOSE(unit=input_iounit,iostat=io_error)
      IF(io_error.NE.0) CALL juDFT_error("IO-Error in Hubbard 1 IO", calledby="write_hubbard1_input_new")

   END SUBROUTINE write_hubbard1_input_new
   

   SUBROUTINE write_hubbard1_input_old(path,i_hia,l,f0,f2,f4,f6,hub1,mu,n)

      USE m_generic_txtio

      IMPLICIT NONE

      CHARACTER(len=*), INTENT(IN)  :: path
      INTEGER,          INTENT(IN)  :: i_hia
      INTEGER,          INTENT(IN)  :: l
      REAL,             INTENT(IN)  :: f0,f2,f4,f6
      TYPE(t_hub1ham),  INTENT(IN)  :: hub1
      REAL,             INTENT(IN)  :: mu
      INTEGER,          INTENT(IN)  :: n
      
      INTEGER :: info, io_error,i_exc
      REAL exc

      !Main input file
      OPEN(unit=input_iounit, file=TRIM(ADJUSTL(path)) // TRIM(ADJUSTL(cfg_file_main)),&
          status="replace", action="write", iostat=io_error)
      IF(io_error.NE.0) CALL juDFT_error("IO-Error in Hubbard 1 IO", calledby="write_hubbard1_input_old")

      CALL header(input_iounit,"Parameters for the atomic Hamiltonian in eV",1)

      CALL comment(input_iounit,"Orbital quantum number",1)
      CALL writeValue(input_iounit,"Lorb",l)

      CALL comment(input_iounit,"Slater Integrals",1)
      CALL writeValue(input_iounit,"Fk",(/f0,f2,f4,f6/))

      CALL comment(input_iounit,"Spin-orbit-coupling parameter",1)
      CALL writeValue(input_iounit,"gfact",hub1%xi(i_hia))

      !calculate the additional exchange splitting
      exc = 0.0
      DO i_exc = 1, hub1%n_exc_given(i_hia)
         exc = exc + hub1%exc(i_hia,i_exc)*hub1%mag_mom(i_hia,i_exc)
      ENDDO

      CALL comment(input_iounit,"External field",1)
      CALL writeValue(input_iounit,"Bz",exc)

      CALL comment(input_iounit,"Inverse temperature",1)
      CALL writeValue(input_iounit,"beta",hub1%beta)

      CALL comment(input_iounit,"Chemical potential",1)
      CALL writeValue(input_iounit,"mu",mu)

      IF(hub1%ccf(i_hia).NE.0.0) THEN
         CALL comment(input_iounit,"Crystal field factor",1)
         CALL writeValue(input_iounit,"ccf",hub1%ccf(i_hia))
      ENDIF

      CALL header(input_iounit,"Parameters for the Solver",1)

      CALL comment(input_iounit,"Minimum and maximum occupation of the orbital",1)
      CALL writeValue(input_iounit,"Nap_min",MAX(0,n-hub1%n_exc))
      CALL writeValue(input_iounit,"Nap_max",MIN(2*(2*l+1),n+hub1%n_exc))

      CALL comment(input_iounit,"Setting the solver to use the power lanczos method",1)
      CALL writeValue(input_iounit, "method_lancz")

      CALL comment(input_iounit,"Number of iterations",1)
      CALL writeValue(input_iounit,"N_lancz_iter",100)

      CALL comment(input_iounit,"Number of eigenstates calculated",1)
      CALL writeValue(input_iounit,"N_lancz_states",35)

      CALL header(input_iounit,"Parameters for the frequency/energy axis",1)

      CALL startSection(input_iounit,"real_freq_axis")
         CALL writeValue(input_iounit, "omegamin", emin)
         CALL writeValue(input_iounit, "omegamax", emax)
         CALL writeValue(input_iounit, "Nomega", ne)
         CALL writeValue(input_iounit, "eps", sigma)
      CALL endSection(input_iounit)

      CALL startSection(input_iounit,"matsub_freq_axis")
         CALL writeValue(input_iounit, "Nmatsub", nmats)
      CALL endSection(input_iounit)

      CLOSE(unit=input_iounit)

      WRITE(*,"(A)") "You are using an old input file format. This does not support the additional arguments"
   END SUBROUTINE write_hubbard1_input_old

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

      INTEGER :: info, io_error,io_unit

      io_unit = 17



      OPEN(unit=io_unit, file=TRIM(ADJUSTL(path)) // TRIM(ADJUSTL(cfg_file_ccf)), status="replace", action="write", iostat=io_error)

      IF(l.EQ.2) THEN
         WRITE(io_unit,"(5f10.6)") ccfmat*hartree_to_ev_const
      ELSE IF(l.EQ.3) THEN
         WRITE(io_unit,"(7f10.6)") ccfmat*hartree_to_ev_const
      ENDIF

      CLOSE(io_unit)

   END SUBROUTINE write_ccfmat

   SUBROUTINE read_ccfmat(path,ccfmat,l)

      CHARACTER(len=*), INTENT(IN)  :: path
      REAL,             INTENT(OUT) :: ccfmat(-l:l,-l:l)
      INTEGER,          INTENT(IN)  :: l

      INTEGER :: info, io_error,io_unit

      ccfmat = 0.0
      OPEN(unit=io_unit, file=TRIM(ADJUSTL(path)) // TRIM(ADJUSTL(cfg_file_ccf)), status="old", action="read", iostat=io_error)
      IF(io_error.NE.0) CALL juDFT_error("IO-error in Hubbard1-IO",calledby="read_ccfmat")

      !The crystal field is assumed to be in eV 
      IF(l.EQ.2) THEN
         READ(io_unit,"(5f10.6)") ccfmat
      ELSE IF(l.EQ.3) THEN
         READ(io_unit,"(7f10.6)") ccfmat
      ENDIF

      !convert to htr (in writing the input file its converted back)
      ccfmat = ccfmat/hartree_to_ev_const

      CLOSE(unit=io_unit)

   END SUBROUTINE read_ccfmat

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
               READ(io_unit,*) selfen(1:matsize,m,i)
            ENDDO
         ENDDO
      ELSE
         OPEN(unit=io_unit, file=TRIM(ADJUSTL(path)) // "se.atom", status="old", action="read", iostat=io_error)

         IF(io_error.NE.0) CALL juDFT_error("IO-Error in reading the self-energy", calledby="read_selfen")
         
         DO i = 1, ne
            READ(io_unit,9010) 
            READ(io_unit,9020) ((tmp(m,n), m= 1, matsize), n= 1, matsize)
            selfen(1:matsize,1:matsize,i) = tmp(1:matsize,1:matsize)/hartree_to_ev_const 
            READ(io_unit,9020) ((tmp(m,n), m= 1, matsize), n= 1, matsize)
            selfen(1:matsize,1:matsize,i) = selfen(1:matsize,1:matsize,i) + ImagUnit * tmp(1:matsize,1:matsize)/hartree_to_ev_const 
         ENDDO
      ENDIF


      CLOSE(io_unit)

9010  FORMAT(f10.5)
9020  FORMAT(7f11.5)
   END SUBROUTINE read_selfen


END MODULE m_hubbard1_io

