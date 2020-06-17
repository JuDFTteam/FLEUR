MODULE m_greensfPostProcess

   USE m_juDFT
   USE m_constants
   USE m_types
   USE m_greensfCalcRealPart
   USE m_greensf_io
   USE m_occmtx
   USE m_excSplitting
   USE m_crystalfield
   USE m_genMTBasis

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE greensfPostProcess(greensFunction,greensfImagPart,atoms,gfinp,input,sym,noco,mpi,&
                                 nococonv,vTot,enpara,hub1inp,hub1data,results)

      !contains all the modules for calculating properties from the greens function

      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_gfinp),             INTENT(IN)     :: gfinp
      TYPE(t_input),             INTENT(IN)     :: input
      TYPE(t_sym),               INTENT(IN)     :: sym
      TYPE(t_noco),              INTENT(IN)     :: noco
      TYPE(t_mpi),               INTENT(IN)     :: mpi
      TYPE(t_nococonv),          INTENT(IN)     :: nococonv
      TYPE(t_hub1inp),           INTENT(IN)     :: hub1inp
      TYPE(t_results),           INTENT(IN)     :: results
      TYPE(t_potden),            INTENT(IN)     :: vTot
      TYPE(t_enpara),            INTENT(IN)     :: enpara
      TYPE(t_hub1data),          INTENT(INOUT)  :: hub1data
      TYPE(t_greensfImagPart),   INTENT(INOUT)  :: greensfImagPart
      TYPE(t_greensf),           INTENT(INOUT)  :: greensFunction(:)

      INTEGER  i_gf,nType,l,lp,atomType,atomTypep,i_elem,indUnique,jspin,ierr
      COMPLEX  mmpmat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,gfinp%n,3)

      REAL, ALLOCATABLE :: u(:,:,:,:,:,:),udot(:,:,:,:,:,:)
      REAL, ALLOCATABLE :: f(:,:,:),g(:,:,:), flo(:,:,:)

      TYPE(t_usdus) :: usdus

#ifdef CPP_HDF
      INTEGER(HID_T) :: greensf_fileID
#endif

      !--------------------------------------------------------------------------------
      ! Obtain the real part of the Green's Function via the Kramers Kronig Integration
      !--------------------------------------------------------------------------------
      CALL timestart("Green's Function: Real Part")
      CALL greensfCalcRealPart(atoms,gfinp,input,sym,noco,mpi,results%ef,greensfImagPart,greensFunction)
      CALL timestop("Green's Function: Real Part")

      IF(mpi%irank==0) THEN
         CALL timestart("Green's Function: Postprocess")
         !-------------------------------------------------------------
         ! Calculate various properties from the greens function
         !-------------------------------------------------------------
         !calculate the crystal field contribution to the local hamiltonian in LDA+Hubbard 1
         IF(atoms%n_hia.GT.0 .AND. ANY(ABS(hub1inp%ccf(:)).GT.1e-12)) THEN
           CALL crystal_field(atoms,gfinp,hub1inp,input,nococonv,greensfImagPart,vTot,results%ef,hub1data)
         ENDIF

         CALL excSplitting(gfinp,input,greensfImagPart,results%ef)
         CALL timestart("Green's Function: Occupation")
         DO i_gf = 1, gfinp%n
            !Occupation matrix
            CALL occmtx(greensFunction(i_gf),gfinp,input,mmpmat(:,:,i_gf,:),l_write=.TRUE.,check=.TRUE.)
         ENDDO
         CALL timestop("Green's Function: Occupation")


         IF(.NOT.gfinp%l_sphavg) THEN
            CALL timestart("Green's Function: Radial Functions")

            !Intializations
            ALLOCATE(u(atoms%jmtd,2,2,2,input%jspins,gfinp%n),source=0.0)
            ALLOCATE(udot(atoms%jmtd,2,2,2,input%jspins,gfinp%n),source=0.0)

            ALLOCATE(f(atoms%jmtd,2,atoms%lmaxd),source=0.0)
            ALLOCATE(g(atoms%jmtd,2,atoms%lmaxd),source=0.0)
            ALLOCATE(flo(atoms%jmtd,2,atoms%nlod),source=0.0)

            CALL usdus%init(atoms,input%jspins)

            DO i_gf = 1, gfinp%n
               l  = gfinp%elem(i_gf)%l
               lp = gfinp%elem(i_gf)%lp
               atomType  = gfinp%elem(i_gf)%atomType
               atomTypep = gfinp%elem(i_gf)%atomTypep

               i_elem = uniqueElements_gfinp(gfinp,ind=i_gf,indUnique=indUnique)

               IF(i_gf/=indUnique) THEN
                  u(:,:,:,:,:,i_gf) = u(:,:,:,:,:,indUnique)
                  udot(:,:,:,:,:,i_gf) = udot(:,:,:,:,:,indUnique)
               ELSE
                  DO jspin = 1, input%jspins
                     CALL genMTBasis(atoms,enpara,vTot,mpi,atomType,jspin,usdus,f,g,flo,hub1inp%l_dftspinpol)

                     u(:,:,1,1,jspin,i_gf) = f(:,:,l)
                     u(:,:,2,1,jspin,i_gf) = f(:,:,lp)

                     udot(:,:,1,1,jspin,i_gf) = g(:,:,l)
                     udot(:,:,2,1,jspin,i_gf) = g(:,:,lp)

                     CALL genMTBasis(atoms,enpara,vTot,mpi,atomTypep,jspin,usdus,f,g,flo,hub1inp%l_dftspinpol)

                     u(:,:,1,2,jspin,i_gf) = f(:,:,l)
                     u(:,:,2,2,jspin,i_gf) = f(:,:,lp)

                     udot(:,:,1,2,jspin,i_gf) = g(:,:,l)
                     udot(:,:,2,2,jspin,i_gf) = g(:,:,lp)

                  ENDDO
               ENDIF
            ENDDO
            CALL timestop("Green's Function: Radial Functions")
         ENDIF


#ifdef CPP_HDF
         CALL timestart("Green's Function: IO/Write")
         CALL openGreensFFile(greensf_fileID, input, gfinp, atoms)
         CALL writeGreensFData(greensf_fileID, input, gfinp, atoms, &
                               GREENSF_GENERAL_CONST, greensFunction, mmpmat,&
                               u=u,udot=udot)
         CALL closeGreensFFile(greensf_fileID)
         CALL timestop("Green's Function: IO/Write")
#endif
         CALL timestop("Green's Function: Postprocess")
      ENDIF

#ifdef CPP_MPI
      CALL MPI_BARRIER(mpi%mpi_comm,ierr)
#endif

   END SUBROUTINE greensfPostProcess

END MODULE m_greensfPostProcess
