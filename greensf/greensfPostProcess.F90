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
   USE m_radovlp

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
      REAL, ALLOCATABLE :: uun21(:,:),udn21(:,:),dun21(:,:),ddn21(:,:)
      REAL, ALLOCATABLE :: f(:,:,:,:),g(:,:,:,:), flo(:,:,:,:)

      TYPE(t_usdus) :: usdus

#ifdef CPP_HDF
      INTEGER(HID_T) :: greensf_fileID
#endif

      !--------------------------------------------------------------------------------
      ! Obtain the real part of the Green's Function via the Kramers Kronig Integration
      !--------------------------------------------------------------------------------
      CALL timestart("Green's Function: Real Part")
      CALL greensfCalcRealPart(atoms,gfinp,input,sym,noco,vTot,enpara,mpi,hub1inp,results%ef,&
                               greensfImagPart,greensFunction)
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
         IF(.NOT.gfinp%l_sphavg) THEN

            ALLOCATE (f(atoms%jmtd,2,0:atoms%lmaxd,input%jspins))
            ALLOCATE (g(atoms%jmtd,2,0:atoms%lmaxd,input%jspins))
            ALLOCATE (flo(atoms%jmtd,2,atoms%nlod,input%jspins))

            ! Initializations
            CALL usdus%init(atoms,input%jspins)
            !Generate the scalar products we need
            DO i_gf = 1, gfinp%n
               atomType = greensFunction(i_gf)%elem%atomType
               DO jspin = 1, input%jspins
                  CALL genMTBasis(atoms,enpara,vTot,mpi,atomType,jspin,usdus,f,g,flo,hub1inp%l_dftspinpol)
               ENDDO
            ENDDO
            DEALLOCATE(f,g,flo)
            !Offdiagonal scalar products
            IF(gfinp%l_mperp) THEN
               !Calculate overlap integrals
               ALLOCATE(uun21(0:atoms%lmaxd,atoms%ntype),source=0.0)
               ALLOCATE(dun21(0:atoms%lmaxd,atoms%ntype),source=0.0)
               ALLOCATE(udn21(0:atoms%lmaxd,atoms%ntype),source=0.0)
               ALLOCATE(ddn21(0:atoms%lmaxd,atoms%ntype),source=0.0)
               CALL rad_ovlp(atoms,usdus,input,hub1inp,vTot%mt,enpara%el0, uun21,udn21,dun21,ddn21)
            ENDIF
         ENDIF
         DO i_gf = 1, gfinp%n
            !Occupation matrix
            l = greensFunction(i_gf)%elem%l
            lp = greensFunction(i_gf)%elem%lp
            atomType = greensFunction(i_gf)%elem%atomType
            atomTypep = greensFunction(i_gf)%elem%atomTypep
            IF(l.NE.lp) CYCLE
            IF(atomType.NE.atomTypep) CYCLE
            IF(gfinp%l_sphavg) THEN
               CALL occmtx(greensFunction(i_gf),gfinp,input,mmpmat(:,:,i_gf,:),l_write=.TRUE.,check=.TRUE.)
            ELSE
               CALL occmtx(greensFunction(i_gf),gfinp,input,mmpmat(:,:,i_gf,:),&
                           ddn=usdus%ddn(l,atomType,:),uun21=uun21(l,atomType),&
                           udn21=udn21(l,atomType),dun21=dun21(l,atomType),&
                           ddn21=ddn21(l,atomType),l_write=.TRUE.,check=.TRUE.)
            ENDIF
         ENDDO
         CALL timestop("Green's Function: Occupation")


         IF(.NOT.gfinp%l_sphavg) THEN
            CALL timestart("Green's Function: Radial Functions")

            !Intializations
            ALLOCATE(u(atoms%jmtd,2,2,2,input%jspins,gfinp%n),source=0.0)
            ALLOCATE(udot(atoms%jmtd,2,2,2,input%jspins,gfinp%n),source=0.0)

            ALLOCATE(f(atoms%jmtd,2,atoms%lmaxd,input%jspins),source=0.0)
            ALLOCATE(g(atoms%jmtd,2,atoms%lmaxd,input%jspins),source=0.0)
            ALLOCATE(flo(atoms%jmtd,2,atoms%nlod,input%jspins),source=0.0)

            DO i_gf = 1, gfinp%n
               l  = gfinp%elem(i_gf)%l
               lp = gfinp%elem(i_gf)%lp
               atomType  = gfinp%elem(i_gf)%atomType
               atomTypep = gfinp%elem(i_gf)%atomTypep

               i_elem = gfinp%uniqueElements(ind=i_gf,indUnique=indUnique)

               IF(i_gf/=indUnique) THEN
                  u(:,:,:,:,:,i_gf) = u(:,:,:,:,:,indUnique)
                  udot(:,:,:,:,:,i_gf) = udot(:,:,:,:,:,indUnique)
               ELSE
                  DO jspin = 1, input%jspins
                     CALL genMTBasis(atoms,enpara,vTot,mpi,atomType,jspin,usdus,&
                                     f(:,:,:,jspin),g(:,:,:,jspin),flo(:,:,:,jspin),hub1inp%l_dftspinpol)

                     u(:,:,1,1,jspin,i_gf) = f(:,:,l,jspin)
                     u(:,:,2,1,jspin,i_gf) = f(:,:,lp,jspin)

                     udot(:,:,1,1,jspin,i_gf) = g(:,:,l,jspin)
                     udot(:,:,2,1,jspin,i_gf) = g(:,:,lp,jspin)

                     CALL genMTBasis(atoms,enpara,vTot,mpi,atomTypep,jspin,usdus,&
                                     f(:,:,:,jspin),g(:,:,:,jspin),flo(:,:,:,jspin),hub1inp%l_dftspinpol)

                     u(:,:,1,2,jspin,i_gf) = f(:,:,l,jspin)
                     u(:,:,2,2,jspin,i_gf) = f(:,:,lp,jspin)

                     udot(:,:,1,2,jspin,i_gf) = g(:,:,l,jspin)
                     udot(:,:,2,2,jspin,i_gf) = g(:,:,lp,jspin)

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
