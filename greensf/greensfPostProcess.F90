MODULE m_greensfPostProcess

   USE m_juDFT
   USE m_constants
   USE m_types
   USE m_greensfCalcRealPart
   USE m_greensf_io
   USE m_occmtx
   USE m_greensfTorgue
   USE m_excSplitting
   USE m_crystalfield
   USE m_genMTBasis
   USE m_radovlp

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE greensfPostProcess(greensFunction,greensfImagPart,atoms,gfinp,input,sym,noco,mpi,&
                                 nococonv,vTot,enpara,hub1inp,sphhar,hub1data,results)

      !contains all the modules for calculating properties from the greens function

      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_gfinp),             INTENT(IN)     :: gfinp
      TYPE(t_input),             INTENT(IN)     :: input
      TYPE(t_sym),               INTENT(IN)     :: sym
      TYPE(t_noco),              INTENT(IN)     :: noco
      TYPE(t_mpi),               INTENT(IN)     :: mpi
      TYPE(t_nococonv),          INTENT(IN)     :: nococonv
      TYPE(t_hub1inp),           INTENT(IN)     :: hub1inp
      TYPE(t_sphhar),            INTENT(IN)     :: sphhar
      TYPE(t_results),           INTENT(IN)     :: results
      TYPE(t_potden),            INTENT(IN)     :: vTot
      TYPE(t_enpara),            INTENT(IN)     :: enpara
      TYPE(t_hub1data),          INTENT(INOUT)  :: hub1data
      TYPE(t_greensfImagPart),   INTENT(INOUT)  :: greensfImagPart
      TYPE(t_greensf),           INTENT(INOUT)  :: greensFunction(:)

      INTEGER  i_gf,nType,l,lp,atomType,atomTypep,i_elem,indUnique,jspin,ierr
      COMPLEX  mmpmat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,gfinp%n,3)
      LOGICAL  l_sphavg,l_check

      REAL :: torgue(3)
      REAL, ALLOCATABLE :: u(:,:,:,:,:,:),udot(:,:,:,:,:,:)
      REAL, ALLOCATABLE :: f(:,:,:,:,:),g(:,:,:,:,:), flo(:,:,:,:,:)

      TYPE(t_usdus)            :: usdus
      TYPE(t_denCoeffsOffDiag) :: denCoeffsOffDiag
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

         CALL excSplitting(gfinp,atoms,input,greensfImagPart,results%ef)

         IF(gfinp%checkRadial()) THEN
            CALL timestart("Green's Function: Radial Functions")
            ALLOCATE (f(atoms%jmtd,2,0:atoms%lmaxd,input%jspins,atoms%nType),source=0.0)
            ALLOCATE (g(atoms%jmtd,2,0:atoms%lmaxd,input%jspins,atoms%nType),source=0.0)
            ALLOCATE (flo(atoms%jmtd,2,atoms%nlod,input%jspins,atoms%nType),source=0.0)

            ! Initializations
            CALL usdus%init(atoms,input%jspins)
            CALL denCoeffsOffDiag%init(atoms,noco,sphhar,.FALSE.,.FALSE.)


            ALLOCATE(u(atoms%jmtd,2,2,2,input%jspins,gfinp%n),source=0.0)
            ALLOCATE(udot(atoms%jmtd,2,2,2,input%jspins,gfinp%n),source=0.0)
            !Generate the scalar products we need
            DO i_gf = 1, gfinp%n
               l  = gfinp%elem(i_gf)%l
               lp = gfinp%elem(i_gf)%lp
               atomType  = gfinp%elem(i_gf)%atomType
               atomTypep = gfinp%elem(i_gf)%atomTypep
               l_sphavg  = gfinp%elem(i_gf)%l_sphavg
               IF(l_sphavg) CYCLE

               i_elem = gfinp%uniqueElements(atoms,ind=i_gf,l_sphavg=l_sphavg,indUnique=indUnique)

               IF(i_gf/=indUnique) THEN
                  u(:,:,:,:,:,i_gf) = u(:,:,:,:,:,indUnique)
                  udot(:,:,:,:,:,i_gf) = udot(:,:,:,:,:,indUnique)
               ELSE
                  DO jspin = 1, input%jspins
                     CALL genMTBasis(atoms,enpara,vTot,mpi,atomType,jspin,usdus,&
                                     f(:,:,:,jspin,atomType),g(:,:,:,jspin,atomType),flo(:,:,:,jspin,atomType),&
                                     hub1inp%l_dftspinpol,l_writeArg=.FALSE.)

                     u(:,:,1,1,jspin,i_gf) = f(:,:,l,jspin,atomType)
                     u(:,:,2,1,jspin,i_gf) = f(:,:,lp,jspin,atomType)

                     udot(:,:,1,1,jspin,i_gf) = g(:,:,l,jspin,atomType)
                     udot(:,:,2,1,jspin,i_gf) = g(:,:,lp,jspin,atomType)

                     CALL genMTBasis(atoms,enpara,vTot,mpi,atomTypep,jspin,usdus,&
                                     f(:,:,:,jspin,atomTypep),g(:,:,:,jspin,atomTypep),flo(:,:,:,jspin,atomTypep),&
                                     hub1inp%l_dftspinpol,l_writeArg=.FALSE.)

                     u(:,:,1,2,jspin,i_gf) = f(:,:,l,jspin,atomTypep)
                     u(:,:,2,2,jspin,i_gf) = f(:,:,lp,jspin,atomTypep)

                     udot(:,:,1,2,jspin,i_gf) = g(:,:,l,jspin,atomTypep)
                     udot(:,:,2,2,jspin,i_gf) = g(:,:,lp,jspin,atomTypep)
                  ENDDO
                  IF(gfinp%l_mperp) THEN
                     CALL denCoeffsOffDiag%addRadFunScalarProducts(atoms,f(:,:,:,:,atomType),g(:,:,:,:,atomType),&
                                                                   flo(:,:,:,:,atomType),atomType)
                     IF(atomType/=atomTypep) THEN
                        CALL denCoeffsOffDiag%addRadFunScalarProducts(atoms,f(:,:,:,:,atomTypep),g(:,:,:,:,atomTypep),&
                                                                      flo(:,:,:,:,atomTypep),atomTypep)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
            CALL timestop("Green's Function: Radial Functions")
         ENDIF
         CALL timestart("Green's Function: Occupation")
         mmpmat = cmplx_0
         DO i_gf = 1, gfinp%n
            !Occupation matrix
            l = greensFunction(i_gf)%elem%l
            lp = greensFunction(i_gf)%elem%lp
            atomType = greensFunction(i_gf)%elem%atomType
            atomTypep = greensFunction(i_gf)%elem%atomTypep
            l_sphavg  = gfinp%elem(i_gf)%l_sphavg
            IF(l.NE.lp) CYCLE
            IF(atomType.NE.atomTypep) CYCLE
            l_check = gfinp%elem(i_gf)%countLOs(atoms)==0 !If there are SCLOs present the occupations can get bigger than 1
            IF(l_sphavg) THEN
               CALL occmtx(greensFunction(i_gf),gfinp,input,atoms,mmpmat(:,:,i_gf,:),l_write=.TRUE.,check=l_check)
            ELSE IF(.NOT.gfinp%l_mperp) THEN
               CALL occmtx(greensFunction(i_gf),gfinp,input,atoms,mmpmat(:,:,i_gf,:),&
                           usdus=usdus,l_write=.TRUE.,check=l_check)
            ELSE
               CALL occmtx(greensFunction(i_gf),gfinp,input,atoms,mmpmat(:,:,i_gf,:),&
                           usdus=usdus,denCoeffsOffDiag=denCoeffsOffDiag,&
                           l_write=.TRUE.,check=l_check)
            ENDIF
         ENDDO
         CALL timestop("Green's Function: Occupation")

         IF(ANY(gfinp%numTorgueElems>0)) THEN
            CALL timestart("Green's Function: Torgue")
            CALL openXMLElementNoAttributes('torgueCalculation')
            WRITE(oUnit,'(/,A)') 'Torgue Calculation:'
            WRITE(oUnit,'(/,A)') '------------------------'
            DO atomType = 1, atoms%nType
               IF(gfinp%numTorgueElems(atomType)==0) CYCLE
               CALL greensfTorgue(greensFunction(gfinp%torgueElem(atomType,:gfinp%numTorgueElems(atomType))),vTot,&
                                  sphhar,atoms,sym,noco,nococonv,input,enpara,hub1inp,mpi,atomType,torgue)
            ENDDO
            CALL closeXMLElement('torgueCalculation')
            CALL timestop("Green's Function: Torgue")
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
