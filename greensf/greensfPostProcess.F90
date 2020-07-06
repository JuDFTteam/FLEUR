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

      REAL :: torgue(3),atomDiff(3)
      REAL, ALLOCATABLE :: f(:,:,:,:,:),g(:,:,:,:,:), flo(:,:,:,:,:)

      TYPE(t_usdus)            :: usdus
      TYPE(t_denCoeffsOffDiag) :: denCoeffsOffDiag
#ifdef CPP_HDF
      INTEGER(HID_T) :: greensf_fileID
#endif

      IF(mpi%irank==0) THEN
         IF(gfinp%checkRadial()) THEN
            CALL timestart("Green's Function: Radial Functions")
            ALLOCATE (f(atoms%jmtd,2,0:atoms%lmaxd,input%jspins,atoms%nType),source=0.0)
            ALLOCATE (g(atoms%jmtd,2,0:atoms%lmaxd,input%jspins,atoms%nType),source=0.0)
            ALLOCATE (flo(atoms%jmtd,2,atoms%nlod,input%jspins,atoms%nType),source=0.0)

            ! Initializations
            CALL usdus%init(atoms,input%jspins)
            CALL denCoeffsOffDiag%init(atoms,noco,sphhar,.FALSE.,.FALSE.)

            !Generate the scalar products we need
            DO i_gf = 1, gfinp%n
               l  = gfinp%elem(i_gf)%l
               lp = gfinp%elem(i_gf)%lp
               atomType  = gfinp%elem(i_gf)%atomType
               atomTypep = gfinp%elem(i_gf)%atomTypep

               l_sphavg  = gfinp%elem(i_gf)%l_sphavg
               IF(l_sphavg) CYCLE

               i_elem = gfinp%uniqueElements(atoms,ind=i_gf,l_sphavg=l_sphavg,indUnique=indUnique)

               IF(i_gf/=indUnique) CYCLE
               DO jspin = 1, input%jspins
                  CALL genMTBasis(atoms,enpara,vTot,mpi,atomType,jspin,usdus,&
                                  f(:,:,:,jspin,atomType),g(:,:,:,jspin,atomType),flo(:,:,:,jspin,atomType),&
                                  hub1inp%l_dftspinpol,l_writeArg=.FALSE.)
                  IF(atomType/=atomTypep) THEN
                     CALL genMTBasis(atoms,enpara,vTot,mpi,atomTypep,jspin,usdus,&
                                     f(:,:,:,jspin,atomTypep),g(:,:,:,jspin,atomTypep),flo(:,:,:,jspin,atomTypep),&
                                     hub1inp%l_dftspinpol,l_writeArg=.FALSE.)
                  ENDIF
               ENDDO
               IF(gfinp%l_mperp) THEN
                  CALL denCoeffsOffDiag%addRadFunScalarProducts(atoms,f(:,:,:,:,atomType),g(:,:,:,:,atomType),&
                                                                flo(:,:,:,:,atomType),atomType)
                  IF(atomType/=atomTypep) THEN
                     CALL denCoeffsOffDiag%addRadFunScalarProducts(atoms,f(:,:,:,:,atomTypep),g(:,:,:,:,atomTypep),&
                                                                   flo(:,:,:,:,atomTypep),atomTypep)
                  ENDIF
               ENDIF
            ENDDO
            CALL timestop("Green's Function: Radial Functions")
         ENDIF
      ENDIF
      !--------------------------------------------------------------------------------
      ! Obtain the real part of the Green's Function via the Kramers Kronig Integration
      !--------------------------------------------------------------------------------
      CALL timestart("Green's Function: Real Part")
      CALL greensfCalcRealPart(atoms,gfinp,input,sym,noco,usdus,denCoeffsOffDiag,mpi,results%ef,&
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

         CALL timestart("Green's Function: Occupation")
         mmpmat = cmplx_0
         DO i_gf = 1, gfinp%n
            !Occupation matrix
            l = greensFunction(i_gf)%elem%l
            lp = greensFunction(i_gf)%elem%lp
            atomType = greensFunction(i_gf)%elem%atomType
            atomTypep = greensFunction(i_gf)%elem%atomTypep
            atomDiff(:) = gfinp%elem(i_gf)%atomDiff(:)
            l_sphavg  = gfinp%elem(i_gf)%l_sphavg
            IF(l.NE.lp) CYCLE
            IF(atomType.NE.atomTypep) CYCLE
            IF(ANY(ABS(atomDiff(:)).GT.1e-12)) CYCLE
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
            CALL openXMLElementNoAttributes('noncollinearTorgue')
            WRITE(oUnit,'(/,A)') 'Torgue Calculation (noco):'
            WRITE(oUnit,'(/,A)') '---------------------------'
            DO atomType = 1, atoms%nType
               IF(gfinp%numTorgueElems(atomType)==0) CYCLE
               CALL greensfTorgue(greensFunction(gfinp%torgueElem(atomType,:gfinp%numTorgueElems(atomType))),&
                                  sphhar,atoms,sym,noco,nococonv,input,f,g,flo,atomType,torgue,vTot)
            ENDDO
            CALL closeXMLElement('noncollinearTorgue')
            IF(noco%l_soc.AND..FALSE.) THEN
               CALL openXMLElementNoAttributes('spinorbitTorgue')
               WRITE(oUnit,'(/,A)') 'Torgue Calculation (spin-orbit):'
               WRITE(oUnit,'(/,A)') '---------------------------------'
               DO atomType = 1, atoms%nType
                  IF(gfinp%numTorgueElems(atomType)==0) CYCLE
                  !CALL greensfSOTorgue(greensFunction(gfinp%torgueElem(atomType,:gfinp%numTorgueElems(atomType))),&
                  !                     sphhar,atoms,sym,noco,nococonv,input,f,g,flo,atomType,torgue,vso)
               ENDDO
               CALL closeXMLElement('spinorbitTorgue')
            ENDIF
            CALL timestop("Green's Function: Torgue")
         ENDIF
#ifdef CPP_HDF
         CALL timestart("Green's Function: IO/Write")
         CALL openGreensFFile(greensf_fileID, input, gfinp, atoms)
         CALL writeGreensFData(greensf_fileID, input, gfinp, atoms, &
                               GREENSF_GENERAL_CONST, greensFunction, mmpmat,&
                               u=f,udot=g,ulo=flo)
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
