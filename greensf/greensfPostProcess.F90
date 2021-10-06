MODULE m_greensfPostProcess

   USE m_juDFT
   USE m_constants
   USE m_types
   USE m_greensfCalcRealPart
   USE m_greensfCalcScalarProducts
   USE m_greensf_io
   USE m_greensfTorque
   USE m_excSplitting
   USE m_crystalfield
   USE m_genMTBasis
   USE m_sointg

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE greensfPostProcess(greensFunction,greensfImagPart,atoms,kpts,cell,gfinp,input,sym,noco,mpi,&
                                 nococonv,vTot,enpara,hub1inp,sphhar,hub1data,results)

      !contains all the modules for calculating properties from the greens function

      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_kpts),              INTENT(IN)     :: kpts
      TYPE(t_cell),              INTENT(IN)     :: cell
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

      INTEGER  i_gf,l,lp,atomType,atomTypep,i_elem,jspin,ierr,i,indUnique
      COMPLEX  mmpmat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,gfinp%n,3)
      LOGICAL  l_sphavg,l_check
      REAL, ALLOCATABLE :: f(:,:,:,:,:),g(:,:,:,:,:), flo(:,:,:,:,:)
      TYPE(t_scalarGF), ALLOCATABLE :: scalarGF(:)

#ifdef CPP_HDF
      INTEGER(HID_T) :: greensf_fileID
#endif

      CALL greensfCalcScalarProducts(gfinp,atoms,input,enpara,noco,sphhar,vTot,mpi,hub1data,&
                                     greensFunctions=greensFunction,scalarProducts=scalarGF,&
                                     fout=f,gout=g,floout=flo)
      !--------------------------------------------------------------------------------
      ! Obtain the real part of the Green's Function via the Kramers Kronig Integration
      !--------------------------------------------------------------------------------
      CALL timestart("Green's Function: Real Part")
      CALL greensfCalcRealPart(atoms,gfinp,sym,input,noco,kpts,mpi,results%ef,&
                               greensfImagPart,greensFunction)
      CALL timestop("Green's Function: Real Part")

      !----------------------------------------------
      ! Torque Calculations
      !----------------------------------------------
      IF(ANY(gfinp%numTorqueElems>0)) THEN
         CALL timestart("Green's Function: Torque")
         CALL greensfTorque(greensFunction,gfinp,mpi,sphhar,atoms,sym,noco,nococonv,input,&
                            f,g,flo,vTot)
         CALL timestop("Green's Function: Torque")
      ENDIF

      IF(mpi%irank==0) THEN
         CALL timestart("Green's Function: Postprocess")
         !-------------------------------------------------------------
         ! Calculate various properties from the greens function
         !-------------------------------------------------------------
         !calculate the crystal field contribution to the local hamiltonian in LDA+Hubbard 1
         IF(atoms%n_hia.GT.0 .AND. ANY(ABS(hub1inp%ccf(:)).GT.1e-12)) THEN
           !CALL crystal_field(atoms,gfinp,input,noco,nococonv,greensfImagPart,vTot,results%ef,hub1data)
         ENDIF

         !Onsite exchange Splitting from difference of center of mass of the bands
         CALL excSplitting(gfinp,atoms,input,scalarGF,greensfImagPart,results%ef)

         !-------------------------------------------------------
         ! Occupation matrix (only for diagonal onsite elements)
         !-------------------------------------------------------
         CALL timestart("Green's Function: Occupation")
         mmpmat = cmplx_0
         DO i_gf = 1, gfinp%n
            IF(gfinp%elem(i_gf)%l_kresolved) CYCLE
            !If there are SCLOs present the occupations can get bigger than 1
            l_check = gfinp%elem(i_gf)%countLOs(atoms)==0 .AND..NOT.gfinp%elem(i_gf)%isOffDiag()
            mmpmat(:,:,i_gf,:) = greensFunction(i_gf)%occupationMatrix(gfinp,input,atoms,noco,nococonv,&
                                                                       l_write=.TRUE.,check=l_check)
         ENDDO
         CALL timestop("Green's Function: Occupation")

#ifdef CPP_HDF
         CALL timestart("Green's Function: IO/Write")
         CALL openGreensFFile(greensf_fileID, input, gfinp, atoms, cell, kpts)
         CALL writeGreensFData(greensf_fileID, input, gfinp, atoms, cell,&
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
