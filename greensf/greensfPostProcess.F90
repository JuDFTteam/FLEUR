MODULE m_greensfPostProcess

   !Contains the main routines called from cdnval and cdngen
   USE m_juDFT
   USE m_constants
   USE m_types
   USE m_greensfCalcRealPart
   USE m_greensf_io
   USE m_greensfUtils

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE greensfPostProcess(greensFunction,greensfImagPart,atoms,gfinp,input,sym,noco,mpi,&
                                 nococonv,vTot,hub1inp,hub1data,results)

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
      TYPE(t_hub1data),          INTENT(INOUT)  :: hub1data
      TYPE(t_greensfImagPart),   INTENT(INOUT)  :: greensfImagPart
      TYPE(t_greensf),           INTENT(INOUT)  :: greensFunction(:)

      INTEGER  i_gf,l,nType
      COMPLEX  mmpmat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,gfinp%n,3)
      LOGICAL  err

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
         IF(atoms%n_hia.GT.0.AND.ANY(ABS(hub1inp%ccf(:)).GT.1e-12)) THEN
           CALL crystal_field(atoms,gfinp,hub1inp,input,nococonv,greensfImagPart,vTot,results%ef,hub1data)
         ENDIF
         IF(input%jspins.EQ.2) THEN
            !CALL eff_excinteraction(greensFunction,gfinp,input,results%ef,greensfImagPart)
         ENDIF
         CALL timestart("Green's Function: Occupation/DOS")
         DO i_gf = 1, gfinp%n
            !IF(l.NE.gfinp%elem(i_gf)%lp) CYCLE
            !IF(nType.NE.gfinp%elem(i_gf)%atomTypep) CYCLE
            !Density of states from Greens function
            !CALL gfDOS(greensFunction,l,nType,i_gf,gfinp,input,results%ef)
            !Occupation matrix
            CALL occmtx(greensFunction(i_gf),i_gf,gfinp,input,mmpmat(:,:,i_gf,:),err,l_write=.TRUE.,check=.TRUE.)
            !Hybridization function
            !CALL hybridization(greensFunction(i_gf),i_gf,gfinp,input,results%ef)
         ENDDO
         CALL timestop("Green's Function: Occupation/DOS")

#ifdef CPP_HDF
         CALL timestart("Green's Function: IO/Write")
         CALL openGreensFFile(greensf_fileID, input, gfinp, atoms)
         CALL writeGreensFData(greensf_fileID, input, gfinp, atoms, &
                              GREENSF_GENERAL_CONST, greensFunction, mmpmat)
         CALL closeGreensFFile(greensf_fileID)
         CALL timestop("Green's Function: IO/Write")
#endif
         CALL timestop("Green's Function: Postprocess")
      ENDIF

   END SUBROUTINE greensfPostProcess

END MODULE m_greensfPostProcess
