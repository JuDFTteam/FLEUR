MODULE m_gfcalc

   !Contains the main routines called from cdnval and cdngen
   USE m_juDFT
   USE m_constants
   USE m_types
   !These dependencies are here so that all these subroutines can be called by simply adding USE m_gfcalc to the file
   USE m_j0
   USE m_gfDOS
   USE m_occmtx
   USE m_hybridization
   USE m_crystalfield

   IMPLICIT NONE


   CONTAINS

   SUBROUTINE bzIntegrationGF(atoms,gfinp,sym,input,ispin,nbands,dosWeights,resWeights,indBound,&
                              wtkpt,ef,eig,denCoeffsOffdiag,usdus,eigVecCoeffs,greensfCoeffs,l21)

      USE m_greensfImag
      USE m_greensfImag21

      TYPE(t_atoms),             INTENT(IN)    :: atoms
      TYPE(t_gfinp),             INTENT(IN)    :: gfinp
      TYPE(t_sym),               INTENT(IN)    :: sym
      TYPE(t_input),             INTENT(IN)    :: input
      TYPE(t_eigVecCoeffs),      INTENT(IN)    :: eigVecCoeffs
      TYPE(t_usdus),             INTENT(IN)    :: usdus
      TYPE(t_denCoeffsOffDiag),  INTENT(IN)    :: denCoeffsOffdiag
      TYPE(t_greensfCoeffs),     INTENT(INOUT) :: greensfCoeffs
      INTEGER,                   INTENT(IN)    :: ispin  !Current spin index
      INTEGER,                   INTENT(IN)    :: nbands !Number of bands to be considered
      REAL,                      INTENT(IN)    :: wtkpt  !Weight of the current k-point
      REAL,                      INTENT(IN)    :: ef
      LOGICAL,                   INTENT(IN)    :: l21    !Calculate spin off-diagonal part ?
      REAL,    ALLOCATABLE,      INTENT(IN)    :: resWeights(:,:)
      REAL,    ALLOCATABLE,      INTENT(IN)    :: dosWeights(:,:) !Precalculated tetrahedron weights for the current k-point
      INTEGER, ALLOCATABLE,      INTENT(IN)    :: indBound(:,:)   !Gives the range where the tetrahedron weights are non-zero
      REAL,                      INTENT(IN)    :: eig(:)          !Eigenvalues for the current k-point

      CALL timestart("Greens Function: Imaginary Part")
      CALL greensfImag(atoms,gfinp,sym,input,ispin,nbands,dosWeights,resWeights,indBound,&
                       wtkpt,ef,eig,usdus,eigVecCoeffs,greensfCoeffs)
      IF(gfinp%l_mperp.AND.l21) THEN
         CALL greensfImag21(atoms,gfinp,sym,input,nbands,dosWeights,resWeights,indBound,&
                            wtkpt,ef,eig,denCoeffsOffdiag,eigVecCoeffs,greensfCoeffs)
      ENDIF
      CALL timestop("Greens Function: Imaginary Part")

   END SUBROUTINE bzIntegrationGF



   SUBROUTINE postProcessGF(greensf,greensfCoeffs,atoms,gfinp,input,sym,noco,nococonv,vTot,hub1inp,hub1data,results)

      !contains all the modules for calculating properties from the greens function
      USE m_onsite
      USE m_rot_gf
      USE m_greensf_io

      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_gfinp),             INTENT(IN)     :: gfinp
      TYPE(t_input),             INTENT(IN)     :: input
      TYPE(t_sym),               INTENT(IN)     :: sym
      TYPE(t_noco),              INTENT(IN)     :: noco
      TYPE(t_nococonv),          INTENT(IN)     :: nococonv
      TYPE(t_hub1inp),           INTENT(IN)     :: hub1inp
      TYPE(t_results),           INTENT(IN)     :: results
      TYPE(t_potden),            INTENT(IN)     :: vTot
      TYPE(t_hub1data),          INTENT(INOUT)  :: hub1data
      TYPE(t_greensfCoeffs),     INTENT(INOUT)  :: greensfCoeffs
      TYPE(t_greensf),           INTENT(INOUT)  :: greensf

      INTEGER  i_gf,l,nType
      COMPLEX  mmpmat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,gfinp%n,3)
      LOGICAL  err

#ifdef CPP_HDF
      INTEGER(HID_T) :: greensf_fileID
#endif


      CALL timestart("Green's Function: Postprocess")
      CALL rot_projDOS(sym,atoms,gfinp,input,greensfCoeffs)
      !Perform the Kramer-Kronigs-Integration if we only have the imaginary part at this point
      CALL calc_onsite(atoms,gfinp,input,sym,noco,results%ef,greensfCoeffs,greensf)
      !-------------------------------------------------------------
      ! Calculate various properties from the greens function
      !-------------------------------------------------------------
      !calculate the crystal field contribution to the local hamiltonian in LDA+Hubbard 1
      IF(atoms%n_hia.GT.0.AND.ANY(ABS(hub1inp%ccf(:)).GT.1e-12)) THEN
        CALL crystal_field(atoms,gfinp,hub1inp,input,nococonv,greensfCoeffs,vTot,results%ef,hub1data)
      ENDIF
      IF(input%jspins.EQ.2) THEN
         CALL eff_excinteraction(greensf,gfinp,input,results%ef,greensfCoeffs)
      ENDIF
      CALL timestart("Green's Function: Occupation/DOS")
      DO i_gf = 1, gfinp%n
         l = gfinp%elem(i_gf)%l
         nType = gfinp%elem(i_gf)%atomType
         IF(l.NE.gfinp%elem(i_gf)%lp) CYCLE
         IF(nType.NE.gfinp%elem(i_gf)%atomTypep) CYCLE
         !Density of states from Greens function
         CALL gfDOS(greensf,l,nType,i_gf,gfinp,input,results%ef)
         !Occupation matrix
         CALL occmtx(greensf,l,nType,gfinp,input,mmpmat(:,:,i_gf,:),err,l_write=.TRUE.,check=.TRUE.)
         !Hybridization function
         !CALL hybridization(greensf,l,nType,gfinp,input,results%ef)
      ENDDO
      CALL timestop("Green's Function: Occupation/DOS")
      CALL timestop("Green's Function: Postprocess")

#ifdef CPP_HDF
      CALL timestart("Green's Function: IO/Write")
      CALL openGreensFFile(greensf_fileID, input, gfinp, atoms, greensf)
      CALL writeGreensFData(greensf_fileID, input, gfinp, greensf)
      CALL closeGreensFFile(greensf_fileID)
      CALL timestop("Green's Function: IO/Write")
#endif

   END SUBROUTINE postProcessGF

END MODULE m_gfcalc
