MODULE m_gfcalc

   !Contains the main routines called from cdnval and cdngen
   USE m_juDFT
   USE m_constants
   !These dependencies are here so that all these subroutines can be called by simply adding USE m_gfcalc to the file
   USE m_j0
   USE m_gfDOS
   USE m_occmtx
   USE m_hybridization
   USE m_crystalfield

   CONTAINS

   SUBROUTINE bzIntegrationGF(atoms,sym,input,angle,ispin,nbands,dosWeights,resWeights,indBound,wtkpt,eig,denCoeffsOffdiag,&
                              usdus,eigVecCoeffs,greensf,greensfCoeffs,l21)

      USE m_greensfImag
      USE m_greensfImag21
      !USE m_greensfRes
      !USE m_greensfRes21

      IMPLICIT NONE

      !-Type Arguments
      TYPE(t_atoms),             INTENT(IN)    :: atoms
      TYPE(t_sym),               INTENT(IN)    :: sym
      TYPE(t_input),             INTENT(IN)    :: input 
      TYPE(t_eigVecCoeffs),      INTENT(IN)    :: eigVecCoeffs
      TYPE(t_usdus),             INTENT(IN)    :: usdus
      TYPE(t_denCoeffsOffDiag),  INTENT(IN)    :: denCoeffsOffdiag
      TYPE(t_greensfCoeffs),     INTENT(INOUT) :: greensfCoeffs
      TYPE(t_greensf),           INTENT(INOUT) :: greensf


      !-Scalar Arguments 
      INTEGER,                   INTENT(IN)    :: ispin  !Current spin index
      INTEGER,                   INTENT(IN)    :: nbands !Number of bands to be considered
      REAL,                      INTENT(IN)    :: wtkpt  !Weight of the current k-point
      LOGICAL,                   INTENT(IN)    :: l21    !Calculate spin off-diagonal part ?

      !-Array Arguments
      REAL,                      INTENT(IN)    :: resWeights(greensfCoeffs%ne,nbands)
      REAL,                      INTENT(IN)    :: dosWeights(greensfCoeffs%ne,nbands) !Precalculated tetrahedron weights for the current k-point
      INTEGER,                   INTENT(IN)    :: indBound(nbands,2)                  !Gives the range where the tetrahedron weights are non-zero
      REAL,                      INTENT(IN)    :: eig(nbands)                         !Eigenvalues for the current k-point
      REAL,                      INTENT(IN)    :: angle(sym%nop)                      !Phases for spin-offdiagonal part
      
      !IF(input%l_resolvent) THEN
         !Calculate greens function directly
      !   CALL timestart("Greens Function: Resolvent")
      !   CALL greensfRes(atoms,sym,input,ispin,nbands,resWeights,indBound,wtkpt,eig,usdus,eigVecCoeffs,greensf)
      !   IF(input%l_gfmperp.AND.l21) THEN
      !      CALL greensfRes21(atoms,sym,angle,input,nbands,resWeights,indBound,wtkpt,eig,denCoeffsOffdiag,eigVecCoeffs,greensf)
      !   ENDIF
      !   CALL timestop("Greens Function: Resolvent")
      !ENDIF
      CALL timestart("Greens Function: Imaginary Part")
      CALL greensfImag(atoms,sym,input,ispin,nbands,dosWeights,resWeights,indBound,wtkpt,eig,usdus,eigVecCoeffs,greensfCoeffs)
      IF(input%l_gfmperp.AND.l21) THEN
         CALL greensfImag21(atoms,sym,angle,input,nbands,dosWeights,resWeights,indBound,wtkpt,eig,denCoeffsOffdiag,eigVecCoeffs,greensfCoeffs)
      ENDIF
      CALL timestop("Greens Function: Imaginary Part")


   END SUBROUTINE bzIntegrationGF



   SUBROUTINE postProcessGF(greensf,greensfCoeffs,atoms,input,sym,noco,vTot,hub1,results,angle)
      
      !contains all the modules for calculating properties from the greens function
      USE m_onsite

      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_input),             INTENT(IN)     :: input
      TYPE(t_sym),               INTENT(IN)     :: sym
      TYPE(t_noco),              INTENT(IN)     :: noco
      TYPE(t_greensfCoeffs),     INTENT(INOUT)  :: greensfCoeffs
      TYPE(t_greensf),           INTENT(INOUT)  :: greensf
      TYPE(t_hub1ham),           INTENT(INOUT)  :: hub1
      TYPE(t_results),           INTENT(IN)     :: results
      TYPE(t_potden),            INTENT(IN)     :: vTot
      REAL,                   INTENT(IN)     :: angle(sym%nop)

      INTEGER  i_gf,l,nType
      COMPLEX  mmpmat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_gf,3)


      CALL timestart("Green's Function: Postprocess")
      !Perform the Kramer-Kronigs-Integration if we only have he imaginary part at this point
      CALL calc_onsite(atoms,input,sym,noco,angle,greensfCoeffs,greensf)
      !-------------------------------------------------------------
      ! Calculate various properties from the greens function 
      !-------------------------------------------------------------
      !calculate the crystal field contribution to the local hamiltonian in LDA+Hubbard 1
      IF(atoms%n_hia.GT.0.AND.ANY(hub1%ccf(:).NE.0.0)) THEN
        CALL crystal_field(atoms,input,greensfCoeffs,hub1,vTot)
      ENDIF
      IF(input%jspins.EQ.2) THEN
         CALL eff_excinteraction(greensf,atoms,input,results%ef,greensfCoeffs)
      ENDIF
      DO i_gf = 1, atoms%n_gf
         l = atoms%gfelem(i_gf)%l
         nType = atoms%gfelem(i_gf)%atomType
         IF(l.NE.atoms%gfelem(i_gf)%lp) CYCLE
         IF(nType.NE.atoms%gfelem(i_gf)%atomTypep) CYCLE
         !Occupation matrix
         CALL occmtx(greensf,l,nType,atoms,sym,input,mmpmat(:,:,i_gf,:),l_write=.TRUE.,check=.TRUE.)
         !Hybridization function
         CALL hybridization(greensf,l,nType,atoms,input,results%ef)
         !Density of states from Greens function
         CALL gfDOS(greensf,l,nType,i_gf,atoms,input,results%ef)
      ENDDO
      CALL timestop("Green's Function: Postprocess")

   END SUBROUTINE postProcessGF

END MODULE m_gfcalc