MODULE m_greensfSpinDiag

   USE m_juDFT
   USE m_types
   USE m_constants

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE greensfSpinDiag(nBands,l,lp,natom,natomp,atomType,atomTypep,spin,&
                              l_sphavg,atoms,usdus,eigVecCoeffs,im)

      INTEGER,                   INTENT(IN)     :: nBands !Bands handled on this rank
      INTEGER,                   INTENT(IN)     :: l,lp,natom,natomp,atomType,atomTypep,spin !Information about the current element
      LOGICAL,                   INTENT(IN)     :: l_sphavg
      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_usdus),             INTENT(IN)     :: usdus
      TYPE(t_eigVecCoeffs),      INTENT(IN)     :: eigVecCoeffs
      COMPLEX,                   INTENT(INOUT)  :: im(-lmaxU_const:,-lmaxU_const:,:,:)

      INTEGER :: m,mp,lm,lmp,ilo,ilop,nLO,nLOp

      im = cmplx_0
      !$OMP parallel do default(none) collapse(2) &
      !$OMP shared(eigVecCoeffs,im,usdus,atoms) &
      !$OMP shared(l,lp,natom,natomp,nBands,atomType,atomTypep,spin,l_sphavg) &
      !$OMP private(m,mp,lm,lmp,ilo,ilop,nLO,nLOp)
      DO m = -l, l
         DO mp = -lp,lp
            nLO  = 0
            nLOp = 0
            lm = l*(l+1)+m
            lmp = lp*(lp+1)+mp

            !-------------------------
            !Contribution from valence states
            !-------------------------
            IF(l_sphavg) THEN
               im(m,mp,:,1) = im(m,mp,:,1) + conjg(eigVecCoeffs%acof(:nBands,lmp,natom,spin))*eigVecCoeffs%acof(:nBands,lm,natom,spin) &
                                           + conjg(eigVecCoeffs%bcof(:nBands,lmp,natom,spin))*eigVecCoeffs%bcof(:nBands,lm,natom,spin) &
                                             * usdus%ddn(l,atomType,spin)
            ELSE
               im(m,mp,:,1) = im(m,mp,:,1) + conjg(eigVecCoeffs%acof(:nBands,lmp,natomp,spin))*eigVecCoeffs%acof(:nBands,lm,natom,spin)
               im(m,mp,:,2) = im(m,mp,:,2) + conjg(eigVecCoeffs%bcof(:nBands,lmp,natomp,spin))*eigVecCoeffs%bcof(:nBands,lm,natom,spin)
               im(m,mp,:,3) = im(m,mp,:,3) + conjg(eigVecCoeffs%acof(:nBands,lmp,natomp,spin))*eigVecCoeffs%bcof(:nBands,lm,natom,spin)
               im(m,mp,:,4) = im(m,mp,:,4) + conjg(eigVecCoeffs%bcof(:nBands,lmp,natomp,spin))*eigVecCoeffs%acof(:nBands,lm,natom,spin)
            END IF

            !------------------------------------------------------------------------------------------------------
            ! add local orbital contribution (not implemented for radial dependence yet and not tested for average)
            !------------------------------------------------------------------------------------------------------
            DO ilo = 1, atoms%nlo(atomType)
               IF(atoms%llo(ilo,atomType).NE.l) CYCLE
               IF(l_sphavg) THEN
                  im(m,mp,:,1) = im(m,mp,:,1) + usdus%uulon(ilo,atomType,spin) &
                                          * ( conjg(eigVecCoeffs%acof(   :nBands,lmp,natom,spin))*eigVecCoeffs%ccof(m,:nBands,ilo,natom,spin) &
                                            + conjg(eigVecCoeffs%ccof(mp,:nBands,ilo,natom,spin))*eigVecCoeffs%acof(  :nBands,lm ,natom,spin) )&
                                          + usdus%dulon(ilo,atomType,spin) &
                                          * ( conjg(eigVecCoeffs%bcof(   :nBands,lmp,natom,spin))*eigVecCoeffs%ccof(m,:nBands,ilo,natom,spin) &
                                            + conjg(eigVecCoeffs%ccof(mp,:nBands,ilo,natom,spin))*eigVecCoeffs%bcof(  :nBands,lm ,natom,spin))
               ELSE
                  nLO = nLO + 1
                  im(m,mp,:,4+(nLO-1)*4+1) = im(m,mp,:,4+(nLO-1)*4+1) + conjg(eigVecCoeffs%acof(   :nBands,lmp,natom,spin))*eigVecCoeffs%ccof(m,:nBands,ilo,natom,spin)
                  im(m,mp,:,4+(nLO-1)*4+2) = im(m,mp,:,4+(nLO-1)*4+2) + conjg(eigVecCoeffs%ccof(mp,:nBands,ilo,natom,spin))*eigVecCoeffs%acof(  :nBands,lm ,natom,spin)
                  im(m,mp,:,4+(nLO-1)*4+3) = im(m,mp,:,4+(nLO-1)*4+3) + conjg(eigVecCoeffs%bcof(   :nBands,lmp,natom,spin))*eigVecCoeffs%ccof(m,:nBands,ilo,natom,spin)
                  im(m,mp,:,4+(nLO-1)*4+4) = im(m,mp,:,4+(nLO-1)*4+4) + conjg(eigVecCoeffs%ccof(mp,:nBands,ilo,natom,spin))*eigVecCoeffs%bcof(  :nBands,lm ,natom,spin)
               ENDIF
            ENDDO
            DO ilo = 1, atoms%nlo(atomType)
               IF(atoms%llo(ilo,atomType).NE.l) CYCLE
               DO ilop = 1, atoms%nlo(atomType)
                  IF (atoms%llo(ilop,atomType).NE.l) CYCLE
                  IF(l_sphavg) THEN
                     im(m,mp,:,1) = im(m,mp,:,1) + usdus%uloulopn(ilo,ilop,atomType,spin) &
                                             * conjg(eigVecCoeffs%ccof(mp,:nBands,ilop,natom,spin))*eigVecCoeffs%ccof(m,:nBands,ilo,natom,spin)
                  ELSE
                     nLOp = nLOp + 1
                     im(m,mp,:,4+nLO*4+nLOp) = im(m,mp,:,4+nLO*4+nLOp) + conjg(eigVecCoeffs%ccof(mp,:nBands,ilop,natom,spin))*eigVecCoeffs%ccof(m,:nBands,ilo,natom,spin)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO!mp
      ENDDO !m
      !$OMP end parallel do

   END SUBROUTINE greensfSpinDiag
END MODULE m_greensfSpinDiag