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

      INTEGER :: iBand,imat
      INTEGER :: m,mp,lm,lmp,ilo,ilop

      im = cmplx_0
      !Loop through bands
      DO iBand = 1, nBands
         DO m = -l, l
            lm = l*(l+1)+m
            DO mp = -lp,lp
               lmp = lp*(lp+1)+mp

               !-------------------------
               !Contribution from valence states
               !-------------------------
               IF(l_sphavg) THEN
                  im(m,mp,iBand,1) = im(m,mp,iBand,1) + conjg(eigVecCoeffs%acof(iBand,lmp,natom,spin))*eigVecCoeffs%acof(iBand,lm,natom,spin) &
                                          + conjg(eigVecCoeffs%bcof(iBand,lmp,natom,spin))*eigVecCoeffs%bcof(iBand,lm,natom,spin) &
                                          * usdus%ddn(l,atomType,spin)
               ELSE
                  im(m,mp,iBand,1) = im(m,mp,iBand,1) + conjg(eigVecCoeffs%acof(iBand,lmp,natomp,spin))*eigVecCoeffs%acof(iBand,lm,natom,spin)
                  im(m,mp,iBand,2) = im(m,mp,iBand,2) + conjg(eigVecCoeffs%bcof(iBand,lmp,natomp,spin))*eigVecCoeffs%bcof(iBand,lm,natom,spin)
                  im(m,mp,iBand,3) = im(m,mp,iBand,3) + conjg(eigVecCoeffs%acof(iBand,lmp,natomp,spin))*eigVecCoeffs%bcof(iBand,lm,natom,spin)
                  im(m,mp,iBand,4) = im(m,mp,iBand,4) + conjg(eigVecCoeffs%bcof(iBand,lmp,natomp,spin))*eigVecCoeffs%acof(iBand,lm,natom,spin)
               END IF

               !------------------------------------------------------------------------------------------------------
               ! add local orbital contribution (not implemented for radial dependence yet and not tested for average)
               !------------------------------------------------------------------------------------------------------
               DO ilo = 1, atoms%nlo(atomType)
                  IF(atoms%llo(ilo,atomType).NE.l) CYCLE
                  IF(l_sphavg) THEN
                     im(m,mp,iBand,1) = im(m,mp,iBand,1) + usdus%uulon(ilo,atomType,spin) &
                                             * ( conjg(eigVecCoeffs%acof(   iBand,lmp,natom,spin))*eigVecCoeffs%ccof(m,iBand,ilo,natom,spin) &
                                               + conjg(eigVecCoeffs%ccof(mp,iBand,ilo,natom,spin))*eigVecCoeffs%acof(  iBand,lm ,natom,spin) )&
                                             + usdus%dulon(ilo,atomType,spin) &
                                             * ( conjg(eigVecCoeffs%bcof(   iBand,lmp,natom,spin))*eigVecCoeffs%ccof(m,iBand,ilo,natom,spin) &
                                               + conjg(eigVecCoeffs%ccof(mp,iBand,ilo,natom,spin))*eigVecCoeffs%bcof(  iBand,lm ,natom,spin))
                  ENDIF
                  DO ilop = 1, atoms%nlo(atomType)
                     IF (atoms%llo(ilop,atomType).NE.l) CYCLE
                     IF(l_sphavg) THEN
                        im(m,mp,iBand,1) = im(m,mp,iBand,1) + usdus%uloulopn(ilo,ilop,atomType,spin) &
                                                * conjg(eigVecCoeffs%ccof(mp,iBand,ilop,natom,spin))*eigVecCoeffs%ccof(m,iBand,ilo,natom,spin)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO!mp
         ENDDO !m
      ENDDO !iBand

   END SUBROUTINE greensfSpinDiag
END MODULE m_greensfSpinDiag