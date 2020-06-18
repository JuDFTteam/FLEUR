MODULE m_greensfSpinOffDiag

   USE m_juDFT
   USE m_types
   USE m_constants

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE greensfSpinOffDiag(nBands,l,lp,natom,natomp,atomType,atomTypep,spin1,spin2,&
                                 l_sphavg,atoms,denCoeffsOffdiag,eigVecCoeffs,im)

      INTEGER,                   INTENT(IN)     :: nBands !Bands handled on this rank
      INTEGER,                   INTENT(IN)     :: l,lp,natom,natomp,atomType,atomTypep,spin1,spin2 !Information about the current element
      LOGICAL,                   INTENT(IN)     :: l_sphavg
      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_denCoeffsOffDiag),  INTENT(IN)     :: denCoeffsOffdiag
      TYPE(t_eigVecCoeffs),      INTENT(IN)     :: eigVecCoeffs
      COMPLEX,                   INTENT(INOUT)  :: im(-lmaxU_const:,-lmaxU_const:,:,:)

      INTEGER :: m,mp,lm,lmp,ilo,ilop

      im = cmplx_0
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP SHARED(eigVecCoeffs,im,denCoeffsOffdiag,atoms) &
      !$OMP SHARED(l,lp,natom,natomp,nBands,atomType,atomTypep,spin1,spin2,l_sphavg) &
      !$OMP PRIVATE(m,mp,lm,lmp,ilo,ilop) COLLAPSE(2)
      DO m = -l, l
         DO mp = -lp,lp
            lm = l*(l+1)+m
            lmp = lp*(lp+1)+mp

            !-------------------------
            !Contribution from valence states
            !-------------------------
            IF(l_sphavg) THEN
               im(m,mp,:,1) = im(m,mp,:,1) + conjg(eigVecCoeffs%acof(:nBands,lmp,natom,spin1))*eigVecCoeffs%acof(:nBands,lm,natom,spin2) * denCoeffsOffdiag%uu21n(l,atomType) &
                                       + conjg(eigVecCoeffs%acof(:nBands,lmp,natom,spin1))*eigVecCoeffs%bcof(:nBands,lm,natom,spin2) * denCoeffsOffdiag%ud21n(l,atomType) &
                                       + conjg(eigVecCoeffs%bcof(:nBands,lmp,natom,spin1))*eigVecCoeffs%acof(:nBands,lm,natom,spin2) * denCoeffsOffdiag%du21n(l,atomType) &
                                       + conjg(eigVecCoeffs%bcof(:nBands,lmp,natom,spin1))*eigVecCoeffs%bcof(:nBands,lm,natom,spin2) * denCoeffsOffdiag%dd21n(l,atomType)
            ELSE
               im(m,mp,:,1) = im(m,mp,:,1) + conjg(eigVecCoeffs%acof(:nBands,lmp,natomp,spin1))*eigVecCoeffs%acof(:nBands,lm,natom,spin2)
               im(m,mp,:,2) = im(m,mp,:,2) + conjg(eigVecCoeffs%bcof(:nBands,lmp,natomp,spin1))*eigVecCoeffs%bcof(:nBands,lm,natom,spin2)
               im(m,mp,:,3) = im(m,mp,:,3) + conjg(eigVecCoeffs%acof(:nBands,lmp,natomp,spin1))*eigVecCoeffs%bcof(:nBands,lm,natom,spin2)
               im(m,mp,:,4) = im(m,mp,:,4) + conjg(eigVecCoeffs%bcof(:nBands,lmp,natomp,spin1))*eigVecCoeffs%acof(:nBands,lm,natom,spin2)
            END IF

            !------------------------------------------------------------------------------------------------------
            ! add local orbital contribution (not implemented for radial dependence yet and not tested for average)
            !------------------------------------------------------------------------------------------------------
            DO ilo = 1, atoms%nlo(atomType)
               IF(atoms%llo(ilo,atomType).NE.l) CYCLE
               IF(l_sphavg) THEN
                  im(m,mp,:,1) = im(m,mp,:,1) + conjg(eigVecCoeffs%acof(   :nBands,lmp,natom,spin1))*eigVecCoeffs%ccof(m,:nBands,ilo,natom,spin2) * denCoeffsOffDiag%uulo21n(ilo,atomType) &
                                          + conjg(eigVecCoeffs%ccof(mp,:nBands,ilo,natom,spin1))*eigVecCoeffs%acof(  :nBands,lm ,natom,spin2) * denCoeffsOffDiag%ulou21n(ilo,atomType) &
                                          + conjg(eigVecCoeffs%bcof(   :nBands,lmp,natom,spin1))*eigVecCoeffs%ccof(m,:nBands,ilo,natom,spin2) * denCoeffsOffDiag%dulo21n(ilo,atomType) &
                                          + conjg(eigVecCoeffs%ccof(mp,:nBands,ilo,natom,spin1))*eigVecCoeffs%bcof(  :nBands,lm ,natom,spin2) * denCoeffsOffDiag%ulod21n(ilo,atomType)
               ENDIF
               DO ilop = 1, atoms%nlo(atomType)
                  IF (atoms%llo(ilop,atomType).NE.l) CYCLE
                  IF(l_sphavg) THEN
                     im(m,mp,:,1) = im(m,mp,:,1) + conjg(eigVecCoeffs%ccof(mp,:nBands,ilop,natom,spin1))*eigVecCoeffs%ccof(m,:nBands,ilo,natom,spin2) * denCoeffsOffDiag%uloulop21n(ilo,ilop,atomType)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO!mp
      ENDDO !m
      !$OMP END PARALLEL DO

   END SUBROUTINE greensfSpinOffDiag
END MODULE m_greensfSpinOffDiag