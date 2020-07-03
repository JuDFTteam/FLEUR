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

      INTEGER :: m,mp,lm,lmp,ilo,ilop,nLO_ind,nLOp_ind,imat

      im = cmplx_0
      IF(l_sphavg.AND.(l.NE.lp.OR.atomType.NE.atomTypep)) RETURN

      !$OMP parallel do default(none) collapse(2)&
      !$OMP shared(eigVecCoeffs,im,denCoeffsOffdiag,atoms) &
      !$OMP shared(l,lp,natom,natomp,nBands,atomType,atomTypep,spin1,spin2,l_sphavg) &
      !$OMP private(m,mp,lm,lmp,ilo,ilop,nLO_ind,nLOp_ind,imat)
      DO m = -l, l
         DO mp = -lp,lp
            nLO_ind  = 0
            nLOp_ind = 0
            lm = l*(l+1)+m
            lmp = lp*(lp+1)+mp

            !-------------------------
            !Contribution from valence states
            !-------------------------
            IF(l_sphavg) THEN
               im(m,mp,:,1) = im(m,mp,:,1) &
                             + denCoeffsOffdiag%uu21n(l,atomType) * conjg(eigVecCoeffs%acof(:nBands,lmp,natom,spin1))&
                                                                  *       eigVecCoeffs%acof(:nBands,lm ,natom,spin2) &
                             + denCoeffsOffdiag%ud21n(l,atomType) * conjg(eigVecCoeffs%acof(:nBands,lmp,natom,spin1))&
                                                                  *       eigVecCoeffs%bcof(:nBands,lm ,natom,spin2) &
                             + denCoeffsOffdiag%du21n(l,atomType) * conjg(eigVecCoeffs%bcof(:nBands,lmp,natom,spin1))&
                                                                  *       eigVecCoeffs%acof(:nBands,lm ,natom,spin2) &
                             + denCoeffsOffdiag%dd21n(l,atomType) * conjg(eigVecCoeffs%bcof(:nBands,lmp,natom,spin1))&
                                                                  *       eigVecCoeffs%bcof(:nBands,lm ,natom,spin2)
            ELSE
               im(m,mp,:,1) = im(m,mp,:,1) + conjg(eigVecCoeffs%acof(:nBands,lmp,natomp,spin1))*eigVecCoeffs%acof(:nBands,lm,natom,spin2)
               im(m,mp,:,2) = im(m,mp,:,2) + conjg(eigVecCoeffs%bcof(:nBands,lmp,natomp,spin1))*eigVecCoeffs%bcof(:nBands,lm,natom,spin2)
               im(m,mp,:,3) = im(m,mp,:,3) + conjg(eigVecCoeffs%acof(:nBands,lmp,natomp,spin1))*eigVecCoeffs%bcof(:nBands,lm,natom,spin2)
               im(m,mp,:,4) = im(m,mp,:,4) + conjg(eigVecCoeffs%bcof(:nBands,lmp,natomp,spin1))*eigVecCoeffs%acof(:nBands,lm,natom,spin2)
            END IF

            !------------------------------------------------------------------------------------------------------
            ! add local orbital contribution (not tested)
            !------------------------------------------------------------------------------------------------------
            DO ilo = 1, atoms%nlo(atomType)
               IF(atoms%llo(ilo,atomType).NE.l) CYCLE
               IF(l_sphavg) THEN
                   im(m,mp,:,1) = im(m,mp,:,1) &
                                 + denCoeffsOffDiag%uulo21n(ilo,atomType) * conjg(eigVecCoeffs%acof(   :nBands,lmp,natom,spin1))&
                                                                          *       eigVecCoeffs%ccof(m ,:nBands,ilo,natom,spin2) &
                                 + denCoeffsOffDiag%ulou21n(ilo,atomType) * conjg(eigVecCoeffs%ccof(mp,:nBands,ilo,natom,spin1))&
                                                                          *       eigVecCoeffs%acof(   :nBands,lm ,natom,spin2) &
                                 + denCoeffsOffDiag%dulo21n(ilo,atomType) * conjg(eigVecCoeffs%bcof(   :nBands,lmp,natom,spin1))&
                                                                          *       eigVecCoeffs%ccof(m ,:nBands,ilo,natom,spin2) &
                                 + denCoeffsOffDiag%ulod21n(ilo,atomType) * conjg(eigVecCoeffs%ccof(mp,:nBands,ilo,natom,spin1))&
                                                                          *       eigVecCoeffs%bcof(   :nBands,lm ,natom,spin2)
               ELSE
                  nLO_ind = nLO_ind + 1
                  imat    = 4+(nLO_ind-1)*2
                  im(m,mp,:,imat+1) = im(m,mp,:,imat+1) + conjg(eigVecCoeffs%acof(  :nBands,lmp,natomp,spin1))&
                                                              * eigVecCoeffs%ccof(m,:nBands,ilo,natom ,spin2)
                  im(m,mp,:,imat+2) = im(m,mp,:,imat+2) + conjg(eigVecCoeffs%bcof(  :nBands,lmp,natomp,spin1))&
                                                              * eigVecCoeffs%ccof(m,:nBands,ilo,natom ,spin2)
               ENDIF
            ENDDO
            IF(.NOT.l_sphavg) THEN
               DO ilo = 1, atoms%nlo(atomTypep)
                  IF(atoms%llo(ilo,atomTypep).NE.lp) CYCLE
                  nLO_ind = nLO_ind + 1
                  imat    = 4+(nLO_ind-1)*2
                  im(m,mp,:,imat+1) = im(m,mp,:,imat+1) + conjg(eigVecCoeffs%ccof(mp,:nBands,ilo,natomp,spin1))&
                                                              * eigVecCoeffs%acof(   :nBands,lm ,natom ,spin2)
                  im(m,mp,:,imat+2) = im(m,mp,:,imat+2) + conjg(eigVecCoeffs%ccof(mp,:nBands,ilo,natomp,spin1))&
                                                              * eigVecCoeffs%bcof(   :nBands,lm ,natom ,spin2)
               ENDDO
            ENDIF
            DO ilo = 1, atoms%nlo(atomType)
               IF(atoms%llo(ilo,atomType).NE.l) CYCLE
               DO ilop = 1, atoms%nlo(atomTypep)
                  IF (atoms%llo(ilop,atomTypep).NE.lp) CYCLE
                  IF(l_sphavg) THEN
                     im(m,mp,:,1) = im(m,mp,:,1) &
                                   + denCoeffsOffDiag%uloulop21n(ilo,ilop,atomType) * conjg(eigVecCoeffs%ccof(mp,:nBands,ilop,natom,spin1))&
                                                                                    *       eigVecCoeffs%ccof(m ,:nBands,ilo ,natom,spin2)
                  ELSE
                     nLOp_ind = nLOp_ind + 1
                     imat = 4+nLO_ind*2+nLOp_ind
                     im(m,mp,:,imat) = im(m,mp,:,imat) + conjg(eigVecCoeffs%ccof(mp,:nBands,ilop,natomp,spin1))&
                                                             * eigVecCoeffs%ccof(m ,:nBands,ilo ,natom ,spin2)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO!mp
      ENDDO !m
      !$OMP end parallel do

   END SUBROUTINE greensfSpinOffDiag
END MODULE m_greensfSpinOffDiag