MODULE m_greensfSpinDiag

   USE m_juDFT
   USE m_types
   USE m_constants

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE greensfSpinDiag(ikpt_i,nBands,i_gf,l,lp,natom,natomp,atomType,atomTypep,spin,&
                              l_sphavg,elementPhase,sym,atoms,usdus,eigVecCoeffs,greensfBZintCoeffs)

      INTEGER,                   INTENT(IN)     :: ikpt_i !current k-point index in cdnvaljob%k_list
      INTEGER,                   INTENT(IN)     :: nBands !Bands handled on this rank
      INTEGER,                   INTENT(IN)     :: i_gf,l,lp,natom,natomp,atomType,atomTypep,spin !Information about the current element
      LOGICAL,                   INTENT(IN)     :: l_sphavg
      COMPLEX,                   INTENT(IN)     :: elementPhase
      TYPE(t_sym),               INTENT(IN)     :: sym
      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_usdus),             INTENT(IN)     :: usdus
      TYPE(t_eigVecCoeffs),      INTENT(IN)     :: eigVecCoeffs
      TYPE(t_greensfBZintCoeffs),INTENT(INOUT)  :: greensfBZintCoeffs

      INTEGER :: iBand,it,is,isi,imat
      INTEGER :: m,mp,lm,lmp,ilo,ilop
      REAL    :: fac
      COMPLEX, ALLOCATABLE :: im(:,:,:)
      COMPLEX, ALLOCATABLE :: im_tmp(:,:,:)


      ALLOCATE(    im(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,MERGE(1,4,l_sphavg)),source=cmplx_0)
      ALLOCATE(im_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,MERGE(1,4,l_sphavg)),source=cmplx_0)

      fac = 1.0/(sym%invarind(natom)*atoms%neq(atomType))

      !Loop through bands
      DO iBand = 1, nBands
         im = cmplx_0
         DO m = -l, l
            lm = l*(l+1)+m
            DO mp = -lp,lp
               lmp = lp*(lp+1)+mp

               !-------------------------
               !Contribution from valence states
               !-------------------------
               IF(l_sphavg) THEN
                  im(m,mp,1) = im(m,mp,1) + conjg(eigVecCoeffs%acof(iBand,lmp,natom,spin))*eigVecCoeffs%acof(iBand,lm,natom,spin) &
                                          + conjg(eigVecCoeffs%bcof(iBand,lmp,natom,spin))*eigVecCoeffs%bcof(iBand,lm,natom,spin) &
                                          * usdus%ddn(l,atomType,spin)
               ELSE
                  im(m,mp,1) = im(m,mp,1) + conjg(eigVecCoeffs%acof(iBand,lmp,natomp,spin))*eigVecCoeffs%acof(iBand,lm,natom,spin)
                  im(m,mp,2) = im(m,mp,2) + conjg(eigVecCoeffs%bcof(iBand,lmp,natomp,spin))*eigVecCoeffs%bcof(iBand,lm,natom,spin)
                  im(m,mp,3) = im(m,mp,3) + conjg(eigVecCoeffs%acof(iBand,lmp,natomp,spin))*eigVecCoeffs%bcof(iBand,lm,natom,spin)
                  im(m,mp,4) = im(m,mp,4) + conjg(eigVecCoeffs%bcof(iBand,lmp,natomp,spin))*eigVecCoeffs%acof(iBand,lm,natom,spin)
               END IF

               !------------------------------------------------------------------------------------------------------
               ! add local orbital contribution (not implemented for radial dependence yet and not tested for average)
               !------------------------------------------------------------------------------------------------------
               DO ilo = 1, atoms%nlo(atomType)
                  IF(atoms%llo(ilo,atomType).NE.l) CYCLE
                  IF(l_sphavg) THEN
                     im(m,mp,1) = im(m,mp,1) + usdus%uulon(ilo,atomType,spin) &
                                             * ( conjg(eigVecCoeffs%acof(   iBand,lmp,natom,spin))*eigVecCoeffs%ccof(m,iBand,ilo,natom,spin) &
                                               + conjg(eigVecCoeffs%ccof(mp,iBand,ilo,natom,spin))*eigVecCoeffs%acof(  iBand,lm ,natom,spin) )&
                                             + usdus%dulon(ilo,atomType,spin) &
                                             * ( conjg(eigVecCoeffs%bcof(   iBand,lmp,natom,spin))*eigVecCoeffs%ccof(m,iBand,ilo,natom,spin) &
                                               + conjg(eigVecCoeffs%ccof(mp,iBand,ilo,natom,spin))*eigVecCoeffs%bcof(  iBand,lm ,natom,spin))
                  ENDIF
                  DO ilop = 1, atoms%nlo(atomType)
                     IF (atoms%llo(ilop,atomType).NE.l) CYCLE
                     IF(l_sphavg) THEN
                        im(m,mp,1) = im(m,mp,1) + usdus%uloulopn(ilo,ilop,atomType,spin) &
                                                * conjg(eigVecCoeffs%ccof(mp,iBand,ilop,natom,spin))*eigVecCoeffs%ccof(m,iBand,ilo,natom,spin)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO!mp
         ENDDO !m
         IF(natom.EQ.natomp.AND.l.EQ.lp) THEN !Rotations only for onsite l diagonal elements
            DO it = 1,sym%invarind(natom)
               DO imat = 1, MERGE(1,4,l_sphavg)
                  is = sym%invarop(natom,it)
                  isi = sym%invtab(is)
                  im_tmp(-l:l,-l:l,imat) = matmul( transpose( conjg(sym%d_wgn(-l:l,-l:l,l,isi)) ) , im(-l:l,-l:l,imat))
                  im_tmp(-l:l,-l:l,imat) = matmul( im_tmp(-l:l,-l:l,imat), sym%d_wgn(-l:l,-l:l,l,isi) )
                  IF(l_sphavg) THEN
                     greensfBZintCoeffs%sphavg(iBand,:,:,ikpt_i,i_gf,spin) = greensfBZintCoeffs%sphavg(iBand,:,:,ikpt_i,i_gf,spin) + CONJG(fac * elementPhase * im_tmp(:,:,imat))
                  ELSE IF(imat.EQ.1) THEN
                     greensfBZintCoeffs%uu(iBand,:,:,ikpt_i,i_gf,spin) = greensfBZintCoeffs%uu(iBand,:,:,ikpt_i,i_gf,spin) + CONJG(fac * elementPhase * im_tmp(:,:,imat))
                  ELSE IF(imat.EQ.2) THEN
                     greensfBZintCoeffs%dd(iBand,:,:,ikpt_i,i_gf,spin) = greensfBZintCoeffs%dd(iBand,:,:,ikpt_i,i_gf,spin) + CONJG(fac * elementPhase * im_tmp(:,:,imat))
                  ELSE IF(imat.EQ.3) THEN
                     greensfBZintCoeffs%ud(iBand,:,:,ikpt_i,i_gf,spin) = greensfBZintCoeffs%ud(iBand,:,:,ikpt_i,i_gf,spin) + CONJG(fac * elementPhase * im_tmp(:,:,imat))
                  ELSE IF(imat.EQ.4) THEN
                     greensfBZintCoeffs%du(iBand,:,:,ikpt_i,i_gf,spin) = greensfBZintCoeffs%du(iBand,:,:,ikpt_i,i_gf,spin) + CONJG(fac * elementPhase * im_tmp(:,:,imat))
                  ENDIF
               ENDDO
            ENDDO!it
         ENDIF
      ENDDO !iBand

   END SUBROUTINE greensfSpinDiag
END MODULE m_greensfSpinDiag