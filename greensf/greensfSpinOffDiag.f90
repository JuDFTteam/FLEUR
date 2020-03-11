MODULE m_greensfSpinOffDiag

   USE m_juDFT
   USE m_types
   USE m_constants

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE greensfSpinOffDiag(ikpt_i,nBands,i_gf,l,lp,natom,natomp,atomType,atomTypep,spin1,spin2,&
                                 l_sphavg,sym,atoms,denCoeffsOffdiag,eigVecCoeffs,greensfBZintCoeffs)

      INTEGER,                   INTENT(IN)     :: ikpt_i !current k-point index in cdnvaljob%k_list
      INTEGER,                   INTENT(IN)     :: nBands !Bands handled on this rank
      INTEGER,                   INTENT(IN)     :: i_gf,l,lp,natom,natomp,atomType,atomTypep,spin1,spin2 !Information about the current element
      LOGICAL,                   INTENT(IN)     :: l_sphavg
      TYPE(t_sym),               INTENT(IN)     :: sym
      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_denCoeffsOffDiag),  INTENT(IN)     :: denCoeffsOffdiag
      TYPE(t_eigVecCoeffs),      INTENT(IN)     :: eigVecCoeffs
      TYPE(t_greensfBZintCoeffs),INTENT(INOUT)  :: greensfBZintCoeffs

      INTEGER :: iBand,it,is,isi,imat
      INTEGER :: m,mp,lm,lmp,ilo,ilop
      REAL    :: fac
      COMPLEX :: phase
      COMPLEX, ALLOCATABLE :: im(:,:,:)
      COMPLEX, ALLOCATABLE :: im_tmp(:,:,:)

      CALL timestart("Green's Function: Spin-OffDiagonal")

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
                  im(m,mp,1) = im(m,mp,1) + conjg(eigVecCoeffs%acof(iBand,lmp,natom,spin1))*eigVecCoeffs%acof(iBand,lm,natom,spin2) * denCoeffsOffdiag%uu21n(l,atomType) &
                                          + conjg(eigVecCoeffs%acof(iBand,lmp,natom,spin1))*eigVecCoeffs%bcof(iBand,lm,natom,spin2) * denCoeffsOffdiag%ud21n(l,atomType) &
                                          + conjg(eigVecCoeffs%bcof(iBand,lmp,natom,spin1))*eigVecCoeffs%acof(iBand,lm,natom,spin2) * denCoeffsOffdiag%du21n(l,atomType) &
                                          + conjg(eigVecCoeffs%bcof(iBand,lmp,natom,spin1))*eigVecCoeffs%bcof(iBand,lm,natom,spin2) * denCoeffsOffdiag%dd21n(l,atomType)
               ELSE
                  im(m,mp,1) = im(m,mp,1) + conjg(eigVecCoeffs%acof(iBand,lmp,natomp,spin1))*eigVecCoeffs%acof(iBand,lm,natom,spin2)
                  im(m,mp,2) = im(m,mp,2) + conjg(eigVecCoeffs%bcof(iBand,lmp,natomp,spin1))*eigVecCoeffs%bcof(iBand,lm,natom,spin2)
                  im(m,mp,3) = im(m,mp,3) + conjg(eigVecCoeffs%acof(iBand,lmp,natomp,spin1))*eigVecCoeffs%bcof(iBand,lm,natom,spin2)
                  im(m,mp,4) = im(m,mp,4) + conjg(eigVecCoeffs%bcof(iBand,lmp,natomp,spin1))*eigVecCoeffs%acof(iBand,lm,natom,spin2)
               END IF

               !------------------------------------------------------------------------------------------------------
               ! add local orbital contribution (not implemented for radial dependence yet and not tested for average)
               !------------------------------------------------------------------------------------------------------
               DO ilo = 1, atoms%nlo(atomType)
                  IF(atoms%llo(ilo,atomType).NE.l) CYCLE
                  IF(l_sphavg) THEN
                     im(m,mp,1) = im(m,mp,1) + conjg(eigVecCoeffs%acof(   iBand,lmp,natom,spin1))*eigVecCoeffs%ccof(m,iBand,ilo,natom,spin2) * denCoeffsOffDiag%uulo21n(ilo,atomType) &
                                             + conjg(eigVecCoeffs%ccof(mp,iBand,ilo,natom,spin1))*eigVecCoeffs%acof(  iBand,lm ,natom,spin2) * denCoeffsOffDiag%ulou21n(ilo,atomType) &
                                             + conjg(eigVecCoeffs%bcof(   iBand,lmp,natom,spin1))*eigVecCoeffs%ccof(m,iBand,ilo,natom,spin2) * denCoeffsOffDiag%dulo21n(ilo,atomType) &
                                             + conjg(eigVecCoeffs%ccof(mp,iBand,ilo,natom,spin1))*eigVecCoeffs%bcof(  iBand,lm ,natom,spin2) * denCoeffsOffDiag%ulod21n(ilo,atomType)
                  ENDIF
                  DO ilop = 1, atoms%nlo(atomType)
                     IF (atoms%llo(ilop,atomType).NE.l) CYCLE
                     IF(l_sphavg) THEN
                        im(m,mp,1) = im(m,mp,1) + conjg(eigVecCoeffs%ccof(mp,iBand,ilop,natom,spin1))*eigVecCoeffs%ccof(m,iBand,ilo,natom,spin2) * denCoeffsOffDiag%uloulop21n(ilo,ilop,atomType)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO!mp
         ENDDO !m

         IF(.FALSE.) THEN !Rotations do not work for offdiagonal elements
            CALL timestart("GF Rotations")
            DO it = 1,sym%invarind(natom)
               DO imat = 1, MERGE(1,4,l_sphavg)
                  is = sym%invarop(natom,it)
                  isi = sym%invtab(is)
                  phase = exp(ImagUnit*sym%phase(isi))
                  im_tmp(:,:,imat) = matmul( transpose( conjg(sym%d_wgn(:,:,l,isi)) ) , im(:,:,imat))
                  im_tmp(:,:,imat) = matmul( im_tmp(:,:,imat), sym%d_wgn(:,:,l,isi) )
                  IF(l_sphavg) THEN
                     greensfBZintCoeffs%sphavg(iBand,:,:,ikpt_i,i_gf,3) = greensfBZintCoeffs%sphavg(iBand,:,:,ikpt_i,i_gf,3) + CONJG(fac * phase * im_tmp(:,:,imat))
                  ELSE IF(imat.EQ.1) THEN
                     greensfBZintCoeffs%uu(iBand,:,:,ikpt_i,i_gf,3) = greensfBZintCoeffs%uu(iBand,:,:,ikpt_i,i_gf,3) + CONJG(fac * phase * im_tmp(:,:,imat))
                  ELSE IF(imat.EQ.2) THEN
                     greensfBZintCoeffs%dd(iBand,:,:,ikpt_i,i_gf,3) = greensfBZintCoeffs%dd(iBand,:,:,ikpt_i,i_gf,3) + CONJG(fac * phase * im_tmp(:,:,imat))
                  ELSE IF(imat.EQ.3) THEN
                     greensfBZintCoeffs%ud(iBand,:,:,ikpt_i,i_gf,3) = greensfBZintCoeffs%ud(iBand,:,:,ikpt_i,i_gf,3) + CONJG(fac * phase * im_tmp(:,:,imat))
                  ELSE IF(imat.EQ.4) THEN
                     greensfBZintCoeffs%du(iBand,:,:,ikpt_i,i_gf,3) = greensfBZintCoeffs%du(iBand,:,:,ikpt_i,i_gf,3) + CONJG(fac * phase * im_tmp(:,:,imat))
                  ENDIF
               ENDDO
            ENDDO!it
            CALL timestop("GF Rotations")
         ELSE
            DO imat = 1, MERGE(1,4,l_sphavg)
               IF(l_sphavg) THEN
                  greensfBZintCoeffs%sphavg(iBand,:,:,ikpt_i,i_gf,3) = greensfBZintCoeffs%sphavg(iBand,:,:,ikpt_i,i_gf,3) + im(:,:,imat)
               ELSE IF(imat.EQ.1) THEN
                  greensfBZintCoeffs%uu(iBand,:,:,ikpt_i,i_gf,3) = greensfBZintCoeffs%uu(iBand,:,:,ikpt_i,i_gf,3) + im(:,:,imat)
               ELSE IF(imat.EQ.2) THEN
                  greensfBZintCoeffs%dd(iBand,:,:,ikpt_i,i_gf,3) = greensfBZintCoeffs%dd(iBand,:,:,ikpt_i,i_gf,3) + im(:,:,imat)
               ELSE IF(imat.EQ.3) THEN
                  greensfBZintCoeffs%ud(iBand,:,:,ikpt_i,i_gf,3) = greensfBZintCoeffs%ud(iBand,:,:,ikpt_i,i_gf,3) + im(:,:,imat)
               ELSE IF(imat.EQ.4) THEN
                  greensfBZintCoeffs%du(iBand,:,:,ikpt_i,i_gf,3) = greensfBZintCoeffs%du(iBand,:,:,ikpt_i,i_gf,3) + im(:,:,imat)
               ENDIF
            ENDDO
         ENDIF
      ENDDO !iBand
      CALL timestop("Green's Function: Spin-OffDiagonal")

   END SUBROUTINE greensfSpinOffDiag
END MODULE m_greensfSpinOffDiag