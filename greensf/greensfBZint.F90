MODULE m_greensfBZint

   USE m_types
   USE m_juDFT
   USE m_constants
   USE m_greensfEigVecCoeffs
   USE m_symMMPmat
   USE m_rotMMPmat

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE greensfBZint(ikpt_i,ikpt,nBands,jspin,gfinp,sym,atoms,noco,nococonv,input,kpts,&
                           scalarGF,eigVecCoeffs,greensfBZintCoeffs)

      INTEGER,                   INTENT(IN)     :: ikpt_i,ikpt        !current k-point index in cdnvaljob%k_list and current k-point
      INTEGER,                   INTENT(IN)     :: nBands             !Bands handled on this rank
      INTEGER,                   INTENT(IN)     :: jspin              !spin index
      TYPE(t_gfinp),             INTENT(IN)     :: gfinp
      TYPE(t_sym),               INTENT(IN)     :: sym
      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_noco),              INTENT(IN)     :: noco
      TYPE(t_nococonv),          INTENT(IN)     :: nococonv
      TYPE(t_input),             INTENT(IN)     :: input
      TYPE(t_kpts),              INTENT(IN)     :: kpts
      TYPE(t_scalarGF),          INTENT(IN)     :: scalarGF(:)
      TYPE(t_eigVecCoeffs),      INTENT(IN)     :: eigVecCoeffs
      TYPE(t_greensfBZintCoeffs),INTENT(INOUT)  :: greensfBZintCoeffs

      INTEGER :: i_gf,l,lp,atomType,atomTypep
      INTEGER :: natom,natomp,natomp_start,natomp_end,natom_start,natom_end
      INTEGER :: i_elem,i_elemLO,nLO,imatSize,imat,iBand,iop
      INTEGER :: spin1,spin2,ispin,spin_start,spin_end,atom,atomp
      COMPLEX :: phase
      REAL    :: atomDiff(3)
      LOGICAL :: l_sphavg,l_intersite
      COMPLEX, ALLOCATABLE :: im(:,:,:,:)
      COMPLEX, ALLOCATABLE :: imSym(:,:)
      TYPE(t_eigVecCoeffs) :: eigVecCoeffs_rot
      TYPE(t_gfelementtype) :: rep_elem

      spin_start = MERGE(1,jspin,gfinp%l_mperp)
      spin_end   = MERGE(3,jspin,gfinp%l_mperp)

      spin_start = MERGE(1           ,spin_start,noco%l_mperp.AND..NOT.gfinp%l_mperp)
      spin_end   = MERGE(input%jspins,spin_end  ,noco%l_mperp.AND..NOT.gfinp%l_mperp)


      eigVecCoeffs_rot = eigVecCoeffs%rotate_to_rep_atom(atoms,sym,lmaxU_const)

      CALL timestart("Green's Function: Brillouin-Zone-Integration")
      DO i_gf = 1, gfinp%n

         !Get the information about the current element
         l  = gfinp%elem(i_gf)%l
         lp = gfinp%elem(i_gf)%lp
         atomType  = gfinp%elem(i_gf)%atomType
         atomTypep = gfinp%elem(i_gf)%atomTypep
         l_sphavg  = gfinp%elem(i_gf)%l_sphavg
         l_intersite = gfinp%elem(i_gf)%isIntersite()
         atomDiff = gfinp%elem(i_gf)%atomDiff
         atom = gfinp%elem(i_gf)%atom
         atomp = gfinp%elem(i_gf)%atomp

         IF(.NOT.gfinp%isUnique(i_gf, distinct_symmetry_equivalent_diffs=.TRUE.)) CYCLE
         IF(gfinp%elem(i_gf)%representative_elem>0) THEN
            rep_elem = gfinp%elem(gfinp%elem(i_gf)%representative_elem)
            IF(rep_elem%atom==atom.AND.rep_elem%atomp==atomp) CYCLE
         ENDIF

         i_elem   = gfinp%uniqueElements(atoms,max_index=i_gf,l_sphavg=l_sphavg)
         i_elemLO = gfinp%uniqueElements(atoms,max_index=i_gf,lo=.TRUE.,l_sphavg=l_sphavg)

         nLO = 0
         imatSize = 1
         IF(.NOT.l_sphavg) THEN
            imatSize = 4
            nLO = gfinp%elem(i_gf)%countLOs(atoms)
            IF(nLO/=0) THEN
               imatSize = 4+4*nLO+nLO**2
            ENDIF
         ENDIF

         ALLOCATE(im(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,nBands,&
                     imatSize),source=cmplx_0)

         natom_start = MERGE(SUM(atoms%neq(:atomType-1)) + 1,atom,.NOT.l_intersite)
         natom_end   = MERGE(SUM(atoms%neq(:atomType))      ,atom,.NOT.l_intersite)
         !Loop over equivalent atoms
         DO natom = natom_start , natom_end

            !Only perform the second atom loop if we calculate intersite elements
            natomp_start = MERGE(natom,atomp,.NOT.l_intersite)
            natomp_end   = MERGE(natom,atomp,.NOT.l_intersite)

            DO natomp = natomp_start, natomp_end

               DO iop = 1, MERGE(1, sym%nop, .NOT.l_intersite)

                  ! natom_rot = sym%mapped_atom(iop,natom)
                  ! natomp_rot = sym%mapped_atom(iop,natomp)

                  IF(l_intersite) THEN
                     i_elem = gfinp%find_symmetry_rotated_bzcoeffs(atoms,sym,i_gf,iop,l_sphavg)
                     i_elemLO = gfinp%find_symmetry_rotated_bzcoeffs(atoms,sym,i_gf,iop,l_sphavg,lo=.TRUE.)
                  ENDIF

                  DO ispin = spin_start, spin_end
                     IF(ispin==3) THEN
                        spin1 = 2
                        spin2 = 1
                     ELSE
                        spin1 = ispin
                        spin2 = ispin
                     ENDIF
                     !which scalar products for intersite and l offdiagonal(IF l_sphavg)
                     !Spin diagonal elements

                     CALL greensfEigVecCoeffs(nBands,l,lp,natom,natomp,atomType,atomTypep,spin1,spin2,&
                                              l_sphavg,atoms,scalarGF(i_gf),eigVecCoeffs_rot,im)


                     !The eigenvector coefficients already contain part of the interstitial phase
                     !but not necessarily the right one
                     IF(natom/=natomp) THEN
                        im = im * exp(-tpi_const*ImagUnit*dot_product(kpts%bk(:,ikpt),  atoms%taual(:,natom) &
                                                                                      - atoms%taual(:,natomp)))
                     ENDIF

                     IF(ispin<3) THEN
                        im = conjg(im)
                     ELSE
                        im = -im
                     ENDIF

#ifndef CPP_NOTYPEPROCINOMP
                     !$omp parallel default(none) &
                     !$omp shared(sym,kpts,atoms,greensfBZintCoeffs,im,nococonv,noco,atomType,imatSize,nBands) &
                     !$omp shared(l_intersite,l,lp,natom,ispin,iop,atomDiff,i_elem,i_elemLO,ikpt,ikpt_i,nLO,l_sphavg) &
                     !$omp private(imat,iBand,imSym)
#endif
                     ALLOCATE(imSym(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const),source=cmplx_0)
#ifndef CPP_NOTYPEPROCINOMP
                     !$omp do collapse(2)
#endif
                     DO imat = 1, imatSize
                        DO iBand = 1, nBands

                           IF(l_intersite) THEN
                              imSym = ImagUnit**(l-lp)/REAL(sym%nop) * symMMPmat(im(:,:,iBand,imat),sym,natom,l,lp=lp,phase=(ispin.EQ.3),&
                                                                                 bk=kpts%bk(:,ikpt),atomDiff=atomDiff,sym_op_list=[iop])
                           ELSE
                              imSym = ImagUnit**(l-lp)/atoms%neq(atomType) * symMMPmat(im(:,:,iBand,imat),sym,natom,l,lp=lp,phase=(ispin.EQ.3))
                           ENDIF

                           !Rotate into the local real frame
                           IF(noco%l_noco) THEN
                              IF (.NOT.l_intersite) THEN
                                 imSym = rotMMPmat(imSym,nococonv%alph(atomType),nococonv%beta(atomType),0.0,l,lp=lp,inverse=.TRUE.)
                              ENDIF
                           ELSE IF(noco%l_soc) THEN
                              imSym = rotMMPmat(imSym,nococonv%phi,nococonv%theta,0.0,l,lp=lp,inverse=.TRUE.)
                           ENDIF

#ifndef CPP_NOTYPEPROCINOMP
                           !$omp critical
                           CALL greensfBZintCoeffs%add_contribution(i_elem, i_elemLO, ikpt_i, iBand, ispin, nLO, imat, l_sphavg, imSym)
                           !$omp end critical
#else
                           CALL greensfBZintCoeffs%add_contribution(i_elem, i_elemLO, ikpt_i, iBand, ispin, nLO, imat, l_sphavg, imSym)
#endif
                        ENDDO
                     ENDDO
#ifndef CPP_NOTYPEPROCINOMP
                     !$omp end do
                     DEALLOCATE(imSym)
                     !$omp end parallel
#else
                     DEALLOCATE(imSym)
#endif

                  ENDDO
               ENDDO
            ENDDO !natomp
         ENDDO !natom
         DEALLOCATE(im)
      ENDDO !i_gf
      CALL timestop("Green's Function: Brillouin-Zone-Integration")


   END SUBROUTINE greensfBZint

END MODULE m_greensfBZint
