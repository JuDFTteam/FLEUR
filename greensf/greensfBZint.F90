MODULE m_greensfBZint

   USE m_types
   USE m_juDFT
   USE m_constants
   USE m_greensfEigVecCoeffs
   USE m_symMMPmat
   USE m_rotMMPmat

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE greensfBZint(ikpt,nBands,spin_ind,gfinp,sym,atoms,noco,nococonv,input,kpts,&
                           scalarGF,eigVecCoeffs,greensfBZintCoeffs)

      INTEGER,                   INTENT(IN)     :: ikpt        ! current k-point
      INTEGER,                   INTENT(IN)     :: nBands      !Bands handled on this rank
      INTEGER,                   INTENT(IN)     :: spin_ind    !spin index
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

      INTEGER :: i_gf,l,lp,atomType,atomTypep,i_gf_rot
      INTEGER :: natom,natomp,natom_start,natom_end
      INTEGER :: i_elem,i_elemLO,nLO,imatSize,imat,iBand,iop
      INTEGER :: spin1,spin2,atom,atomp
      COMPLEX :: phase
      REAL    :: atomDiff(3)
      LOGICAL :: l_sphavg,l_intersite
      COMPLEX, ALLOCATABLE :: im(:,:,:,:)
      COMPLEX, ALLOCATABLE :: imSym(:,:)
      TYPE(t_eigVecCoeffs) :: eigVecCoeffs_rot
      TYPE(t_gfelementtype) :: rep_elem

      IF(spin_ind==3) THEN
         spin1 = 2
         spin2 = 1
      ELSE
         spin1 = spin_ind
         spin2 = spin_ind
      ENDIF

      call greensfBZintCoeffs%reset()
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

            if (l_intersite) then
               natomp = atomp
            else
               natomp = natom
            endif

            DO iop = 1, MERGE(1, sym%nop, .NOT.l_intersite)

               IF(l_intersite) THEN
                  i_gf_rot = gfinp%find_symmetry_rotated_greensf(atoms,sym,i_gf,iop,distinct_kresolved_int=.false.)
                  i_elem   = gfinp%uniqueElements(atoms,max_index=i_gf_rot,l_sphavg=l_sphavg)
                  i_elemLO = gfinp%uniqueElements(atoms,max_index=i_gf_rot,lo=.TRUE.,l_sphavg=l_sphavg)
               ENDIF

               CALL greensfEigVecCoeffs(nBands,l,lp,natom,natomp,atomType,atomTypep,spin1,spin2,&
                                          l_sphavg,atoms,scalarGF(i_gf),eigVecCoeffs_rot,im)


               !The eigenvector coefficients already contains part of the interstitial phase
               !but not necessarily the right one
               IF(natom/=natomp) THEN
                  im = im * exp(-tpi_const*ImagUnit*dot_product(kpts%bk(:,ikpt),  atoms%taual(:,natom) &
                                                                                 - atoms%taual(:,natomp)))
               ENDIF

               IF(spin_ind<3) THEN
                  im = conjg(im)
               ELSE
                  im = -im
               ENDIF

#ifndef CPP_NOTYPEPROCINOMP
               !$omp parallel default(none) &
               !$omp shared(sym,kpts,atoms,greensfBZintCoeffs,im,nococonv,noco,atomType,imatSize,nBands) &
               !$omp shared(l_intersite,l,lp,natom,natomp,spin_ind,iop,atomDiff,i_elem,i_elemLO,ikpt,nLO,l_sphavg) &
               !$omp private(imat,iBand,imSym)
#endif
               ALLOCATE(imSym(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const),source=cmplx_0)
#ifndef CPP_NOTYPEPROCINOMP
               !$omp do collapse(2)
#endif
               DO imat = 1, imatSize
                  DO iBand = 1, nBands

                     IF(l_intersite) THEN
                        imSym = ImagUnit**(l-lp)/REAL(sym%nop) * symMMPmat(im(:,:,iBand,imat),sym,natom,l,lp=lp,phase=(spin_ind.EQ.3),&
                                                                           bk=kpts%bk(:,ikpt),atomDiff=atomDiff,&
                                                                           taualdiff=atoms%taual(:,natom) - atoms%taual(:,natomp),sym_op_list=[iop])
                     ELSE
                        imSym = ImagUnit**(l-lp)/atoms%neq(atomType) * symMMPmat(im(:,:,iBand,imat),sym,natom,l,lp=lp,phase=(spin_ind.EQ.3))
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
                     CALL greensfBZintCoeffs%add_contribution(i_elem, i_elemLO, iBand, nLO, imat, l_sphavg, imSym)
                     !$omp end critical
#else
                     CALL greensfBZintCoeffs%add_contribution(i_elem, i_elemLO, iBand, nLO, imat, l_sphavg, imSym)
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
         ENDDO !natom
         DEALLOCATE(im)
      ENDDO !i_gf
      CALL timestop("Green's Function: Brillouin-Zone-Integration")


   END SUBROUTINE greensfBZint

END MODULE m_greensfBZint
