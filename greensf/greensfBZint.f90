MODULE m_greensfBZint

   USE m_types
   USE m_juDFT
   USE m_constants
   USE m_greensfSpinDiag
   USE m_greensfSpinOffDiag
   USE m_greensfSym

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE greensfBZint(ikpt_i,ikpt,nBands,jspin,gfinp,sym,atoms,noco,input,kpts,&
                           usdus,denCoeffsOffDiag,eigVecCoeffs,greensfBZintCoeffs)

      INTEGER,                   INTENT(IN)     :: ikpt_i,ikpt        !current k-point index in cdnvaljob%k_list and current k-point
      INTEGER,                   INTENT(IN)     :: nBands             !Bands handled on this rank
      INTEGER,                   INTENT(IN)     :: jspin              !spin index
      TYPE(t_gfinp),             INTENT(IN)     :: gfinp
      TYPE(t_sym),               INTENT(IN)     :: sym
      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_noco),              INTENT(IN)     :: noco
      TYPE(t_input),             INTENT(IN)     :: input
      TYPE(t_kpts),              INTENT(IN)     :: kpts
      TYPE(t_usdus),             INTENT(IN)     :: usdus
      TYPE(t_denCoeffsOffDiag),  INTENT(IN)     :: denCoeffsOffdiag
      TYPE(t_eigVecCoeffs),      INTENT(IN)     :: eigVecCoeffs
      TYPE(t_greensfBZintCoeffs),INTENT(INOUT)  :: greensfBZintCoeffs

      INTEGER :: i_gf,l,lp,atomType,atomTypep,indUnique
      INTEGER :: natom,natomp,natomp_start,natomp_end,natom_start,natom_end
      INTEGER :: i_elem,i_elemLO,nLO,imatSize
      INTEGER :: spin1,spin2,ispin,spin_start,spin_end
      COMPLEX :: phase
      REAL    :: atomFactor,atomDiff(3)
      LOGICAL :: l_sphavg
      COMPLEX, ALLOCATABLE :: im(:,:,:,:,:)

      spin_start = MERGE(1,jspin,gfinp%l_mperp)
      spin_end   = MERGE(3,jspin,gfinp%l_mperp)

      spin_start = MERGE(1           ,spin_start,noco%l_mperp.AND..NOT.gfinp%l_mperp)
      spin_end   = MERGE(input%jspins,spin_end  ,noco%l_mperp.AND..NOT.gfinp%l_mperp)


      CALL timestart("Green's Function: Brillouin-Zone-Integration")
      DO i_gf = 1, gfinp%n

         !Get the information about the current element
         l  = gfinp%elem(i_gf)%l
         lp = gfinp%elem(i_gf)%lp
         atomType  = gfinp%elem(i_gf)%atomType
         atomTypep = gfinp%elem(i_gf)%atomTypep
         l_sphavg  = gfinp%elem(i_gf)%l_sphavg
         atomDiff(:) = gfinp%elem(i_gf)%atomDiff(:)
         atomFactor = MERGE(1.0,1.0/atoms%neq(atomType),l.NE.lp)
         atomFactor = MERGE(1.0,atomFactor,atomType.NE.atomTypep)

         i_elem   = gfinp%uniqueElements(atoms,ind=i_gf,l_sphavg=l_sphavg,indUnique=indUnique)
         i_elemLO = gfinp%uniqueElements(atoms,ind=i_gf,lo=.TRUE.,l_sphavg=l_sphavg,indUnique=indUnique)

         IF(i_gf/=indUnique) CYCLE

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
                     imatSize,spin_start:spin_end),source=cmplx_0)

         natom_start = SUM(atoms%neq(:atomType-1)) + 1
         natom_end   = MERGE(SUM(atoms%neq(:atomType-1)) + 1,SUM(atoms%neq(:atomType)),l.NE.lp)
         !Loop over equivalent atoms
         DO natom = natom_start , natom_end

            !Only perform the second atom loop if we calculate intersite elements
            !natomp_start = MERGE(natom,SUM(atoms%neq(:atomTypep-1)) + 1,atomType==atomTypep)
            !natomp_end   = MERGE(natom,SUM(atoms%neq(:atomTypep))      ,atomType==atomTypep)

            !Deactivate this loop (notice natomp_end) (only calculate intersite between representative atoms)
            natomp_start = MERGE(natom,SUM(atoms%neq(:atomTypep-1)) + 1,atomType==atomTypep)
            natomp_end   = MERGE(natom,SUM(atoms%neq(:atomTypep-1)) + 1,atomType==atomTypep)

            DO natomp = natomp_start, natomp_end

               !Phase factor for intersite elements (Does nothing atm)
               IF(ANY(ABS(atomDiff).GT.1e-12)) THEN
                  phase = exp(ImagUnit*dot_product(kpts%bk(:,ikpt),atomDiff(:)))
               ELSE
                  phase = cmplx_1
               ENDIF

               !l-offdiagonal phase
               phase = phase * ImagUnit**(l-lp)

               DO ispin = spin_start, spin_end
                  IF(ispin==3) THEN
                     spin1 = 2
                     spin2 = 1
                  ELSE
                     spin1 = ispin
                     spin2 = ispin
                  ENDIF
                  !which scalar products for intersite and l offdiagonal(IF l_sphavg)
                  !Can these be unified ?
                  !Spin diagonal elements
                  IF(spin1==spin2) THEN
                     CALL greensfSpinDiag(nBands,l,lp,natom,natomp,atomType,atomTypep,spin1,&
                                          l_sphavg,atoms,usdus,eigVecCoeffs,im(:,:,:,:,ispin))
                  ELSE
                     !Spin offdiagonal elements
                     CALL greensfSpinOffDiag(nBands,l,lp,natom,natomp,atomType,atomTypep,spin1,spin2,&
                                             l_sphavg,atoms,denCoeffsOffdiag,eigVecCoeffs,im(:,:,:,:,ispin))
                  ENDIF

                  CALL greensfSym(ikpt_i,i_elem,i_elemLO,nLO,natom,l,natom.EQ.natomp.AND.l.EQ.lp.AND.ALL(ABS(atomDiff).LT.1e-12),&
                                  l_sphavg,ispin,sym,atomFactor,phase,im(:,:,:,:,ispin),greensfBZintCoeffs)
               ENDDO

            ENDDO !natomp
         ENDDO !natom
         DEALLOCATE(im)
      ENDDO !i_gf
      CALL timestop("Green's Function: Brillouin-Zone-Integration")


   END SUBROUTINE greensfBZint

END MODULE m_greensfBZint