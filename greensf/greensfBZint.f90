MODULE m_greensfBZint

   USE m_types
   USE m_juDFT
   USE m_constants
   USE m_greensfEigVecCoeffs
   USE m_greensfSym

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

      INTEGER :: i_gf,i_gf_p,l,lp,atomType,atomTypep,n_op
      INTEGER :: natom,natomp,natomp_start,natomp_end,natom_start,natom_end
      INTEGER :: i_elem,i_elemLO,nLO,imatSize
      INTEGER :: spin1,spin2,ispin,spin_start,spin_end
      COMPLEX :: phase
      REAL    :: atomFactor,atomDiff(3)
      LOGICAL :: l_sphavg,l_intersite
      COMPLEX, ALLOCATABLE :: im(:,:,:,:,:)
      INTEGER :: repr_ops(gfinp%n)

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
         atomFactor = 1.0/atoms%neq(atomType)

         IF(.NOT.gfinp%isUnique(i_gf)) CYCLE

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

         IF(ANY(ABS(atomDiff).GT.1e-12)) THEN
            n_op = 1
            repr_ops(1) = 1
            ! DO i_gf_p = 1, gfinp%n
            !    IF(gfinp%elem(i_gf_p)%representative_elem==i_gf) THEN
            !       n_op = n_op + 1
            !       repr_ops(n_op) = gfinp%elem(i_gf_p)%representative_op
            !    ENDIF
            ! ENDDO
         ENDIF

         ALLOCATE(im(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,nBands,&
                     imatSize,spin_start:spin_end),source=cmplx_0)

         natom_start = SUM(atoms%neq(:atomType-1)) + 1
         natom_end   = SUM(atoms%neq(:atomType))
         !Loop over equivalent atoms
         DO natom = natom_start , natom_end

            !Only perform the second atom loop if we calculate intersite elements
            natomp_start = MERGE(natom,SUM(atoms%neq(:atomTypep-1)) + 1,atomType==atomTypep.AND.ALL(ABS(atomDiff).LT.1e-12))
            natomp_end   = MERGE(natom,SUM(atoms%neq(:atomTypep-1)) + 1,atomType==atomTypep.AND.ALL(ABS(atomDiff).LT.1e-12))

            DO natomp = natomp_start, natomp_end

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
                                           l_sphavg,atoms,scalarGF(i_gf),eigVecCoeffs,im(:,:,:,:,ispin))

                  !The eigenvector coefficients already contain part of the interstitial phase
                  !but not necessarily the right one
                  im(:,:,:,:,ispin) = im(:,:,:,:,ispin) &
                                    * exp(-tpi_const*ImagUnit*dot_product(kpts%bk(:,ikpt),  atoms%taual(:,natom) &
                                                                                          - atoms%taual(:,natomp)))

                  IF(ispin==3) THEN
                     im(:,:,:,:,ispin) = CMPLX(-REAL(im(:,:,:,:,ispin)), AIMAG(im(:,:,:,:,ispin)))
                  ENDIF

                  !l-offdiagonal phase
                  phase = ImagUnit**(l-lp)

                  CALL greensfSym(ikpt_i,ikpt,i_elem,i_elemLO,nLO,atomType,natom,l,lp,ANY(ABS(atomDiff).GT.1e-12),l_sphavg,ispin,&
                                  sym,kpts,atomFactor,atomDiff,phase,repr_ops(:n_op),noco,nococonv,im(:,:,:,:,ispin),greensfBZintCoeffs)

               ENDDO

            ENDDO !natomp
         ENDDO !natom
         DEALLOCATE(im)
      ENDDO !i_gf
      CALL timestop("Green's Function: Brillouin-Zone-Integration")


   END SUBROUTINE greensfBZint

END MODULE m_greensfBZint
