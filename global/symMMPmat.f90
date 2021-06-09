MODULE m_symMMPmat

   USE m_rotMMPmat
   USE m_types
   USE m_constants
   USE m_juDFT

   IMPLICIT NONE

   PRIVATE
   PUBLIC symMMPmat

   INTERFACE symMMPmat
      PROCEDURE :: symMMPmatFull, symMMPmatoneSpin
   END INTERFACE


   CONTAINS

   FUNCTION symMMPmatFull(mmpmat,sym,natom,l,lp,phase,atomDiff,kpt_indices,sym_op_list,kpts) Result(mmpmatSym)

      COMPLEX,                INTENT(IN)  :: mmpmat(-lmaxU_const:,-lmaxU_const:,:)
      TYPE(t_sym),            INTENT(IN)  :: sym
      INTEGER,                INTENT(IN)  :: natom
      INTEGER,                INTENT(IN)  :: l
      INTEGER,OPTIONAL,       INTENT(IN)  :: lp
      LOGICAL,OPTIONAL,       INTENT(IN)  :: phase !multiply spin-offdiagonal phase
                                                   !(if the full matrix is not given)
      REAL   ,OPTIONAL,       INTENT(IN)  :: atomDiff(:)
      INTEGER,OPTIONAL,       INTENT(IN)  :: kpt_indices(:)
      INTEGER,OPTIONAL,       INTENT(IN)  :: sym_op_list(:)
      TYPE(t_kpts), OPTIONAL, INTENT(IN)  :: kpts

      COMPLEX :: mmpmatSym(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,SIZE(mmpmat,3))
      COMPLEX :: mmpmat_kpt(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,SIZE(mmpmat,3))
      INTEGER :: i_op,n_ops,is,isi,lpArg,nkpts,ikpt,kpt,sym_kpt
      INTEGER, ALLOCATABLE :: sym_ops(:)
      COMPLEX :: symPhase, intersite_phase
      INTEGER :: rrot(3,3)
      COMPLEX :: rrot_dwgn_l(-lmaxU_const:lmaxU_const, -lmaxU_const:lmaxU_const)
      COMPLEX :: rrot_dwgn_lp(-lmaxU_const:lmaxU_const, -lmaxU_const:lmaxU_const)
      REAL    :: symFac,rotbk(3),kpt_parent(3)

      mmpmatSym = cmplx_0

      lpArg=l
      IF(PRESENT(lp)) lpArg = lp

      IF(PRESENT(sym_op_list)) THEN
         sym_ops = sym_op_list
         n_ops = SIZE(sym_op_list)
      ELSE
         n_ops = sym%invarind(natom)
         sym_ops = sym%invarop(natom,:)
      ENDIF

      IF(PRESENT(kpt_indices)) THEN
         IF(.NOT.PRESENT(kpts).OR..NOT.PRESENT(atomDiff)) THEN
            CALL juDFT_error('Not all arguments available for intersite phases',&
                             hint='This is BUG in FLEUR, please report')
         ENDIF
         nkpts = SIZE(kpt_indices)
      ELSE
         nkpts = 1
      ENDIF

      symFac = 1.0/(REAL(n_ops)*REAL(nkpts))

      DO i_op = 1, n_ops
         is  = sym_ops(i_op)
         isi = sym%invtab(is)

         symPhase = cmplx_1
         IF(PRESENT(phase)) THEN
            IF(phase) symPhase = exp(ImagUnit*sym%phase(isi))
         ENDIF

         DO ikpt = 1, nkpts

            intersite_phase = cmplx_1
            IF(PRESENT(kpt_indices)) THEN
               kpt = kpt_indices(ikpt)
               kpt_parent = kpts%bk(:,kpts%bkp(kpt))
               sym_kpt = sym%invtab(kpts%bksym(kpt))

               IF(sym_kpt.LE.sym%nop) THEN
                  rrot = transpose(sym%mrot(:,:,sym_kpt))
                  rrot_dwgn_l = transpose(sym%d_wgn(:,:,l,sym_kpt))
                  rrot_dwgn_lp = transpose(sym%d_wgn(:,:,lpArg,sym_kpt))
               ELSE
                  rrot = -transpose(sym%mrot(:,:,sym%invtab(sym_kpt-sym%nop)))
                  rrot_dwgn_l = -transpose(sym%d_wgn(:,:,l,sym%invtab(sym_kpt-sym%nop)))
                  rrot_dwgn_lp = -transpose(sym%d_wgn(:,:,lpArg,sym%invtab(sym_kpt-sym%nop)))
               ENDIF

               rotbk = matmul(rrot,kpt_parent)
               mmpmat_kpt = rotMMPmat(mmpmat,dwgn =rrot_dwgn_lp,&
                                             dwgnp=rrot_dwgn_l)

               intersite_phase = exp(-tpi_const*ImagUnit*dot_product(rotbk,matmul(sym%mrot(:,:,isi),atomDiff)))
            ELSE
               mmpmat_kpt = mmpmat
            ENDIF

            !The complex conjugation is taken from n_mat
            !It seems there is an inconsistency here that should be resolved at aome point
            mmpmatSym = mmpmatSym + symFac * symPhase * intersite_phase &
                                   * conjg(rotMMPmat(mmpmat_kpt,dwgn =sym%d_wgn(:,:,lpArg,isi),&
                                                                dwgnp=sym%d_wgn(:,:,l    ,isi)))
         ENDDO

      ENDDO

   END FUNCTION symMMPmatFull

   FUNCTION symMMPmatoneSpin(mmpmat,sym,natom,l,lp,phase,atomDiff,kpt_indices,sym_op_list,kpts) Result(mmpmatSym)

      COMPLEX,                INTENT(IN)  :: mmpmat(-lmaxU_const:,-lmaxU_const:)
      TYPE(t_sym),            INTENT(IN)  :: sym
      INTEGER,                INTENT(IN)  :: natom
      INTEGER,                INTENT(IN)  :: l
      INTEGER,OPTIONAL,       INTENT(IN)  :: lp
      REAL   ,OPTIONAL,       INTENT(IN)  :: atomDiff(:)
      INTEGER,OPTIONAL,       INTENT(IN)  :: kpt_indices(:)
      LOGICAL,OPTIONAL,       INTENT(IN)  :: phase
      INTEGER,OPTIONAL,       INTENT(IN)  :: sym_op_list(:)
      TYPE(t_kpts), OPTIONAL, INTENT(IN)  :: kpts

      INTEGER :: ilow(2),iup(2)
      COMPLEX :: mmpmatSym(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      COMPLEX, ALLOCATABLE :: mmpmatOut2(:,:,:),mmpmatIn(:,:,:)

      !Add "extra spin dimension"
      ilow = LBOUND(mmpmat)
      iup  = UBOUND(mmpmat)
      ALLOCATE(mmpmatIn(ilow(1):iup(1),ilow(2):iup(2),1),source=cmplx_0)
      mmpmatIn(:,:,1) = mmpmat

      ALLOCATE(mmpmatOut2,mold=mmpMatIn)
      mmpmatOut2 = symMMPmatFull(mmpmatIn,sym,natom,l,lp=lp,atomDiff=atomDiff,&
                                 kpt_indices=kpt_indices,phase=phase,sym_op_list=sym_op_list,&
                                 kpts=kpts)

      mmpmatSym = mmpmatOut2(:,:,1)

   END FUNCTION symMMPmatoneSpin

END MODULE m_symMMPmat