MODULE m_symMMPmat

   USE m_rotMMPmat
   USE m_types
   USE m_constants

   IMPLICIT NONE

   PRIVATE
   PUBLIC symMMPmat

   INTERFACE symMMPmat
      PROCEDURE :: symMMPmatFull, symMMPmatoneSpin
   END INTERFACE


   CONTAINS

   PURE FUNCTION symMMPmatFull(mmpmat,sym,natom,l,lp,atomDiff,bk,phase, sym_op_list) Result(mmpmatSym)

      COMPLEX,          INTENT(IN)  :: mmpmat(-lmaxU_const:,-lmaxU_const:,:)
      TYPE(t_sym),      INTENT(IN)  :: sym
      INTEGER,          INTENT(IN)  :: natom
      INTEGER,          INTENT(IN)  :: l
      INTEGER,OPTIONAL, INTENT(IN)  :: lp
      REAL   ,OPTIONAL, INTENT(IN)  :: atomDiff(:)
      REAL   ,OPTIONAL, INTENT(IN)  :: bk(:)
      LOGICAL,OPTIONAL, INTENT(IN)  :: phase !multiply spin-offdiagonal phase
                                             !(if the full matrix is not given)
      INTEGER,OPTIONAL, INTENT(IN)  :: sym_op_list(:)

      COMPLEX :: mmpmatSym(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,SIZE(mmpmat,3))
      INTEGER :: i_op, n_ops,is,isi,lpArg
      INTEGER, ALLOCATABLE :: operations(:)
      COMPLEX :: symPhase

      mmpmatSym = cmplx_0

      lpArg=l
      IF(PRESENT(lp)) lpArg = lp

      IF(PRESENT(sym_op_list)) THEN
         operations = sym_op_list
         n_ops = SIZE(sym_op_list)
      ELSE
         n_ops = sym%invarind(natom)
         operations = sym%invarop(natom,:)
      ENDIF

      DO i_op = 1, n_ops
         is  = operations(i_op)
         isi = sym%invtab(is)

         symPhase = cmplx_1
         IF(PRESENT(phase)) THEN
            IF(phase) symPhase = exp(ImagUnit*sym%phase(isi))
         ENDIF

         IF(PRESENT(atomDiff).AND.PRESENT(bk)) THEN
            symPhase = symPhase * exp(-tpi_const*ImagUnit*dot_product(bk,matmul(sym%mrot(:,:,isi),atomDiff)))
         ENDIF

         !The complex conjugation is taken from n_mat
         !It seems there is an inconsistency here that should be resolved at aome point
         mmpmatSym = mmpmatSym + 1.0/REAL(n_ops) * symPhase * conjg(rotMMPmat(mmpmat,dwgn =sym%d_wgn(:,:,lpArg,isi),&
                                                                                     dwgnp=sym%d_wgn(:,:,l    ,isi)))

      ENDDO

   END FUNCTION symMMPmatFull

   PURE FUNCTION symMMPmatoneSpin(mmpmat,sym,natom,l,lp,atomDiff,bk,phase, sym_op_list) Result(mmpmatSym)

      COMPLEX,          INTENT(IN)  :: mmpmat(-lmaxU_const:,-lmaxU_const:)
      TYPE(t_sym),      INTENT(IN)  :: sym
      INTEGER,          INTENT(IN)  :: natom
      INTEGER,          INTENT(IN)  :: l
      INTEGER,OPTIONAL, INTENT(IN)  :: lp
      REAL   ,OPTIONAL, INTENT(IN)  :: atomDiff(:)
      REAL   ,OPTIONAL, INTENT(IN)  :: bk(:)
      LOGICAL,OPTIONAL, INTENT(IN)  :: phase
      INTEGER,OPTIONAL, INTENT(IN)  :: sym_op_list(:)

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
                                 bk=bk,phase=phase,sym_op_list=sym_op_list)

      mmpmatSym = mmpmatOut2(:,:,1)

   END FUNCTION symMMPmatoneSpin

END MODULE m_symMMPmat