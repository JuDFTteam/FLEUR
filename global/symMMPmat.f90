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

   FUNCTION symMMPmatFull(mmpmat,sym,natom,l,lp,phase,bk,atomDiff,sym_op_list) Result(mmpmatSym)

      COMPLEX,                INTENT(IN)  :: mmpmat(-lmaxU_const:,-lmaxU_const:,:)
      TYPE(t_sym),            INTENT(IN)  :: sym
      INTEGER,                INTENT(IN)  :: natom
      INTEGER,                INTENT(IN)  :: l
      INTEGER,OPTIONAL,       INTENT(IN)  :: lp
      LOGICAL,OPTIONAL,       INTENT(IN)  :: phase !multiply spin-offdiagonal phase
                                                   !(if the full matrix is not given)
      REAL   ,OPTIONAL,       INTENT(IN)  :: bk(:)
      REAL   ,OPTIONAL,       INTENT(IN)  :: atomDiff(:)
      INTEGER,OPTIONAL,       INTENT(IN)  :: sym_op_list(:)

      COMPLEX :: mmpmatSym(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,SIZE(mmpmat,3))
      INTEGER :: i_op,n_ops,is,isi,lpArg
      INTEGER, ALLOCATABLE :: sym_ops(:)
      COMPLEX :: symPhase, intersite_phase
      REAL    :: symFac,rotbk(3),rotdiff(3),rrot_rec(3,3)

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

      IF(PRESENT(atomDiff).OR.PRESENT(bk)) THEN
         IF(.NOT.(PRESENT(atomDiff).AND.PRESENT(bk))) THEN
            CALL juDFT_error('Not all arguments available for intersite phases',&
                             hint='This is BUG in FLEUR, please report')
         ENDIF
      ENDIF

      symFac = 1.0/REAL(n_ops)

      DO i_op = 1, n_ops
         is  = sym_ops(i_op)

         symPhase = cmplx_1
         IF(PRESENT(phase)) THEN
            IF(phase) symPhase = exp(ImagUnit*sym%phase(sym%invtab(is)))
         ENDIF

         intersite_phase = cmplx_1
         IF(PRESENT(atomDiff)) THEN
            IF(is <= sym%nop) THEN
               rrot_rec = transpose(sym%mrot(:,:,sym%invtab(is)))
            ELSE
               rrot_rec = -transpose(sym%mrot(:,:,sym%invtab(is-sym%nop)))
            ENDIF

            !TODO: Add phases from non-symorphic symmetries
            rotbk = matmul(rrot_rec,bk)
            !TODO: Add phases from backfolding the kpoints
            rotbk = rotbk - CEILING(rotbk-[0.5,0.5,0.5])
            intersite_phase = exp(tpi_const*ImagUnit*dot_product(rotbk,atomDiff))
         ENDIF

         mmpmatSym = mmpmatSym + symFac * symPhase * intersite_phase &
                                   * rotMMPmat(mmpmat,sym,is,l,lp=lp,reciprocal=.TRUE.)

      ENDDO

   END FUNCTION symMMPmatFull

   FUNCTION symMMPmatoneSpin(mmpmat,sym,natom,l,lp,phase,bk,atomDiff,sym_op_list) Result(mmpmatSym)

      COMPLEX,                INTENT(IN)  :: mmpmat(-lmaxU_const:,-lmaxU_const:)
      TYPE(t_sym),            INTENT(IN)  :: sym
      INTEGER,                INTENT(IN)  :: natom
      INTEGER,                INTENT(IN)  :: l
      INTEGER,OPTIONAL,       INTENT(IN)  :: lp
      REAL   ,OPTIONAL,       INTENT(IN)  :: atomDiff(:)
      REAL   ,OPTIONAL,       INTENT(IN)  :: bk(:)
      LOGICAL,OPTIONAL,       INTENT(IN)  :: phase
      INTEGER,OPTIONAL,       INTENT(IN)  :: sym_op_list(:)

      INTEGER :: ilow(2),iup(2)
      COMPLEX :: mmpmatSym(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      COMPLEX, ALLOCATABLE :: mmpmatOut2(:,:,:),mmpmatIn(:,:,:)

      !Add "extra spin dimension"
      ilow = LBOUND(mmpmat)
      iup  = UBOUND(mmpmat)
      ALLOCATE(mmpmatIn(ilow(1):iup(1),ilow(2):iup(2),1),source=cmplx_0)
      mmpmatIn(:,:,1) = mmpmat

      ALLOCATE(mmpmatOut2,mold=mmpMatIn)
      mmpmatOut2 = symMMPmatFull(mmpmatIn,sym,natom,l,lp=lp,bk=bk,atomDiff=atomDiff,&
                                 phase=phase,sym_op_list=sym_op_list)

      mmpmatSym = mmpmatOut2(:,:,1)

   END FUNCTION symMMPmatoneSpin

END MODULE m_symMMPmat