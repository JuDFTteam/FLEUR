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

   PURE FUNCTION symMMPmatFull(mmpmat,sym,natom,l,lp,atomDiff,bk,phase) Result(mmpmatSym)

      COMPLEX,          INTENT(IN)  :: mmpmat(-lmaxU_const:,-lmaxU_const:,:)
      TYPE(t_sym),      INTENT(IN)  :: sym
      INTEGER,          INTENT(IN)  :: natom
      INTEGER,          INTENT(IN)  :: l
      INTEGER,OPTIONAL, INTENT(IN)  :: lp
      REAL   ,OPTIONAL, INTENT(IN)  :: atomDiff(:)
      REAL   ,OPTIONAL, INTENT(IN)  :: bk(:)
      LOGICAL,OPTIONAL, INTENT(IN)  :: phase !multiply spin-offdiagonal phase
                                             !(if the full matrix is not given)

      COMPLEX :: mmpmatSym(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,SIZE(mmpmat,3))
      REAL    :: symFac
      INTEGER :: it,is,isi,lpArg
      COMPLEX :: symPhase

      mmpmatSym = cmplx_0

      lpArg=l
      IF(PRESENT(lp)) lpArg = lp

      symFac = 1.0/sym%invarind(natom)

      DO it = 1, sym%invarind(natom)
         is  = sym%invarop(natom,it)
         isi = sym%invtab(is)

         symPhase = cmplx_1
         IF(PRESENT(phase)) THEN
            IF(phase) symPhase = exp(ImagUnit*sym%phase(isi))
         ENDIF

         IF(PRESENT(atomDiff).AND.PRESENT(bk)) THEN
            symPhase = symPhase * exp(ImagUnit*dot_product(bk,matmul(TRANSPOSE(sym%mrot(:,:,isi)),atomDiff)))
         ENDIF

         mmpmatSym = mmpmatSym + symFac * symPhase * rotMMPmat(mmpmat,dwgn =sym%d_wgn(:,:,l.   ,isi),&
                                                                      dwgnp=sym%d_wgn(:,:,lpArg,isi))

      ENDDO

   END FUNCTION symMMPmatFull

   PURE FUNCTION symMMPmatoneSpin(mmpmat,sym,natom,l,lp,atomDiff,bk,phase) Result(mmpmatSym)

      COMPLEX,          INTENT(IN)  :: mmpmat(-lmaxU_const:,-lmaxU_const:)
      TYPE(t_sym),      INTENT(IN)  :: sym
      INTEGER,          INTENT(IN)  :: natom
      INTEGER,          INTENT(IN)  :: l
      INTEGER,OPTIONAL, INTENT(IN)  :: lp
      REAL   ,OPTIONAL, INTENT(IN)  :: atomDiff(:)
      REAL   ,OPTIONAL, INTENT(IN)  :: bk(:)
      LOGICAL,OPTIONAL, INTENT(IN)  :: phase

      INTEGER :: ilow(2),iup(2)
      COMPLEX :: mmpmatSym(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      COMPLEX, ALLOCATABLE :: mmpmatOut2(:,:,:),mmpmatIn(:,:,:)

      !Add "extra spin dimension"
      ilow = LBOUND(mmpmat)
      iup  = UBOUND(mmpmat)
      ALLOCATE(mmpmatIn(ilow(1):iup(1),ilow(2):iup(2),1),source=cmplx_0)
      mmpmatIn(:,:,1) = mmpmat

      ALLOCATE(mmpmatOut2,mold=mmpMatIn)
      mmpmatOut2 = symMMPmatFull(mmpmatIn,sym,natom,l,lp=lp,atomDiff=atomDiff,bk=bk,phase=phase)

      mmpmatSym = mmpmatOut2(:,:,1)

   END FUNCTION symMMPmatoneSpin

END MODULE m_symMMPmat