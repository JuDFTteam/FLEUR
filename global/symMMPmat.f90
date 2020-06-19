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

   PURE FUNCTION symMMPmatFull(mmpmat,sym,natom,l,phase) Result(mmpmatSym)

      COMPLEX,          INTENT(IN)  :: mmpmat(-lmaxU_const:,-lmaxU_const:,:)
      TYPE(t_sym),      INTENT(IN)  :: sym
      INTEGER,          INTENT(IN)  :: natom
      INTEGER,          INTENT(IN)  :: l
      LOGICAL,OPTIONAL, INTENT(IN)  :: phase !multiply spin-offdiagonal phase
                                             !(if the full matrix is not given)

      COMPLEX, ALLOCATABLE :: mmpmatSym(:,:,:)
      REAL    :: symFac
      INTEGER :: it,is,isi
      COMPLEX :: offdPhase


      IF(.NOT.ALLOCATED(mmpmatSym)) ALLOCATE(mmpmatSym,mold=mmpmat)
      mmpmatSym = cmplx_0

      symFac = 1.0/sym%invarind(natom)

      DO it = 1, sym%invarind(natom)
         is  = sym%invarop(natom,it)
         isi = sym%invtab(is)

         offdPhase = cmplx_1
         IF(PRESENT(phase)) THEN
            IF(phase) offdPhase = exp(ImagUnit*sym%phase(isi))
         ENDIF

         mmpmatSym = mmpmatSym + symFac * offdPhase * rotMMPmat(mmpmat,dwgn=sym%d_wgn(:,:,l,isi))

      ENDDO

   END FUNCTION symMMPmatFull

   PURE FUNCTION symMMPmatoneSpin(mmpmat,sym,natom,l,phase) Result(mmpmatSym)

      COMPLEX,          INTENT(IN)  :: mmpmat(-lmaxU_const:,-lmaxU_const:)
      TYPE(t_sym),      INTENT(IN)  :: sym
      INTEGER,          INTENT(IN)  :: natom
      INTEGER,          INTENT(IN)  :: l
      LOGICAL,OPTIONAL, INTENT(IN)  :: phase

      INTEGER :: ilow(2),iup(2)
      COMPLEX, ALLOCATABLE :: mmpmatSym(:,:)
      COMPLEX, ALLOCATABLE :: mmpmatOut2(:,:,:),mmpmatIn(:,:,:)

      !Add "extra spin dimension"
      ilow = LBOUND(mmpmat)
      iup  = UBOUND(mmpmat)
      ALLOCATE(mmpmatIn(ilow(1):iup(1),ilow(2):iup(2),1),source=cmplx_0)
      mmpmatIn(:,:,1) = mmpmat

      ALLOCATE(mmpmatOut2,mold=mmpMatIn)
      mmpmatOut2 = symMMPmatFull(mmpmatIn,sym,natom,l,phase=phase)

      IF(.NOT.ALLOCATED(mmpmatSym)) ALLOCATE(mmpmatSym,mold=mmpmat)
      mmpmatSym = mmpmatOut2(:,:,1)

   END FUNCTION symMMPmatoneSpin

END MODULE m_symMMPmat