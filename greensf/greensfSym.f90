MODULE m_greensfSym

   USE m_constants
   USE m_types
   USE m_symMMPmat

   CONTAINS

   SUBROUTINE greensfSym(ikpt_i,i_elem,natom,l,l_onsite,l_sphavg,spin_start,spin_end,&
                         sym,atomFactor,phase,im,greensfBZintCoeffs)

      INTEGER,                      INTENT(IN)     :: ikpt_i
      INTEGER,                      INTENT(IN)     :: i_elem
      INTEGER,                      INTENT(IN)     :: natom
      INTEGER,                      INTENT(IN)     :: l
      LOGICAL,                      INTENT(IN)     :: l_onsite
      LOGICAL,                      INTENT(IN)     :: l_sphavg
      INTEGER,                      INTENT(IN)     :: spin_start
      INTEGER,                      INTENT(IN)     :: spin_end
      TYPE(t_sym),                  INTENT(IN)     :: sym
      REAL,                         INTENT(IN)     :: atomFactor
      COMPLEX,                      INTENT(IN)     :: phase
      COMPLEX,                      INTENT(IN)     :: im(-lmaxU_const:,-lmaxU_const:,:,:,:)
      TYPE(t_greensfBZintCoeffs),   INTENT(INOUT)  :: greensfBZintCoeffs

      INTEGER imat,iBand
      COMPLEX, ALLOCATABLE :: imSym(:,:,:)

      !$OMP parallel do default(none) &
      !$OMP shared(ikpt_i,i_elem,natom,l,l_onsite,l_sphavg)&
      !$OMP shared(spin_start,spin_end,sym,atomFactor,phase,im,greensfBZintCoeffs)&
      !$OMP private(imat,iBand,imSym) collapse(2)
      DO imat = 1, SIZE(im,4)
         DO iBand = 1, SIZE(im,3)
            IF(l_onsite) THEN !These rotations are only available for the onsite elements
               imSym = symMMPmat(im(:,:,iBand,imat,:),sym,natom,l)
            ELSE
               imSym = im(:,:,iBand,imat,:)
            ENDIF
            IF(l_sphavg) THEN
               greensfBZintCoeffs%sphavg(iBand,:,:,ikpt_i,i_elem,spin_start:spin_end) = &
                  greensfBZintCoeffs%sphavg(iBand,:,:,ikpt_i,i_elem,spin_start:spin_end) + atomFactor * phase * imSym
            ELSE IF(imat.EQ.1) THEN
               greensfBZintCoeffs%uu(iBand,:,:,ikpt_i,i_elem,spin_start:spin_end) = &
                  greensfBZintCoeffs%uu(iBand,:,:,ikpt_i,i_elem,spin_start:spin_end) + atomFactor * phase * imSym
            ELSE IF(imat.EQ.2) THEN
               greensfBZintCoeffs%dd(iBand,:,:,ikpt_i,i_elem,spin_start:spin_end) = &
                  greensfBZintCoeffs%dd(iBand,:,:,ikpt_i,i_elem,spin_start:spin_end) + atomFactor * phase * imSym
            ELSE IF(imat.EQ.3) THEN
               greensfBZintCoeffs%ud(iBand,:,:,ikpt_i,i_elem,spin_start:spin_end) = &
                  greensfBZintCoeffs%ud(iBand,:,:,ikpt_i,i_elem,spin_start:spin_end) + atomFactor * phase * imSym
            ELSE IF(imat.EQ.4) THEN
               greensfBZintCoeffs%du(iBand,:,:,ikpt_i,i_elem,spin_start:spin_end) = &
                  greensfBZintCoeffs%du(iBand,:,:,ikpt_i,i_elem,spin_start:spin_end) + atomFactor * phase * imSym
            ENDIF
         ENDDO
      ENDDO
      !$OMP end parallel do

   END SUBROUTINE greensfSym

END MODULE m_greensfSym