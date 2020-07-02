MODULE m_greensfSym

   USE m_constants
   USE m_types
   USE m_symMMPmat

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE greensfSym(ikpt_i,i_elem,natom,l,l_onsite,l_sphavg,ispin,&
                         sym,atomFactor,addPhase,im,greensfBZintCoeffs)

      INTEGER,                      INTENT(IN)     :: ikpt_i
      INTEGER,                      INTENT(IN)     :: i_elem
      INTEGER,                      INTENT(IN)     :: natom
      INTEGER,                      INTENT(IN)     :: l
      LOGICAL,                      INTENT(IN)     :: l_onsite
      LOGICAL,                      INTENT(IN)     :: l_sphavg
      INTEGER,                      INTENT(IN)     :: ispin
      TYPE(t_sym),                  INTENT(IN)     :: sym
      REAL,                         INTENT(IN)     :: atomFactor
      COMPLEX,                      INTENT(IN)     :: addPhase
      COMPLEX,                      INTENT(IN)     :: im(-lmaxU_const:,-lmaxU_const:,:,:)
      TYPE(t_greensfBZintCoeffs),   INTENT(INOUT)  :: greensfBZintCoeffs

      INTEGER imat,iBand
      COMPLEX, ALLOCATABLE :: imSym(:,:)

      !$OMP parallel default(none) &
      !$OMP shared(ikpt_i,i_elem,natom,l,l_onsite,l_sphavg)&
      !$OMP shared(ispin,sym,atomFactor,addPhase,im,greensfBZintCoeffs)&
      !$OMP private(imat,iBand,imSym)
      ALLOCATE(imSym(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const),source=cmplx_0)
      !$OMP do collapse(2)
      DO imat = 1, SIZE(im,4)
         DO iBand = 1, SIZE(im,3)
            IF(l_onsite) THEN !These rotations are only available for the onsite elements
               imSym = symMMPmat(im(:,:,iBand,imat),sym,natom,l,phase=(ispin.EQ.3))
            ELSE
               imSym = conjg(im(:,:,iBand,imat))
            ENDIF
            IF(l_sphavg) THEN
               greensfBZintCoeffs%sphavg(iBand,:,:,i_elem,ikpt_i,ispin) = &
                  greensfBZintCoeffs%sphavg(iBand,:,:,i_elem,ikpt_i,ispin) + atomFactor * addPhase * imSym
            ELSE IF(imat.EQ.1) THEN
               greensfBZintCoeffs%uu(iBand,:,:,i_elem,ikpt_i,ispin) = &
                  greensfBZintCoeffs%uu(iBand,:,:,i_elem,ikpt_i,ispin) + atomFactor * addPhase * imSym
            ELSE IF(imat.EQ.2) THEN
               greensfBZintCoeffs%dd(iBand,:,:,i_elem,ikpt_i,ispin) = &
                  greensfBZintCoeffs%dd(iBand,:,:,i_elem,ikpt_i,ispin) + atomFactor * addPhase * imSym
            ELSE IF(imat.EQ.3) THEN
               greensfBZintCoeffs%ud(iBand,:,:,i_elem,ikpt_i,ispin) = &
                  greensfBZintCoeffs%ud(iBand,:,:,i_elem,ikpt_i,ispin) + atomFactor * addPhase * imSym
            ELSE IF(imat.EQ.4) THEN
               greensfBZintCoeffs%du(iBand,:,:,i_elem,ikpt_i,ispin) = &
                  greensfBZintCoeffs%du(iBand,:,:,i_elem,ikpt_i,ispin) + atomFactor * addPhase * imSym
            ENDIF
         ENDDO
      ENDDO
      !$OMP end do
      !$OMP end parallel

   END SUBROUTINE greensfSym

END MODULE m_greensfSym