MODULE m_greensfSym

   USE m_constants
   USE m_types
   USE m_symMMPmat

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE greensfSym(ikpt_i,i_elem,i_elemLO,nLO,natom,l,l_diagonal,l_intersite,l_sphavg,ispin,&
                         sym,atomFactor,atomDiff,bk,addPhase,im,greensfBZintCoeffs)

      INTEGER,                      INTENT(IN)     :: ikpt_i
      INTEGER,                      INTENT(IN)     :: i_elem
      INTEGER,                      INTENT(IN)     :: i_elemLO
      INTEGER,                      INTENT(IN)     :: nLO
      INTEGER,                      INTENT(IN)     :: natom
      INTEGER,                      INTENT(IN)     :: l
      LOGICAL,                      INTENT(IN)     :: l_diagonal
      LOGICAL,                      INTENT(IN)     :: l_intersite
      LOGICAL,                      INTENT(IN)     :: l_sphavg
      INTEGER,                      INTENT(IN)     :: ispin
      TYPE(t_sym),                  INTENT(IN)     :: sym
      REAL,                         INTENT(IN)     :: atomFactor
      REAL,                         INTENT(IN)     :: atomDiff(:)
      REAL,                         INTENT(IN)     :: bk(:)
      COMPLEX,                      INTENT(IN)     :: addPhase
      COMPLEX,                      INTENT(IN)     :: im(-lmaxU_const:,-lmaxU_const:,:,:)
      TYPE(t_greensfBZintCoeffs),   INTENT(INOUT)  :: greensfBZintCoeffs

      INTEGER imat,iBand,iLO
      COMPLEX, ALLOCATABLE :: imSym(:,:)

      !$OMP parallel default(none) &
      !$OMP shared(ikpt_i,i_elem,i_elemLO,nLO,natom,l,l_diagonal,l_intersite,l_sphavg)&
      !$OMP shared(ispin,sym,atomFactor,addPhase,bk,atomDiff,im,greensfBZintCoeffs)&
      !$OMP private(imat,iBand,imSym,iLO)
      ALLOCATE(imSym(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const),source=cmplx_0)
      !$OMP do collapse(2)
      DO imat = 1, SIZE(im,4)
         DO iBand = 1, SIZE(im,3)
            IF(l_diagonal.AND.l_intersite) THEN !These rotations are only available for the onsite elements
               imSym = symMMPmat(im(:,:,iBand,imat),sym,natom,l,bk=bk,atomDiff=atomDiff,phase=(ispin.EQ.3))
            ELSE IF (l_diagonal) THEN
               imSym = symMMPmat(im(:,:,iBand,imat),sym,natom,l,phase=(ispin.EQ.3))
            ELSE
               imSym = addPhase * conjg(im(:,:,iBand,imat))
            ENDIF
            IF(l_sphavg) THEN
               !Spherically averaged (already multiplied with scalar products)
               greensfBZintCoeffs%sphavg(iBand,:,:,i_elem,ikpt_i,ispin) = &
                  greensfBZintCoeffs%sphavg(iBand,:,:,i_elem,ikpt_i,ispin) + atomFactor * imSym
            ELSE IF(imat.EQ.1) THEN
               !imat 1-4: coefficients for Valence-Valence contribution
               greensfBZintCoeffs%uu(iBand,:,:,i_elem,ikpt_i,ispin) = &
                  greensfBZintCoeffs%uu(iBand,:,:,i_elem,ikpt_i,ispin) + atomFactor * imSym
            ELSE IF(imat.EQ.2) THEN
               greensfBZintCoeffs%dd(iBand,:,:,i_elem,ikpt_i,ispin) = &
                  greensfBZintCoeffs%dd(iBand,:,:,i_elem,ikpt_i,ispin) + atomFactor * imSym
            ELSE IF(imat.EQ.3) THEN
               greensfBZintCoeffs%ud(iBand,:,:,i_elem,ikpt_i,ispin) = &
                  greensfBZintCoeffs%ud(iBand,:,:,i_elem,ikpt_i,ispin) + atomFactor * imSym
            ELSE IF(imat.EQ.4) THEN
               greensfBZintCoeffs%du(iBand,:,:,i_elem,ikpt_i,ispin) = &
                  greensfBZintCoeffs%du(iBand,:,:,i_elem,ikpt_i,ispin) + atomFactor * imSym
            ELSE IF((imat-4.0)/2.0<=nLO) THEN
               !imat 5 - 4+2*numberofLOs: coefficients for Valence-LO contribution
               iLO = CEILING(REAL(imat-4.0)/2.0)
               IF(MOD(imat-4,2)==1) THEN
                  greensfBZintCoeffs%uulo(iBand,:,:,iLO,i_elemLO,ikpt_i,ispin) = &
                     greensfBZintCoeffs%uulo(iBand,:,:,iLO,i_elemLO,ikpt_i,ispin) + atomFactor * imSym
               ELSE IF(MOD(imat-4,2)==0) THEN
                  greensfBZintCoeffs%dulo(iBand,:,:,iLO,i_elemLO,ikpt_i,ispin) = &
                     greensfBZintCoeffs%dulo(iBand,:,:,iLO,i_elemLO,ikpt_i,ispin) + atomFactor * imSym
               ENDIF
            ELSE IF((imat-4.0)/2.0<=2.0*nLO) THEN
               !imat 4+2*numberofLOs+1 - 4+4*numberofLOs: coefficients for LO-Valence contribution
               iLO = CEILING(REAL(imat-4.0-2*nLO)/2.0)
               IF(MOD(imat-4-2*nLO,2)==1) THEN
                  greensfBZintCoeffs%ulou(iBand,:,:,iLO,i_elemLO,ikpt_i,ispin) = &
                     greensfBZintCoeffs%ulou(iBand,:,:,iLO,i_elemLO,ikpt_i,ispin) + atomFactor * imSym
               ELSE IF(MOD(imat-4-2*nLO,2)==0) THEN
                  greensfBZintCoeffs%ulod(iBand,:,:,iLO,i_elemLO,ikpt_i,ispin) = &
                     greensfBZintCoeffs%ulod(iBand,:,:,iLO,i_elemLO,ikpt_i,ispin) + atomFactor * imSym
               ENDIF
            ELSE
               !imat 4+4*numberofLOs+1 - 4+4*numberofLOs+numberofLOs**2: coefficients for LO-LO contribution
               iLO = imat - 4 - 4*nLO
               greensfBZintCoeffs%uloulop(iBand,:,:,iLO,i_elemLO,ikpt_i,ispin) = &
                     greensfBZintCoeffs%uloulop(iBand,:,:,iLO,i_elemLO,ikpt_i,ispin) + atomFactor * imSym
            ENDIF
         ENDDO
      ENDDO
      !$OMP end do
      DEALLOCATE(imSym)
      !$OMP end parallel

   END SUBROUTINE greensfSym

END MODULE m_greensfSym