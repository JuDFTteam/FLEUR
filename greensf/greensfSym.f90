MODULE m_greensfSym

   USE m_constants
   USE m_types
   USE m_symMMPmat
   USE m_rotMMPmat

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE greensfSym(ikpt_i,i_elem,i_elemLO,nLO,atomType,natom,l,lp,l_intersite,l_sphavg,ispin,&
                         sym,atomFactor,atomDiff,bk,addPhase,noco,nococonv,im,greensfBZintCoeffs)

      TYPE(t_noco),                 INTENT(IN)     :: noco
      TYPE(t_nococonv),             INTENT(IN)     :: nococonv
      INTEGER,                      INTENT(IN)     :: ikpt_i
      INTEGER,                      INTENT(IN)     :: i_elem
      INTEGER,                      INTENT(IN)     :: i_elemLO
      INTEGER,                      INTENT(IN)     :: nLO
      INTEGER,                      INTENT(IN)     :: atomType
      INTEGER,                      INTENT(IN)     :: natom
      INTEGER,                      INTENT(IN)     :: l
      INTEGER,                      INTENT(IN)     :: lp
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
      COMPLEX, ALLOCATABLE :: imSym_tmp(:,:,:)

      !$OMP parallel default(none) &
      !$OMP shared(ikpt_i,i_elem,i_elemLO,nLO,atomType,natom,l,lp,l_intersite,l_sphavg)&
      !$OMP shared(ispin,sym,atomFactor,addPhase,bk,atomDiff,im,greensfBZintCoeffs,noco,nococonv)&
      !$OMP private(imat,iBand,imSym,imSym_tmp,iLO)
      ALLOCATE(imSym(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const),source=cmplx_0)
      ALLOCATE(imSym_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,1),source=cmplx_0)
      !$OMP do collapse(2)
      DO imat = 1, SIZE(im,4)
         DO iBand = 1, SIZE(im,3)
            IF(l_intersite) THEN
               imSym = atomFactor * addPhase * conjg(im(:,:,iBand,imat)) &
                      *exp(-tpi_const*ImagUnit*dot_product(bk,atomDiff))
            ELSE
               imSym = atomFactor * addPhase * symMMPmat(im(:,:,iBand,imat),sym,natom,l,lp=lp,phase=(ispin.EQ.3))
            ENDIF

            !Rotate into the local real frame
            IF(noco%l_noco) THEN
               imSym_tmp(:,:,1) = imSym
               imSym_tmp = rotMMPmat(imSym_tmp,0.0,-nococonv%beta(atomType),-nococonv%alph(atomType),l)
               imSym = imSym_tmp(:,:,1)
            ELSE IF(noco%l_soc) THEN
               imSym_tmp(:,:,1) = imSym
               imSym_tmp = rotMMPmat(imSym_tmp,0.0,-nococonv%theta,-nococonv%phi,l)
               imSym = imSym_tmp(:,:,1)
            ENDIF

            IF(l_sphavg) THEN
               !Spherically averaged (already multiplied with scalar products)
               greensfBZintCoeffs%sphavg(iBand,:,:,i_elem,ikpt_i,ispin) = &
                  greensfBZintCoeffs%sphavg(iBand,:,:,i_elem,ikpt_i,ispin) + imSym
            ELSE IF(imat.EQ.1) THEN
               !imat 1-4: coefficients for Valence-Valence contribution
               greensfBZintCoeffs%uu(iBand,:,:,i_elem,ikpt_i,ispin) = &
                  greensfBZintCoeffs%uu(iBand,:,:,i_elem,ikpt_i,ispin) + imSym
            ELSE IF(imat.EQ.2) THEN
               greensfBZintCoeffs%dd(iBand,:,:,i_elem,ikpt_i,ispin) = &
                  greensfBZintCoeffs%dd(iBand,:,:,i_elem,ikpt_i,ispin) + imSym
            ELSE IF(imat.EQ.3) THEN
               greensfBZintCoeffs%ud(iBand,:,:,i_elem,ikpt_i,ispin) = &
                  greensfBZintCoeffs%ud(iBand,:,:,i_elem,ikpt_i,ispin) + imSym
            ELSE IF(imat.EQ.4) THEN
               greensfBZintCoeffs%du(iBand,:,:,i_elem,ikpt_i,ispin) = &
                  greensfBZintCoeffs%du(iBand,:,:,i_elem,ikpt_i,ispin) + imSym
            ELSE IF((imat-4.0)/2.0<=nLO) THEN
               !imat 5 - 4+2*numberofLOs: coefficients for Valence-LO contribution
               iLO = CEILING(REAL(imat-4.0)/2.0)
               IF(MOD(imat-4,2)==1) THEN
                  greensfBZintCoeffs%uulo(iBand,:,:,iLO,i_elemLO,ikpt_i,ispin) = &
                     greensfBZintCoeffs%uulo(iBand,:,:,iLO,i_elemLO,ikpt_i,ispin) + imSym
               ELSE IF(MOD(imat-4,2)==0) THEN
                  greensfBZintCoeffs%dulo(iBand,:,:,iLO,i_elemLO,ikpt_i,ispin) = &
                     greensfBZintCoeffs%dulo(iBand,:,:,iLO,i_elemLO,ikpt_i,ispin) + imSym
               ENDIF
            ELSE IF((imat-4.0)/2.0<=2.0*nLO) THEN
               !imat 4+2*numberofLOs+1 - 4+4*numberofLOs: coefficients for LO-Valence contribution
               iLO = CEILING(REAL(imat-4.0-2*nLO)/2.0)
               IF(MOD(imat-4-2*nLO,2)==1) THEN
                  greensfBZintCoeffs%ulou(iBand,:,:,iLO,i_elemLO,ikpt_i,ispin) = &
                     greensfBZintCoeffs%ulou(iBand,:,:,iLO,i_elemLO,ikpt_i,ispin) + imSym
               ELSE IF(MOD(imat-4-2*nLO,2)==0) THEN
                  greensfBZintCoeffs%ulod(iBand,:,:,iLO,i_elemLO,ikpt_i,ispin) = &
                     greensfBZintCoeffs%ulod(iBand,:,:,iLO,i_elemLO,ikpt_i,ispin) + imSym
               ENDIF
            ELSE
               !imat 4+4*numberofLOs+1 - 4+4*numberofLOs+numberofLOs**2: coefficients for LO-LO contribution
               iLO = imat - 4 - 4*nLO
               greensfBZintCoeffs%uloulop(iBand,:,:,iLO,i_elemLO,ikpt_i,ispin) = &
                     greensfBZintCoeffs%uloulop(iBand,:,:,iLO,i_elemLO,ikpt_i,ispin) + imSym
            ENDIF
         ENDDO
      ENDDO
      !$OMP end do
      DEALLOCATE(imSym)
      !$OMP end parallel

   END SUBROUTINE greensfSym

END MODULE m_greensfSym