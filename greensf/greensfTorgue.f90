MODULE m_greensfTorgue

   USE m_types
   USE m_juDFT
   USE m_constants
   USE m_intgr
   USE m_gaunt
   USE m_xmlOutput
   USE m_lattHarmsSphHarmsConv

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE greensfTorgue(greensFunction,sphhar,atoms,sym,noco,nococonv,input,f,g,flo,atomType,torgue,vTot)

      !--------------------------------------------------------------------------
      ! This Subroutine implements the formula:
      !   alpha     1                  ->    Ef        i ->    alpha     ->->
      !  J      = - -  Im Tr      int dx  int   dE    B (x) sigma    G  (x,x,E)
      !   i         pi      sigma           -infinity  xc             ii
      !
      ! For the evaluation of the torgue at site i (argument atomType)
      !--------------------------------------------------------------------------

      TYPE(t_greensf),        INTENT(IN)     :: greensFunction(:)
      TYPE(t_sphhar),         INTENT(IN)     :: sphhar
      TYPE(t_atoms),          INTENT(IN)     :: atoms
      TYPE(t_sym),            INTENT(IN)     :: sym
      TYPE(t_noco),           INTENT(IN)     :: noco
      TYPE(t_nococonv),       INTENT(IN)     :: nococonv
      TYPE(t_input),          INTENT(IN)     :: input
      REAL,                   INTENT(IN)     :: f(:,:,0:,:,:)
      REAL,                   INTENT(IN)     :: g(:,:,0:,:,:)
      REAL,                   INTENT(IN)     :: flo(:,:,:,:,:)
      INTEGER,                INTENT(IN)     :: atomType
      REAL,                   INTENT(INOUT)  :: torgue(:)
      TYPE(t_potden),         INTENT(IN)     :: vTot

      INTEGER :: i_gf,l,lp,iContour,iGrid,ispin
      INTEGER :: lh,mu,m,mp,iz,ipm,jr,alpha,lhmu
      COMPLEX :: phaseFactor, weight
      REAL    :: realIntegral, imagIntegral
      COMPLEX :: sigma(2,2,3),torgue_cmplx(3),g_Spin(2,2)
      CHARACTER(LEN=20) :: attributes(5)

      COMPLEX, ALLOCATABLE :: bxc(:,:)
      COMPLEX, ALLOCATABLE :: integrand(:)
      COMPLEX, ALLOCATABLE :: g_ii(:,:,:,:)
      COMPLEX, ALLOCATABLE :: vlm(:,:,:)

      CALL timestart("Green's Function Torgue: init")
      !Get Bxc from the total potential (local frame)
      !TODO: FFN components
      ALLOCATE(vlm(atoms%jmtd,atoms%lmaxd*(atoms%lmaxd+2)+1,input%jspins),source=cmplx_0)
      vlm = cmplx_0
      DO ispin = 1, input%jspins
         CALL lattHarmsRepToSphHarms(sym, atoms, sphhar, atomType, vTot%mt(:,0:,atomType,ispin), vlm(:,:,ispin))
      ENDDO
      !Get the Bxc part of the potential
      ALLOCATE(bxc(SIZE(vlm,1),SIZE(vlm,2)))
      bxc = (vlm(:,:,1) - vlm(:,:,2))/2.0
      DEALLOCATE(vlm)

      !L=0 of potential has an additional rescaling of r/sqrt(4pi)
      bxc(:,1) = bxc(:,1) * sfp_const/atoms%rmsh(:atoms%jri(atomType),atomType)

      ! sigma are the Pauli matrices
      sigma=cmplx_0
      sigma(1,2,1)=CMPLX(1.0,0.0)
      sigma(2,1,1)=CMPLX(1.0,0.0)
      sigma(1,2,2)=CMPLX(0.0,-1.0)
      sigma(2,1,2)=CMPLX(0.0,1.0)
      sigma(1,1,3)=CMPLX(1.0,0.0)
      sigma(2,2,3)=CMPLX(-1.0,0.0)

      CALL timestop("Green's Function Torgue: init")
      CALL timestart("Green's Function Torgue: Integration")
      torgue_cmplx = cmplx_0
      iContour = -1
      DO i_gf = 1, SIZE(greensFunction)

         IF(greensFunction(i_gf)%elem%atomType.NE.atomType.OR.&
            greensFunction(i_gf)%elem%atomTypep.NE.atomType) THEN
            CALL juDFT_error("Provided greensFunction for wrong atomType", calledby="greensFunctionTorgue")
         ENDIF

         l  = greensFunction(i_gf)%elem%l
         lp = greensFunction(i_gf)%elem%lp

         IF(iContour == -1) THEN
            iContour = greensFunction(i_gf)%elem%iContour
         ELSE IF(greensFunction(i_gf)%elem%iContour/=iContour) THEN
            CALL juDFT_error("Provided different energy contours", calledby="greensFunctionTorgue")
         ENDIF

!         !$OMP parallel default(none) &
!         !$OMP shared(sphhar,atoms,greensFunction,f,g,flo,sigma,bxc) &
!         !$OMP shared(l,lp,i_gf,atomType,torgue_cmplx) &
!         !$OMP private(lh,m,mu,mp,lhmu,phaseFactor,ipm,iz,alpha,jr) &
!         !$OMP private(realIntegral,imagIntegral,integrand,g_ii,g_Spin)
         ALLOCATE(integrand(atoms%jmtd),source=cmplx_0)
         ALLOCATE(g_ii(2,2,atoms%jmtd,greensFunction(i_gf)%contour%nz),source=cmplx_0)
!         !$OMP do collapse(2) reduction(+:torgue_cmplx)
         DO lh = 0, atoms%lmaxd
            DO m = -l, l
               IF(MOD(lh+l+lp,2) .NE. 0) CYCLE
               IF(lh.GT.l+lp) CYCLE
               IF(lh.LT.abs(l-lp)) CYCLE
               DO mu = -lh, lh
                  lhmu = lh * (lh+1) + mu + 1
                  mp = m - mu
                  IF(ABS(mp).GT.lp) CYCLE
                  phaseFactor = gaunt1(lp,lh,l,mp,mu,m,atoms%lmaxd)
                  IF(ABS(phaseFactor).LT.1e-12) CYCLE
                  DO ipm = 1, 2
                     CALL greensFunction(i_gf)%getRadialSpin(atoms,m,mp,ipm==2,f,g,flo,g_ii)
                     DO iz = 1, SIZE(g_ii,4)
                        weight = greensFunction(i_gf)%contour%de(iz) * phaseFactor
                        DO alpha = 1, 3 !(x,y,z)
                           DO jr = 1, atoms%jri(atomType)
                              IF(ipm == 1) THEN
                                 g_Spin = matmul(sigma(:,:,alpha),g_ii(:,:,jr,iz))
                                 integrand(jr) = (g_Spin(1,1) + g_Spin(2,2)) * bxc(jr,lhmu)
                              ELSE
                                 g_Spin = matmul(conjg(sigma(:,:,alpha)),g_ii(:,:,jr,iz))
                                 integrand(jr) = (g_Spin(1,1) + g_Spin(2,2)) * conjg(bxc(jr,lhmu))
                              ENDIF
                           ENDDO
                           CALL intgr3(REAL(integrand(:)),atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),realIntegral)
                           CALL intgr3(AIMAG(integrand),atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),imagIntegral)
                           torgue_cmplx(alpha) = torgue_cmplx(alpha) - 1/(2*ImagUnit*pi_const) * (-1)**(ipm-1) * (realIntegral+ImagUnit*imagIntegral) &
                                                * MERGE(weight,conjg(weight),ipm.EQ.1)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
!         !$OMP end do
         DEALLOCATE(integrand,g_ii)
!         !$OMP end parallel

      ENDDO
      torgue = REAL(torgue_cmplx)
      CALL timestop("Green's Function Torgue: Integration")

      WRITE(oUnit,'(A,I4,A,3f14.8,A)') '  atom: ', atomType, '   torgue: ', torgue * hartree_to_ev_const * 1000, ' meV'

      attributes = ''
      WRITE(attributes(1),'(i0)') atomType
      WRITE(attributes(2),'(f14.8)') torgue(1) * hartree_to_ev_const * 1000
      WRITE(attributes(3),'(f14.8)') torgue(2) * hartree_to_ev_const * 1000
      WRITE(attributes(4),'(f14.8)') torgue(3) * hartree_to_ev_const * 1000
      WRITE(attributes(5),'(a3)') 'meV'
      CALL writeXMLElementForm('torgue',['atomType','sigma_x ','sigma_y ','sigma_z ','unit    '],&
                               attributes,reshape([8,7,7,7,4,6,14,14,14,3],[5,2]))


   END SUBROUTINE greensfTorgue

   SUBROUTINE greensfSOTorgue(greensFunction,sphhar,atoms,sym,noco,nococonv,input,enpara,mpi,f,g,flo,atomType,torgue,vTot)

      USE m_sointg
      USE m_spnorb
      USE m_fourProduct

      !--------------------------------------------------------------------------
      ! This Subroutine implements the formula:
      !   alpha     1                  ->    Ef        i ->    alpha     ->->
      !  J      = - -  Im Tr      int dx  int   dE    V (x) sigma    G  (x,x,E)
      !   i         pi      sigma           -infinity  SO             ii
      !
      ! For the evaluation of the torgue at site i (argument atomType) caused by Spin-orbit-coupling
      !--------------------------------------------------------------------------

      TYPE(t_greensf),        INTENT(IN)     :: greensFunction(:)
      TYPE(t_atoms),          INTENT(IN)     :: atoms
      TYPE(t_sphhar),         INTENT(IN)     :: sphhar
      TYPE(t_sym),            INTENT(IN)     :: sym
      TYPE(t_noco),           INTENT(IN)     :: noco
      TYPE(t_nococonv),       INTENT(IN)     :: nococonv
      TYPE(t_input),          INTENT(IN)     :: input
      TYPE(t_enpara),         INTENT(IN)     :: enpara
      TYPE(t_mpi),            INTENT(IN)     :: mpi
      REAL,                   INTENT(IN)     :: f(:,:,0:,:,:)
      REAL,                   INTENT(IN)     :: g(:,:,0:,:,:)
      REAL,                   INTENT(IN)     :: flo(:,:,:,:,:)
      INTEGER,                INTENT(IN)     :: atomType
      REAL,                   INTENT(INOUT)  :: torgue(:)
      TYPE(t_potden),         INTENT(IN)     :: vTot

      INTEGER :: jspin,na,nsym,nh,i_gf,l,lp,spin,iContour,i
      INTEGER :: lh,mh,mhp,m,mp,iz,ipm,jr,alpha,jspin1,jspin2
      COMPLEX :: phaseFactor
      REAL    :: realIntegral, imagIntegral, e
      COMPLEX :: sigma(2,2,3),chi(2,2),torgue_cmplx(3),g_Spin(2,2)
      CHARACTER(LEN=20) :: attributes(5)

      REAL :: v0(atoms%jmtd),vso_tmp(atoms%jmtd,2)
      COMPLEX, ALLOCATABLE :: integrand(:),g_iiSpin(:,:,:,:)
      COMPLEX,ALLOCATABLE :: soangl(:,:,:,:,:,:)
      COMPLEX,ALLOCATABLE :: vso(:,:,:,:,:,:)


      !
      !---> common spin-orbit integrant V   (average spin directions)
      !                                  SO
      ALLOCATE(soangl(atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,2,atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,2))
      ALLOCATE(vso(atoms%jmtd,2,2,-atoms%lmaxd:atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,0:atoms%lmaxd))
      IF(.NOT.noco%l_noco) THEN
         CALL spnorb_angles(atoms,mpi,nococonv%beta(atomType),nococonv%alph(atomType),soangl)
      ELSE
         CALL spnorb_angles(atoms,mpi,nococonv%theta,nococonv%phi,soangl)
      ENDIF

      DO l = 0, atoms%lmaxd
         v0(:) = 0.0
         DO i = 1,atoms%jri(atomType)
            v0(i) = (vtot%mt(i,0,atomType,1)+vtot%mt(i,0,atomType,input%jspins))/2.
         END DO
         e = (enpara%el0(l,atomType,1)+enpara%el0(l,atomType,input%jspins))/2.

         CALL sointg(atomType,e,vtot%mt(:,0,atomType,:),v0,atoms,input,vso_tmp)
         IF (.TRUE.) THEN
            DO i= 1,atoms%jmtd
               vso_tmp(i,1)= (vso_tmp(i,1)+vso_tmp(i,2))/2.
               vso_tmp(i,2)= vso_tmp(i,1)
            ENDDO
         ENDIF
         DO m = -l, l
            DO mp = -l,l
               DO jspin1 = 1,2
                  DO jspin2 = 1,2
                     vso(:,jspin1,jspin2,m,mp,l) = vso_tmp(:,1) * soangl(l,m,jspin1,l,mp,jspin2)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO



      na=SUM(atoms%neq(:atomType-1))+1

      ! sigma are the Pauli matrices
      sigma=cmplx_0
      sigma(1,2,1)=CMPLX(1.0,0.0)
      sigma(2,1,1)=CMPLX(1.0,0.0)
      sigma(1,2,2)=CMPLX(0.0,-1.0)
      sigma(2,1,2)=CMPLX(0.0,1.0)
      sigma(1,1,3)=CMPLX(1.0,0.0)
      sigma(2,2,3)=CMPLX(-1.0,0.0)

      CALL timestart("Green's Function SOTorgue: Integration")
      torgue_cmplx = cmplx_0
      iContour = -1
      ALLOCATE(integrand(atoms%jmtd),source=cmplx_0)
      DO i_gf = 1, SIZE(greensFunction)

         IF(greensFunction(i_gf)%elem%atomType.NE.atomType.OR.&
            greensFunction(i_gf)%elem%atomTypep.NE.atomType) THEN
            CALL juDFT_error("Provided greensFunction for wrong atomType", calledby="greensFunctionTorgue")
         ENDIF

         l  = greensFunction(i_gf)%elem%l
         lp = greensFunction(i_gf)%elem%lp

         IF(iContour == -1) THEN
            iContour = greensFunction(i_gf)%elem%iContour
         ELSE IF(greensFunction(i_gf)%elem%iContour/=iContour) THEN
            CALL juDFT_error("Provided different energy contours", calledby="greensFunctionTorgue")
         ENDIF

         DO lh = 0, SIZE(vso,6)-1
            DO mh = -lh,lh
               DO mhp = -lh, lh
                  IF(MAXVAL(ABS(vso(:,:,:,mh,mhp,lh))).LT.1e-12) CYCLE
                  DO m = -l, l
                     DO mp = -lp, lp
                        phaseFactor = fourProduct(lp,lh,l,lh,mp,mh,m,mhp,atoms%lmaxd)
                        IF(ABS(phaseFactor).LT.1e-12) CYCLE !Naive approach just skip all elements with zero gaunt coefficient
                        DO ipm = 1, 2
                           CALL greensFunction(i_gf)%getRadialSpin(atoms,m,mp,ipm==2,f,g,flo,g_iiSpin)
                           DO iz = 1, SIZE(g_iiSpin,4)
                              DO alpha = 1, 3 !(x,y,z)
                                 DO jr = 1, atoms%jri(atomType)
                                    IF(ipm==1) THEN
                                       g_Spin = matmul(sigma(:,:,alpha),g_iiSpin(:,:,jr,iz))
                                       g_Spin = matmul(vso(jr,:,:,mh,mhp,lh),g_Spin)
                                    ELSE
                                       g_Spin = matmul(conjg(sigma(:,:,alpha)),g_iiSpin(:,:,jr,iz))
                                       g_Spin = matmul(conjg(vso(jr,:,:,mh,mhp,lh)),g_Spin)
                                    ENDIF
                                    integrand(jr) = g_Spin(1,1) + g_Spin(2,2)
                                 ENDDO
                                 CALL intgr3(REAL(integrand(:)),atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),realIntegral)
                                 CALL intgr3(AIMAG(integrand(:)),atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),imagIntegral)
                                 torgue_cmplx(alpha) = torgue_cmplx(alpha) - 1/(2*ImagUnit*pi_const) * (-1)**(ipm-1) * (realIntegral+ImagUnit*imagIntegral) &
                                                      * MERGE(phaseFactor*greensFunction(i_gf)%contour%de(iz),conjg(phaseFactor*greensFunction(i_gf)%contour%de(iz)),ipm.EQ.1)
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

      ENDDO
      torgue = REAL(torgue_cmplx)
      CALL timestop("Green's Function SOTorgue: Integration")

      WRITE(oUnit,'(A,I4,A,3f14.8,A)') '  atom: ', atomType, '   torgue: ', torgue * hartree_to_ev_const * 1000, ' meV'

      attributes = ''
      WRITE(attributes(1),'(i0)') atomType
      WRITE(attributes(2),'(f14.8)') torgue(1) * hartree_to_ev_const * 1000
      WRITE(attributes(3),'(f14.8)') torgue(2) * hartree_to_ev_const * 1000
      WRITE(attributes(4),'(f14.8)') torgue(3) * hartree_to_ev_const * 1000
      WRITE(attributes(5),'(a3)') 'meV'
      CALL writeXMLElementForm('torgue',['atomType','sigma_x ','sigma_y ','sigma_z ','unit    '],&
                               attributes,reshape([8,7,7,7,4,6,14,14,14,3],[5,2]))


   END SUBROUTINE greensfSOTorgue
END MODULE m_greensfTorgue
