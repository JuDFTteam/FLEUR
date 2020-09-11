MODULE m_greensfTorgue

   USE m_types
   USE m_juDFT
   USE m_constants
   USE m_intgr
   USE m_gaunt
   USE m_xmlOutput
   USE m_mt_tofrom_grid

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

      INTEGER :: na,nsym,nh,i_gf,l,lp,iContour,iGrid,ispin,lhmu
      INTEGER :: lh,mem,mu,m,mp,iz,ipm,lamda,jr,alpha
      COMPLEX :: phaseFactor, weight
      REAL    :: realIntegral, imagIntegral
      COMPLEX :: sigma(2,2,3),torgue_cmplx(3),g_Spin(2,2)
      CHARACTER(LEN=20) :: attributes(5)

      COMPLEX, ALLOCATABLE :: bxc(:,:)
      COMPLEX, ALLOCATABLE :: integrand(:),g_iiSpin(:,:,:,:)

      TYPE(t_potden) :: vTotProcess
      TYPE(t_gradients) :: grad
      COMPLEX, ALLOCATABLE :: vlm(:,:,:)
      REAL, ALLOCATABLE :: vTotch(:,:)


      IF(input%jspins.NE.2) CALL juDFT_error("Torgue calculation only for magnetic systems", calledby="greensFunctionTorgue")
      IF(sym%nop>1) CALL juDFT_warn("Torgue calculation only without symmetries", calledby="greensFunctionTorgue")

      na=SUM(atoms%neq(:atomType-1))+1
      nsym = sym%ntypsy(na)
      nh = sphhar%nlh(nsym)

      CALL timestart("Green's Function Torgue: init")
      !Get Bxc from the total potential (local frame)
      vTotProcess = vTot
      ALLOCATE(vlm(atoms%jmtd,0:MAXVAL(sphhar%llh)*(MAXVAL(sphhar%llh)+2),input%jspins),source=cmplx_0)
      CALL init_mt_grid(input%jspins, atoms, sphhar, .FALSE., sym, l_mdependency=.TRUE.)
      !                          sigma
      !Decompose potential into V(r)
      !                          lm
      DO ispin =1, input%jspins
         DO iGrid=1,atoms%jri(atomType)
            vTotProcess%mt(iGrid,:,atomType,ispin)=vTotProcess%mt(iGrid,:,atomType,ispin)*atoms%rmsh(iGrid,atomType)**2
         END DO
      ENDDO
      ALLOCATE(vTotch(atoms%nsp()*atoms%jri(atomType),input%jspins))
      CALL mt_to_grid(.FALSE., input%jspins, atoms,sym,sphhar,.True.,vTotProcess%mt(:,0:,atomType,:),atomType,noco,grad,vTotch)
      !modified mt_from_grid with lm index
      vlm = cmplx_0
      CALL mt_from_gridlm(atoms, sym, sphhar, atomType, input%jspins, vTotch, vlm)
      CALL finish_mt_grid()
      !Get the Bxc part of the potential
      ALLOCATE(bxc(SIZE(vlm,1),0:SIZE(vlm,2)-1))
      bxc = vlm(:,:,1) - vlm(:,:,2)

      DEALLOCATE(vTotch,vlm)

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
!         !$OMP shared(nh,nsym,l,lp,i_gf,atomType,torgue_cmplx) &
!         !$OMP private(lh,m,lamda,mem,mu,mp,phaseFactor,ipm,iz,alpha,jr) &
!         !$OMP private(realIntegral,imagIntegral,g_ii,g_iiSpin,g_Spin)
         ALLOCATE(integrand(atoms%jmtd),source=cmplx_0)
         ALLOCATE(g_iiSpin(2,2,atoms%jmtd,greensFunction(i_gf)%contour%nz),source=cmplx_0)
!         !$OMP do collapse(2) reduction(+:torgue_cmplx)
         DO lh = 0, nh
            DO m = -l, l
               lamda = sphhar%llh(lh,nsym)
               IF(MOD(lamda+l+lp,2) .NE. 0) CYCLE
               IF(lamda.GT.l+lp) CYCLE
               IF(lamda.LT.abs(l-lp)) CYCLE
               DO mem = 1,sphhar%nmem(lh,nsym)
                  mu = sphhar%mlh(mem,lh,nsym)
                  lhmu = lh * (lh+1) + mu
                  mp = m - mu
                  IF(ABS(mp).GT.lp) CYCLE
                  phaseFactor = gaunt1(lp,lamda,l,mp,mu,m,atoms%lmaxd)
                  IF(ABS(phaseFactor).LT.1e-12) CYCLE !Naive approach just skip all elements with zero gaunt coefficient
                  DO ipm = 1, 2
                     CALL greensFunction(i_gf)%getRadialSpin(atoms,m,mp,ipm==2,f,g,flo,g_iiSpin)
                     DO iz = 1, SIZE(g_ii,2)
                        weight = greensFunction(i_gf)%contour%de(iz)
                        DO alpha = 1, 3 !(x,y,z)
                           DO jr = 1, atoms%jri(atomType)
                              IF(ipm == 1) THEN
                                 g_Spin = matmul(sigma(:,:,alpha),g_iiSpin(:,:,jr,iz))
                                 integrand(jr) = (g_Spin(1,1) + g_Spin(2,2)) * bxc(jr,lhmu)
                              ELSE
                                 g_Spin = matmul(conjg(sigma(:,:,alpha)),g_iiSpin(:,:,jr,iz))
                                 integrand(jr) = (g_Spin(1,1) + g_Spin(2,2)) * conjg(bxc(jr,lhmu))
                              ENDIF
                           ENDDO
                           CALL intgr3(REAL(integrand(:)),atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),realIntegral)
                           CALL intgr3(AIMAG(integrand(:)),atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),imagIntegral)
                           torgue_cmplx(alpha) = torgue_cmplx(alpha) - 1/(2*ImagUnit*pi_const) * (-1)**(ipm-1) * (realIntegral+ImagUnit*imagIntegral) &
                                                * phaseFactor * MERGE(weight,conjg(weight),ipm.EQ.1)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
!         !$OMP end do
         DEALLOCATE(g_ii,g_iiSpin)
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

   SUBROUTINE greensfSOTorgue(greensFunction,sphhar,atoms,sym,noco,nococonv,input,f,g,flo,atomType,torgue,vso)

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
      REAL,                   INTENT(IN)     :: f(:,:,0:,:,:)
      REAL,                   INTENT(IN)     :: g(:,:,0:,:,:)
      REAL,                   INTENT(IN)     :: flo(:,:,:,:,:)
      INTEGER,                INTENT(IN)     :: atomType
      REAL,                   INTENT(INOUT)  :: torgue(:)
      REAL,                   INTENT(IN)     :: vso(:,0:)

      INTEGER :: jspin,na,nsym,nh,i_gf,l,lp,spin,iContour
      INTEGER :: lh,mems,mem,mh,m,mp,iz,ipm,lamda,jr,alpha
      COMPLEX :: phaseFactor
      REAL    :: realIntegral, imagIntegral
      COMPLEX :: sigma(2,2,3),chi(2,2),torgue_cmplx(3,2),g_Spin(2,2)
      CHARACTER(LEN=20) :: attributes(5)

      COMPLEX, ALLOCATABLE :: g_ii(:,:),g_iiSpin(:,:,:,:)

      TYPE(t_usdus) :: usdus

      IF(input%jspins.NE.2) CALL juDFT_error("Torgue calculation only for magnetic systems", calledby="greensFunctionTorgue")
      IF(sym%nop>1) CALL juDFT_warn("Torgue calculation only without symmetries", calledby="greensFunctionTorgue")

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
      ALLOCATE(g_ii(atoms%jmtd,greensFunction(1)%contour%nz),source=cmplx_0)
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

         DO lh = 0, SIZE(vso,2)-1
            DO mh = -lh,lh
               DO m = -l, l
                  DO mp = -lp, lp
                     phaseFactor = gaunt1(lp,lh,l,mp,mh,m,atoms%lmaxd)
                     IF(ABS(phaseFactor).LT.1e-12) CYCLE !Naive approach just skip all elements with zero gaunt coefficient
                     DO ipm = 1, 2
                        CALL greensFunction(i_gf)%getRadialSpin(atoms,m,mp,ipm==2,f,g,flo,g_iiSpin)
                        DO iz = 1, SIZE(g_ii,2)
                           DO alpha = 1, 3 !(x,y,z)
                              DO jr = 1, atoms%jri(atomType)
                                 IF(ipm==1) THEN
                                    g_Spin = matmul(sigma(:,:,alpha),g_iiSpin(:,:,jr,iz))
                                 ELSE
                                    g_Spin = matmul(conjg(sigma(:,:,alpha)),g_iiSpin(:,:,jr,iz))
                                 ENDIF
                                 g_ii(jr,iz) = g_Spin(1,1) + g_Spin(2,2)
                              ENDDO
                              CALL intgr3(REAL(g_ii(:,iz)*vso(:,lh)),atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),realIntegral)
                              CALL intgr3(AIMAG(g_ii(:,iz)*vso(:,lh)),atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),imagIntegral)
                              torgue_cmplx(alpha,ipm) = torgue_cmplx(alpha,ipm) - 1/(2*ImagUnit*pi_const) * (-1)**(ipm-1) * (realIntegral+ImagUnit*imagIntegral) &
                                                   * MERGE(phaseFactor*greensFunction(i_gf)%contour%de(iz),conjg(phaseFactor*greensFunction(i_gf)%contour%de(iz)),ipm.EQ.1)
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

      ENDDO
      WRITE(*,*) torgue_cmplx(:,1) * hartree_to_ev_const * 1000
      WRITE(*,*) torgue_cmplx(:,2) * hartree_to_ev_const * 1000
      torgue = REAL(torgue_cmplx(:,1)+torgue_cmplx(:,2))
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
