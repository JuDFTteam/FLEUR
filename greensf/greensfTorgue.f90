MODULE m_greensfTorgue

   USE m_types
   USE m_juDFT
   USE m_constants
   USE m_genMTBasis
   USE m_intgr
   USE m_gaunt

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE greensfTorgue(greensFunction,vTot,sphhar,atoms,sym,noco,nococonv,input,enpara,&
                            hub1inp,fmpi,atomType,torgue)

      !--------------------------------------------------------------------------
      ! This Subroutine implements the formula:
      !   alpha     1                  ->    Ef        i ->    alpha     ->->
      !  J      = - -  Im Tr      int dx  int   dE    B (x) sigma    G  (x,x,E)
      !   i         pi      sigma           -infinity  xc             ii
      !
      ! For the evaluation of the torgue at site i (argument atomType)
      !--------------------------------------------------------------------------

      TYPE(t_greensf),        INTENT(IN)     :: greensFunction(:)
      TYPE(t_potden),         INTENT(IN)     :: vTot
      TYPE(t_sphhar),         INTENT(IN)     :: sphhar
      TYPE(t_atoms),          INTENT(IN)     :: atoms
      TYPE(t_sym),            INTENT(IN)     :: sym
      TYPE(t_noco),           INTENT(IN)     :: noco
      TYPE(t_nococonv),       INTENT(IN)     :: nococonv
      TYPE(t_input),          INTENT(IN)     :: input
      TYPE(t_enpara),         INTENT(IN)     :: enpara
      TYPE(t_hub1inp),        INTENT(IN)     :: hub1inp
      TYPE(t_mpi),            INTENT(IN)     :: fmpi
      INTEGER,                INTENT(IN)     :: atomType
      REAL,                   INTENT(INOUT)  :: torgue(:)

      INTEGER :: jspin,na,nsym,nh,i_gf,l,lp,spin,iContour
      INTEGER :: lh,mems,mem,mu,m,mp,iz,ipm,lamda,jr,alpha
      COMPLEX :: phaseFactor
      REAL    :: realIntegral, imagIntegral
      COMPLEX :: sigma(2,2,3),chi(2,2),torgue_cmplx(3)
      REAL,    ALLOCATABLE :: bxc(:,:)
      COMPLEX, ALLOCATABLE :: g_ii(:,:),g_iiSpin(:,:,:,:)
      REAL,    ALLOCATABLE :: f(:,:,:,:), g(:,:,:,:),flo(:,:,:,:)

      TYPE(t_usdus) :: usdus

      IF(input%jspins.NE.2) CALL juDFT_error("Torgue calculation only for magnetic systems", calledby="greensFunctionTorgue")
      IF(sym%nop>1) CALL juDFT_error("Torgue calculation only without symmetries", calledby="greensFunctionTorgue")


      CALL timestart("Green's Function Torgue: init")
      !Get Bxc from the total potential (local frame)
      ALLOCATE(bxc(SIZE(vTot%mt,1),0:SIZE(vTot%mt,2)-1))
      bxc = (vTot%mt(:,:,atomType,2) - vTot%mt(:,:,atomType,1))/2.0

      ALLOCATE (f(atoms%jmtd,2,0:atoms%lmaxd,input%jspins),source=0.0)
      ALLOCATE (g(atoms%jmtd,2,0:atoms%lmaxd,input%jspins),source=0.0)
      ALLOCATE (flo(atoms%jmtd,2,atoms%nlod,input%jspins),source=0.0)
      CALL usdus%init(atoms,input%jspins)
      DO jspin = 1, input%jspins
         CALL genMTBasis(atoms,enpara,vTot,fmpi,atomType,jspin,usdus,&
                         f(:,:,:,jspin),g(:,:,:,jspin),flo(:,:,:,jspin),hub1inp%l_dftspinpol)
      ENDDO
      na=SUM(atoms%neq(:atomType-1))+1
      nsym = sym%ntypsy(na)
      nh = sphhar%nlh(nsym)

      ! sigma are the Pauli matrices
      sigma=cmplx_0
      sigma(1,2,1)=CMPLX(1.0,0.0)
      sigma(2,1,1)=CMPLX(1.0,0.0)
      sigma(1,2,2)=CMPLX(0.0,-1.0)
      sigma(2,1,2)=CMPLX(0.0,1.0)
      sigma(1,1,3)=CMPLX(1.0,0.0)
      sigma(2,2,3)=CMPLX(-1.0,0.0)

      chi(1,1) =  exp(ImagUnit*nococonv%alph(atomType)/2)*cos(nococonv%beta(atomType)/2)
      chi(1,2) = -EXP(ImagUnit*nococonv%alph(atomType)/2)*SIN(nococonv%beta(atomType)/2)
      chi(2,1) =  EXP(-ImagUnit*nococonv%alph(atomType)/2)*SIN(nococonv%beta(atomType)/2)
      chi(2,2) =  EXP(-ImagUnit*nococonv%alph(atomType)/2)*COS(nococonv%beta(atomType)/2)

      sigma(:,:,1)=MATMUL(CONJG(TRANSPOSE(chi)), MATMUL(sigma(:,:,1),chi))
      sigma(:,:,2)=MATMUL(CONJG(TRANSPOSE(chi)), MATMUL(sigma(:,:,2),chi))
      sigma(:,:,3)=MATMUL(CONJG(TRANSPOSE(chi)), MATMUL(sigma(:,:,3),chi))

      CALL timestop("Green's Function Torgue: init")
      CALL timestart("Green's Function Torgue: Integration")
      torgue_cmplx = cmplx_0
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

         DO alpha = 1, 3 !(x,y,z)
            DO lh = 0, nh
               lamda = sphhar%llh(lh,nsym)
               mems = sphhar%nmem(lh,nsym)
               DO mem = 1,mems
                  mu = sphhar%mlh(mem,lh,nsym)
                  DO m = -l, l
                     DO mp = -lp, lp
                        phaseFactor = (sphhar%clnu(mem,lh,nsym))*gaunt1(lp,lamda,l,mp,mu,m,atoms%lmaxd)
                        DO ipm = 1, 2
                           CALL greensFunction(i_gf)%getRadialSpin(m,mp,ipm==2,f,g,g_iiSpin)
                           ALLOCATE(g_ii(SIZE(g_iiSpin,1),SIZE(g_iiSpin,4)),source=cmplx_0)
                           DO iz = 1, SIZE(g_ii,2)
                              DO jr = 1, atoms%jri(atomType)
                                 g_iiSpin(jr,:,:,iz) = matmul(sigma(:,:,alpha),g_iiSpin(jr,:,:,iz))
                                 g_ii(jr,iz) = g_iiSpin(jr,1,1,iz) + g_iiSpin(jr,2,2,iz)
                              ENDDO
                              CALL intgr3(REAL(g_ii(:,iz)*bxc(:,lh)),atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),realIntegral)
                              CALL intgr3(AIMAG(g_ii(:,iz)*bxc(:,lh)),atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),imagIntegral)
                              torgue_cmplx(alpha) = torgue_cmplx(alpha) - 1/(2*ImagUnit*pi_const) * phaseFactor * (-1)**(ipm-1) * (realIntegral+ImagUnit*imagIntegral) &
                                                   * MERGE(greensFunction(i_gf)%contour%de(iz),conjg(greensFunction(i_gf)%contour%de(iz)),ipm.EQ.1)
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

      ENDDO
      torgue = REAL(torgue_cmplx)
      CALL timestop("Green's Function Torgue: Integration")

   END SUBROUTINE greensfTorgue


END MODULE m_greensfTorgue