MODULE m_crystalfieldCoeffs

   USE m_types
   USE m_constants
   USE m_juDFT
   USE m_genMTBasis
   USE m_intgr
   USE m_gaunt,ONLY:gaunt1
   USE m_mt_tofrom_grid
   USE m_xmlOutput

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE crystalfieldCoeffs(input,atoms,sphhar,sym,noco,vTot,hub1data)

      !Calculates the crystal field coefficients for 4f orbitals according to
      !J. Phys.: Condens. Matter 31, 305901 (2019)

      TYPE(t_input),       INTENT(IN)     :: input
      TYPE(t_atoms),       INTENT(IN)     :: atoms
      TYPE(t_sphhar),      INTENT(IN)     :: sphhar
      TYPE(t_sym),         INTENT(IN)     :: sym
      TYPE(t_noco),        INTENT(IN)     :: noco
      TYPE(t_potden),      INTENT(IN)     :: vTot
      TYPE(t_hub1data),    INTENT(INOUT)  :: hub1data

      INTEGER, PARAMETER :: lmax = 6
      INTEGER, PARAMETER :: lcf = 3 !For which l do we calculate the coefficients (make adjustable)

      INTEGER :: iType,ispin,l,m,lm,i_hia
      REAL    :: n_0Norm
      CHARACTER(LEN=20) :: attributes(5)

      REAL, ALLOCATABLE :: vlm(:,:,:),vTotch(:,:)
      REAL, ALLOCATABLE :: Blm(:,:,:),Alm(:,:,:)
      REAL :: n_0(atoms%jmtd)
      REAL :: alphalm(0:48)

      TYPE(t_gradients)     :: grad

      !Conversion Factor between Blm and Alm<r^l>
      alphalm(0:)= [1.0,0.0,0.0,0.0,SQRT(6.0)/2.0,0.0,1.0/2.0,0.0,SQRT(6.0)/2.0, &
                    0.0,0.0,0.0,0.0,0.0,0.0,0.0,SQRT(70.0)/8.0,-SQRT(35.0)/2.0,SQRT(10.0)/4.0,&
                    0.0,1.0/8.0,0.0,SQRT(10.0)/4.0,-SQRT(35.0)/2.0,SQRT(70.0)/8.0,0.0,0.0,0.0,0.0,&
                    0.0,0.0,0.0,0.0,0.0,0.0,0.0,SQRT(231.0)/16.0,0.0,3*SQRT(14.0)/16.0,-SQRT(105.0)/8.0,&
                    SQRT(105.0)/16.0,0.0,1./16.0,0.0,SQRT(105.0)/16.0,-SQRT(105.0)/8.0,3*SQRT(14.0)/16.0,&
                    0.0,SQRT(231.0)/16.0]


      CALL timestart('Crystal Field Coefficients')

      !initializations
      ALLOCATE(vlm(atoms%jmtd,0:MAXVAL(sphhar%llh)*(MAXVAL(sphhar%llh)+2),input%jspins),source=0.0)
      ALLOCATE(Blm(0:lmax,-lmax:lmax,input%jspins),source=0.0)
      !Question rotate to global frame ??

      CALL openXMLElementNoAttributes('crystalfieldCoefficients')

      DO iType = 1, atoms%ntype

         !                          sigma
         !Decompose potential into V(r)
         !                          lm
         CALL init_mt_grid(input%jspins, atoms, sphhar, .FALSE., sym, l_mdependency=.TRUE.)
         CALL mt_to_grid(.FALSE., input%jspins, atoms,sym,sphhar,.True.,vTot%mt(:,0:,iType,:),iType,noco,grad,vTotch)
         !modified mt_from_grid with lm index
         vlm = 0.0
         CALL mt_from_gridlm(atoms, sym, sphhar, iType, input%jspins, vTotch, vlm)
         CALL finish_mt_grid()

         !Calculate n_4f^0(r) (normed spherical part of the 4f charge density)
         n_0 = hub1data%cdn_spherical(:,lcf,iType)
         !Norm to int r^2 n_4f(r) dr = 1
         CALL intgr3(n_0,atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),n_0Norm)
         n_0 = n_0/n_0Norm

         Blm = 0.0
         DO ispin = 1, input%jspins
            !Perform integration (spin dependent because of V_xc)
            !------------------------------------------------------
            !   sigma    2l+1          R_MT               sigma
            !  B      = (----)^0.5  int    r^2 n_4f^0(r) V(r)  dr
            !   lm        4pi                             lm
            !------------------------------------------------------
            DO l = 0, lmax
               DO m = -l, l
                  lm = l*(l+1) + m
                  CALL intgr3(vlm(:,lm,ispin)*n_0(:),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),Blm(l,m,ispin))
                  Blm(l,m,ispin) = ((2*l+1)/(4*pi_const))**0.5 * Blm(l,m,ispin)
                  Alm(l,m,ispin) = alphalm(lm) * Blm(l,m,ispin)
               ENDDO
            ENDDO
         ENDDO

         DO ispin = 1, input%jspins
            CALL openXMLElementPoly('ccfCoeffs',['atomType','spin    '],[iType,ispin])
            DO l = 1, lmax
               DO m = -l, l
                  IF(ABS(Blm(l,m,ispin)).LT.1e-12.AND.ABS(Alm(l,m,ispin)).LT.1e-12) CYCLE
                  attributes = ''
                  WRITE(attributes(1),'(i0)') l
                  WRITE(attributes(2),'(i0)') m
                  WRITE(attributes(3),'(f12.7)') Blm(l,m,ispin)*hartree_to_ev_const*1000.0
                  WRITE(attributes(4),'(f12.7)') Alm(l,m,ispin)*hartree_to_ev_const*1000.0
                  WRITE(attributes(5),'(a3)') 'meV'
                  CALL writeXMLElementForm('ccf',['l   ','m   ','Blm ','Alm ','unit'],&
                                           attributes,reshape([1,1,3,3,4,1,1,12,12,3],[5,2]))
               ENDDO
            ENDDO
            CALL closeXMLElement('ccfCoeffs')
         ENDDO

         IF(atoms%n_hia > 0) THEN
            !Construct the crystal field potential in the correlated subshell
            DO i_hia = atoms%n_u + 1, atoms%n_u + atoms%n_hia
               IF(atoms%lda_u(i_hia)%atomType.NE.iType.OR.atoms%lda_u(i_hia)%l.NE.lcf) CYCLE
               CALL crystalfieldPot(input,lmax,lcf,Blm,hub1data%ccfmat(i_hia-atoms%n_u,:,:))
            ENDDO
         ENDIF

      ENDDO

      CALL closeXMLElement('crystalfieldCoefficients')

      CALL timestop('Crystal Field Coefficients')
   END SUBROUTINE crystalfieldCoeffs

   SUBROUTINE crystalfieldPot(input,lmax,lcf,Blm,ccf)

      !For Hubbard 1:
      !Construct the crystalfield potential in the correlated subshell
      !linear combinations with Gaunt coefficients
      !---------------------------------------------------
      !   l sigma                             lh  l  l'*
      !  V        = sum    int dOmega B      Y   Y  Y
      !   m m'       lh mh             lh mh  mh  m  m'
      !---------------------------------------------------
      !Here we average over spins

      TYPE(t_input),    INTENT(IN)     :: input
      INTEGER,          INTENT(IN)     :: lmax,lcf
      REAL,             INTENT(IN)     :: Blm(0:,-lmax:,:)
      REAL,             INTENT(INOUT)  :: ccf(-lmaxU_const:,-lmaxU_const:)

      INTEGER :: lcoeff,mcoeff,m,mp
      REAL    :: BlmAvg, coef

      ccf = 0.0
      DO lcoeff = 0, lmax
         DO mcoeff = -lcoeff, lcoeff
            IF(ALL(ABS(Blm(lcoeff,mcoeff,:)).LT.1e-12)) CYCLE
            BlmAvg = (Blm(lcoeff,mcoeff,1) + Blm(lcoeff,mcoeff,input%jspins))/2.0
            DO m = -lcf,lcf
               DO mp = -lcf,lcf
                  coef = gaunt1(lcf,lcoeff,lcf,m,mcoeff,mp,lmax)
                  IF(ABS(coef).LT.1e-12) CYCLE
                  ccf(m,mp) = ccf(m,mp) + coef * BlmAvg
               ENDDO
            ENDDO
         ENDDO
      ENDDO


   END SUBROUTINE crystalfieldPot

END MODULE m_crystalfieldCoeffs