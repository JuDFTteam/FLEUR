MODULE m_greensfCalcScalarProducts

   USE m_types
   USE m_constants
   USE m_juDFT
   USE m_genMTBasis

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE greensfCalcScalarProducts(gfinp,atoms,input,enpara,noco,sphhar,vTot,fmpi,hub1data,scalarProducts,greensFunctions,fout,gout,floout)
      !-----------------------------------------------------------
      ! Calculate the needed radial functions and scalar products
      !-----------------------------------------------------------
      TYPE(t_gfinp),             INTENT(IN)     :: gfinp
      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_input),             INTENT(IN)     :: input
      TYPE(t_mpi),               INTENT(IN)     :: fmpi
      TYPE(t_potden),            INTENT(IN)     :: vTot
      TYPE(t_enpara),            INTENT(IN)     :: enpara
      TYPE(t_noco),              INTENT(IN)     :: noco
      TYPE(t_sphhar),            INTENT(IN)     :: sphhar
      TYPE(t_hub1data), OPTIONAL,INTENT(IN)     :: hub1data

      TYPE(t_scalarGF), OPTIONAL, ALLOCATABLE, INTENT(OUT) :: scalarProducts(:)
      TYPE(t_greensf), OPTIONAL, INTENT(INOUT) :: greensFunctions(:)
      REAL, OPTIONAL, ALLOCATABLE :: fout(:,:,:,:,:),gout(:,:,:,:,:), floout(:,:,:,:,:)

      INTEGER  i_gf,l,lp,atomType,atomTypep,jspin,i,indUnique,ispin
      LOGICAL  l_sphavg

      REAL, ALLOCATABLE :: f(:,:,:,:,:),g(:,:,:,:,:), flo(:,:,:,:,:)

      TYPE(t_usdus)            :: usdus
      TYPE(t_denCoeffsOffDiag) :: denCoeffsOffDiag
      TYPE(t_scalarGF), ALLOCATABLE :: scalarGF(:)

      IF(.NOT.PRESENT(scalarProducts).AND..NOT.PRESENT(greensFunctions)) THEN
         CALL juDFT_error('Either list of scalarGF or greensf have to be provided', calledby='greensfCalcScalarProducts')
      ENDIF

      CALL timestart("Green's Function: Radial Functions")
      ALLOCATE (f(atoms%jmtd,2,0:atoms%lmaxd,input%jspins,atoms%nType),source=0.0)
      ALLOCATE (g(atoms%jmtd,2,0:atoms%lmaxd,input%jspins,atoms%nType),source=0.0)
      ALLOCATE (flo(atoms%jmtd,2,atoms%nlod,input%jspins,atoms%nType),source=0.0)

      ! Initializations
      CALL usdus%init(atoms,input%jspins)
      CALL denCoeffsOffDiag%init(atoms,noco,sphhar,.FALSE.,.FALSE.)

      ALLOCATE(scalarGF(gfinp%n))
      !Generate the scalar products we need
      DO i_gf = 1, gfinp%n
         l  = gfinp%elem(i_gf)%l
         lp = gfinp%elem(i_gf)%lp
         atomType  = gfinp%elem(i_gf)%atomType
         atomTypep = gfinp%elem(i_gf)%atomTypep
         CALL scalarGF(i_gf)%init(atoms,input)

         IF(.NOT.gfinp%isUnique(i_gf)) THEN
            indUnique = gfinp%getuniqueElement(i_gf)
            scalarGF(i_gf) = scalarGF(indUnique)
            CYCLE
         ENDIF
         DO jspin = 1, input%jspins
            CALL genMTBasis(atoms,enpara,vTot,fmpi,atomType,jspin,usdus,&
                            f(:,:,:,jspin,atomType),g(:,:,:,jspin,atomType),flo(:,:,:,jspin,atomType),&
                            hub1data=hub1data,l_writeArg=.FALSE.)
            IF(atomType/=atomTypep) THEN
               CALL genMTBasis(atoms,enpara,vTot,fmpi,atomTypep,jspin,usdus,&
                               f(:,:,:,jspin,atomTypep),g(:,:,:,jspin,atomTypep),flo(:,:,:,jspin,atomTypep),&
                               hub1data=hub1data,l_writeArg=.FALSE.)
            ENDIF
         ENDDO
         IF(gfinp%elem(i_gf)%isOffDiag()) THEN
            CALL scalarGF(i_gf)%addOffdScalarProduct(l,lp,atomType,atomTypep,gfinp%elem(i_gf)%isIntersite(),&
                                                     gfinp%l_mperp,atoms,input,f,g,flo)
         ELSE
            DO ispin = 1, input%jspins
               scalarGF(i_gf)%uun(ispin,ispin) = 1.0
               scalarGF(i_gf)%dun(ispin,ispin) = 0.0
               scalarGF(i_gf)%udn(ispin,ispin) = 0.0
               scalarGF(i_gf)%ddn(ispin,ispin) = usdus%ddn(l,atomType,ispin)

               scalarGF(i_gf)%uulon(:,ispin,ispin) = usdus%uulon(:,atomType,ispin)
               scalarGF(i_gf)%uloun(:,ispin,ispin) = usdus%uulon(:,atomType,ispin)
               scalarGF(i_gf)%dulon(:,ispin,ispin) = usdus%dulon(:,atomType,ispin)
               scalarGF(i_gf)%ulodn(:,ispin,ispin) = usdus%dulon(:,atomType,ispin)

               scalarGF(i_gf)%uloulopn(:,:,ispin,ispin) = usdus%uloulopn(:,:,atomType,ispin)

            ENDDO
            IF(gfinp%l_mperp) THEN
               CALL denCoeffsOffDiag%addRadFunScalarProducts(atoms,f(:,:,:,:,atomType),g(:,:,:,:,atomType),&
                                                          flo(:,:,:,:,atomType),atomType)
               IF(atomType/=atomTypep) THEN
                  CALL denCoeffsOffDiag%addRadFunScalarProducts(atoms,f(:,:,:,:,atomTypep),g(:,:,:,:,atomTypep),&
                                                                flo(:,:,:,:,atomTypep),atomTypep)
               ENDIF
               scalarGF(i_gf)%uun(1,2) = denCoeffsOffDiag%uu21n(l,atomType)
               scalarGF(i_gf)%uun(2,1) = denCoeffsOffDiag%uu21n(l,atomType)
               scalarGF(i_gf)%dun(1,2) = denCoeffsOffDiag%du21n(l,atomType)
               scalarGF(i_gf)%dun(2,1) = denCoeffsOffDiag%du21n(l,atomType)
               scalarGF(i_gf)%udn(1,2) = denCoeffsOffDiag%ud21n(l,atomType)
               scalarGF(i_gf)%udn(2,1) = denCoeffsOffDiag%ud21n(l,atomType)
               scalarGF(i_gf)%ddn(1,2) = denCoeffsOffDiag%dd21n(l,atomType)
               scalarGF(i_gf)%ddn(2,1) = denCoeffsOffDiag%dd21n(l,atomType)

               scalarGF(i_gf)%uulon(:,1,2) = denCoeffsOffDiag%uulo21n(:,atomType)
               scalarGF(i_gf)%uulon(:,2,1) = denCoeffsOffDiag%uulo21n(:,atomType)
               scalarGF(i_gf)%uloun(:,1,2) = denCoeffsOffDiag%ulou21n(:,atomType)
               scalarGF(i_gf)%uloun(:,2,1) = denCoeffsOffDiag%ulou21n(:,atomType)
               scalarGF(i_gf)%dulon(:,1,2) = denCoeffsOffDiag%dulo21n(:,atomType)
               scalarGF(i_gf)%dulon(:,2,1) = denCoeffsOffDiag%dulo21n(:,atomType)
               scalarGF(i_gf)%ulodn(:,1,2) = denCoeffsOffDiag%ulod21n(:,atomType)
               scalarGF(i_gf)%ulodn(:,2,1) = denCoeffsOffDiag%ulod21n(:,atomType)

               scalarGF(i_gf)%uloulopn(:,:,1,2) = denCoeffsOffDiag%uloulop21n(:,:,atomType)
               scalarGF(i_gf)%uloulopn(:,:,2,1) = denCoeffsOffDiag%uloulop21n(:,:,atomType)
            ENDIF
         ENDIF
      ENDDO

      IF(PRESENT(greensFunctions)) THEN
         DO i_gf = 1, gfinp%n
            greensFunctions(i_gf)%scalarProducts = scalarGF(i_gf)
         ENDDO
      ENDIF
      IF(PRESENT(scalarProducts)) THEN
         CALL move_alloc(scalarGF, scalarProducts)
      ENDIF

      IF(PRESENT(fout)) CALL move_alloc(f,fout)
      IF(PRESENT(gout)) CALL move_alloc(g,gout)
      IF(PRESENT(floout)) CALL move_alloc(flo,floout)

      CALL timestop("Green's Function: Radial Functions")

   END SUBROUTINE greensfCalcScalarProducts
END MODULE m_greensfCalcScalarProducts