!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rinpXML
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!
   !!! The routine r_inpXML reads in the inp.xml file
   !!!
   !!!                               GM'16
   !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
   SUBROUTINE r_inpXML(&
      atoms,vacuum,input,stars,sliceplot,banddos,forcetheo,field,&
      cell,sym,xcpot,noco,oneD,hybrid,kpts,enpara,coreSpecInput,wann,&
      noel,namex,relcor,a1,a2,a3,dtild)

      USE iso_c_binding
      USE m_juDFT
      USE m_types
      USE m_types_forcetheo_extended
      USE m_symdata , ONLY : nammap, ord2, l_c2
      !USE m_rwsymfile
      !USE m_xmlIntWrapFort
      USE m_inv3
      !USE m_spg2set
      !USE m_closure, ONLY : check_close
      !USE m_symproperties
      USE m_calculator
      USE m_constants
      !USE m_inpeig
      USE m_sort
      USE m_types_xcpot_inbuild
#ifdef CPP_LIBXC
      USE xc_f03_lib_m
#endif
      IMPLICIT NONE

      TYPE(t_input),INTENT(INOUT)   :: input!ok
      TYPE(t_sym),INTENT(INOUT)     :: sym !ok
      TYPE(t_stars),INTENT(INOUT)   :: stars !ok
      TYPE(t_atoms),INTENT(INOUT)   :: atoms
      TYPE(t_vacuum),INTENT(INOUT)   :: vacuum !ok
      TYPE(t_kpts),INTENT(INOUT)     :: kpts!ok
      TYPE(t_oneD),INTENT(INOUT)     :: oneD !ok
      TYPE(t_hybrid),INTENT(INOUT)   :: hybrid !ok
      TYPE(t_cell),INTENT(INOUT)     :: cell!ok
      TYPE(t_banddos),INTENT(INOUT)  :: banddos!ok
      TYPE(t_sliceplot),INTENT(INOUT):: sliceplot !ok
      CLASS(t_xcpot),INTENT(INOUT),ALLOCATABLE :: xcpot!ok
      TYPE(t_noco),INTENT(INOUT)     :: noco!ok
      TYPE(t_enpara)   ,INTENT(OUT)  :: enpara !ok
      TYPE(t_field), INTENT(INOUT)   :: field !ok
      CLASS(t_forcetheo),ALLOCATABLE,INTENT(OUT):: forcetheo!ok
      TYPE(t_coreSpecInput),INTENT(OUT) :: coreSpecInput !ok
      TYPE(t_wann)   ,INTENT(INOUT)  :: wann !ok

      !+odim
      INTEGER MM,vM,m_cyl
      LOGICAL invs1,zrfs1
      INTEGER chi,rot
      LOGICAL d1,band
      NAMELIST /odim/ d1,MM,vM,m_cyl,chi,rot,invs1,zrfs1
      !-odim
      ! ..
      ! ..  Local Variables
      REAL     :: scpos  ,zc
      INTEGER ieq,i,k,na,n,ii,vxc_id_c,vxc_id_x,exc_id_c,exc_id_x
      REAL s3,ah,a,hs2,rest,thetaj
      LOGICAL l_hyb,l_sym,ldum
      INTEGER :: ierr
      ! ..
      !...  Local Arrays
      !   CHARACTER :: helpchar(atoms%ntype)
      CHARACTER(len=  4) :: chntype
      CHARACTER(len= 41) :: chform
      CHARACTER(len=100) :: line

      CHARACTER(len=20) :: tempNumberString
      CHARACTER(len=150) :: format
      CHARACTER(len=20) :: mixingScheme
      CHARACTER(len=10) :: loType
      LOGICAL :: kptGamma, l_relcor,ldummy
      INTEGER :: iAtomType
      CHARACTER(len=100) :: xPosString, yPosString, zPosString
      CHARACTER(len=20),ALLOCATABLE :: speciesName(:)
      
      !   REAL :: tempTaual(3,atoms%nat)
      TYPE(t_econfig)  :: econf

      INTEGER            :: iType, iLO, iSpecies, lNumCount, nNumCount, iLLO, jsp, j, l, absSum, numTokens
      INTEGER            :: numberNodes, nodeSum, numSpecies, n2spg, n1, n2, ikpt, iqpt
      INTEGER            :: atomicNumber,  gridPoints, lmax, lnonsphr, lmaxAPW
      INTEGER            :: latticeDef, symmetryDef, nop48, firstAtomOfType, errorStatus
      INTEGER            :: loEDeriv, ntp1, ios, ntst, jrc, minNeigd
      INTEGER            :: nv, nv2, kq1, kq2, kq3, nprncTemp, kappaTemp, tempInt
      INTEGER            :: ldau_l(4), numVac, numU
      INTEGER            :: speciesEParams(0:3)
      INTEGER            :: mrotTemp(3,3,48)
      REAL               :: tauTemp(3,48)
      REAL               :: bk(3)
      LOGICAL            :: flipSpin, l_eV, invSym, l_qfix, relaxX, relaxY, relaxZ
      LOGICAL            :: coreConfigPresent, l_enpara, l_orbcomp, tempBool, l_nocoinp
      REAL               :: magMom, radius, logIncrement, qsc(3), latticeScale, dr
      REAL               :: aTemp, zp, rmtmax, sumWeight, ldau_u(4), ldau_j(4), tempReal
      REAL               :: ldau_phi(4),ldau_theta(4)
      REAL               :: weightScale, eParamUp, eParamDown
      LOGICAL            :: l_amf(4)
      INTEGER            :: lcutm,lcutwf,hybSelect(4)
      REAL               :: evac0Temp(2,2)

      CHARACTER(LEN=200,KIND=c_char) :: schemaFilename, docFilename
      CHARACTER(LEN=255) :: valueString, lString, nString, token
      CHARACTER(LEN=255) :: xPathA, xPathB, xPathC, xPathD, xPathE
      CHARACTER(LEN=11)  :: latticeType
      CHARACTER(LEN=50)  :: versionString
      CHARACTER(LEN=150) :: kPointsPrefix

      INTEGER            :: altKPointSetIndex,  altKPointSetIndices(2)
      LOGICAL            :: ldaSpecies
      REAL               :: socscaleSpecies,b_field_mtspecies,vcaspecies

      INTEGER, ALLOCATABLE :: lNumbers(:), nNumbers(:), speciesLLO(:)
      INTEGER, ALLOCATABLE :: loOrderList(:)
      INTEGER, ALLOCATABLE :: speciesNLO(:)
      INTEGER, ALLOCATABLE :: multtab(:,:), invOps(:), optype(:)
      INTEGER, ALLOCATABLE :: lmx1(:), nq1(:), nlhtp1(:)
      INTEGER, ALLOCATABLE :: speciesLOEDeriv(:)
      REAL,    ALLOCATABLE :: speciesLOeParams(:), speciesLLOReal(:)
      LOGICAL, ALLOCATABLE :: wannAtomList(:)

      ! Variables for MT radius testing:

      REAL                 :: dtild1,kmax1,dvac1
      LOGICAL              :: l_test
      INTEGER, ALLOCATABLE :: jri1(:), lmax1(:)
      REAL, ALLOCATABLE    :: rmt1(:), dx1(:)

      EXTERNAL prp_xcfft_box

      INTERFACE
         FUNCTION dropInputSchema() BIND(C, name="dropInputSchema")
            USE iso_c_binding
            INTEGER(c_int) dropInputSchema
         END FUNCTION dropInputSchema
      END INTERFACE

      errorStatus = 0
      errorStatus = dropInputSchema()
      IF(errorStatus.NE.0) THEN
         CALL juDFT_error('Error: Cannot print out FleurInputSchema.xsd')
      END IF

      schemaFilename = "FleurInputSchema.xsd"//C_NULL_CHAR
      docFilename = "inp.xml"//C_NULL_CHAR

    
      !WRITE(*,*) 'Start reading of inp.xml file'
      CALL xmlInitInterface()
      CALL xmlParseSchema(schemaFilename)
      CALL xmlParseDoc(docFilename)
      CALL xmlValidateDoc()
      CALL xmlInitXPath()

      ! Check version of inp.xml
      versionString = xml%GetAttributeValue('/fleurInput/@fleurInputVersion')
      IF((TRIM(ADJUSTL(versionString)).NE.'0.30')) THEN
         CALL juDFT_error('version number of inp.xml file is not compatible with this fleur version')
      END IF


      !Some types have their own readers
      CALL kpts%read_xml(xml)
      CALL sym%read_sym(xml)

      

      ! Get number of atoms, atom types, and atom species

      numberNodes = xmlGetNumberOfNodes('/fleurInput/atomGroups/atomGroup/relPos')
      numberNodes = numberNodes + xmlGetNumberOfNodes('/fleurInput/atomGroups/atomGroup/absPos')
      numberNodes = numberNodes + xmlGetNumberOfNodes('/fleurInput/atomGroups/atomGroup/filmPos')

   
      numSpecies = xmlGetNumberOfNodes('/fleurInput/atomSpecies/species')

    

      ! Read in constants

      xPathA = '/fleurInput/constants/constant'
      numberNodes = xmlGetNumberOfNodes(xPathA)
      DO i = 1, numberNodes
         WRITE(xPathB,*) TRIM(ADJUSTL(xPathA)), '[',i,']'
         tempReal = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@value'))
         valueString = xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@name')
         CALL ASSIGN_var(valueString,tempReal)
      END DO


      ! Read general cutoff parameters

     

      

   

     
  IF (ABS(input%fixed_moment)>1E-8.AND.(input%jspins==1.OR.noco%l_noco)) CALL judft_error("Fixed moment only in collinear calculations with two spins")

 
      ! Check for alternative k point sets for the chosen FLEUR mode

      xPathA = '/fleurInput/output'
      numberNodes = xmlGetNumberOfNodes(xPathA)
     
      altKPointSetIndices(:) = -1
      xPathA = '/fleurInput/calculationSetup/bzIntegration/altKPointSet'
      numberNodes = xmlGetNumberOfNodes(xPathA)
      IF(numberNodes.NE.0) THEN
         DO i = 1, numberNodes
            WRITE(xPathA,*) '/fleurInput/calculationSetup/bzIntegration/altKPointSet[',i,']/@purpose'
            valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA)))))
            IF((altKPointSetIndices(2).EQ.-1).AND.(TRIM(ADJUSTL(valueString)).EQ.'GW')) THEN
               altKPointSetIndices(2) = i
            ELSE IF((altKPointSetIndices(1).EQ.-1).AND.(TRIM(ADJUSTL(valueString)).EQ.'bands')) THEN
               altKPointSetIndices(1) = i
            END IF
         END DO
      END IF

   
      

     
    
      IF (altKPointSetIndex.EQ.-1) THEN
         WRITE(kPointsPrefix,*) '/fleurInput/calculationSetup/bzIntegration'
      END IF

      call judft_error("BUG reading of kpts")
      ! Option kPointDensity
 
         numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(kPointsPrefix))//'/kPointCount/specialPoint')
         IF(numberNodes.EQ.1) THEN
            CALL juDFT_error('Error: Single special k point provided. This does not make sense!')
         END IF
         kpts%numSpecialPoints = numberNodes
         IF(kpts%numSpecialPoints.GE.2) THEN
            DEALLOCATE(kpts%specialPoints)
            ALLOCATE(kpts%specialPoints(3,kpts%numSpecialPoints))
            ALLOCATE(kpts%specialPointNames(kpts%numSpecialPoints))
            DO i = 1, kpts%numSpecialPoints
               WRITE(xPathA,*) TRIM(ADJUSTL(kPointsPrefix))//'/kPointCount/specialPoint[',i,']'
               valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA)))))
               kpts%specialPoints(1,i) = evaluatefirst(valueString)
               kpts%specialPoints(2,i) = evaluatefirst(valueString)
               kpts%specialPoints(3,i) = evaluatefirst(valueString)
               kpts%specialPointNames(i) = xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@name')
            END DO
         END IF

      ! Option kPointList
      numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(kPointsPrefix))//'/kPointList')
      IF (numberNodes.EQ.1) THEN
         l_kpts = .TRUE.
         numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(kPointsPrefix))//'/kPointList/kPoint')
         kpts%nkpt = numberNodes
         kpts%l_gamma = .FALSE.
         ALLOCATE(kpts%bk(3,kpts%nkpt))
         ALLOCATE(kpts%wtkpt(kpts%nkpt))
         kpts%bk = 0.0
         kpts%wtkpt = 0.0
 
          weightScale = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(kPointsPrefix))//'/kPointList/@weightScale'))

         DO i = 1, kpts%nkpt
            WRITE(xPathA,*) TRIM(ADJUSTL(kPointsPrefix))//'/kPointList/kPoint[',i,']'
            valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA)))))
            READ(valueString,*) kpts%bk(1,i), kpts%bk(2,i), kpts%bk(3,i)
            kpts%bk(:,i)=kpts%bk(:,i)
            kpts%wtkpt(i) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@weight'))
            kpts%wtkpt(i) = kpts%wtkpt(i) / weightScale
         END DO
      END IF

      ! Option kPointListFile
      xPathA = TRIM(ADJUSTL(kPointsPrefix))//'/kPointListFile'
      numberNodes = xmlGetNumberOfNodes(xPathA)
      IF (numberNodes.EQ.1) THEN
         valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@filename')))
         OPEN (41,file=TRIM(ADJUSTL(valueString)),form='formatted',status='old')
            READ (41,*) kpts%nkpt
         CLOSE (41)
         ALLOCATE(kpts%bk(3,kpts%nkpt))
         ALLOCATE(kpts%wtkpt(kpts%nkpt))
         kpts%bk = 0.0
         kpts%wtkpt = 0.0
         kpts%l_gamma = .FALSE.
         l_kpts = .TRUE.
     
         

         !CALL inpeig(atoms,cell,input,.FALSE.,kpts,kptsFilename=TRIM(ADJUSTL(valueString)))
      END IF



      

      

     
    
     
   

     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End of calculationSetup section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Start of cell section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   
      ! Construction of amat requires additional information about the lattice
      ! and is done later (scroll down)!

         IF(input%film) THEN
            cell%z1 = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@dVac'))
            dtild = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@dTilda'))
            
            a3(3) = dtild
            evac0Temp = eVac0Default_const
            xPathB = TRIM(ADJUSTL(xPathA))//'/vacuumEnergyParameters'
            numberNodes = xmlGetNumberOfNodes(xPathB)
            IF(numberNodes.GE.1) THEN
               DO i = 1, numberNodes
                  xPathC = ''
                  WRITE(xPathC,'(a,i0,a)') TRIM(ADJUSTL(xPathB))//'[',i,']'
                  numVac = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathC))//'/@vacuum'))
                  eParamUp = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathC))//'/@spinUp'))
                  eParamDown = eParamUp
                  IF (xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathC))//'/@spinDown').GE.1) THEN
                     eParamDown = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathC))//'/@spinDown'))
                  END IF
                  evac0Temp(numVac,1) = eParamUp
                  IF(input%jspins.GT.1) evac0Temp(numVac,2) = eParamDown
                  IF(i == 1) THEN
                     evac0Temp(3-numVac,1) = eParamUp
                     IF(input%jspins.GT.1) evac0Temp(3-numVac,2) = eParamDown
                  END IF
               END DO
            END IF
      

      ! Calculate missing symmetry and cell properties and check consistency of parameters.


      IF (latticeScale.EQ.0.0) latticeScale = 1.0
      IF (.NOT.input%film) vacuum%dvac = a3(3)
      vacuum%dvac = latticeScale*vacuum%dvac
      dtild = latticeScale*dtild



      ! Construction of missing symmetry information
      IF ((symmetryDef.EQ.2).OR.(symmetryDef.EQ.3)) THEN
         CALL sym%init(cell,input%film)
      END IF
 
    
      !some settings for film calculations
     
      IF (sym%zrfs.OR.sym%invs) vacuum%nvac = 1
      IF (oneD%odd%d1) vacuum%nvac = 1
    
      IF (oneD%odd%d1) vacuum%delz = 20.0 / vacuum%nmz
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End of cell section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Start of XC functional section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Read in libxc parameters if present
      xPathA = '/fleurInput/xcFunctional/LibXCID'
      xPathB = '/fleurInput/xcFunctional/LibXCName'

      IF(xmlGetNumberOfNodes(xPathA) == 1 .AND. xmlGetNumberOfNodes(xPathB) == 1) THEN
         CALL judft_error("LibXC is given both by Name and ID and is therefore overdetermined", calledby="r_inpXML")
      ENDIF


      ! LibXCID 
      IF (xmlGetNumberOfNodes(xPathA) == 1) THEN
#ifdef CPP_LIBXC
         vxc_id_x=evaluateFirstOnly(xmlGetAttributeValue(xPathA // '/@exchange'))
         vxc_id_c=evaluateFirstOnly(xmlGetAttributeValue(xPathA // '/@correlation'))

         IF(xmlGetNumberOfNodes(TRIM(xPathA) // '/@etot_exchange') == 1) THEN
            exc_id_x = evaluateFirstOnly(xmlGetAttributeValue(xPathA // '/@etot_exchange'))
         ELSE
            exc_id_x = vxc_id_x
         ENDIF
         
         IF(xmlGetNumberOfNodes(TRIM(xPathA) // '/@exc_correlation') == 1) THEN
            exc_id_c = evaluateFirstOnly(xmlGetAttributeValue(xPathA // '/@exc_correlation'))
         ELSE
            exc_id_c = vxc_id_c
         ENDIF
#else
         CALL judft_error("To use libxc functionals you have to compile with libXC support")
#endif
      ! LibXCName 
      ELSEIF (xmlGetNumberOfNodes(TRIM(xPathB)) == 1) THEN
#ifdef CPP_LIBXC
         valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(xPathB) // '/@exchange')))
         vxc_id_x =  xc_f03_functional_get_number(TRIM(valueString))

         valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(xPathB) // '/@correlation')))
         vxc_id_c =  xc_f03_functional_get_number(TRIM(valueString))
         
         IF(xmlGetNumberOfNodes(TRIM(xPathB) // '/@etot_exchange') == 1) THEN
            valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(xPathB) // '/@etot_exchange')))
            exc_id_x =  xc_f03_functional_get_number(TRIM(valueString))
         ELSE
            exc_id_x = vxc_id_x
         ENDIF
         
         IF(xmlGetNumberOfNodes(TRIM(xPathB) // '/@etot_correlation') == 1) THEN
            valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(xPathB) // '/@etot_correlation')))
            exc_id_c =  xc_f03_functional_get_number(TRIM(valueString))
         ELSE
            exc_id_c = vxc_id_c
         ENDIF
#else
         CALL judft_error("To use libxc functionals you have to compile with libXC support")
#endif
      ELSE
         vxc_id_x=0; vxc_id_c=0;
         exc_id_x=0; exc_id_c=0;
      ENDIF

      ! Read in xc functional parameters
      valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL('/fleurInput/xcFunctional/@name')))))
      namex(1:4) = valueString(1:4)
      l_relcor = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/xcFunctional/@relativisticCorrections'))

      relcor = 'non-relativi'
      IF (l_relcor) THEN
         relcor = 'relativistic'
      END IF

      !now initialize the xcpot variable
      CALL setXCParameters(atoms,valueString,l_relcor,input%jspins,vxc_id_x,vxc_id_c,exc_id_x, exc_id_c, xcpot)

      xPathA = '/fleurInput/calculationSetup/cutoffs/@GmaxXC'
      numberNodes = xmlGetNumberOfNodes(xPathA)
      xcpot%gmaxxc = stars%gmax
      IF(numberNodes.EQ.1) THEN
         xcpot%gmaxxc = evaluateFirstOnly(xmlGetAttributeValue(xPathA))
      END IF
      

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End of XC functional section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Start of species section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ALLOCATE (speciesNLO(numSpecies))
      
      ALLOCATE(speciesName(numSpecies))
    
   
      DO iSpecies = 1, numSpecies
         ! Attributes of species
         WRITE(xPathA,*) '/fleurInput/atomSpecies/species[',iSpecies,']'
         speciesName(iSpecies) = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@name')))
         atomicNumber = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@atomicNumber'))
         magMom = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@magMom'))
         flipSpin = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@flipSpin'))

         ! Attributes of mtSphere element of species
         radius = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/mtSphere/@radius'))
         gridPoints = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/mtSphere/@gridPoints'))
         logIncrement = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/mtSphere/@logIncrement'))

         ! Attributes of atomicCutoffs element of species
         lmax = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/atomicCutoffs/@lmax'))
         lnonsphr = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/atomicCutoffs/@lnonsphr'))
         lmaxAPW = -1
         numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/atomicCutoffs/@lmaxAPW')
         IF (numberNodes.EQ.1) THEN
            lmaxAPW = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/atomicCutoffs/@lmaxAPW'))
         END IF

         numU = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/ldaU')
        IF (numU.GT.4) CALL juDFT_error("Too many U parameters provided for a certain species (maximum is 4).",calledby ="r_inpXML")
         ldau_l = -1
         ldau_u = 0.0
         ldau_j = 0.0
         l_amf = .FALSE.
         DO i = 1, numU
            WRITE(xPathB,*) i
            ldau_l(i) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/ldaU['//TRIM(ADJUSTL(xPathB))//']/@l'))
            ldau_u(i) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/ldaU['//TRIM(ADJUSTL(xPathB))//']/@U'))
            ldau_j(i) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/ldaU['//TRIM(ADJUSTL(xPathB))//']/@J'))
            ldau_phi(i) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/ldaU['//TRIM(ADJUSTL(xPathB))//']/@phi'))
            ldau_theta(i) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/ldaU['//TRIM(ADJUSTL(xPathB))//']/@theta'))
            l_amf(i) = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/ldaU['//TRIM(ADJUSTL(xPathB))//']/@l_amf'))
         END DO

         speciesNLO(iSpecies) = 0
         WRITE(xPathA,*) '/fleurInput/atomSpecies/species[',iSpecies,']/lo'
         numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA)))
         DO iLO = 1, numberNodes
            WRITE(xPathB,*) TRIM(ADJUSTL(xPathA)),'[',iLO,']/@l'
            WRITE(xPathC,*) TRIM(ADJUSTL(xPathA)),'[',iLO,']/@n'
            lString = xmlGetAttributeValue(TRIM(ADJUSTL(xPathB)))
            nString = xmlGetAttributeValue(TRIM(ADJUSTL(xPathC)))
            CALL getIntegerSequenceFromString(TRIM(ADJUSTL(lString)), lNumbers, lNumCount)
            CALL getIntegerSequenceFromString(TRIM(ADJUSTL(nString)), nNumbers, nNumCount)
            IF(lNumCount.NE.nNumCount) THEN
               CALL judft_error('Error in LO input: l quantum number count does not equal n quantum number count')
            END IF
            speciesNLO(iSpecies) = speciesNLO(iSpecies) + lNumCount
            DEALLOCATE (lNumbers, nNumbers)
         END DO
         ! Special switches for species
         vcaspecies=0.0
         WRITE(xPathA,*) '/fleurInput/atomSpecies/species[',iSpecies,']/special'
         numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA)))
         IF (numberNodes==1) THEN
            vcaSpecies   = evaluateFirstOnly(TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@vca_charge'))))
         ENDIF
       
 
    
 
      CALL enpara%init(atoms,input%jspins)
      enpara%evac0(:,:) = evac0Temp(:,:)

            
               ! Explicit xc functional
               SELECT TYPE(xcpot)
               TYPE IS(t_xcpot_inbuild)
                  xcpot%lda_atom(iType)=ldaSpecies
               END SELECT
               noco%socscale(iType)=socscaleSpecies
               IF (field%l_b_field) THEN
                  IF (.NOT.ALLOCATED(field%b_field_mt)) THEN
                     ALLOCATE(field%b_field_mt(atoms%ntype))
                     field%b_field_mt=0.0
                  ENDIF
                  field%b_field_mt(itype)=b_field_mtSpecies
               ENDIF
            END IF
         END DO
         DEALLOCATE(loOrderList)
         DEALLOCATE(speciesLLO,speciesLOeParams,speciesLOEDeriv)
      END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End of species section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Start of atomGroup section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      atoms%l_geo = .FALSE.
      atoms%relax = 0
      na = 0
      firstAtomOfType = 1
      DO iType = 1, atoms%ntype
         WRITE(xPathA,*) '/fleurInput/atomGroups/atomGroup[',iType,']'

         ! Read in force parameters
         xPathB = TRIM(ADJUSTL(xPathA))//'/force'
         numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathB)))
         IF (numberNodes.GE.1) THEN
            atoms%l_geo(iType) = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@calculate'))
            valueString = xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@relaxXYZ')
            READ(valueString,'(3l1)') relaxX, relaxY, relaxZ
            IF (relaxX) atoms%relax(1,iType) = 1
            IF (relaxY) atoms%relax(2,iType) = 1
            IF (relaxZ) atoms%relax(3,iType) = 1
         END IF

         ! Obtain number of equivalent atoms
         numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/relPos')
         numberNodes = numberNodes + xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/absPos')
         numberNodes = numberNodes + xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/filmPos')
         atoms%neq(iType) = numberNodes

         IF (iType.GE.2) THEN
            firstAtomOfType = firstAtomOfType + atoms%neq(iType-1)
         END IF

         ! Read in atom positions
         numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/relPos')
         DO i = 1, numberNodes
            na = na + 1
            WRITE(xPathB,*) TRIM(ADJUSTL(xPathA)),'/relPos[',i,']'
            IF(xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathB))//'/@label').NE.0) THEN
               atoms%label(na) = xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@label')
            ELSE
               WRITE(atoms%label(na),'(i0)') na
            END IF
            valueString = xmlGetAttributeValue(TRIM(ADJUSTL(xPathB)))
            atoms%taual(1,na) = evaluatefirst(valueString)
            atoms%taual(2,na) = evaluatefirst(valueString)
            atoms%taual(3,na) = evaluatefirst(valueString)
            atoms%pos(:,na) = MATMUL(cell%amat,atoms%taual(:,na))
            l_orbcomp = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@orbcomp'))
            IF(l_orbcomp) THEN
               IF(banddos%l_orb) THEN
                  CALL juDFT_error("Multiple orbcomp flags set.", calledby = "r_inpXML")
               END IF
               banddos%l_orb = .TRUE.
               banddos%orbCompAtom = na
            END IF
            wannAtomList(na) = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@wannier'))
         END DO

         numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/absPos')
         DO i = 1, numberNodes
            na = na + 1
            CALL judft_error('absPos not yet implemented!')
         END DO

         numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/filmPos')
         DO i = 1, numberNodes
            na = na + 1
            WRITE(xPathB,*) TRIM(ADJUSTL(xPathA)),'/filmPos[',i,']'
            IF(xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathB))//'/@label').NE.0) THEN
               atoms%label(na) = xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@label')
            ELSE
               WRITE(atoms%label(na),'(i0)') na
            END IF
            valueString = xmlGetAttributeValue(TRIM(ADJUSTL(xPathB)))
            atoms%taual(1,na) = evaluatefirst(valueString)
            atoms%taual(2,na) = evaluatefirst(valueString)
            atoms%taual(3,na) = evaluatefirst(valueString) / cell%amat(3,3)
            atoms%pos(:,na) = MATMUL(cell%amat,atoms%taual(:,na))
            wannAtomList(na) = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@wannier'))
         END DO

        
      END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End of atomGroup section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Start of force-theorem section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      xPathA = '/fleurInput/forceTheorem'
      numberNodes = xmlGetNumberOfNodes(xPathA)
      IF (numberNodes.EQ.1) THEN
         !Magnetic anisotropy...
         xPathA = '/fleurInput/forceTheorem/MAE'
         numberNodes = xmlGetNumberOfNodes(xPathA)
         IF (numberNodes.EQ.1) THEN
         ENDIF
         !spin-spiral dispersion
         xPathA = '/fleurInput/forceTheorem/spinSpiralDispersion'
         numberNodes = xmlGetNumberOfNodes(xPathA)
         IF (numberNodes.EQ.1) THEN
         ENDIF
         !dmi
         
         !jij
         xPathA = '/fleurInput/forceTheorem/Jij'
         numberNodes = xmlGetNumberOfNodes(xPathA)
         IF (numberNodes.EQ.1) THEN
         ENDIF

      ELSE
         ALLOCATE(t_forcetheo::forcetheo) !default no forcetheorem type
      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End of force-theorem section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         
         

         
      

    
      END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End of output section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Start of non-XML input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Read in enpara file iff available

      l_enpara = .FALSE.
      INQUIRE (file ='enpara',exist= l_enpara)
      IF (l_enpara) THEN
         CALL enpara%READ(atoms,input%jspins,input%film,.FALSE.)
      END IF
      hybrid%l_hybrid=xcpot%is_hybrid()

    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End of non-XML input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !CALL xmlFreeResources()
      
      !WRITE(*,*) 'Reading of inp.xml file finished'

      DEALLOCATE(speciesNLO)

   END SUBROUTINE r_inpXML

   SUBROUTINE setXCParameters(atoms,namex,relcor,jspins,vxc_id_x,vxc_id_c,exc_id_x,exc_id_c,xcpot)
      USE m_juDFT
      USE m_types
      USE m_types_xcpot_inbuild
      USE m_types_xcpot_libxc

      IMPLICIT NONE
      TYPE(t_atoms),INTENT(IN)          :: atoms
      CHARACTER(LEN=*),     INTENT(IN)  :: namex
      LOGICAL,              INTENT(IN)  :: relcor
      INTEGER,              INTENT(IN)  :: jspins,vxc_id_c,vxc_id_x,exc_id_x,exc_id_c
      CLASS(t_xcpot),INTENT(OUT),ALLOCATABLE      :: xcpot

      IF (namex(1:5)=='LibXC') THEN
         ALLOCATE(t_xcpot_libxc::xcpot)
      ELSE
         ALLOCATE(t_xcpot_inbuild::xcpot)
      ENDIF

      SELECT TYPE(xcpot)
      TYPE IS(t_xcpot_inbuild)
         CALL xcpot%init(namex(1:4),relcor,atoms%ntype)
      TYPE IS(t_xcpot_libxc)
         CALL xcpot%init(jspins,vxc_id_x,vxc_id_c,exc_id_x,exc_id_c)

      END SELECT

      CALL set_xcpot_usage(xcpot)
   END SUBROUTINE setXCParameters

   SUBROUTINE set_xcpot_usage(xcpot)
      use m_judft_usage
      USE m_types
      USE m_types_xcpot_inbuild
      USE m_types_xcpot_libxc
      implicit none
      class(t_xcpot), intent(in)    :: xcpot

      ! give some information about XC functional to usage.json
      ! 1 -> LDA
      ! 2 -> GGA
      ! 3 -> MetaGGA
      ! 4 -> Hybrid functional
      if(xcpot%vxc_is_lda()) then
         call add_usage_data("XC-treatment", 1)
         return
      endif

      if(xcpot%exc_is_MetaGGA()) then
         call add_usage_data("XC-treatment", 3)
         return
      endif

      if(xcpot%vxc_is_GGA()) then
         call add_usage_data("XC-treatment", 2)
         return
      endif

      if(xcpot%is_hybrid()) then
         call add_usage_data("XC-treatment", 4)
         return
      endif

   END SUBROUTINE set_xcpot_usage

   SUBROUTINE getIntegerSequenceFromString(string, sequence, count)
      use m_juDFT_stop
      IMPLICIT NONE

      CHARACTER(*),         INTENT(IN)  :: string
      INTEGER, ALLOCATABLE, INTENT(OUT) :: sequence(:)
      INTEGER,              INTENT(OUT) :: count

      INTEGER :: i, length, start, lastNumber, currentNumber, index
      LOGICAL singleNumber, comma, dash

      ! 3 cases: 1. a single number
      !          2. number - number
      !          3. comma separated numbers

      length = LEN(string)
      count = 0
      start = 1
      singleNumber = .TRUE.
      comma = .FALSE.
      dash = .FALSE.
      lastNumber = 0
      count = 0

      ! 1. Determine number count

      DO i = 1, length
         SELECT CASE (string(i:i))
         CASE ('0':'9')
         CASE (',')
            IF ((start.EQ.i).OR.(dash)) THEN
               CALL judft_error('String has wrong syntax (in getIntegerSequenceFromString)')
            END IF
            singleNumber = .FALSE.
            comma = .TRUE.
            READ(string(start:i-1),*) lastNumber
            count = count + 1
            start = i+1
         CASE ('-')
            IF ((start.EQ.i).OR.(dash).OR.(comma)) THEN
               CALL judft_error('String has wrong syntax (in getIntegerSequenceFromString)')
            END IF
            singleNumber = .FALSE.
            dash = .TRUE.
            READ(string(start:i-1),*) lastNumber
            start = i+1
         CASE DEFAULT
            CALL judft_error('String has wrong syntax (in getIntegerSequenceFromString)')
         END SELECT
      END DO
      IF(start.GT.length) THEN
         CALL judft_error('String has wrong syntax (in getIntegerSequenceFromString)')
      END IF
      READ(string(start:length),*) currentNumber
      IF (dash) THEN
         count = currentNumber - lastNumber + 1
      ELSE
         count = count + 1
      END IF

      IF (ALLOCATED(sequence)) THEN
         DEALLOCATE(sequence)
      END IF
      ALLOCATE(sequence(count))

      ! 2. Read in numbers iff comma separation ...and store numbers in any case

      IF (singleNumber) THEN
         sequence(1) = currentNumber
      ELSE IF (dash) THEN
         DO i = 1, count
            sequence(i) = lastNumber + i - 1
         END DO
      ELSE
         index = 1
         start = 1
         DO i = 1, length
            SELECT CASE (string(i:i))
            CASE (',')
               comma = .TRUE.
               READ(string(start:i-1),*) lastNumber
               start = i+1
               SEQUENCE(index) = lastNumber
               index = index + 1
            END SELECT
         END DO
         sequence(index) = currentNumber
      END IF

   END SUBROUTINE getIntegerSequenceFromString

  
END MODULE m_rinpXML
