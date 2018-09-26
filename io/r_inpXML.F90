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
      atoms,obsolete,vacuum,input,stars,sliceplot,banddos,DIMENSION,forcetheo,&
      cell,sym,xcpot,noco,oneD,hybrid,kpts,enpara,coreSpecInput,wann,&
      noel,namex,relcor,a1,a2,a3,dtild,xmlElectronStates,&
      xmlPrintCoreStates,xmlCoreOccs,atomTypeSpecies,speciesRepAtomType,&
      l_kpts)

      USE iso_c_binding
      USE m_juDFT
      USE m_types
      USE m_types_forcetheo_extended
      USE m_symdata , ONLY : nammap, ord2, l_c2
      USE m_rwsymfile
      USE m_xmlIntWrapFort
      USE m_inv3
      USE m_spg2set
      USE m_closure, ONLY : check_close
      USE m_symproperties
      USE m_calculator
      USE m_constants
      USE m_inpeig
      USE m_sort
      USE m_types_xcpot_inbuild
#ifdef CPP_LIBXC
      USE xc_f03_lib_m
#endif
      IMPLICIT NONE

      TYPE(t_input),INTENT(INOUT)   :: input
      TYPE(t_sym),INTENT(INOUT)     :: sym
      TYPE(t_stars),INTENT(INOUT)   :: stars
      TYPE(t_atoms),INTENT(INOUT)   :: atoms
      TYPE(t_vacuum),INTENT(INOUT)   :: vacuum
      TYPE(t_obsolete),INTENT(INOUT) :: obsolete
      TYPE(t_kpts),INTENT(INOUT)     :: kpts
      TYPE(t_oneD),INTENT(INOUT)     :: oneD
      TYPE(t_hybrid),INTENT(INOUT)   :: hybrid
      TYPE(t_cell),INTENT(INOUT)     :: cell
      TYPE(t_banddos),INTENT(INOUT)  :: banddos
      TYPE(t_sliceplot),INTENT(INOUT):: sliceplot
      CLASS(t_xcpot),INTENT(INOUT),ALLOCATABLE :: xcpot
      TYPE(t_noco),INTENT(INOUT)     :: noco
      TYPE(t_dimension),INTENT(OUT)  :: DIMENSION
      TYPE(t_enpara)   ,INTENT(OUT)  :: enpara
      CLASS(t_forcetheo),ALLOCATABLE,INTENT(OUT):: forcetheo
      TYPE(t_coreSpecInput),INTENT(OUT) :: coreSpecInput
      TYPE(t_wann)   ,INTENT(INOUT)  :: wann
      LOGICAL, INTENT(OUT)           :: l_kpts
      INTEGER,          ALLOCATABLE, INTENT(INOUT) :: xmlElectronStates(:,:)
      INTEGER,          ALLOCATABLE, INTENT(INOUT) :: atomTypeSpecies(:)
      INTEGER,          ALLOCATABLE, INTENT(INOUT) :: speciesRepAtomType(:)
      REAL,             ALLOCATABLE, INTENT(INOUT) :: xmlCoreOccs(:,:,:)
      LOGICAL,          ALLOCATABLE, INTENT(INOUT) :: xmlPrintCoreStates(:,:)
      CHARACTER(len=3), ALLOCATABLE, INTENT(INOUT) :: noel(:)
      CHARACTER(len=4), INTENT(OUT)  :: namex
      CHARACTER(len=12), INTENT(OUT) :: relcor
      REAL, INTENT(OUT)              :: a1(3),a2(3),a3(3)
      REAL, INTENT(OUT)              :: dtild

      CHARACTER(len=8) :: name(10)

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
      CHARACTER(len=150) :: FORMAT
      CHARACTER(len=20) :: mixingScheme
      CHARACTER(len=10) :: loType
      LOGICAL :: kptGamma, l_relcor,ldummy
      INTEGER :: iAtomType, startCoreStates, endCoreStates
      CHARACTER(len=100) :: xPosString, yPosString, zPosString
      CHARACTER(len=200) :: coreStatesString
      !   REAL :: tempTaual(3,atoms%nat)
      REAL             :: coreStateOccs(29,2)
      INTEGER          :: coreStateNprnc(29), coreStateKappa(29)
      INTEGER          :: speciesXMLElectronStates(29)
      REAL             :: speciesXMLCoreOccs(2,29)
      LOGICAL          :: speciesXMLPrintCoreStates(29)

      INTEGER            :: iType, iLO, iSpecies, lNumCount, nNumCount, iLLO, jsp, j, l, absSum, numTokens
      INTEGER            :: numberNodes, nodeSum, numSpecies, n2spg, n1, n2, ikpt, iqpt
      INTEGER            :: atomicNumber, coreStates, gridPoints, lmax, lnonsphr, lmaxAPW
      INTEGER            :: latticeDef, symmetryDef, nop48, firstAtomOfType, errorStatus
      INTEGER            :: loEDeriv, ntp1, ios, ntst, jrc, minNeigd, providedCoreStates, providedStates
      INTEGER            :: nv, nv2, kq1, kq2, kq3, nprncTemp, kappaTemp, tempInt
      INTEGER            :: ldau_l(4), numVac, numU
      INTEGER            :: speciesEParams(0:3)
      INTEGER            :: mrotTemp(3,3,48)
      REAL               :: tauTemp(3,48)
      REAL               :: bk(3)
      LOGICAL            :: flipSpin, l_eV, invSym, l_qfix, relaxX, relaxY, relaxZ
      LOGICAL            :: l_vca, coreConfigPresent, l_enpara, l_orbcomp, tempBool
      REAL               :: magMom, radius, logIncrement, qsc(3), latticeScale, dr
      REAL               :: aTemp, zp, rmtmax, sumWeight, ldau_u(4), ldau_j(4), tempReal
      REAL               :: weightScale, eParamUp, eParamDown
      LOGICAL            :: l_amf(4)
      REAL, PARAMETER    :: boltzmannConst = 3.1668114e-6 ! value is given in Hartree/Kelvin
      INTEGER            :: lcutm,lcutwf,hybSelect(4)
      REAL               :: evac0Temp(2,2)

      CHARACTER(LEN=200,KIND=c_char) :: schemaFilename, docFilename
      CHARACTER(LEN=255) :: valueString, lString, nString, token
      CHARACTER(LEN=255) :: xPathA, xPathB, xPathC, xPathD, xPathE
      CHARACTER(LEN=11)  :: latticeType
      CHARACTER(LEN=50)  :: versionString

      LOGICAL            :: ldaSpecies
      REAL               :: socscaleSpecies

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

      !TODO! these switches should be in the inp-file
      input%l_core_confpot=.TRUE. !former CPP_CORE
      input%l_useapw=.FALSE.   !former CPP_APW
      !WRITE(*,*) 'Start reading of inp.xml file'
      CALL xmlInitInterface()
      CALL xmlParseSchema(schemaFilename)
      CALL xmlParseDoc(docFilename)
      CALL xmlValidateDoc()
      CALL xmlInitXPath()

      ! Check version of inp.xml
      versionString = xmlGetAttributeValue('/fleurInput/@fleurInputVersion')
      IF((TRIM(ADJUSTL(versionString)).NE.'0.27').AND.(TRIM(ADJUSTL(versionString)).NE.'0.28').AND.&
         (TRIM(ADJUSTL(versionString)).NE.'0.29')) THEN
         CALL juDFT_error('version number of inp.xml file is not compatible with this fleur version')
      END IF

      ! Get number of atoms, atom types, and atom species

      numberNodes = xmlGetNumberOfNodes('/fleurInput/atomGroups/atomGroup/relPos')
      numberNodes = numberNodes + xmlGetNumberOfNodes('/fleurInput/atomGroups/atomGroup/absPos')
      numberNodes = numberNodes + xmlGetNumberOfNodes('/fleurInput/atomGroups/atomGroup/filmPos')

      atoms%nat = numberNodes

      numberNodes = xmlGetNumberOfNodes('/fleurInput/atomGroups/atomGroup')

      atoms%ntype = numberNodes

      numSpecies = xmlGetNumberOfNodes('/fleurInput/atomSpecies/species')

      ALLOCATE(atoms%nz(atoms%ntype))     !nz and zatom have the same content!
      ALLOCATE(atoms%zatom(atoms%ntype))  !nz and zatom have the same content!
      ALLOCATE(atoms%jri(atoms%ntype))
      ALLOCATE(atoms%dx(atoms%ntype))
      ALLOCATE(atoms%lmax(atoms%ntype))
      ALLOCATE(atoms%nlo(atoms%ntype))
      ALLOCATE(atoms%ncst(atoms%ntype))
      ALLOCATE(atoms%lnonsph(atoms%ntype))
      ALLOCATE(atoms%nflip(atoms%ntype))
      ALLOCATE(atoms%l_geo(atoms%ntype))
      ALLOCATE(atoms%lda_u(4*atoms%ntype))
      ALLOCATE(atoms%bmu(atoms%ntype))
      ALLOCATE(atoms%relax(3,atoms%ntype))
      ALLOCATE(atoms%neq(atoms%ntype))
      ALLOCATE(atoms%taual(3,atoms%nat))
      ALLOCATE(atoms%label(atoms%nat))
      ALLOCATE(atoms%pos(3,atoms%nat))
      ALLOCATE(atoms%rmt(atoms%ntype))
      ALLOCATE(atoms%numStatesProvided(atoms%ntype))
      ALLOCATE(atoms%namex(atoms%ntype))
      ALLOCATE(atoms%icorr(atoms%ntype))
      ALLOCATE(atoms%igrd(atoms%ntype))
      ALLOCATE(atoms%krla(atoms%ntype))
      ALLOCATE(atoms%relcor(atoms%ntype))

      atoms%namex = ''
      atoms%icorr = -99

      ALLOCATE(atoms%ncv(atoms%ntype)) ! For what is this?
      ALLOCATE(atoms%ngopr(atoms%nat)) ! For what is this?
      ALLOCATE(atoms%lapw_l(atoms%ntype)) ! Where do I put this?
      ALLOCATE(atoms%invsat(atoms%nat)) ! Where do I put this?

      ALLOCATE(noco%l_relax(atoms%ntype),noco%b_con(2,atoms%ntype))
      ALLOCATE(noco%alphInit(atoms%ntype),noco%alph(atoms%ntype),noco%beta(atoms%ntype))
      ALLOCATE(noco%socscale(atoms%ntype))

      DEALLOCATE(atomTypeSpecies,speciesRepAtomType)
      ALLOCATE(atomTypeSpecies(atoms%ntype))
      ALLOCATE(speciesRepAtomType(numSpecies))
      atomTypeSpecies = -1
      speciesRepAtomType = -1

      DEALLOCATE(xmlElectronStates,xmlPrintCoreStates,xmlCoreOccs)
      ALLOCATE(xmlElectronStates(29,atoms%ntype))
      ALLOCATE(xmlPrintCoreStates(29,atoms%ntype))
      ALLOCATE(xmlCoreOccs(2,29,atoms%ntype))
      xmlElectronStates = noState_const
      xmlPrintCoreStates = .FALSE.
      xmlCoreOccs = 0.0

      ALLOCATE (kpts%ntetra(4,kpts%ntet),kpts%voltet(kpts%ntet))

      ALLOCATE (wannAtomList(atoms%nat))

      ! Read in constants

      xPathA = '/fleurInput/constants/constant'
      numberNodes = xmlGetNumberOfNodes(xPathA)
      DO i = 1, numberNodes
         WRITE(xPathB,*) TRIM(ADJUSTL(xPathA)), '[',i,']'
         tempReal = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@value'))
         valueString = xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@name')
         CALL ASSIGN_var(valueString,tempReal)
      END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Comment section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      input%comment = '        '
      xPathA = '/fleurInput/comment'
      valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA)))))
      DO i = 1, LEN(TRIM(ADJUSTL(valueString)))
         IF (valueString(i:i) == ACHAR(10)) valueString(i:i) = ' ' !remove line breaks
      END DO
      valueString = TRIM(ADJUSTL(valueString))
      DO i = 1, 10
         j = (i-1) * 8 + 1
         input%comment(i) = valueString(j:j+7)
      END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Start of calculationSetup section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Read general cutoff parameters

      input%rkmax = evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/cutoffs/@Kmax'))
      stars%gmax = evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/cutoffs/@Gmax'))

      stars%gmaxInit = stars%gmax

      xPathA = '/fleurInput/calculationSetup/cutoffs/@numbands'
      numberNodes = xmlGetNumberOfNodes(xPathA)
      DIMENSION%neigd = 0
      IF(numberNodes.EQ.1) THEN
         valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA)))))
         IF(TRIM(ADJUSTL(valueString)).EQ.'all') THEN
            CALL juDFT_error('Feature to calculate all eigenfunctions not yet implemented.')
         ELSE
            READ(valueString,*) DIMENSION%neigd
         END IF
      END IF

      ! Read SCF loop parametrization

      input%itmax = evaluateFirstIntOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@itmax'))
      input%minDistance = evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@minDistance'))
      input%maxiter = evaluateFirstIntOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@maxIterBroyd'))

      valueString = TRIM(ADJUSTL(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@imix')))
      SELECT CASE (valueString)
      CASE ('straight')
         input%imix = 0
      CASE ('Broyden1')
         input%imix = 3
      CASE ('Broyden2')
         input%imix = 5
      CASE ('Anderson')
         input%imix = 7
      CASE DEFAULT
         CALL juDFT_error('Error: unknown mixing scheme selected!')
      END SELECT

      input%alpha = evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@alpha'))
input%preconditioning_param = evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@preconditioning_param'))
      input%spinf = evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@spinf'))

      ! Get parameters for core electrons

      input%ctail = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/coreElectrons/@ctail'))
      IF((TRIM(ADJUSTL(versionString)).EQ.'0.27')) THEN
         input%coretail_lmax = 99
      ELSE
       input%coretail_lmax = evaluateFirstIntOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/coreElectrons/@coretail_lmax'))
      END IF
      input%frcor = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/coreElectrons/@frcor'))
      input%kcrel = evaluateFirstIntOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/coreElectrons/@kcrel'))

      ! Read in magnetism parameters

      input%jspins = evaluateFirstIntOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/magnetism/@jspins'))
      noco%l_noco = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/magnetism/@l_noco'))
      input%swsp = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/magnetism/@swsp'))
      input%lflip = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/magnetism/@lflip'))
      input%fixed_moment=evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/magnetism/@fixed_moment'))

      IF (ABS(input%fixed_moment)>1E-8.AND.(input%jspins==1.OR.noco%l_noco)) CALL juDFT_error("Fixed moment only in collinear calculations with two spins")
      DIMENSION%jspd = input%jspins

      ! Read in Brillouin zone integration parameters

      kpts%nkpt3 = 0
      l_kpts = .FALSE.

      valueString = TRIM(ADJUSTL(xmlGetAttributeValue('/fleurInput/calculationSetup/bzIntegration/@mode')))
      SELECT CASE (valueString)
      CASE ('hist')
         input%gauss = .FALSE.
         input%tria = .FALSE.
      CASE ('gauss')
         input%gauss = .TRUE.
         input%tria = .FALSE.
      CASE ('tria')
         input%gauss = .FALSE.
         input%tria = .TRUE.
      CASE DEFAULT
         CALL juDFT_error('Invalid bzIntegration mode selected!')
      END SELECT

      nodeSum = 0
      xPathA = '/fleurInput/calculationSetup/bzIntegration/@fermiSmearingEnergy'
      numberNodes = xmlGetNumberOfNodes(xPathA)
      nodeSum = nodeSum + numberNodes
      IF (numberNodes.EQ.1) THEN
         input%tkb = evaluateFirstOnly(xmlGetAttributeValue(xPathA))
      END IF
      xPathA = '/fleurInput/calculationSetup/bzIntegration/@fermiSmearingTemp'
      numberNodes = xmlGetNumberOfNodes(xPathA)
      nodeSum = nodeSum + numberNodes
      IF (numberNodes.EQ.1) THEN
         input%tkb = evaluateFirstOnly(xmlGetAttributeValue(xPathA))
         input%tkb = boltzmannConst * input%tkb
      END IF
      IF(nodeSum.GE.2) THEN
         CALL juDFT_error('Error: Multiple fermi Smearing parameters provided in input file!')
      END IF

      xPathA = '/fleurInput/calculationSetup/bzIntegration/@valenceElectrons'
      numberNodes = xmlGetNumberOfNodes(xPathA)
      IF (numberNodes.EQ.1) THEN
         input%zelec = evaluateFirstOnly(xmlGetAttributeValue(xPathA))
      ELSE
         CALL juDFT_error('Error: Optionality of valence electrons in input file not yet implemented!')
      END IF

      ! Option kPointDensity
      kpts%kPointDensity(:) = 0.0
      xPathA = '/fleurInput/calculationSetup/bzIntegration/kPointDensity'
      numberNodes = xmlGetNumberOfNodes(xPathA)
      IF (numberNodes.EQ.1) THEN
         l_kpts = .FALSE.
         kpts%kPointDensity(1) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@denX'))
         kpts%kPointDensity(2) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@denY'))
         kpts%kPointDensity(3) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@denZ'))
         kpts%l_gamma = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@gamma'))
         kpts%specificationType = 4
      END IF

      ! Option kPointMesh
      xPathA = '/fleurInput/calculationSetup/bzIntegration/kPointMesh'
      numberNodes = xmlGetNumberOfNodes(xPathA)
      IF (numberNodes.EQ.1) THEN
         l_kpts = .FALSE.
         kpts%nkpt3(1) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nx'))
         kpts%nkpt3(2) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@ny'))
         kpts%nkpt3(3) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nz'))
         kpts%l_gamma = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@gamma'))
         kpts%nkpt = kpts%nkpt3(1) * kpts%nkpt3(2) * kpts%nkpt3(3)
         kpts%specificationType = 2
      END IF

      ! Option kPointCount
      xPathA = '/fleurInput/calculationSetup/bzIntegration/kPointCount'
      numberNodes = xmlGetNumberOfNodes(xPathA)
      IF (numberNodes.EQ.1) THEN
         l_kpts = .FALSE.
         kpts%nkpt = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@count'))
         kpts%l_gamma = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@gamma'))
         kpts%nkpt = kpts%nkpt
         kpts%specificationType = 1

         ALLOCATE(kpts%bk(3,kpts%nkpt))
         ALLOCATE(kpts%wtkpt(kpts%nkpt))
         kpts%bk = 0.0
         kpts%wtkpt = 0.0
         kpts%posScale = 1.0

         numberNodes = xmlGetNumberOfNodes('/fleurInput/calculationSetup/bzIntegration/kPointCount/specialPoint')
         IF(numberNodes.EQ.1) THEN
            CALL juDFT_error('Error: Single special k point provided. This does not make sense!')
         END IF
         kpts%numSpecialPoints = numberNodes
         IF(kpts%numSpecialPoints.GE.2) THEN
            DEALLOCATE(kpts%specialPoints)
            ALLOCATE(kpts%specialPoints(3,kpts%numSpecialPoints))
            ALLOCATE(kpts%specialPointNames(kpts%numSpecialPoints))
            DO i = 1, kpts%numSpecialPoints
               WRITE(xPathA,*) '/fleurInput/calculationSetup/bzIntegration/kPointCount/specialPoint[',i,']'
               valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA)))))
               kpts%specialPoints(1,i) = evaluatefirst(valueString)
               kpts%specialPoints(2,i) = evaluatefirst(valueString)
               kpts%specialPoints(3,i) = evaluatefirst(valueString)
               kpts%specialPointNames(i) = xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@name')
            END DO
         END IF
      ELSE
         DEALLOCATE(kpts%specialPoints)
         ALLOCATE(kpts%specialPoints(3,kpts%numSpecialPoints))
         ALLOCATE(kpts%specialPointNames(kpts%numSpecialPoints))
      END IF

      ! Option kPointList
      numberNodes = xmlGetNumberOfNodes('/fleurInput/calculationSetup/bzIntegration/kPointList')
      IF (numberNodes.EQ.1) THEN
         l_kpts = .TRUE.
         numberNodes = xmlGetNumberOfNodes('/fleurInput/calculationSetup/bzIntegration/kPointList/kPoint')
         kpts%nkpt = numberNodes
         kpts%l_gamma = .FALSE.
         ALLOCATE(kpts%bk(3,kpts%nkpt))
         ALLOCATE(kpts%wtkpt(kpts%nkpt))
         kpts%bk = 0.0
         kpts%wtkpt = 0.0
         kpts%specificationType = 3

         kpts%posScale = evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/bzIntegration/kPointList/@posScale'))
         weightScale = evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/bzIntegration/kPointList/@weightScale'))

         DO i = 1, kpts%nkpt
            WRITE(xPathA,*) '/fleurInput/calculationSetup/bzIntegration/kPointList/kPoint[',i,']'
            valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA)))))
            READ(valueString,*) kpts%bk(1,i), kpts%bk(2,i), kpts%bk(3,i)
            kpts%bk(:,i)=kpts%bk(:,i)/kpts%posScale
            kpts%wtkpt(i) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@weight'))
            kpts%wtkpt(i) = kpts%wtkpt(i) / weightScale
         END DO
         kpts%posScale=1.0
      END IF

      ! Read in optional SOC parameters if present

      xPathA = '/fleurInput/calculationSetup/soc'
      numberNodes = xmlGetNumberOfNodes(xPathA)

      noco%l_soc = .FALSE.
      noco%theta = 0.0
      noco%phi = 0.0

      IF (numberNodes.EQ.1) THEN
         noco%theta=evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@theta'))
         noco%phi=evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@phi'))
         noco%l_soc = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_soc'))
         noco%l_spav = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@spav'))
      END IF

      ! Read in optional noco parameters if present

      xPathA = '/fleurInput/calculationSetup/nocoParams'
      numberNodes = xmlGetNumberOfNodes(xPathA)

      noco%l_ss = .FALSE.
      noco%l_mperp = .FALSE.
      noco%l_constr = .FALSE.
      noco%mix_b = 0.0
      noco%qss = 0.0

      noco%l_relax(:) = .FALSE.
      noco%alphInit(:) = 0.0
      noco%alph(:) = 0.0
      noco%beta(:) = 0.0
      noco%b_con(:,:) = 0.0

      IF ((noco%l_noco).AND.(numberNodes.EQ.0)) THEN
         CALL juDFT_error('Error: l_noco is true but no noco parameters set in xml input file!')
      END IF

      IF (numberNodes.EQ.1) THEN
         noco%l_ss = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_ss'))
         noco%l_mperp = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_mperp'))
         noco%l_constr = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_constr'))

         noco%mix_b = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@mix_b'))

         valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/qss')))
         READ(valueString,*) noco%qss(1), noco%qss(2), noco%qss(3)

         !WRITE(*,*) 'Note: TODO: Calculation of q points!'

         numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/qsc')
         IF (numberNodes.EQ.1) THEN
            valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/qsc')))
            READ(valueString,*) qsc(1), qsc(2), qsc(3)
            DO i = 1, 3
               noco%qss(i) = noco%qss(i) / qsc(i)
            END DO
            !WRITE(*,*) 'Note: TODO: Integrate qsc directly into qss in input file!'
            !WRITE(*,*) '(no problem for users)'
         END IF
      END IF

      ! Read in optional 1D parameters if present

      xPathA = '/fleurInput/calculationSetup/oneDParams'
      numberNodes = xmlGetNumberOfNodes(xPathA)

      oneD%odd%d1 = .FALSE.

      IF (numberNodes.EQ.1) THEN
         oneD%odd%d1 = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@d1'))
         oneD%odd%M = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@MM'))
         oneD%odd%mb = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@vM'))
         oneD%odd%m_cyl = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@m_cyl'))
         oneD%odd%chi = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@chi'))
         oneD%odd%rot = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@rot'))
         oneD%odd%invs = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@invs1'))
         oneD%odd%zrfs = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@zrfs1'))
      END IF

      ! Read in optional expert modes switches

      xPathA = '/fleurInput/calculationSetup/expertModes'
      numberNodes = xmlGetNumberOfNodes(xPathA)

      input%gw = 0
      input%isec1 = 999999
      input%secvar = .FALSE.

      IF (numberNodes.EQ.1) THEN
         input%gw = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@gw'))
         input%isec1 = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@isec1'))
         input%secvar = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@secvar'))
      END IF

      ! Read in optional geometry optimization parameters

      xPathA = '/fleurInput/calculationSetup/geometryOptimization'
      numberNodes = xmlGetNumberOfNodes(xPathA)

      input%l_f = .FALSE.
      input%qfix = 0

      IF (numberNodes.EQ.1) THEN
         input%l_f = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_f'))
         input%xa = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@xa'))
         input%thetad = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@thetad'))
         input%epsdisp = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@epsdisp'))
         input%epsforce = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@epsforce'))

         numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@qfix')
         IF (numberNodes.EQ.1) THEN
            input%qfix = 1
            l_qfix = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@qfix'))
            IF (l_qfix) THEN
               input%qfix = 2
            END IF
         END IF
      END IF

      ! Read in optional general LDA+U parameters

      xPathA = '/fleurInput/calculationSetup/ldaU'
      numberNodes = xmlGetNumberOfNodes(xPathA)
      IF (numberNodes.EQ.1) THEN
         input%ldauLinMix = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_linMix'))
         input%ldauMixParam = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@mixParam'))
         input%ldauSpinf = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@spinf'))
      END IF

      ! Read in optional q point mesh for spin spirals

      xPathA = '/fleurInput/calculationSetup/spinSpiralQPointMesh'
      numberNodes = xmlGetNumberOfNodes(xPathA)

      !   IF ((noco%l_ss).AND.(numberNodes.EQ.0)) THEN
      !      call juDFT_error('Error: l_ss is true but no q point mesh set in xml input file!')
      !   END IF

      ! Read in optional E-Field parameters

      xPathA = '/fleurInput/calculationSetup/eField'
      numberNodes = xmlGetNumberOfNodes(xPathA)

      IF (numberNodes.EQ.1) THEN
         !input%efield%zsigma = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@zsigma'))
         !input%efield%sig_b(1) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@sig_b_1'))
         !input%efield%sig_b(2) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@sig_b_2'))
         !input%efield%plot_charge = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@plot_charge'))
         !input%efield%plot_rho = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@plot_rho'))
         !input%efield%autocomp = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@autocomp'))
         !input%efield%dirichlet = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@dirichlet'))
         !l_eV = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@eV'))

         CALL juDFT_error('Error: Reading input for E-Fields not yet implemented completely!')
         !      ALLOCATE(input%efield%sigEF(3*k1d, 3*k2d, nvac))
         !      input%efield%sigEF = 0.0
         !IF (l_eV) THEN
         !   input%efield%sig_b(:) = input%efield%sig_b/hartree_to_ev_const
         !         input%efield%sigEF(:,:,:) = input%efield%sigEF/hartree_to_ev_const
         !END IF
      END IF

      ! Read in optional energy parameter limits

      xPathA = '/fleurInput/calculationSetup/energyParameterLimits'
      numberNodes = xmlGetNumberOfNodes(xPathA)

      IF (numberNodes.EQ.1) THEN
         input%ellow = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@ellow'))
         input%elup = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@elup'))
      END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End of calculationSetup section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Start of cell section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Read in lattice parameters

      a1 = 0.0
      a2 = 0.0
      a3 = 0.0
      cell%z1 = 0.0
      dtild = 0.0
      input%film = .FALSE.
      latticeType = 'bulkLattice'
      latticeDef = 0
      symmetryDef = 0
      cell%latnam = 'any'

      numberNodes = xmlGetNumberOfNodes('/fleurInput/cell/filmLattice')

      IF (numberNodes.EQ.1) THEN
         input%film = .TRUE.
         latticeType = 'filmLattice'
      END IF

      xPathA = '/fleurInput/cell/'//latticeType
      numberNodes = xmlGetNumberOfNodes(xPathA)

      IF (numberNodes.EQ.1) THEN
         latticeScale = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@scale'))
         input%scaleCell = latticeScale
         valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@latnam')))
         READ(valueString,*) cell%latnam

         IF(input%film) THEN
            cell%z1 = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@dVac'))
            dtild = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@dTilda'))
            vacuum%dvac = cell%z1
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
         END IF

         numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/a1')
         IF (numberNodes.EQ.1) THEN
            latticeDef = 1
            input%scaleA1 = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/a1/@scale'))
            a1(1) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/a1'))
            a1(1) = a1(1) * input%scaleA1
            numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/a2')
            IF (numberNodes.EQ.1) THEN
               latticeDef = 2
               input%scaleA2 = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/a2/@scale'))
               a2(2) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/a2'))
               a2(2) = a2(2) * input%scaleA2
            END IF
            IF(.NOT.input%film) THEN
               input%scaleC = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/c/@scale'))
               a3(3) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/c'))
               a3(3) = a3(3) * input%scaleC
            END IF
         END IF

         numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/row-1')
         IF (numberNodes.EQ.1) THEN
            latticeDef = 3
            valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/row-1')))
            a1(1) = evaluateFirst(valueString)
            a1(2) = evaluateFirst(valueString)
            valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/row-2')))
            a2(1) = evaluateFirst(valueString)
            a2(2) = evaluateFirst(valueString)
            IF(.NOT.input%film) THEN
               input%scaleC = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/c/@scale'))
               a3(3) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/c'))
               a3(3) = a3(3) * input%scaleC
            END IF
         END IF

         numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/bravaisMatrix')
         IF (numberNodes.EQ.1) THEN
            latticeDef = 4
            valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/bravaisMatrix/row-1')))
            a1(1) = evaluateFirst(valueString)
            a1(2) = evaluateFirst(valueString)
            IF(.NOT.input%film) THEN
               a1(3) = evaluateFirst(valueString)
            END IF
            valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/bravaisMatrix/row-2')))
            a2(1) = evaluateFirst(valueString)
            a2(2) = evaluateFirst(valueString)
            IF(.NOT.input%film) THEN
               a2(3) = evaluateFirst(valueString)
            END IF
            valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/bravaisMatrix/row-3')))
            IF(.NOT.input%film) THEN
               a3(1) = evaluateFirst(valueString)
               a3(2) = evaluateFirst(valueString)
               a3(3) = evaluateFirst(valueString)
            ELSE
               !WRITE(*,*) 'Note: For film calculations only the upper left 2x2 part of the Bravais matrix is considered.'
            END IF
         END IF
      END IF ! Note: Further ways to define lattices might be added later. (1D lattice,...)

      ! Construction of amat requires additional information about the lattice
      ! and is done later (scroll down)!

      ! Read in symmetry parameters

      sym%namgrp = 'any'

      xPathA = '/fleurInput/cell/symmetry'
      numberNodes = xmlGetNumberOfNodes('/fleurInput/cell/symmetry')

      IF (numberNodes.EQ.1) THEN
         sym%symSpecType = 2
         symmetryDef = 1
         valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@spgrp')))
         READ(valueString,*) sym%namgrp
         sym%invs = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@invs'))
         sym%zrfs = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@zrfs'))
         sym%invs2 = sym%invs.AND.sym%zrfs

         IF (sym%namgrp.EQ.'any ') THEN
            sym%nop = 48
            ! Read in sym.out file if sym%namgrp='any' set.
            CALL rw_symfile('r',94,'sym.out',48,cell%bmat,&
               &                        mrotTemp,tauTemp,sym%nop,sym%nop2,sym%symor)
            IF (ALLOCATED(sym%mrot)) THEN
               DEALLOCATE(sym%mrot)
            END IF
            ALLOCATE(sym%mrot(3,3,sym%nop))
            IF (ALLOCATED(sym%tau)) THEN
               DEALLOCATE(sym%tau)
            END IF
            ALLOCATE(sym%tau(3,sym%nop))

            DO k = 1, sym%nop
               DO i = 1, 3
                  DO j = 1, 3
                     sym%mrot(j,i,k) = mrotTemp(j,i,k)
                  END DO
                  sym%tau(i,k) = tauTemp(i,k)
               END DO
            END DO
         ELSE
            n2spg = 0
            DO i = 1, 20
               IF (sym%namgrp.EQ.nammap(i)) n2spg = i
            END DO
            IF (n2spg == 0 ) THEN
               WRITE (*,*) 'Spacegroup ',sym%namgrp,' not known! Choose one of:'
               WRITE (*,'(20(a4,1x))') (nammap(i),i=1,20)
               CALL juDFT_error("Could not determine spacegroup!", calledby = "r_inpXML")
            END IF
            IF ((n2spg.GE.13).AND.(n2spg.LE.17)) THEN
               IF (.NOT.((cell%latnam.EQ.'hx3').OR.(cell%latnam.EQ.'hex'))) THEN
                  CALL juDFT_error("Use only hex or hx3 with p3, p3m1, p31m, p6 or p6m!", calledby ="r_inpXML")
               END IF
            END IF
            sym%nop = ord2(n2spg)
            IF (sym%invs) THEN
               sym%nop = 2*sym%nop
               IF (sym%zrfs.AND.(.NOT.l_c2(n2spg))) sym%nop = 2*sym%nop
            ELSE
               IF (sym%zrfs) sym%nop = 2*sym%nop
            END IF
            IF (ALLOCATED(sym%mrot)) THEN
               DEALLOCATE(sym%mrot)
            END IF
            ALLOCATE(sym%mrot(3,3,sym%nop))
            IF (ALLOCATED(sym%tau)) THEN
               DEALLOCATE(sym%tau)
            END IF
            ALLOCATE(sym%tau(3,sym%nop))
            CALL spg2set(sym%nop,sym%zrfs,sym%invs,sym%namgrp,cell%latnam,&
               &                     sym%mrot,sym%tau,sym%nop2,sym%symor)
         END IF
      END IF

      xPathA = '/fleurInput/cell/symmetryFile'
      numberNodes = xmlGetNumberOfNodes(xPathA)

      IF (numberNodes.EQ.1) THEN
         symmetryDef = 2
         sym%symSpecType = 1
         sym%nop = 48
         valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@filename')))

         CALL rw_symfile('r',94,TRIM(ADJUSTL(valueString)),48,cell%bmat,&
            &                     mrotTemp,tauTemp,sym%nop,sym%nop2,sym%symor)

         IF (ALLOCATED(sym%mrot)) THEN
            DEALLOCATE(sym%mrot)
         END IF
         ALLOCATE(sym%mrot(3,3,sym%nop))
         IF (ALLOCATED(sym%tau)) THEN
            DEALLOCATE(sym%tau)
         END IF
         ALLOCATE(sym%tau(3,sym%nop))

         sym%invs = .FALSE.
         sym%zrfs = .FALSE.

         DO k = 1, sym%nop
            absSum = 0
            DO i = 1, 3
               DO j = 1, 3
                  sym%mrot(j,i,k) = mrotTemp(j,i,k)
                  absSum = absSum + ABS(sym%mrot(j,i,k))
               END DO
               sym%tau(i,k) = tauTemp(i,k)
            END DO
            IF (absSum.EQ.3) THEN
               IF (ALL(sym%tau(:,k).EQ.0.0)) THEN
                  IF ((sym%mrot(1,1,k).EQ.-1).AND.(sym%mrot(2,2,k).EQ.-1).AND.(sym%mrot(3,3,k).EQ.-1)) sym%invs = .TRUE.
                  IF ((sym%mrot(1,1,k).EQ.1).AND.(sym%mrot(2,2,k).EQ.1).AND.(sym%mrot(3,3,k).EQ.-1)) sym%zrfs = .TRUE.
               END IF
            END IF
         END DO

         sym%invs2 = sym%invs.AND.sym%zrfs
      END IF

      xPathA = '/fleurInput/cell/symmetryOperations'
      numberNodes = xmlGetNumberOfNodes(xPathA)

      IF (numberNodes.EQ.1) THEN
         sym%symSpecType = 3
         symmetryDef = 3

         numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/symOp')
         sym%nop = numberNodes

         IF (ALLOCATED(sym%mrot)) DEALLOCATE(sym%mrot)
         ALLOCATE(sym%mrot(3,3,sym%nop))
         IF (ALLOCATED(sym%tau)) DEALLOCATE(sym%tau)
         ALLOCATE(sym%tau(3,sym%nop))
         sym%symor = .TRUE.
         DO i = 1, sym%nop
            WRITE(xPathB,*) TRIM(ADJUSTL(xPathA))//'/symOp[',i,']/row-1'
            valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB)))))
            READ(valueString,*) sym%mrot(1,1,i), sym%mrot(1,2,i), sym%mrot(1,3,i), sym%tau(1,i)

            WRITE(xPathB,*) TRIM(ADJUSTL(xPathA))//'/symOp[',i,']/row-2'
            valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB)))))
            READ(valueString,*) sym%mrot(2,1,i), sym%mrot(2,2,i), sym%mrot(2,3,i), sym%tau(2,i)

            WRITE(xPathB,*) TRIM(ADJUSTL(xPathA))//'/symOp[',i,']/row-3'
            valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB)))))
            READ(valueString,*) sym%mrot(3,1,i), sym%mrot(3,2,i), sym%mrot(3,3,i), sym%tau(3,i)

            IF ((sym%tau(1,i)**2 + sym%tau(2,i)**2 + sym%tau(3,i)**2).GT.1.e-8) THEN
               sym%symor = .FALSE.
            END IF
            DO j = 1,3
               IF (ABS(sym%tau(j,i)-0.33333) < 0.00001) THEN
                  sym%tau(j,i) = 1./3.
               ENDIF
               IF (ABS(sym%tau(j,i)+0.33333) < 0.00001) THEN
                  sym%tau(j,i) = -1./3.
               ENDIF
               IF (ABS(sym%tau(j,i)-0.66667) < 0.00001) THEN
                  sym%tau(j,i) = 2./3.
               ENDIF
               IF (ABS(sym%tau(j,i)+0.66667) < 0.00001) THEN
                  sym%tau(j,i) = -2./3.
               ENDIF
            ENDDO
         END DO
      END IF

      ! Calculate missing symmetry and cell properties and check consistency of parameters.

      ! Construction of amat
      SELECT CASE (latticeDef)
      CASE (1)
         IF (cell%latnam.EQ.'squ') THEN
            a2(2) = a1(1)
         ELSE IF (cell%latnam.EQ.'c-b') THEN
            aTemp = a1(1)
            a1(1) = aTemp*0.5*SQRT(2.0)
            a1(2) = -aTemp*0.5
            a2(1) = aTemp*0.5*SQRT(2.0)
            a2(2) = aTemp*0.5
         ELSE IF (cell%latnam.EQ.'hex') THEN
            aTemp = 0.5*a1(1)
            a1(1) = aTemp*SQRT(3.0)
            a1(2) = -aTemp
            a2(1) = a1(1)
            a2(2) = aTemp
         ELSE IF (cell%latnam.EQ.'hx3') THEN
            aTemp = 0.5*a1(1)
            a1(1) = aTemp
            a1(2) = -aTemp*SQRT(3.0)
            a2(1) = a1(1)
            a2(2) = -a1(2)
         ELSE IF (cell%latnam.EQ.'fcc') THEN
            aTemp = a1(1)
            a1(1) =       0.0 ; a1(2) = 0.5*aTemp ; a1(3) = 0.5*aTemp
            a2(1) = 0.5*aTemp ; a2(2) =       0.0 ; a2(3) = 0.5*aTemp
            a3(1) = 0.5*aTemp ; a3(2) = 0.5*aTemp ; a3(3) =       0.0
         ELSE IF (cell%latnam.EQ.'bcc') THEN
            aTemp = a1(1)
            a1(1) =-0.5*aTemp ; a1(2) = 0.5*aTemp ; a1(3) = 0.5*aTemp
            a2(1) = 0.5*aTemp ; a2(2) =-0.5*aTemp ; a2(3) = 0.5*aTemp
            a3(1) = 0.5*aTemp ; a3(2) = 0.5*aTemp ; a3(3) =-0.5*aTemp
         ELSE
            CALL juDFT_error("latnam is incompatible to parametrization of lattice (1)", calledby ="r_inpXML")
         END IF
      CASE (2)
         IF ((cell%latnam.EQ.'c-r').OR.(cell%latnam.EQ.'p-r')) THEN
            IF (cell%latnam.EQ.'c-r') THEN
               a1(2) = -a2(2)
               a2(1) =  a1(1)
            END IF
         ELSE
            CALL juDFT_error("latnam is incompatible to parametrization of lattice (2)", calledby ="r_inpXML")
         END IF
      CASE (3)
         IF (.NOT.(cell%latnam.EQ.'obl')) THEN
            CALL juDFT_error("latnam is incompatible to parametrization of lattice (3)", calledby ="r_inpXML")
         END IF
      CASE (4)
         IF (.NOT.(cell%latnam.EQ.'any')) THEN
            CALL juDFT_error("latnam is incompatible to parametrization of lattice (4)", calledby ="r_inpXML")
         END IF
      CASE DEFAULT
         CALL juDFT_error("Illegal lattice definition", calledby ="r_inpXML")
      END SELECT

      IF (latticeScale.EQ.0.0) latticeScale = 1.0
      IF (.NOT.input%film) vacuum%dvac = a3(3)
      vacuum%dvac = latticeScale*vacuum%dvac
      dtild = latticeScale*dtild

      cell%amat(:,1) = a1(:) * latticeScale
      cell%amat(:,2) = a2(:) * latticeScale
      cell%amat(:,3) = a3(:) * latticeScale

      CALL inv3(cell%amat,cell%bmat,cell%omtil)
      cell%bmat(:,:) = tpi_const*cell%bmat(:,:)
      cell%omtil = ABS(cell%omtil)

      IF (input%film.AND..NOT.oneD%odd%d1) THEN
         cell%vol = (cell%omtil/cell%amat(3,3))*vacuum%dvac
         cell%area = cell%omtil/cell%amat(3,3)
         !-odim
      ELSE IF (oneD%odd%d1) THEN
         cell%area = tpi_const*cell%amat(3,3)
         cell%vol = pi_const*(vacuum%dvac**2)*cell%amat(3,3)/4.0
         !+odim
      ELSE
         cell%vol = cell%omtil
         cell%area = cell%amat(1,1)*cell%amat(2,2)-cell%amat(1,2)*cell%amat(2,1)
         IF (cell%area < 1.0e-7) THEN
            IF (cell%latnam.EQ.'any') THEN
               cell%area = 1.
            ELSE
               CALL juDFT_error("area = 0",calledby ="r_inpXML")
            END IF
         END IF
      END IF

      ! Construction of missing symmetry information
      IF ((symmetryDef.EQ.2).OR.(symmetryDef.EQ.3)) THEN
         nop48 = 48
         ALLOCATE (invOps(sym%nop),multtab(sym%nop,sym%nop),optype(nop48))
         CALL check_close(sym%nop,sym%mrot,sym%tau,&
            &                      multtab,invOps,optype)

         CALL symproperties(nop48,optype,input%film,sym%nop,multtab,cell%amat,&
            &                        sym%symor,sym%mrot,sym%tau,&
            &                        invSym,sym%invs,sym%zrfs,sym%invs2,sym%nop,sym%nop2)
         DEALLOCATE(invOps,multtab,optype)
         IF (.NOT.input%film) sym%nop2=sym%nop
         IF (input%film) THEN
            DO n = 1, sym%nop
               DO i = 1, 3
                  IF (ABS(sym%tau(i,n)) > 0.00001) THEN
                     CALL juDFT_error("nonsymmorphic symmetries not yet implemented for films!",calledby ="r_inpXML")
                  ENDIF
               END DO
            END DO
         END IF
      END IF
      sym%invs2 = sym%invs.AND.sym%zrfs

      ALLOCATE (sym%invarop(atoms%nat,sym%nop),sym%invarind(atoms%nat))
      ALLOCATE (sym%multab(sym%nop,sym%nop),sym%invtab(sym%nop))
      ALLOCATE (sym%invsatnr(atoms%nat),sym%d_wgn(-3:3,-3:3,3,sym%nop))

      !some settings for film calculations
      vacuum%nmzd = 250
      vacuum%nmzxyd = 100
      vacuum%nvac = 2
      IF (sym%zrfs.OR.sym%invs) vacuum%nvac = 1
      IF (oneD%odd%d1) vacuum%nvac = 1
      cell%z1 = vacuum%dvac/2
      vacuum%nmz = vacuum%nmzd
      vacuum%delz = 25.0/vacuum%nmz
      IF (oneD%odd%d1) vacuum%delz = 20.0 / vacuum%nmz
      IF (vacuum%nmz.GT.vacuum%nmzd) CALL juDFT_error("nmzd",calledby ="inped")
      vacuum%nmzxy = vacuum%nmzxyd
      IF (vacuum%nmzxy.GT.vacuum%nmzxyd) CALL juDFT_error("nmzxyd",calledby ="inped")

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
            write (*,*) "read exc_id_x", exc_id_x
         ELSE
            exc_id_x = vxc_id_x
            write (*,*) "ignore exc_id_x", exc_id_x
         ENDIF
         
         IF(xmlGetNumberOfNodes(TRIM(xPathA) // '/@exc_correlation') == 1) THEN
            exc_id_c = evaluateFirstOnly(xmlGetAttributeValue(xPathA // '/@exc_correlation'))
            write (*,*) "read exc_id_c", exc_id_x
         ELSE
            exc_id_c = vxc_id_c
            write (*,*) "ignore exc_id_c", exc_id_x
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
            write (*,*) "read exc_id_x"
         ELSE
            exc_id_x = vxc_id_x
            write (*,*) "ignore exc_id_x"
         ENDIF
         
         IF(xmlGetNumberOfNodes(TRIM(xPathB) // '/@etot_correlation') == 1) THEN
            valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(xPathB) // '/@etot_correlation')))
            exc_id_c =  xc_f03_functional_get_number(TRIM(valueString))
            write (*,*) "read exc_id_c"
         ELSE
            exc_id_c = vxc_id_c
            write (*,*) "ignore exc_id_c"
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
      hybrid%l_hybrid=xcpot%is_hybrid()

      ALLOCATE(hybrid%lcutm1(atoms%ntype),hybrid%lcutwf(atoms%ntype),hybrid%select1(4,atoms%ntype))

      obsolete%lwb=.FALSE.
      IF (xcpot%vxc_is_gga()) THEN
         obsolete%ndvgrd=6
         obsolete%chng=-0.1e-11
      END IF

      IF (xcpot%vxc_is_gga()) THEN
         obsolete%ndvgrd = MAX(obsolete%ndvgrd,3)
      END IF

      hybrid%gcutm1 = input%rkmax - 0.5
      hybrid%tolerance1 = 1.0e-4
      hybrid%ewaldlambda = 3
      hybrid%lexp = 16
      hybrid%bands1 = DIMENSION%neigd

      numberNodes = xmlGetNumberOfNodes('/fleurInput/calculationSetup/prodBasis')
      IF (numberNodes==0) THEN
         IF (hybrid%l_hybrid) CALL judft_error("Mixed product basis input missing in inp.xml")
      ELSE
         hybrid%gcutm1=evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/prodBasis/@gcutm'))
         hybrid%tolerance1=evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/prodBasis/@tolerance'))
         hybrid%ewaldlambda=evaluateFirstIntOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/prodBasis/@ewaldlambda'))
         hybrid%lexp=evaluateFirstIntOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/prodBasis/@lexp'))
         hybrid%bands1=evaluateFirstIntOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/prodBasis/@bands'))
      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End of XC functional section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Start of species section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ALLOCATE (speciesNLO(numSpecies))
      ALLOCATE(atoms%speciesName(numSpecies))

      atoms%numStatesProvided = 0
      atoms%lapw_l(:) = -1
      atoms%n_u = 0

      DEALLOCATE(noel)
      ALLOCATE(noel(atoms%ntype))

      DO iSpecies = 1, numSpecies
         ! Attributes of species
         WRITE(xPathA,*) '/fleurInput/atomSpecies/species[',iSpecies,']'
         atoms%speciesName(iSpecies) = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@name')))
         atomicNumber = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@atomicNumber'))
         coreStates = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@coreStates'))
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

         DO iType = 1, atoms%ntype
            WRITE(xPathA,*) '/fleurInput/atomGroups/atomGroup[',iType,']/@species'
            valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA)))))
            IF(TRIM(ADJUSTL(atoms%speciesName(iSpecies))).EQ.TRIM(ADJUSTL(valueString))) THEN
               atoms%nz(iType) = atomicNumber
               IF (atoms%nz(iType).EQ.0) THEN
                  WRITE(*,*) 'Note: Replacing atomic number 0 by 1.0e-10 on atom type ', iType
                  atoms%zatom(iType) = 1.0e-10
               END IF
               noel(iType) = namat_const(atoms%nz(iType))
               atoms%zatom(iType) = atoms%nz(iType)
               atoms%rmt(iType) = radius
               atoms%jri(iType) = gridPoints
               atoms%dx(iType) = logIncrement
               atoms%lmax(iType) = lmax
               atoms%nlo(iType) = speciesNLO(iSpecies)
               atoms%ncst(iType) = coreStates
               atoms%lnonsph(iType) = lnonsphr
               atoms%lapw_l(iType) = lmaxAPW
               IF (flipSpin) THEN
                  atoms%nflip(iType) = 1
               ELSE
                  atoms%nflip(iType) = 0
               ENDIF
               atoms%bmu(iType) = magMom
               DO i = 1, numU
                  atoms%n_u = atoms%n_u + 1
                  atoms%lda_u(atoms%n_u)%l = ldau_l(i)
                  atoms%lda_u(atoms%n_u)%u = ldau_u(i)
                  atoms%lda_u(atoms%n_u)%j = ldau_j(i)
                  atoms%lda_u(atoms%n_u)%l_amf = l_amf(i)
                  atoms%lda_u(atoms%n_u)%atomType = iType
               END DO
               atomTypeSpecies(iType) = iSpecies
               IF(speciesRepAtomType(iSpecies).EQ.-1) speciesRepAtomType(iSpecies) = iType
            END IF
         END DO
      END DO

      atoms%lmaxd = MAXVAL(atoms%lmax(:))
      atoms%llod  = 0
      atoms%nlod = 0
      DO iType = 1, atoms%ntype
         atoms%nlod = MAX(atoms%nlod,atoms%nlo(iType))
      END DO
      atoms%nlod = MAX(atoms%nlod,2) ! for chkmt
      ALLOCATE(atoms%llo(atoms%nlod,atoms%ntype)); atoms%llo=-1
      ALLOCATE(atoms%ulo_der(atoms%nlod,atoms%ntype))
      ALLOCATE(atoms%l_dulo(atoms%nlod,atoms%ntype)) ! For what is this?

      DIMENSION%nstd = 29

      ALLOCATE(atoms%coreStateOccs(DIMENSION%nstd,2,atoms%ntype)); atoms%coreStateOccs=0.0
      ALLOCATE(atoms%coreStateNprnc(DIMENSION%nstd,atoms%ntype))
      ALLOCATE(atoms%coreStateKappa(DIMENSION%nstd,atoms%ntype))

      CALL enpara%init(atoms,input%jspins)
      enpara%evac0(:,:) = evac0Temp(:,:)

      DO iSpecies = 1, numSpecies
         ALLOCATE(speciesLLO(speciesNLO(iSpecies)))
         ALLOCATE(speciesLOeParams(speciesNLO(iSpecies)))
         ALLOCATE(speciesLOEDeriv(speciesNLO(iSpecies)))

         ! Attributes of energyParameters element of species
         WRITE(xPathA,*) '/fleurInput/atomSpecies/species[',iSpecies,']'
         speciesEParams(0) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/energyParameters/@s'))
         speciesEParams(1) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/energyParameters/@p'))
         speciesEParams(2) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/energyParameters/@d'))
         speciesEParams(3) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/energyParameters/@f'))

         ! Parameters for hybrid functionals
         IF (hybrid%l_hybrid) THEN
            WRITE(xPathA,*) '/fleurInput/atomSpecies/species[',iSpecies,']/prodBasis'
            numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA)))
            IF (numberNodes.NE.1) CALL judft_error("Parameters for mixed basis are missing for some specified")
            lcutm =evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@lcutm'))
            lcutwf=evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@lcutwf'))
            xPathA=xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@select')
            hybSelect(1) = NINT(evaluateFirst(xPathA))
            hybSelect(2) = NINT(evaluateFirst(xPathA))
            hybSelect(3) = NINT(evaluateFirst(xPathA))
            hybSelect(4) = NINT(evaluateFirst(xPathA))
         ENDIF

         ! Special switches for species
         ldaspecies=.FALSE.
         socscalespecies=1.0
         WRITE(xPathA,*) '/fleurInput/atomSpecies/species[',iSpecies,']/special'
         numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA)))
         IF (numberNodes==1) THEN
            ldaSpecies = evaluateFirstBoolOnly(TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@lda'))))
            socscaleSpecies   = evaluateFirstOnly(TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@socscale'))))
         ENDIF
         ! Explicitely provided core configurations

         coreConfigPresent = .FALSE.
         providedCoreStates = 0
         providedStates = 0
         coreStateOccs = 0.0
         speciesXMLElectronStates = noState_const
         speciesXMLCoreOccs = -1.0
         speciesXMLPrintCoreStates = .FALSE.
         WRITE(xPathA,*) '/fleurInput/atomSpecies/species[',iSpecies,']/electronConfig'
         numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA)))
         IF (numberNodes.EQ.1) THEN
            coreConfigPresent = .TRUE.
            valueString = xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/coreConfig')
            token = popFirstStringToken(valueString)
            DO WHILE (token.NE.' ')
               IF (token(1:1).EQ.'[') THEN
                  DO i = 1, 6
                     IF (TRIM(ADJUSTL(token)).EQ.nobleGasConfigList_const(i)) THEN
                        IF (providedCoreStates+nobleGasNumStatesList_const(i).GT.29) THEN
                           CALL judft_error('Error: Too many core states provided in xml input file!')
                        END IF
                        DO j = providedCoreStates+1, providedCoreStates+nobleGasNumStatesList_const(i)
                           coreStateOccs(j-providedCoreStates,:) = coreStateNumElecsList_const(j)
                           coreStateNprnc(j-providedCoreStates) = coreStateNprncList_const(j)
                           coreStateKappa(j-providedCoreStates) = coreStateKappaList_const(j)
                           speciesXMLElectronStates(j) = coreState_const
                        END DO
                        providedCoreStates = providedCoreStates + nobleGasNumStatesList_const(i)
                     END IF
                  END DO
               ELSE
                  DO i = 1, 29
                     IF (TRIM(ADJUSTL(token)).EQ.coreStateList_const(i)) THEN
                        providedCoreStates = providedCoreStates + 1
                        IF (providedCoreStates.GT.29) THEN
                           CALL judft_error('Error: Too many core states provided in xml input file!')
                        END IF
                        coreStateOccs(providedCoreStates,:) = coreStateNumElecsList_const(i)
                        coreStateNprnc(providedCoreStates) = coreStateNprncList_const(i)
                        coreStateKappa(providedCoreStates) = coreStateKappaList_const(i)
                        speciesXMLElectronStates(i) = coreState_const
                     END IF
                  END DO
               END IF
               token = popFirstStringToken(valueString)
            END DO
            numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/valenceConfig')
            providedStates = providedCoreStates
            IF(numberNodes.EQ.1) THEN
               valueString = xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/valenceConfig')
               token = popFirstStringToken(valueString)
               DO WHILE (token.NE.' ')
                  DO i = 1, 29
                     IF (TRIM(ADJUSTL(token)).EQ.coreStateList_const(i)) THEN
                        providedStates = providedStates + 1
                        IF (providedStates.GT.29) THEN
                           CALL judft_error('Error: Too many valence states provided in xml input file!')
                        END IF
                        coreStateOccs(providedStates,:) = coreStateNumElecsList_const(i)
                        coreStateNprnc(providedStates) = coreStateNprncList_const(i)
                        coreStateKappa(providedStates) = coreStateKappaList_const(i)
                        speciesXMLElectronStates(i) = valenceState_const
                     END IF
                  END DO
                  token = popFirstStringToken(valueString)
               END DO
            END IF
         END IF

         ! Explicitely provided core occupations

         WRITE(xPathA,*) '/fleurInput/atomSpecies/species[',iSpecies,']/electronConfig/stateOccupation'
         numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA)))
         IF (numberNodes.GE.1) THEN
            IF (.NOT.coreConfigPresent) THEN
               WRITE(*,*) 'Note: This just has to be implemented:'
               CALL judft_error('Error: Core occupation given while core config not set!')
            END IF
            DO i = 1, numberNodes
               WRITE(xPathB,*) TRIM(ADJUSTL(xPathA)),'[',i,']'
               valueString = xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@state')
               nprncTemp = 0
               kappaTemp = 0
               DO j = 1, 29
                  IF (TRIM(ADJUSTL(valueString)).EQ.coreStateList_const(j)) THEN
                     nprncTemp = coreStateNprncList_const(j)
                     kappaTemp = coreStateKappaList_const(j)
                     speciesXMLPrintCoreStates(j) = .TRUE.
                     speciesXMLCoreOccs(1,j) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@spinUp'))
                     speciesXMLCoreOccs(2,j) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@spinDown'))
                  END IF
               END DO
               DO j = 1, providedStates
                  IF ((nprncTemp.EQ.coreStateNprnc(j)).AND.(kappaTemp.EQ.coreStateKappa(j))) THEN
                     coreStateOccs(j,1) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@spinUp'))
                     coreStateOccs(j,2) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@spinDown'))
                  END IF
               END DO
            END DO
         END IF

         ! local orbitals

         WRITE(xPathA,*) '/fleurInput/atomSpecies/species[',iSpecies,']/lo'
         numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA)))
         iLLO = 1
         DO iLO = 1, numberNodes
            WRITE(xPathB,*) TRIM(ADJUSTL(xPathA)),'[',iLO,']/@l'
            WRITE(xPathC,*) TRIM(ADJUSTL(xPathA)),'[',iLO,']/@n'
            WRITE(xPathD,*) TRIM(ADJUSTL(xPathA)),'[',iLO,']/@type'
            WRITE(xPathE,*) TRIM(ADJUSTL(xPathA)),'[',iLO,']/@eDeriv'
            valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathD)))))
            lString = xmlGetAttributeValue(TRIM(ADJUSTL(xPathB)))
            nString = xmlGetAttributeValue(TRIM(ADJUSTL(xPathC)))
            CALL getIntegerSequenceFromString(TRIM(ADJUSTL(lString)), lNumbers, lNumCount)
            CALL getIntegerSequenceFromString(TRIM(ADJUSTL(nString)), nNumbers, nNumCount)
            IF(lNumCount.NE.nNumCount) THEN
               CALL judft_error('Error in LO input: l quantum number count does not equal n quantum number count')
            END IF
            loEDeriv = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathE))))
            DO i = 1, lNumCount
               speciesLLO(iLLO) = lNumbers(i)
               speciesLOeParams(iLLO) = nNumbers(i)
               IF(TRIM(ADJUSTL(valueString)).EQ.'HELO') THEN
                  speciesLOeParams(iLLO) = -speciesLOeParams(iLLO)
               END IF
               speciesLOEDeriv(iLLO) = loEDeriv
               iLLO = iLLO + 1
            END DO
            DEALLOCATE (lNumbers, nNumbers)
         END DO

         ! sort LOs according to l quantum number

         ALLOCATE (loOrderList(speciesNLO(iSpecies)),speciesLLOReal(speciesNLO(iSpecies)))
         DO iLLO = 1, speciesNLO(iSpecies)
            speciesLLOReal(iLLO) = speciesLLO(iLLO)
         END DO
         CALL sort(speciesNLO(iSpecies),speciesLLOReal,loOrderList)
         DEALLOCATE(speciesLLOReal)

         ! apply species parameters to atom groups

         DO iType = 1, atoms%ntype
            WRITE(xPathA,*) '/fleurInput/atomGroups/atomGroup[',iType,']/@species'
            valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA)))))
            IF(TRIM(ADJUSTL(atoms%speciesName(iSpecies))).EQ.TRIM(ADJUSTL(valueString))) THEN
               atoms%numStatesProvided(iType) = providedStates
               IF (coreConfigPresent) THEN
                  IF (providedCoreStates.NE.atoms%ncst(iType)) THEN
                     WRITE(6,*) " providedCoreStates:",providedCoreStates
                     WRITE(6,*) "atoms%ncst(iType)  :",atoms%ncst(iType)
                     CALL judft_error('Wrong number of core states provided!')
                  END IF
                  DO k = 1, providedStates !atoms%ncst(iType)
                     atoms%coreStateOccs(k,1,iType) = coreStateOccs(k,1)
                     atoms%coreStateOccs(k,2,iType) = coreStateOccs(k,2)
                     atoms%coreStateNprnc(k,iType) = coreStateNprnc(k)
                     atoms%coreStateKappa(k,iType) = coreStateKappa(k)
                     xmlElectronStates(k,iType) = speciesXMLElectronStates(k)
                     xmlPrintCoreStates(k,iType) = speciesXMLPrintCoreStates(k)
                     xmlCoreOccs (1,k,iType) = speciesXMLCoreOccs(1,k)
                     xmlCoreOccs (2,k,iType) = speciesXMLCoreOccs(2,k)
                  END DO
               END IF
               DO iLLO = 1, speciesNLO(iSpecies)
                  atoms%llo(iLLO,iType) = speciesLLO(loOrderList(iLLO))
                  atoms%ulo_der(iLLO,iType) = speciesLOEDeriv(loOrderList(iLLO))
                  atoms%llod = MAX(ABS(atoms%llo(iLLO,iType)),atoms%llod)
                  DO jsp = 1, input%jspins
                     enpara%ello0(iLLO,iType,jsp) = speciesLOeParams(loOrderList(iLLO))
                     IF (enpara%ello0(iLLO,iType,jsp)==NINT(enpara%ello0(iLLO,iType,jsp))) THEN
                        enpara%qn_ello(iLLO,iType,jsp)=NINT(enpara%ello0(iLLO,iType,jsp))
                        enpara%ello0(iLLO,iType,jsp)=0
                     ELSE
                        enpara%qn_ello(iLLO,iType,jsp)=0
                     ENDIF
                     enpara%skiplo(iType,jsp)=enpara%skiplo(iType,jsp)+(2*atoms%llo(iLLO,itype)+1)
                  END DO
               END DO
               ! Energy parameters
               DO jsp = 1, input%jspins
                  DO l = 0, 3
                     enpara%el0(l,iType,jsp) = speciesEParams(l)
                     IF (enpara%el0(l,iType,jsp)==NINT(enpara%el0(l,iType,jsp))) THEN
                        enpara%qn_el(l,iType,jsp)=NINT(enpara%el0(l,iType,jsp))
                        enpara%el0(l,iType,jsp)=0
                     ELSE
                        enpara%qn_el(l,iType,jsp)=0
                     ENDIF
                  END DO
                  DO l = 4,atoms%lmax(iType)
                     enpara%el0(l,iType,jsp) = enpara%el0(3,iType,jsp)
                  END DO
               END DO
               !Hybrid functional stuff
               hybrid%lcutm1(iType) = 4
               hybrid%lcutwf(iType) = atoms%lmax(iType) - atoms%lmax(iType) / 10
               hybrid%select1(:,iType) = (/4, 0, 4, 2 /)
               IF (hybrid%l_hybrid) THEN
                  hybrid%lcutm1(iType)=lcutm
                  hybrid%lcutwf(iType)=lcutwf
                  hybrid%select1(:,iType)=hybSelect
               ENDIF
               ! Explicit xc functional
               SELECT TYPE(xcpot)
               TYPE IS(t_xcpot_inbuild)
                  xcpot%lda_atom(iType)=ldaSpecies
               END SELECT
               noco%socscale(iType)=socscaleSpecies
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

      banddos%l_orb = .FALSE.
      banddos%orbCompAtom = 0
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

         !Read in atom group specific noco parameters
         xPathB = TRIM(ADJUSTL(xPathA))//'/nocoParams'
         numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathB)))
         IF (numberNodes.GE.1) THEN
            noco%l_relax(iType) = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@l_relax'))
            noco%alphInit(iType) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@alpha'))
            noco%alph(iType) = noco%alphInit(iType)
            noco%beta(iType) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@beta'))
            noco%b_con(1,iType) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@b_cons_x'))
            noco%b_con(2,iType) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@b_cons_y'))
         END IF
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
            lString=xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@theta')
            nString=xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@phi')
            ALLOCATE(t_forcetheo_mae::forcetheo)
            SELECT TYPE(forcetheo)
            TYPE IS(t_forcetheo_mae) !this is ok, we just allocated the type...
               CALL forcetheo%init(lString,nString)
            END SELECT
         ENDIF
         !spin-spiral dispersion
         xPathA = '/fleurInput/forceTheorem/spinSpiralDispersion'
         numberNodes = xmlGetNumberOfNodes(xPathA)
         IF (numberNodes.EQ.1) THEN
            ALLOCATE(t_forcetheo_ssdisp::forcetheo)
            SELECT TYPE(forcetheo)
            TYPE IS(t_forcetheo_ssdisp) !this is ok, we just allocated the type...
               CALL forcetheo%init(priv_read_q_list(xPathA))
            END SELECT
         ENDIF
         !dmi
         xPathA = '/fleurInput/forceTheorem/DMI'
         numberNodes = xmlGetNumberOfNodes(xPathA)
         IF (numberNodes.EQ.1) THEN
            lString=xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@theta')
            nString=xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@phi')
            ALLOCATE(t_forcetheo_dmi::forcetheo)
            SELECT TYPE(forcetheo)
            TYPE IS(t_forcetheo_dmi) !this is ok, we just allocated the type...
               CALL forcetheo%init(priv_read_q_list(TRIM(ADJUSTL(xPathA))//'/qVectors'),lstring,nstring)
            END SELECT
         ENDIF
         !jij
         xPathA = '/fleurInput/forceTheorem/Jij'
         numberNodes = xmlGetNumberOfNodes(xPathA)
         IF (numberNodes.EQ.1) THEN
            thetaj=evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@thetaj'))
            ALLOCATE(t_forcetheo_jij::forcetheo)
            SELECT TYPE(forcetheo)
            TYPE IS(t_forcetheo_jij) !this is ok, we just allocated the type...
               CALL forcetheo%init(priv_read_q_list(TRIM(ADJUSTL(xPathA))//'/qVectors'),thetaj,atoms)
            END SELECT
         ENDIF

      ELSE
         ALLOCATE(t_forcetheo::forcetheo) !default no forcetheorem type
      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End of force-theorem section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Start of output section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      banddos%dos = .FALSE.
      banddos%band = .FALSE.
      banddos%vacdos = .FALSE.
      sliceplot%slice = .FALSE.
      input%l_coreSpec = .FALSE.
      input%l_wann = .FALSE.

      input%vchk = .FALSE.
      input%cdinf = .FALSE.

      sliceplot%iplot = .FALSE.
      input%score = .FALSE.
      sliceplot%plpot = .FALSE.

      input%eonly = .FALSE.
      input%l_bmt = .FALSE.

      xPathA = '/fleurInput/output'
      numberNodes = xmlGetNumberOfNodes(xPathA)

      IF (numberNodes.EQ.1) THEN

         ! Read in general output switches
         banddos%dos = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@dos'))
         banddos%band = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@band'))
         banddos%vacdos = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@vacdos'))
         sliceplot%slice = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@slice'))
         input%l_coreSpec = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@coreSpec'))
         input%l_wann = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@wannier'))
         banddos%l_mcd = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@mcd'))

         ! Read in optional switches for checks

         xPathA = '/fleurInput/output/checks'
         numberNodes = xmlGetNumberOfNodes(xPathA)

         IF (numberNodes.EQ.1) THEN
            input%vchk = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@vchk'))
            input%cdinf = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@cdinf'))

         END IF

         ! Read in optional plotting parameters

         xPathA = '/fleurInput/output/plotting'
         numberNodes = xmlGetNumberOfNodes(xPathA)

         IF (numberNodes.EQ.1) THEN
            sliceplot%iplot = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@iplot'))
            input%score = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@score'))
            sliceplot%plpot = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@plplot'))
         END IF

         ! Read in optional specialOutput switches

         xPathA = '/fleurInput/output/specialOutput'
         numberNodes = xmlGetNumberOfNodes(xPathA)

         IF (numberNodes.EQ.1) THEN
            input%eonly = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@eonly'))
            input%l_bmt = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@bmt'))
         END IF

         ! Read in optional densityOfStates output parameters

         xPathA = '/fleurInput/output/densityOfStates'
         numberNodes = xmlGetNumberOfNodes(xPathA)

         IF ((banddos%dos).AND.(numberNodes.EQ.0)) THEN
            CALL juDFT_error("dos is true but densityOfStates parameters are not set!", calledby = "r_inpXML")
         END IF

         IF (numberNodes.EQ.1) THEN
            banddos%ndir = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@ndir'))
            banddos%e2_dos = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@minEnergy'))
            banddos%e1_dos = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@maxEnergy'))
            banddos%sig_dos = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@sigma'))
         END IF

         ! Read in optional vacuumDOS parameters

         xPathA = '/fleurInput/output/vacuumDOS'
         numberNodes = xmlGetNumberOfNodes(xPathA)

         IF ((banddos%vacdos).AND.(numberNodes.EQ.0)) THEN
            CALL juDFT_error("vacdos is true but vacDOS parameters are not set!", calledby = "r_inpXML")
         END IF

         vacuum%layers = 1
         input%integ = .FALSE.
         vacuum%starcoeff = .FALSE.
         vacuum%nstars = 0
         vacuum%locx = 0.0
         vacuum%locy = 0.0
         vacuum%nstm = 0
         vacuum%tworkf = 0.0
         IF ((banddos%vacdos).AND.(numberNodes.EQ.1)) THEN
            vacuum%layers = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@layers'))
            input%integ = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@integ'))
            vacuum%starcoeff = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@star'))
            vacuum%nstars = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nstars'))
            vacuum%locx(1) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locx1'))
            vacuum%locx(2) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locx2'))
            vacuum%locy(1) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locy1'))
            vacuum%locy(2) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locy2'))
            vacuum%nstm = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nstm'))
            vacuum%tworkf = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@tworkf'))
         END IF
         vacuum%layerd = vacuum%layers
         ALLOCATE(vacuum%izlay(vacuum%layerd,2))

         ! Read in optional chargeDensitySlicing parameters

         xPathA = '/fleurInput/output/chargeDensitySlicing'
         numberNodes = xmlGetNumberOfNodes(xPathA)

         IF ((sliceplot%slice).AND.(numberNodes.EQ.0)) THEN
            CALL juDFT_error("slice is true but chargeDensitySlicing parameters are not set!", calledby = "r_inpXML")
         END IF

         IF (numberNodes.EQ.1) THEN
            sliceplot%kk = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@numkpt'))
            sliceplot%e1s = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@minEigenval'))
            sliceplot%e2s = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@maxEigenval'))
            sliceplot%nnne = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nnne'))
            input%pallst = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@pallst'))
         END IF

         IF (banddos%band) THEN
            banddos%dos=.TRUE.
            banddos%ndir = -4
            WRITE(*,*) 'band="T" --> Overriding "dos" and "ndir"!'
         ENDIF

         ! Read in optional core spectrum (EELS) input parameters

         xPathA = '/fleurInput/output/coreSpectrum'
         numberNodes = xmlGetNumberOfNodes(xPathA)

         IF ((input%l_coreSpec).AND.(numberNodes.EQ.0)) THEN
            CALL juDFT_error("coreSpec is true but coreSpectrum parameters are not set!", calledby = "r_inpXML")
         END IF

         IF (numberNodes.EQ.1) THEN
            coreSpecInput%verb = 0
            tempBool = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@verbose'))
            IF(tempBool) coreSpecInput%verb = 1
            coreSpecInput%ek0 = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@eKin'))
            coreSpecInput%atomType = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@atomType'))
            coreSpecInput%lx = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@lmax'))
            coreSpecInput%edge = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@edgeType')))
            coreSpecInput%emn = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@eMin'))
            coreSpecInput%emx = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@eMax'))
            tempInt = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@numPoints'))
            coreSpecInput%ein = (coreSpecInput%emx - coreSpecInput%emn) / (tempInt - 1.0)
            xPathB = TRIM(ADJUSTL(xPathA))//'/edgeIndices'
            xPathB = TRIM(ADJUSTL(xPathB))//'/text()'
            valueString = xmlGetAttributeValue(TRIM(ADJUSTL(xPathB)))
            numTokens = countStringTokens(valueString)
            coreSpecInput%edgeidx(:) = 0
            IF(numTokens.GT.SIZE(coreSpecInput%edgeidx)) THEN
               CALL juDFT_error('More EELS edge indices provided than allowed.',calledby='r_inpXML')
            END IF
            DO i = 1, MAX(numTokens,SIZE(coreSpecInput%edgeidx))
               coreSpecInput%edgeidx(i) = evaluateFirstIntOnly(popFirstStringToken(valueString))
            END DO
         END IF

         ! Read in optional Wannier functions parameters

         xPathA = '/fleurInput/output/wannier'
         numberNodes = xmlGetNumberOfNodes(xPathA)

         IF ((input%l_wann).AND.(numberNodes.EQ.0)) THEN
            CALL juDFT_error("wannier is true but Wannier parameters are not set!", calledby = "r_inpXML")
         END IF

         IF (numberNodes.EQ.1) THEN
            wann%l_ms = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@ms'))
            wann%l_sgwf = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@sgwf'))
            wann%l_socgwf = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@socgwf'))
            wann%l_bs_comf = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@bsComf'))
            wann%l_atomlist = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@atomList'))
         END IF

         xPathA = '/fleurInput/output/wannier/bandSelection'
         numberNodes = xmlGetNumberOfNodes(xPathA)

         IF (numberNodes.EQ.1) THEN
            wann%l_byindex=.TRUE.
            wann%band_min(1) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@minSpinUp'))
            wann%band_max(1) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@maxSpinUp'))
            xPathA = '/fleurInput/output/wannier/bandSelection/@minSpinDown'
            numberNodes = xmlGetNumberOfNodes(xPathA)
            IF (numberNodes.EQ.1) THEN
               wann%band_min(2) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))))
            ELSE
               wann%band_min(2) = wann%band_min(1)
            END IF
            xPathA = '/fleurInput/output/wannier/bandSelection/@maxSpinDown'
            numberNodes = xmlGetNumberOfNodes(xPathA)
            IF (numberNodes.EQ.1) THEN
               wann%band_max(2) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))))
            ELSE
               wann%band_max(2) = wann%band_max(1)
            END IF
            wann%l_byindex = .TRUE.
            IF(input%l_wann) THEN
               IF (DIMENSION%neigd.LT.MAX(wann%band_max(1),wann%band_max(2))) THEN
                  DIMENSION%neigd = MAX(wann%band_max(1),wann%band_max(2))
               END IF
            END IF
         END IF

         xPathA = '/fleurInput/output/wannier/jobList'
         numberNodes = xmlGetNumberOfNodes(xPathA)

         IF (numberNodes.EQ.1) THEN
            xPathA = '/fleurInput/output/wannier/jobList/text()'

            ! Note: At the moment only 255 characters for the text in this node. Maybe this is not enough.
            valueString = xmlGetAttributeValue(TRIM(ADJUSTL(xPathA)))
            numTokens = countStringTokens(valueString)
            ALLOCATE(wann%jobList(numTokens))
            DO i = 1, numTokens
               wann%jobList(i) = popFirstStringToken(valueString)
            END DO
         END IF

         ! Read in optional magnetic circular dichroism parameters
         xPathA = '/fleurInput/output/magneticCircularDichroism'
         numberNodes = xmlGetNumberOfNodes(xPathA)

         IF ((banddos%l_mcd).AND.(numberNodes.EQ.0)) THEN
            CALL juDFT_error("mcd is true but magneticCircularDichroism parameters are not set!", calledby = "r_inpXML")
         END IF

         IF (numberNodes.EQ.1) THEN
            banddos%e_mcd_lo = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@energyLo'))
            banddos%e_mcd_up = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@energyUp'))
         END IF

         ! Read in optional parameter for unfolding bandstructure of supercell
         xPathA = '/fleurInput/output/unfoldingBand'
         numberNodes = xmlGetNumberOfNodes(xPathA)

         IF (numberNodes.EQ.1) THEN
            banddos%unfoldband = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@unfoldband'))
            banddos%s_cell_x = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@supercellX'))
            banddos%s_cell_y = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@supercellY'))
            banddos%s_cell_z = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@supercellZ'))
         END IF

      END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End of output section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Generate / fill wann%atomlist(:) array
      IF (wann%l_atomlist) THEN
         absSum = 0
         DO i = 1, atoms%nat
            IF (wannAtomList(i)) absSum = absSum + 1
         END DO
         wann%atomlist_num = absSum
         ALLOCATE(wann%atomlist(wann%atomlist_num))
         j = 1
         DO i = 1, atoms%nat
            IF (wannAtomList(i)) THEN
               wann%atomlist(j) = i
               j = j + 1
            END IF
         END DO
      ELSE
         wann%atomlist_num = atoms%nat
         ALLOCATE(wann%atomlist(wann%atomlist_num))
         DO i = 1, atoms%nat
            wann%atomlist(i) = i
         END DO
      END IF

      DEALLOCATE(wannAtomList)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Start of non-XML input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Read in enpara file iff available

      l_enpara = .FALSE.
      INQUIRE (file ='enpara',exist= l_enpara)
      IF (l_enpara) THEN
         CALL enpara%READ(atoms,input%jspins,input%film,.FALSE.)
      END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End of non-XML input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CALL xmlFreeResources()

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

   END SUBROUTINE setXCParameters

   SUBROUTINE getIntegerSequenceFromString(string, SEQUENCE, count)
      USE m_juDFT

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

   FUNCTION popFirstStringToken(line) RESULT(firstToken)

      IMPLICIT NONE

      CHARACTER(*), INTENT(INOUT) :: line
      CHARACTER(LEN = LEN(line)) :: firstToken

      INTEGER separatorIndex

      separatorIndex = 0
      line = TRIM(ADJUSTL(line))

      separatorIndex = INDEX(line,' ')
      IF (separatorIndex.LE.1) THEN
         firstToken = ' '
      ELSE
         firstToken = line(1:separatorIndex-1)
         line = line(separatorIndex+1:)
      END IF

   END FUNCTION popFirstStringToken

   FUNCTION countStringTokens(line) RESULT(tokenCount)

      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: line
      INTEGER                  :: tokenCount

      CHARACTER(LEN=LEN(line)) :: tempLine
      INTEGER separatorIndex

      tokenCount = 0

      tempLine = TRIM(ADJUSTL(line))

      DO WHILE (tempLine.NE.'')
         separatorIndex = 0
         separatorIndex = INDEX(tempLine,' ')
         IF (separatorIndex.EQ.0) THEN
            tempLine = ''
         ELSE
            tempLine = TRIM(ADJUSTL(tempLine(separatorIndex+1:)))
         END IF
         tokenCount = tokenCount + 1
      END DO

   END FUNCTION countStringTokens

   FUNCTION priv_read_q_list(path)RESULT(q)
      USE m_xmlIntWrapFort
      USE m_calculator
      IMPLICIT NONE
      CHARACTER(len=*),INTENT(in):: path
      REAL,ALLOCATABLE::q(:,:)

      INTEGER:: n,i
      CHARACTER(len=256):: xpatha,valueString

      n=xmlGetNumberOfNodes(TRIM(ADJUSTL(path))//'/q')
      ALLOCATE(q(3,n))
      DO i = 1, n
         PRINT *, path,'/q[',i,']'
         WRITE(xPathA,"(a,a,i0,a)") TRIM(ADJUSTL(path)),'/q[',i,']'
         PRINT *,xpatha
         valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA)))))
         PRINT *,"Q:",valueString
         READ(valueString,*) q(1,i),q(2,i),q(3,i)
      END DO
   END FUNCTION priv_read_q_list
END MODULE m_rinpXML
