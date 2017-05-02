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
  &                   atoms,obsolete,vacuum,input,stars,sliceplot,banddos,dimension,&
  &                   cell,sym,xcpot,noco,jij,oneD,hybrid,kpts,enpara,sphhar,l_opti,&
       &                   noel,namex,relcor,a1,a2,a3,scale,dtild,xmlElectronStates,&
       &                   xmlPrintCoreStates,xmlCoreOccs,atomTypeSpecies,speciesRepAtomType)

  USE iso_c_binding
  USE m_juDFT
  USE m_types
  USE m_symdata , ONLY : nammap, ord2, l_c2
  USE m_rwsymfile
  USE m_xmlIntWrapFort
  USE m_inv3
  USE m_spg2set
  USE m_closure, ONLY : check_close
  USE m_symproperties
  USE m_calculator
  USE m_icorrkeys
  USE m_constants
  USE m_hybridmix, ONLY : aMix_VHSE, omega_VHSE
  USE m_julia
  USE m_kptgen_hybrid
  USE m_od_kptsgen
  USE m_strgndim
  USE m_strgn
  USE m_od_strgn1
  USE m_localsym
  USE m_od_chisym
  USE m_ylm
  USE m_convndim
  USE m_dwigner
  USE m_mapatom
  USE m_od_mapatom
  USE m_inpeig
  USE m_prpqfft
  USE m_prpxcfft
  USE m_stepf
  USE m_cdn_io
  USE m_convn
  USE m_efield
  USE m_writegw
  USE m_apwsdim
  USE m_sort
  USE m_nocoInputCheck
  USE m_enpara,    ONLY : r_enpara

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
  TYPE(t_Jij),INTENT(INOUT)      :: Jij
  TYPE(t_cell),INTENT(INOUT)     :: cell
  TYPE(t_banddos),INTENT(INOUT)  :: banddos
  TYPE(t_sliceplot),INTENT(INOUT):: sliceplot
  TYPE(t_xcpot),INTENT(INOUT)    :: xcpot
  TYPE(t_noco),INTENT(INOUT)     :: noco
  TYPE(t_dimension),INTENT(OUT)  :: dimension
  TYPE(t_enpara)   ,INTENT(OUT)  :: enpara
  TYPE(t_sphhar)   ,INTENT(OUT)  :: sphhar
  LOGICAL, INTENT(OUT)           :: l_opti
  INTEGER,          ALLOCATABLE, INTENT(INOUT) :: xmlElectronStates(:,:)
  INTEGER,          ALLOCATABLE, INTENT(INOUT) :: atomTypeSpecies(:)
  INTEGER,          ALLOCATABLE, INTENT(INOUT) :: speciesRepAtomType(:)
  REAL,             ALLOCATABLE, INTENT(INOUT) :: xmlCoreOccs(:,:,:)
  LOGICAL,          ALLOCATABLE, INTENT(INOUT) :: xmlPrintCoreStates(:,:)
  CHARACTER(len=3), ALLOCATABLE, INTENT(INOUT) :: noel(:)
  CHARACTER(len=4), INTENT(OUT)  :: namex
  CHARACTER(len=12), INTENT(OUT) :: relcor
  REAL, INTENT(OUT)              :: a1(3),a2(3),a3(3)
  REAL, INTENT(OUT)              :: scale, dtild
  
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
  INTEGER ieq,i,k,na,n,ii
  REAL s3,ah,a,hs2,rest
  LOGICAL l_hyb,l_sym,ldum
  INTEGER :: ierr
  ! ..
  !...  Local Arrays
  !   CHARACTER :: helpchar(atoms%ntype)
  CHARACTER(len=  4) :: chntype
  CHARACTER(len= 41) :: chform
  CHARACTER(len=100) :: line
  
  !     added for HF and hybrid functionals
  REAL                  ::  aMix,omega
  INTEGER               :: idum
  CHARACTER (len=1)     ::  check

  CHARACTER(len=20) :: tempNumberString, speciesName
  CHARACTER(len=150) :: format
  CHARACTER(len=20) :: mixingScheme
  CHARACTER(len=10) :: loType
  LOGICAL :: kptGamma, l_relcor
  INTEGER :: iAtomType, startCoreStates, endCoreStates
  CHARACTER(len=100) :: xPosString, yPosString, zPosString
  CHARACTER(len=200) :: coreStatesString
  !   REAL :: tempTaual(3,atoms%nat)
  CHARACTER(len=7) :: coreStateList(29) !'(1s1/2)'
  CHARACTER(len=4) :: nobleGasConfigList(6) !'[He]'
  INTEGER          :: nobleGasNumStatesList(6)
  REAL             :: coreStateNumElecsList(29) !(per spin)
  INTEGER          :: coreStateNprncList(29)
  INTEGER          :: coreStateKappaList(29)
  REAL             :: coreStateOccs(29,2)
  INTEGER          :: coreStateNprnc(29), coreStateKappa(29)
  INTEGER          :: speciesXMLElectronStates(29)
  REAL             :: speciesXMLCoreOccs(2,29)
  LOGICAL          :: speciesXMLPrintCoreStates(29)

  INTEGER            :: iType, iLO, iSpecies, lNumCount, nNumCount, iLLO, jsp, j, l
  INTEGER            :: numberNodes, nodeSum, numSpecies, n2spg, n1, n2, ikpt, iqpt
  INTEGER            :: atomicNumber, coreStates, gridPoints, lmax, lnonsphr, lmaxAPW
  INTEGER            :: latticeDef, symmetryDef, nop48, firstAtomOfType, errorStatus
  INTEGER            :: loEDeriv, ntp1, ios, ntst, jrc, minNeigd, providedCoreStates, providedStates
  INTEGER            :: nv, nv2, kq1, kq2, kq3, nprncTemp, kappaTemp
  INTEGER            :: ldau_l, numVac
  INTEGER            :: speciesEParams(0:3)
  INTEGER            :: mrotTemp(3,3,48)
  REAL               :: tauTemp(3,48)
  REAL               :: bk(3)
  LOGICAL            :: flipSpin, l_eV, invSym, l_qfix, relaxX, relaxY, relaxZ, l_gga, l_kpts
  LOGICAL            :: l_vca, coreConfigPresent, l_enpara
  REAL               :: magMom, radius, logIncrement, qsc(3), latticeScale, dr
  REAL               :: aTemp, zp, rmtmax, sumWeight, ldau_u, ldau_j, tempReal
  REAL               :: weightScale, eParamUp, eParamDown
  LOGICAL            :: l_amf
  REAL, PARAMETER    :: boltzmannConst = 3.1668114e-6 ! value is given in Hartree/Kelvin
  REAL, PARAMETER    :: htr_eV   = 27.21138386 ! eV



  CHARACTER(LEN=200,KIND=c_char) :: schemaFilename, docFilename
  CHARACTER(LEN=255) :: valueString, lString, nString, token
  CHARACTER(LEN=255) :: xPathA, xPathB, xPathC, xPathD, xPathE
  CHARACTER(LEN=11)  :: latticeType
  CHARACTER(LEN=50)  :: versionString

  INTEGER, ALLOCATABLE :: lNumbers(:), nNumbers(:), speciesLLO(:)
  INTEGER, ALLOCATABLE :: loOrderList(:)
  CHARACTER(LEN=50), ALLOCATABLE :: speciesNames(:)
  INTEGER, ALLOCATABLE :: speciesNLO(:)
  INTEGER, ALLOCATABLE :: multtab(:,:), invOps(:), optype(:)
  INTEGER, ALLOCATABLE :: lmx1(:), nq1(:), nlhtp1(:)
  INTEGER, ALLOCATABLE :: speciesLOEDeriv(:)
  REAL,    ALLOCATABLE :: speciesLOeParams(:), speciesLLOReal(:)

  EXTERNAL prp_xcfft_box

  interface
     function dropInputSchema() bind(C, name="dropInputSchema")
       use iso_c_binding
       INTEGER(c_int) dropInputSchema
     end function dropInputSchema
  end interface

  errorStatus = 0
  errorStatus = dropInputSchema()
  IF(errorStatus.NE.0) THEN
     STOP 'Error: Cannot print out FleurInputSchema.xsd'
  END IF

  schemaFilename = "FleurInputSchema.xsd"//C_NULL_CHAR
  docFilename = "inp.xml"//C_NULL_CHAR

  DATA coreStateList / '(1s1/2)','(2s1/2)','(2p1/2)','(2p3/2)','(3s1/2)',&
       &                       '(3p1/2)','(3p3/2)','(3d3/2)','(3d5/2)','(4s1/2)',&
       &                       '(4p1/2)','(4p3/2)','(5s1/2)','(4d3/2)','(4d5/2)',&
       &                       '(5p1/2)','(5p3/2)','(6s1/2)','(4f5/2)','(4f7/2)',&
       &                       '(5d3/2)','(5d5/2)','(6p1/2)','(6p3/2)','(7s1/2)',&
       &                       '(5f5/2)','(5f7/2)','(6d3/2)','(6d5/2)' /

  DATA nobleGasConfigList / '[He]','[Ne]','[Ar]','[Kr]','[Xe]','[Rn]' /
  DATA nobleGasNumStatesList / 1, 4, 7, 12, 17, 24 /
  DATA coreStateNumElecsList / 1, 1, 1, 2, 1, 1, 2, 2, 3, 1, 1, 2, 1, 2,&
       &                               3, 1, 2, 1, 3, 4, 2, 3, 1, 2, 1, 3, 4, 2, 3 /
  DATA coreStateNprncList    / 1, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 4, 4,&
       &                               5, 5, 6, 4, 4, 5, 5, 6, 6, 7, 5, 5, 6, 6 /
  DATA coreStateKappaList    /-1,-1, 1,-2,-1, 1,-2, 2,-3,-1, 1,-2,-1, 2,-3,&
       &                               1,-2,-1, 3,-4, 2,-3, 1,-2,-1, 3,-4, 2,-3 /


  !TODO! these switches should be in the inp-file
  input%l_core_confpot=.true. !former CPP_CORE
  input%l_useapw=.false.   !former CPP_APW
  WRITE(*,*) 'Start reading of inp.xml file'
  CALL xmlInitInterface()
  CALL xmlParseSchema(schemaFilename)
  CALL xmlParseDoc(docFilename)
  CALL xmlValidateDoc()
  CALL xmlInitXPath()

  ! Check version of inp.xml
  versionString = xmlGetAttributeValue('/fleurInput/@fleurInputVersion')
  IF((TRIM(ADJUSTL(versionString)).NE.'0.27').AND.(TRIM(ADJUSTL(versionString)).NE.'0.28')) THEN
     STOP 'version number of inp.xml file is not compatible with this fleur version'
  END IF

  ! Get number of atoms, atom types, and atom species

  numberNodes = xmlGetNumberOfNodes('/fleurInput/atomGroups/atomGroup/relPos')
  numberNodes = numberNodes + xmlGetNumberOfNodes('/fleurInput/atomGroups/atomGroup/absPos')
  numberNodes = numberNodes + xmlGetNumberOfNodes('/fleurInput/atomGroups/atomGroup/filmPos')

  atoms%nat = numberNodes
  atoms%nat = numberNodes

  numberNodes = xmlGetNumberOfNodes('/fleurInput/atomGroups/atomGroup')

  atoms%ntype = numberNodes
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
  ALLOCATE(atoms%lda_u(atoms%ntype))
  ALLOCATE(atoms%bmu(atoms%ntype))
  ALLOCATE(atoms%relax(3,atoms%ntype))
  ALLOCATE(atoms%neq(atoms%ntype))
  ALLOCATE(atoms%taual(3,atoms%nat))
  ALLOCATE(atoms%pos(3,atoms%nat))
  ALLOCATE(atoms%rmt(atoms%ntype))
  ALLOCATE(atoms%numStatesProvided(atoms%ntype))

  ALLOCATE(atoms%ncv(atoms%ntype)) ! For what is this?
  ALLOCATE(atoms%ngopr(atoms%nat)) ! For what is this?
  ALLOCATE(atoms%lapw_l(atoms%ntype)) ! Where do I put this?
  ALLOCATE(atoms%invsat(atoms%nat)) ! Where do I put this?

  ALLOCATE(noco%soc_opt(atoms%ntype+2),noco%l_relax(atoms%ntype),noco%b_con(2,atoms%ntype))
  ALLOCATE(noco%alphInit(atoms%ntype),noco%alph(atoms%ntype),noco%beta(atoms%ntype))

  ALLOCATE (Jij%alph1(atoms%ntype),Jij%l_magn(atoms%ntype),Jij%M(atoms%ntype))
  ALLOCATE (Jij%magtype(atoms%ntype),Jij%nmagtype(atoms%ntype))

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
     IF (valueString(i:i).EQ.achar(10)) valueString(i:i) = ' ' !remove line breaks
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

  xPathA = '/fleurInput/calculationSetup/cutoffs/@GmaxXC'
  numberNodes = xmlGetNumberOfNodes(xPathA)
  xcpot%gmaxxc = stars%gmax
  IF(numberNodes.EQ.1) THEN
     xcpot%gmaxxc = evaluateFirstOnly(xmlGetAttributeValue(xPathA))
  END IF
  stars%gmaxInit = stars%gmax

  xPathA = '/fleurInput/calculationSetup/cutoffs/@numbands'
  numberNodes = xmlGetNumberOfNodes(xPathA)
  dimension%neigd = 0
  IF(numberNodes.EQ.1) THEN
     valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA)))))
     IF(TRIM(ADJUSTL(valueString)).EQ.'all') THEN
        STOP 'Feature to calculate all eigenfunctions not yet implemented.'
     ELSE
        READ(valueString,*) dimension%neigd
     END IF
  END IF

  ! Read SCF loop parametrization

  input%itmax = evaluateFirstIntOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@itmax'))
  input%minDistance = evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@minDistance'))
  input%maxiter = evaluateFirstIntOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@maxIterBroyd'))

  valueString = TRIM(ADJUSTL(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@imix')))
  SELECT CASE (valueString)
  CASE ('straight')
     input%imix = 1
  CASE ('Broyden1')
     input%imix = 3
  CASE ('Broyden2')
     input%imix = 5
  CASE ('Anderson')
     input%imix = 7
  CASE DEFAULT
     STOP 'Error: unknown mixing scheme selected!'
  END SELECT

  input%alpha = evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@alpha'))
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
  Jij%l_J = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/magnetism/@l_J'))
  input%swsp = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/magnetism/@swsp'))
  input%lflip = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/magnetism/@lflip'))

  dimension%jspd = input%jspins

  ! Read in Brillouin zone integration parameters

  kpts%nkpt3 = 0
  kpts%nmop = 0
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
     STOP 'Invalid bzIntegration mode selected!'
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
     STOP 'Error: Multiple fermi Smearing parameters provided in input file!'
  END IF

  xPathA = '/fleurInput/calculationSetup/bzIntegration/@valenceElectrons'
  numberNodes = xmlGetNumberOfNodes(xPathA)
  IF (numberNodes.EQ.1) THEN
     input%zelec = evaluateFirstOnly(xmlGetAttributeValue(xPathA))
  ELSE
     STOP 'Error: Optionality of valence electrons in input file not yet implemented!'
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
     kpts%nmop(1) = kpts%nkpt3(1)
     kpts%nmop(2) = kpts%nkpt3(2)
     kpts%nmop(3) = kpts%nkpt3(3)
     kpts%nkpt = kpts%nkpt3(1) * kpts%nkpt3(2) * kpts%nkpt3(3)
  END IF

  ! Option kPointCount
  xPathA = '/fleurInput/calculationSetup/bzIntegration/kPointCount'
  numberNodes = xmlGetNumberOfNodes(xPathA)
  IF (numberNodes.EQ.1) THEN
     l_kpts = .FALSE.
     kpts%nkpt = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@count'))
     kpts%l_gamma = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@gamma'))
     kpts%nkpt = kpts%nkpt

     ALLOCATE(kpts%bk(3,kpts%nkpt))
     ALLOCATE(kpts%wtkpt(kpts%nkpt))
     kpts%bk = 0.0
     kpts%wtkpt = 0.0
     kpts%posScale = 1.0

     numberNodes = xmlGetNumberOfNodes('/fleurInput/calculationSetup/bzIntegration/kPointCount/specialPoint')
     IF(numberNodes.EQ.1) THEN
        STOP 'Error: Single special k point provided. This does not make sense!'
     END IF
     kpts%numSpecialPoints = numberNodes
     IF(kpts%numSpecialPoints.GE.2) THEN
        DEALLOCATE(kpts%specialPoints)
        ALLOCATE(kpts%specialPoints(3,kpts%numSpecialPoints))
        ALLOCATE(kpts%specialPointNames(kpts%numSpecialPoints))
        DO i = 1, kpts%numSpecialPoints
           WRITE(xPathA,*) '/fleurInput/calculationSetup/bzIntegration/kPointCount/specialPoint[',i,']'
           valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA)))))
           READ(valueString,*) kpts%specialPoints(1,i), kpts%specialPoints(2,i), kpts%specialPoints(3,i)
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
     kpts%nkpt = numberNodes
     ALLOCATE(kpts%bk(3,kpts%nkpt))
     ALLOCATE(kpts%wtkpt(kpts%nkpt))
     kpts%bk = 0.0
     kpts%wtkpt = 0.0

     kpts%posScale = evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/bzIntegration/kPointList/@posScale'))
     weightScale = evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/bzIntegration/kPointList/@weightScale'))

     DO i = 1, kpts%nkpt
        WRITE(xPathA,*) '/fleurInput/calculationSetup/bzIntegration/kPointList/kPoint[',i,']'
        valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA)))))
        READ(valueString,*) kpts%bk(1,i), kpts%bk(2,i), kpts%bk(3,i)
        kpts%wtkpt(i) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@weight'))
        kpts%wtkpt(i) = kpts%wtkpt(i) / weightScale
     END DO
  END IF

  ! Read in optional SOC parameters if present

  xPathA = '/fleurInput/calculationSetup/soc'
  numberNodes = xmlGetNumberOfNodes(xPathA)

  noco%l_soc = .FALSE.
  noco%theta = 0.0
  noco%phi = 0.0
  noco%soc_opt(atoms%ntype+2) = .FALSE.
  noco%soc_opt(atoms%ntype+1) = .FALSE.

  IF (numberNodes.EQ.1) THEN
     noco%theta = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@theta'))
     noco%phi = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@phi'))
     noco%l_soc = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_soc'))
     noco%soc_opt(atoms%ntype+2) = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@spav'))
     noco%soc_opt(atoms%ntype+1) = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@off'))
  END IF

  ! Read in optional noco parameters if present

  xPathA = '/fleurInput/calculationSetup/nocoParams'
  numberNodes = xmlGetNumberOfNodes(xPathA)

  noco%l_ss = .FALSE.
  noco%l_mperp = .FALSE.
  noco%l_constr = .FALSE.
  Jij%l_disp = .FALSE.
  input%sso_opt = .FALSE.
  noco%mix_b = 0.0
  Jij%thetaJ = 0.0
  Jij%nmagn=1
  Jij%nsh = 0
  noco%qss = 0.0

  noco%l_relax(:) = .FALSE.
  noco%alphInit(:) = 0.0
  noco%alph(:) = 0.0
  noco%beta(:) = 0.0
  noco%b_con(:,:) = 0.0

  Jij%M(:) = 0.0
  Jij%l_magn(:) = .FALSE.
  Jij%l_wr=.TRUE.
  Jij%nqptd=1
  Jij%mtypes=1
  Jij%phnd=1

  IF ((noco%l_noco).AND.(numberNodes.EQ.0)) THEN
     STOP 'Error: l_noco is true but no noco parameters set in xml input file!'
  END IF

  IF (numberNodes.EQ.1) THEN
     noco%l_ss = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_ss'))
     noco%l_mperp = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_mperp'))
     noco%l_constr = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_constr'))
     Jij%l_disp = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_disp'))

     valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@sso_opt')))
     READ(valueString,'(2l1)') input%sso_opt(1),input%sso_opt(2)

     noco%mix_b = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@mix_b'))
     Jij%thetaJ = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@thetaJ'))
     Jij%nsh = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nsh'))

     valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/qss')))
     READ(valueString,*) noco%qss(1), noco%qss(2), noco%qss(3)

     WRITE(*,*) 'Note: TODO: Calculation of q points!'

     numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/qsc')
     IF (numberNodes.EQ.1) THEN
        valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/qsc')))
        READ(valueString,*) qsc(1), qsc(2), qsc(3)
        DO i = 1, 3
           noco%qss(i) = noco%qss(i) / qsc(i)
        END DO
        WRITE(*,*) 'Note: TODO: Integrate qsc directly into qss in input file!'
        WRITE(*,*) '(no problem for users)'
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
  obsolete%pot8 = .FALSE.
  input%isec1 = 999999
  input%secvar = .FALSE.

  IF (numberNodes.EQ.1) THEN
     input%gw = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@gw'))
     obsolete%pot8 = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@pot8'))
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

  ! Read in optional q point mesh for spin spirals

  xPathA = '/fleurInput/calculationSetup/spinSpiralQPointMesh'
  numberNodes = xmlGetNumberOfNodes(xPathA)

  !   IF ((noco%l_ss).AND.(numberNodes.EQ.0)) THEN
  !      STOP 'Error: l_ss is true but no q point mesh set in xml input file!'
  !   END IF

  Jij%nmopq = 0
  IF (numberNodes.EQ.1) THEN
     Jij%nmopq(1) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@qx'))
     Jij%nmopq(2) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@qy'))
     Jij%nmopq(3) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@qz'))
  END IF

  ! Read in optional E-Field parameters

  xPathA = '/fleurInput/calculationSetup/eField'
  numberNodes = xmlGetNumberOfNodes(xPathA)

  IF (numberNodes.EQ.1) THEN
     input%efield%zsigma = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@zsigma'))
     input%efield%sig_b(1) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@sig_b_1'))
     input%efield%sig_b(2) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@sig_b_2'))
     input%efield%plot_charge = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@plot_charge'))
     input%efield%plot_rho = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@plot_rho'))
     input%efield%autocomp = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@autocomp'))
     input%efield%dirichlet = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@dirichlet'))
     l_eV = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@eV'))

     STOP 'Error: Reading input for E-Fields not yet implemented completely!'
     !      ALLOCATE(input%efield%sigEF(3*k1d, 3*k2d, nvac))
     !      input%efield%sigEF = 0.0
     IF (l_eV) THEN
        input%efield%sig_b(:) = input%efield%sig_b/htr_eV
        !         input%efield%sigEF(:,:,:) = input%efield%sigEF/htr_eV
     END IF
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

  ALLOCATE(enpara%evac0(2,input%jspins))
  ALLOCATE(enpara%lchg_v(2,input%jspins),enpara%skiplo(atoms%ntype,input%jspins))
  ALLOCATE(enpara%enmix(input%jspins))

  enpara%lchg_v = .FALSE.

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
     scale = latticeScale
     valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@latnam')))
     READ(valueString,*) cell%latnam

     IF(input%film) THEN
        cell%z1 = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@dVac'))
        dtild = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@dTilda'))
        vacuum%dvac = cell%z1
        a3(3) = dtild
        enpara%evac0 = eVac0Default_const
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
              enpara%evac0(numVac,1) = eParamUp
              IF(input%jspins.GT.1) enpara%evac0(numVac,2) = eParamDown
              IF(i.EQ.1) THEN
                 enpara%evac0(3-numVac,1) = eParamUp
                 IF(input%jspins.GT.1) enpara%evac0(3-numVac,2) = eParamDown
              END IF
           END DO
        END IF
     END IF

     numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/a1')
     IF (numberNodes.EQ.1) THEN
        latticeDef = 1
        a1(1) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/a1'))
        numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/a2')
        IF (numberNodes.EQ.1) THEN
           latticeDef = 2
           a2(2) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/a2'))
        END IF
        IF(.NOT.input%film) THEN
           a3(3) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/c'))
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
           a3(3) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/c'))
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
           WRITE(*,*) 'Note: For film calculations only the upper left 2x2 part of the Bravais matrix is considered.'
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
           IF (.not.((cell%latnam.EQ.'hx3').OR.(cell%latnam.EQ.'hex'))) THEN
              CALL juDFT_error("Use only hex or hx3 with p3, p3m1, p31m, p6 or p6m!", calledby ="r_inpXML")
           END IF
        END IF
        sym%nop = ord2(n2spg)
        IF (sym%invs) THEN
           sym%nop = 2*sym%nop
           IF (sym%zrfs.and.(.not.l_c2(n2spg))) sym%nop = 2*sym%nop
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

     DO k = 1, sym%nop
        DO i = 1, 3
           DO j = 1, 3
              sym%mrot(j,i,k) = mrotTemp(j,i,k)
           END DO
           sym%tau(i,k) = tauTemp(i,k)
        END DO
     END DO
  END IF

  xPathA = '/fleurInput/cell/symmetryOperations'
  numberNodes = xmlGetNumberOfNodes(xPathA)

  IF (numberNodes.EQ.1) THEN
     symmetryDef = 3

     numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/symOp')
     sym%nop = numberNodes

     IF (ALLOCATED(sym%mrot)) THEN
        DEALLOCATE(sym%mrot)
     END IF
     ALLOCATE(sym%mrot(3,3,sym%nop))
     IF (ALLOCATED(sym%tau)) THEN
        DEALLOCATE(sym%tau)
     END IF
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
        a1(1) = aTemp*0.5*sqrt(2.0)
        a1(2) = -aTemp*0.5
        a2(1) = aTemp*0.5*sqrt(2.0)
        a2(2) = aTemp*0.5
     ELSE IF (cell%latnam.EQ.'hex') THEN
        aTemp = 0.5*a1(1)
        a1(1) = aTemp*sqrt(3.0)
        a1(2) = -aTemp
        a2(1) = a1(1)
        a2(2) = aTemp
     ELSE IF (cell%latnam.EQ.'hx3') THEN
        aTemp = 0.5*a1(1)
        a1(1) = aTemp
        a1(2) = -aTemp*sqrt(3.0)
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
  cell%omtil = abs(cell%omtil)

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
     IF (cell%area.lt.1.0e-7) THEN
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
     IF (.not.input%film) sym%nop2=sym%nop
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

  ! Read in xc functional parameters

  namex = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL('/fleurInput/xcFunctional/@name')))))
  l_relcor = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/xcFunctional/@relativisticCorrections'))

  relcor = 'non-relativi'
  IF (l_relcor) THEN
     relcor = 'relativistic'
  END IF

  xcpot%icorr = -99
  !  l91: lsd(igrd=0) with dsprs=1.d-19 in pw91.
  IF (namex.EQ.'exx ') xcpot%icorr = icorr_exx
  IF (namex.EQ.'hf  ') xcpot%icorr = icorr_hf
  IF (namex.EQ.'l91 ') xcpot%icorr = -1
  IF (namex.EQ.'x-a ') xcpot%icorr =  0
  IF (namex.EQ.'wign') xcpot%icorr =  1
  IF (namex.EQ.'mjw')  xcpot%icorr =  2
  IF (namex.EQ.'hl')   xcpot%icorr =  3
  IF (namex.EQ.'bh')   xcpot%icorr =  3
  IF (namex.EQ.'vwn')  xcpot%icorr =  4
  IF (namex.EQ.'pz')   xcpot%icorr =  5
  IF (namex.EQ.'pw91') xcpot%icorr =  6
  !  pbe: easy_pbe [Phys.Rev.Lett. 77, 3865 (1996)]
  !  rpbe: rev_pbe [Phys.Rev.Lett. 80, 890 (1998)]
  !  Rpbe: Rev_pbe [Phys.Rev.B 59, 7413 (1999)]
  IF (namex.eq.'pbe')  xcpot%icorr =  7
  IF (namex.eq.'rpbe') xcpot%icorr =  8
  IF (namex.eq.'Rpbe') xcpot%icorr =  9
  IF (namex.eq.'wc')   xcpot%icorr = 10
  !  wc: Wu & Cohen, [Phys.Rev.B 73, 235116 (2006)]
  IF (namex.eq.'PBEs') xcpot%icorr = 11
  !  PBEs: PBE for solids ( arXiv:0711.0156v2 )
  IF (namex.eq.'pbe0') xcpot%icorr = icorr_pbe0
  !  hse: Heyd, Scuseria, Ernzerhof, JChemPhys 118, 8207 (2003)
  IF (namex.eq.'hse ') xcpot%icorr = icorr_hse
  IF (namex.eq.'vhse') xcpot%icorr = icorr_vhse
  ! local part of HSE
  IF (namex.eq.'lhse') xcpot%icorr = icorr_hseloc

  IF (xcpot%icorr == -99) THEN
     WRITE(6,*) 'Name of XC-potential not recognized. Use one of:'
     WRITE(6,*) 'x-a,wign,mjw,hl,bh,vwn,pz,l91,pw91,pbe,rpbe,Rpbe,'//&
          &                'wc,PBEs,pbe0,hf,hse,lhse'
     CALL juDFT_error("Wrong name of XC-potential!", calledby="r_inpXML")
  END IF
  xcpot%igrd = 0
  IF (xcpot%icorr.GE.6) THEN
     xcpot%igrd = 1
     ! Am I sure about the following 3 lines? They were included in a similar section in rw_inp
     obsolete%lwb=.false.
     obsolete%ndvgrd=6
     obsolete%chng=-0.1e-11
  END IF
  input%krla = 0
  IF (l_relcor) THEN 
     input%krla = 1    
     IF (xcpot%igrd.EQ.1) THEN
        WRITE(6,'(18a,a4)') 'Use XC-potential: ',namex
        WRITE(6,*) 'only without relativistic corrections !'
        CALL juDFT_error("relativistic corrections + GGA not implemented",&
             &                         calledby ="r_inpXML")
     END IF
  END IF

  IF (xcpot%icorr.eq.0) WRITE(6,*) 'WARNING: using X-alpha for XC!'
  IF (xcpot%icorr.eq.1) WRITE(6,*) 'INFO   : using Wigner  for XC!'
  IF ((xcpot%icorr.eq.2).and.(namex.NE.'mjw')) WRITE(6,*) 'CAUTION: using MJW(BH) for XC!'

  IF ((xcpot%icorr.EQ.-1).OR.(xcpot%icorr.GE.6)) THEN
     obsolete%ndvgrd = max(obsolete%ndvgrd,3)
     IF ((xcpot%igrd.NE.0).AND.(xcpot%igrd.NE.1)) THEN 
        WRITE (6,*) 'selecting l91 or pw91 as XC-Potental you should'
        WRITE (6,*) ' have 2 lines like this in your inp-file:'
        WRITE (6,*) 'igrd=1,lwb=F,ndvgrd=4,idsprs=0,chng= 1.000E-16'
        WRITE (6,*) 'iggachk=1,idsprs0=1,idsprsl=1,idsprsi=1,idsprsv=1'
        CALL juDFT_error("igrd =/= 0 or 1",calledby ="inped")
     END IF
  END IF

  WRITE(*,*) 'Note: hybrid functionals input has to be realized at some point!'
  IF (namex.EQ.'vhse') THEN
     ! overwrite if sane input
     IF (aMix > 0 .and. aMix <= 1) THEN
        aMix = aMix_VHSE( aMix )
     ELSE
        aMix = aMix_VHSE()
     END IF
     ! overwrite if sane input
     IF (omega > 0) THEN
        omega = omega_VHSE(omega)
     ELSE
        omega = omega_VHSE()
     END IF
     !       WRITE (6,9041) namex,relcor,aMix,omega
  ELSE
     !       WRITE (6,9040) namex,relcor
  END IF

  l_gga = .FALSE.
  IF (xcpot%icorr.GE.6) l_gga = .TRUE.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End of XC functional section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Start of species section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ALLOCATE (speciesNames(numSpecies), speciesNLO(numSpecies))

  atoms%numStatesProvided = 0

  DO iSpecies = 1, numSpecies
     ! Attributes of species
     WRITE(xPathA,*) '/fleurInput/atomSpecies/species[',iSpecies,']'
     speciesNames(iSpecies) = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@name')))
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

     WRITE(*,*) 'APW+lo cutoffs ignored for the moment'

     numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/ldaU')
     ldau_l = -1
     ldau_u = 0.0
     ldau_j = 0.0
     l_amf = .FALSE.
     DO i = 1, numberNodes
        IF (i.GT.1) THEN
           WRITE (*,*) 'Not yet implemented:'
           STOP 'ERROR: More than 1 U parameter provided for a certain species.'
        END IF
        ldau_l = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/ldaU/@l'))
        ldau_u = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/ldaU/@U'))
        ldau_j = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/ldaU/@J'))
        l_amf = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/ldaU/@l_amf'))
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
           STOP 'Error in LO input: l quantum number count does not equal n quantum number count'
        END IF
        speciesNLO(iSpecies) = speciesNLO(iSpecies) + lNumCount
        DEALLOCATE (lNumbers, nNumbers)
     END DO

     DO iType = 1, atoms%ntype
        WRITE(xPathA,*) '/fleurInput/atomGroups/atomGroup[',iType,']/@species'
        valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA)))))
        IF(TRIM(ADJUSTL(speciesNames(iSpecies))).EQ.TRIM(ADJUSTL(valueString))) THEN
           atoms%nz(iType) = atomicNumber
           IF (atoms%nz(iType).EQ.0) THEN
              WRITE(*,*) 'Note: Replacing atomic number 0 by 1.0e-10 on atom type ', iType
              atoms%zatom(iType) = 1.0e-10
           END IF
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
           atoms%lda_u(iType)%l = ldau_l
           atoms%lda_u(iType)%u = ldau_u
           atoms%lda_u(iType)%j = ldau_j
           atoms%lda_u(iType)%l_amf = l_amf
           atomTypeSpecies(iType) = iSpecies
           IF(speciesRepAtomType(iSpecies).EQ.-1) speciesRepAtomType(iSpecies) = iType
        END IF
     END DO
  END DO

  atoms%lmaxd = maxval(atoms%lmax(:))
  atoms%llod  = 0
  atoms%nlod = 0
  DO iType = 1, atoms%ntype
     atoms%nlod = max(atoms%nlod,atoms%nlo(iType))
  END DO
  atoms%nlod = max(atoms%nlod,2) ! for chkmt
  ALLOCATE(atoms%llo(atoms%nlod,atoms%ntype))
  ALLOCATE(atoms%ulo_der(atoms%nlod,atoms%ntype))
  ALLOCATE(enpara%ello0(atoms%nlod,atoms%ntype,input%jspins))
  ALLOCATE(enpara%llochg(atoms%nlod,atoms%ntype,input%jspins))
  ALLOCATE(enpara%el0(0:atoms%lmaxd,atoms%ntype,input%jspins))
  ALLOCATE(enpara%lchange(0:atoms%lmaxd,atoms%ntype,input%jspins))
  ALLOCATE(atoms%l_dulo(atoms%nlod,atoms%ntype)) ! For what is this?

  enpara%el0 = 0.0
  enpara%ello0 = 0.0
  enpara%lchange = .FALSE.
  dimension%nstd = 29

  ALLOCATE(atoms%coreStateOccs(dimension%nstd,2,atoms%ntype))
  ALLOCATE(atoms%coreStateNprnc(dimension%nstd,atoms%ntype))
  ALLOCATE(atoms%coreStateKappa(dimension%nstd,atoms%ntype))

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
                 IF (TRIM(ADJUSTL(token)).EQ.nobleGasConfigList(i)) THEN
                    IF (providedCoreStates+nobleGasNumStatesList(i).GT.29) THEN
                       STOP 'Error: Too many core states provided in xml input file!'
                    END IF
                    DO j = providedCoreStates+1, providedCoreStates+nobleGasNumStatesList(i)
                       coreStateOccs(j-providedCoreStates,:) = coreStateNumElecsList(j)
                       coreStateNprnc(j-providedCoreStates) = coreStateNprncList(j)
                       coreStateKappa(j-providedCoreStates) = coreStateKappaList(j)
                       speciesXMLElectronStates(j) = coreState_const
                    END DO
                    providedCoreStates = providedCoreStates + nobleGasNumStatesList(i)
                 END IF
              END DO
           ELSE
              DO i = 1, 29
                 IF (TRIM(ADJUSTL(token)).EQ.coreStateList(i)) THEN
                    providedCoreStates = providedCoreStates + 1
                    IF (providedCoreStates.GT.29) THEN
                       STOP 'Error: Too many core states provided in xml input file!'
                    END IF
                    coreStateOccs(providedCoreStates,:) = coreStateNumElecsList(i)
                    coreStateNprnc(providedCoreStates) = coreStateNprncList(i)
                    coreStateKappa(providedCoreStates) = coreStateKappaList(i)
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
                 IF (TRIM(ADJUSTL(token)).EQ.coreStateList(i)) THEN
                    providedStates = providedStates + 1
                    IF (providedStates.GT.29) THEN
                       STOP 'Error: Too many valence states provided in xml input file!'
                    END IF
                    coreStateOccs(providedStates,:) = coreStateNumElecsList(i)
                    coreStateNprnc(providedStates) = coreStateNprncList(i)
                    coreStateKappa(providedStates) = coreStateKappaList(i)
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
           STOP 'Error: Core occupation given while core config not set!'
        END IF
        DO i = 1, numberNodes
           WRITE(xPathB,*) TRIM(ADJUSTL(xPathA)),'[',i,']'
           valueString = xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@state')
           nprncTemp = 0
           kappaTemp = 0
           DO j = 1, 29
              IF (TRIM(ADJUSTL(valueString)).EQ.coreStateList(j)) THEN
                 nprncTemp = coreStateNprncList(j)
                 kappaTemp = coreStateKappaList(j)
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
           STOP 'Error in LO input: l quantum number count does not equal n quantum number count'
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
        IF(TRIM(ADJUSTL(speciesNames(iSpecies))).EQ.TRIM(ADJUSTL(valueString))) THEN
           atoms%numStatesProvided(iType) = providedStates
           IF (coreConfigPresent) THEN
              IF (providedCoreStates.NE.atoms%ncst(iType)) THEN
                 STOP 'Wrong number of core states provided!'
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
              atoms%llod = max(abs(atoms%llo(iLLO,iType)),atoms%llod)
              DO jsp = 1, input%jspins
                 enpara%ello0(iLLO,iType,jsp) = speciesLOeParams(loOrderList(iLLO))
              END DO
           END DO
           ! Energy parameters
           DO jsp = 1, input%jspins
              DO l = 0, 3
                 enpara%el0(l,iType,jsp) = speciesEParams(l)
              END DO
              DO l = 4,atoms%lmax(iType)
                 enpara%el0(l,iType,jsp) = enpara%el0(3,iType,jsp)
              END DO
           END DO
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
        valueString = xmlGetAttributeValue(TRIM(ADJUSTL(xPathB)))
        atoms%taual(1,na) = evaluatefirst(valueString)
        atoms%taual(2,na) = evaluatefirst(valueString)
        atoms%taual(3,na) = evaluatefirst(valueString)
        atoms%pos(:,na) = matmul(cell%amat,atoms%taual(:,na))
     END DO

     numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/absPos')
     DO i = 1, numberNodes
        na = na + 1
        STOP 'absPos not yet implemented!'
     END DO

     numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/filmPos')
     DO i = 1, numberNodes
        na = na + 1
        WRITE(xPathB,*) TRIM(ADJUSTL(xPathA)),'/filmPos[',i,']'
        valueString = xmlGetAttributeValue(TRIM(ADJUSTL(xPathB)))
        atoms%taual(1,na) = evaluatefirst(valueString)
        atoms%taual(2,na) = evaluatefirst(valueString)
        atoms%taual(3,na) = evaluatefirst(valueString) / cell%amat(3,3)
        atoms%pos(:,na) = matmul(cell%amat,atoms%taual(:,na))
     END DO

     !Read in atom group specific noco parameters
     xPathB = TRIM(ADJUSTL(xPathA))//'/nocoParams'
     numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathB)))
     IF (numberNodes.GE.1) THEN
        noco%l_relax(iType) = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@l_relax'))
        Jij%l_magn(iType) = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@l_magn'))
        Jij%M(iType) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@M'))
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
!!! Start of output section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  banddos%dos = .FALSE.
  banddos%band = .FALSE.
  banddos%vacdos = .FALSE.
  sliceplot%slice = .FALSE.

  input%vchk = .FALSE.
  input%cdinf = .FALSE.
  obsolete%disp = .FALSE.

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

     ! Read in optional switches for checks

     xPathA = '/fleurInput/output/checks'
     numberNodes = xmlGetNumberOfNodes(xPathA)

     IF (numberNodes.EQ.1) THEN
        input%vchk = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@vchk'))
        input%cdinf = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@cdinf'))
        obsolete%disp = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@disp'))
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

  END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End of output section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Start of input postprocessing (calculate missing parameters)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Check the LO stuff and call setlomap (from inped):

  ALLOCATE(atoms%lo1l(0:atoms%llod,atoms%ntype))
  ALLOCATE(atoms%nlol(0:atoms%llod,atoms%ntype))

  DO iType = 1, atoms%ntype
     IF (atoms%nlo(iType).GE.1) THEN
        IF (input%secvar) THEN
           CALL juDFT_error("LO + sevcar not implemented",calledby ="r_inpXML")
        END IF
        IF (input%isec1<input%itmax) THEN
           CALL juDFT_error("LO + Wu not implemented",calledby ="r_inpXML")
        END IF
        IF (atoms%nlo(iType).GT.atoms%nlod) THEN
           WRITE (6,*) 'nlo(n) =',atoms%nlo(iType),' > nlod =',atoms%nlod
           CALL juDFT_error("nlo(n)>nlod",calledby ="r_inpXML")
        END IF
        DO j=1,atoms%nlo(iType)
           IF (.NOT.input%l_useapw) THEN
              IF (atoms%llo(j,iType).LT.0) THEN ! CALL juDFT_error("llo<0 ; compile with DCPP_APW!",calledby="inped")
                 WRITE(6,'(A)') 'Info: l_useapw not set.'
                 WRITE(6,'(A,I2,A,I2,A)') '      LO #', j, ' at atom type', iType, ' is an e-derivative.'
              END IF
           ENDIF
           IF ( (atoms%llo(j,iType).GT.atoms%llod).OR.(mod(-atoms%llod,10)-1).GT.atoms%llod ) THEN
              WRITE (6,*) 'llo(j,n) =',atoms%llo(j,iType),' > llod =',atoms%llod
              CALL juDFT_error("llo(j,n)>llod",calledby ="r_inpXML")
           END IF
        END DO

        ! Replace call to setlomap with the following 3 loops (preliminary).
        ! atoms%nlol and atoms%lo1l arrays are strange. This should be solved differently.
        DO l = 0,atoms%llod
           atoms%nlol(l,iType) = 0
           atoms%lo1l(l,iType) = 0
        END DO

        DO ilo = 1,atoms%nlod
           atoms%l_dulo(ilo,iType) = .FALSE.
        END DO

        DO ilo = 1,atoms%nlo(iType)
           if (input%l_useapW) THEN

              IF (atoms%ulo_der(ilo,iType).EQ.1) THEN
                 atoms%l_dulo(ilo,iType) = .TRUE.
              END IF
           endif
           WRITE(6,'(A,I2,A,I2)') 'I use',atoms%ulo_der(ilo,iType),'. derivative of l =',atoms%llo(ilo,iType)
           IF (atoms%llo(ilo,iType)>atoms%llod) CALL juDFT_error(" l > llod!!!",calledby="r_inpXML")
           l = atoms%llo(ilo,iType)
           IF (ilo.EQ.1) THEN
              atoms%lo1l(l,iType) = ilo
           ELSE
              IF (l.NE.atoms%llo(ilo-1,iType)) THEN
                 atoms%lo1l(l,iType) = ilo
              END IF
           END IF
           atoms%nlol(l,iType) = atoms%nlol(l,iType) + 1
        END DO
        WRITE (6,*) 'atoms%lapw_l(n) = ',atoms%lapw_l(iType)
     END IF

     enpara%skiplo(iType,:) = 0
     DO j = 1, atoms%nlo(iType)
        enpara%skiplo(iType,:) = enpara%skiplo(iType,1) + (2*atoms%llo(j,iType)+1)
     END DO
  END DO

  ! Check lda+u stuff (from inped)

  atoms%n_u = 0
  DO iType = 1,atoms%ntype
     IF (atoms%lda_u(iType)%l.GE.0)  THEN
        atoms%n_u = atoms%n_u + 1
        IF (atoms%nlo(iType).GE.1) THEN
           DO iLLO = 1, atoms%nlo(iType)
              IF ((abs(atoms%llo(iLLO,iType)).EQ.atoms%lda_u(iType)%l).AND.&
                   .NOT.atoms%l_dulo(iLLO,iType)) THEN
                 WRITE (*,*) 'LO and LDA+U for same l not implemented'
              END IF
           END DO
        END IF
     END IF
  END DO
  IF (atoms%n_u.GT.0) THEN
     IF (input%secvar) CALL juDFT_error("LDA+U and sevcar not implemented",calledby ="r_inpXML")
     IF (input%isec1<input%itmax) CALL juDFT_error("LDA+U and Wu not implemented",calledby ="r_inpXML")
     IF (noco%l_mperp) CALL juDFT_error("LDA+U and l_mperp not implemented",calledby ="r_inpXML")
  END IF

  ! Check DOS related stuff (from inped)

  IF ((banddos%ndir.LT.0).AND..NOT.banddos%dos) THEN
     CALL juDFT_error('STOP banddos: the inbuild dos-program  <0'//&
          ' can only be used if dos = true',calledby ="r_inpXML")
  END IF

  IF ((banddos%ndir.LT.0).AND.banddos%dos) THEN
     IF (banddos%e1_dos-banddos%e2_dos.LT.1e-3) THEN
        CALL juDFT_error("STOP banddos: no valid energy window for "//&
             "internal dos-program",calledby ="r_inpXML")
     END IF
     IF (banddos%sig_dos.LT.0) THEN
        CALL juDFT_error("STOP DOS: no valid broadening (sig_dos) for "//&
             "internal dos-PROGRAM",calledby ="r_inpXML")
     END IF
  END IF

  IF (banddos%vacdos) THEN
     IF (.NOT.banddos%dos) THEN
        CALL juDFT_error("STOP DOS: only set vacdos = .true. if dos = .true.",calledby ="r_inpXML")
     END IF
     IF (.NOT.vacuum%starcoeff.AND.(vacuum%nstars.NE.1))THEN
        CALL juDFT_error("STOP banddos: if stars = f set vacuum=1",calledby ="r_inpXML")
     END IF
     IF (vacuum%layers.LT.1) THEN
        CALL juDFT_error("STOP DOS: specify layers if vacdos = true",calledby ="r_inpXML")
     END IF
     DO i=1,vacuum%layers
        IF (vacuum%izlay(i,1).LT.1) THEN
           CALL juDFT_error("STOP DOS: all layers must be at z>0",calledby ="r_inpXML")
        END IF
     END DO
  END IF

  ! Check noco stuff and calculate missing noco parameters

  IF (noco%l_noco) THEN
     CALL nocoInputCheck(atoms,input,vacuum,jij,noco)

     IF (.not.jij%l_j.and.noco%l_ss) THEN

        !--->    the angle beta is relative to the spiral in a spin-spiral
        !--->    calculation, i.e. if beta = 0 for all atoms in the unit cell
        !--->    that means that the moments are "in line" with the spin-spiral
        !--->    (beta = qss * taual). note: this means that only atoms within
        !--->    a plane perpendicular to qss can be equivalent!

        na = 1
        DO iType = 1,atoms%ntype
           noco%phi = tpi_const*dot_product(noco%qss,atoms%taual(:,na))
           noco%alph(iType) = noco%alphInit(iType) + noco%phi
           na = na + atoms%neq(iType)
        END DO
     END IF
  END IF

  ! Calculate missing kpts parameters

  IF (.not.l_kpts) THEN
     IF (.NOT.oneD%odd%d1) THEN
        IF (jij%l_J) THEN
           n1=sym%nop
           n2=sym%nop2
           sym%nop=1
           sym%nop2=1
           CALL julia(sym,cell,input,noco,banddos,kpts,.FALSE.,.TRUE.)
           sym%nop=n1
           sym%nop2=n2
        ELSE IF(kpts%l_gamma .and. banddos%ndir .eq. 0) THEN
           STOP 'Error: No kpoint set generation for gamma=T yet!'
           CALL kptgen_hybrid(kpts%nmop(1),kpts%nmop(2),kpts%nmop(3),&
                kpts%nkpt,sym%invs,noco%l_soc,sym%nop,&
                sym%mrot,sym%tau)
        ELSE
           CALL julia(sym,cell,input,noco,banddos,kpts,.FALSE.,.TRUE.)
        END IF
     ELSE
        STOP 'Error: No kpoint set generation for 1D systems yet!'
        CALL od_kptsgen (kpts%nkpt)
     END IF
  END IF
  sumWeight = 0.0
  DO i = 1, kpts%nkpt
     sumWeight = sumWeight + kpts%wtkpt(i)
     kpts%bk(:,i) = kpts%bk(:,i) / kpts%posScale
  END DO
  kpts%posScale = 1.0
  DO i = 1, kpts%nkpt
     kpts%wtkpt(i) = kpts%wtkpt(i) / sumWeight
  END DO

  ! Generate missing general parameters

  minNeigd = NINT(0.75*input%zelec) + 1
  IF (noco%l_soc.and.(.not.noco%l_noco)) minNeigd = 2 * minNeigd
  IF (noco%l_soc.and.noco%l_ss) minNeigd=(3*minNeigd)/2
  IF (dimension%neigd.LT.minNeigd) THEN
     WRITE(*,*) 'numbands is too small. Setting parameter to default value.'
     dimension%neigd = minNeigd
     WRITE(*,*) 'changed numbands (dimension%neigd) to ',dimension%neigd
  END IF

  kpts%nkpt = kpts%nkpt
  dimension%nvd = 0 ; dimension%nv2d = 0
  stars%kq1_fft = 0 ; stars%kq2_fft = 0 ; stars%kq3_fft = 0
  obsolete%l_u2f = .FALSE.
  obsolete%l_f2u = .FALSE.
  !cell%aamat=matmul(transpose(cell%amat),cell%amat)
  cell%bbmat=matmul(cell%bmat,transpose(cell%bmat))
  jij%nqpt=1
  IF (jij%l_J) THEN
     WRITE(*,*) 'jij%nqpt has to be corrected. Not yet done!'
  END IF

  DO iqpt=1,jij%nqpt
     IF(jij%l_J) THEN
        WRITE(*,*) 'select noco%qss here. ...not yet implemented'
     END IF
     DO ikpt = 1,kpts%nkpt
        DO i = 1, 3
           bk(i) = kpts%bk(i,ikpt)
        END DO
        IF (input%film .OR.oneD%odd%d1) THEN
           WRITE(*,*) 'There might be additional work required for the k points here!'
           WRITE(*,*) '...in r_inpXML. See inpeig_dim for comparison!'
        END IF
        CALL apws_dim(bk(:),cell,input,noco,oneD,nv,nv2,kq1,kq2,kq3)
        stars%kq1_fft = max(kq1,stars%kq1_fft)
        stars%kq2_fft = max(kq2,stars%kq2_fft)
        stars%kq3_fft = max(kq3,stars%kq3_fft)
        dimension%nvd = max(dimension%nvd,nv)
        dimension%nv2d = max(dimension%nv2d,nv2)
     END DO ! k-pts
  END DO ! q-pts

  obsolete%lepr = 0

  IF (noco%l_noco) dimension%neigd = 2*dimension%neigd

  ! Generate missing parameters for atoms and calculate volume of the different regions

  cell%volint = cell%vol
  atoms%jmtd = maxval(atoms%jri(:))
  CALL ylmnorm_init(atoms%lmaxd)
  dimension%nspd=(atoms%lmaxd+1+mod(atoms%lmaxd+1,2))*(2*atoms%lmaxd+1)
  rmtmax = maxval(atoms%rmt(:))
  rmtmax = rmtmax*stars%gmax
  CALL convn_dim(rmtmax,dimension%ncvd)
  dimension%msh = 0
  ALLOCATE(atoms%rmsh(atoms%jmtd,atoms%ntype))
  ALLOCATE(atoms%volmts(atoms%ntype))
  ALLOCATE(atoms%vr0(atoms%ntype))  ! This should actually not be in the atoms type!
  atoms%vr0(:) = 0.0
  na = 0
  DEALLOCATE(noel)
  ALLOCATE(noel(atoms%ntype))
  DO iType = 1, atoms%ntype
     l_vca = .FALSE.
     INQUIRE (file="vca.in", exist=l_vca)
     IF (l_vca) THEN
        WRITE(*,*) 'Note: Implementation for virtual crystal approximation should be changed in r_inpXML!'
        WRITE(*,*) 'I am not sure whether the implementation actually makes any sense. It is from inped.'
        WRITE(*,*) 'We have to get rid of the file vca.in!'
        OPEN (17,file='vca.in',form='formatted')
        DO i= 1, iType
           READ (17,*,IOSTAT=ios) ntst,zp
           IF (ios /= 0) EXIT
           IF (ntst == iType) THEN
              atoms%zatom(iType) = atoms%zatom(iType) + zp
           END IF
        END DO
        CLOSE (17)
     END IF

     ! Calculate mesh for valence states
     radius = atoms%rmt(iType)*exp(atoms%dx(iType)*(1-atoms%jri(iType)))
     dr = exp(atoms%dx(iType))
     DO i = 1, atoms%jri(iType)
        atoms%rmsh(i,iType) = radius
        radius = radius*dr
     END DO
     ! Calculate mesh dimension for core states
     radius = atoms%rmt(iType)
     jrc = atoms%jri(iType)
     DO WHILE (radius < atoms%rmt(iType) + 20.0)
        jrc = jrc + 1
        radius = radius*dr
     END DO
     dimension%msh = max(dimension%msh,jrc)

     atoms%volmts(iType) = (fpi_const/3.0)*atoms%rmt(iType)**3
     cell%volint = cell%volint - atoms%volmts(iType)*atoms%neq(iType)

     noel(iType) = namat_const(atoms%nz(iType))
  END DO

  ! Read in enpara file iff available

  l_enpara = .FALSE.
  INQUIRE (file ='enpara',exist= l_enpara)
  IF (l_enpara) THEN
     OPEN (40,file ='enpara',form='formatted',status='old')
     DO jsp = 1,input%jspins
        CALL r_enpara(atoms,input,jsp,enpara)
     END DO !dimension%jspd
     CLOSE (40)
  END IF

  ! Dimensioning of lattice harmonics

  ALLOCATE(atoms%nlhtyp(atoms%ntype),atoms%ntypsy(atoms%nat))
  ALLOCATE(sphhar%clnu(1,1,1),sphhar%nlh(1),sphhar%llh(1,1),sphhar%nmem(1,1),sphhar%mlh(1,1,1))
  sphhar%ntypsd = 0
  IF (.NOT.oneD%odd%d1) THEN
     CALL local_sym(atoms%lmaxd,atoms%lmax,sym%nop,sym%mrot,sym%tau,&
          atoms%nat,atoms%ntype,atoms%neq,cell%amat,cell%bmat,&
          atoms%taual,sphhar%nlhd,sphhar%memd,sphhar%ntypsd,.true.,&
          atoms%nlhtyp,atoms%ntypsy,sphhar%nlh,sphhar%llh,&
          sphhar%nmem,sphhar%mlh,sphhar%clnu)
  ELSE IF (oneD%odd%d1) THEN
     WRITE(*,*) 'Note: I would be surprised if lattice harmonics generation works'
     WRITE(*,*) 'Dimensioning of local arrays seems to be inconsistent with routine local_sym'
     ALLOCATE (nq1(atoms%nat),lmx1(atoms%nat),nlhtp1(atoms%nat))
     ii = 1
     nq1=1
     DO i = 1,atoms%ntype
        DO j = 1,atoms%neq(i)
           lmx1(ii) = atoms%lmax(i)
           ii = ii + 1
        END DO
     END DO
     CALL local_sym(atoms%lmaxd,lmx1,sym%nop,sym%mrot,sym%tau,&
          atoms%nat,atoms%nat,nq1,cell%amat,cell%bmat,atoms%taual,&
          sphhar%nlhd,sphhar%memd,sphhar%ntypsd,.true.,nlhtp1,&
          atoms%ntypsy,sphhar%nlh,sphhar%llh,sphhar%nmem,&
          sphhar%mlh,sphhar%clnu)        
     ii = 1
     DO i = 1,atoms%ntype
        atoms%nlhtyp(i) = nlhtp1(ii)
        ii = ii + atoms%neq(i)
     END DO
     DEALLOCATE (nq1,lmx1,nlhtp1)
  END IF
  DEALLOCATE(sphhar%clnu,sphhar%nlh,sphhar%llh,sphhar%nmem,sphhar%mlh)

  ALLOCATE(sphhar%clnu(sphhar%memd,0:sphhar%nlhd,sphhar%ntypsd))
  ALLOCATE(sphhar%llh(0:sphhar%nlhd,sphhar%ntypsd))
  ALLOCATE(sphhar%mlh(sphhar%memd,0:sphhar%nlhd,sphhar%ntypsd))
  ALLOCATE(sphhar%nlh(sphhar%ntypsd),sphhar%nmem(0:sphhar%nlhd,sphhar%ntypsd))

  ! Dimensioning of stars

  IF (input%film.OR.(sym%namgrp.ne.'any ')) THEN
     CALL strgn1_dim(stars%gmax,cell%bmat,sym%invs,sym%zrfs,sym%mrot,&
          sym%tau,sym%nop,sym%nop2,stars%mx1,stars%mx2,stars%mx3,&
          stars%ng3,stars%ng2,oneD%odd)

  ELSE
     CALL strgn2_dim(stars%gmax,cell%bmat,sym%invs,sym%zrfs,sym%mrot,&
          sym%tau,sym%nop,stars%mx1,stars%mx2,stars%mx3,&
          stars%ng3,stars%ng2)
     oneD%odd%n2d = stars%ng2
     oneD%odd%nq2 = stars%ng2
     oneD%odd%nop = sym%nop
  END IF

  dimension%nn2d = (2*stars%mx1+1)*(2*stars%mx2+1)
  dimension%nn3d = (2*stars%mx1+1)*(2*stars%mx2+1)*(2*stars%mx3+1)
  IF (oneD%odd%d1) THEN
     oneD%odd%k3 = stars%mx3
     oneD%odd%nn2d = (2*(oneD%odd%k3)+1)*(2*(oneD%odd%M)+1)
  ELSE
     oneD%odd%k3 = 0
     oneD%odd%M = 0
     oneD%odd%nn2d = 1
     oneD%odd%mb = 0
  END IF
  ALLOCATE (stars%ig(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,-stars%mx3:stars%mx3))
  ALLOCATE (stars%ig2(stars%ng3))
  ALLOCATE (stars%kv2(2,stars%ng2),stars%kv3(3,stars%ng3))
  ALLOCATE (stars%nstr2(stars%ng2),stars%nstr(stars%ng3))
  ALLOCATE (stars%sk2(stars%ng2),stars%sk3(stars%ng3),stars%phi2(stars%ng2))
  ALLOCATE (stars%igfft(0:dimension%nn3d-1,2),stars%igfft2(0:dimension%nn2d-1,2))
  ALLOCATE (stars%rgphs(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,-stars%mx3:stars%mx3))
  ALLOCATE (stars%pgfft(0:dimension%nn3d-1),stars%pgfft2(0:dimension%nn2d-1))
  ALLOCATE (stars%ufft(0:27*stars%mx1*stars%mx2*stars%mx3-1),stars%ustep(stars%ng3))

  stars%sk2(:) = 0.0
  stars%phi2(:) = 0.0

  ! Initialize xc fft box

  CALL prp_xcfft_box(xcpot%gmaxxc,cell%bmat,stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft)

  ! Initialize missing 1D code arrays

  ALLOCATE (oneD%ig1(-oneD%odd%k3:oneD%odd%k3,-oneD%odd%M:oneD%odd%M))
  ALLOCATE (oneD%kv1(2,oneD%odd%n2d),oneD%nstr1(oneD%odd%n2d))
  ALLOCATE (oneD%ngopr1(atoms%nat),oneD%mrot1(3,3,oneD%odd%nop),oneD%tau1(3,oneD%odd%nop))
  ALLOCATE (oneD%invtab1(oneD%odd%nop),oneD%multab1(oneD%odd%nop,oneD%odd%nop))
  ALLOCATE (oneD%igfft1(0:oneD%odd%nn2d-1,2),oneD%pgfft1(0:oneD%odd%nn2d-1))

  ! Initialize missing hybrid functionals arrays

  ALLOCATE (hybrid%nindx(0:atoms%lmaxd,atoms%ntype))
  ALLOCATE (hybrid%select1(4,atoms%ntype),hybrid%lcutm1(atoms%ntype))
  ALLOCATE (hybrid%select2(4,atoms%ntype),hybrid%lcutm2(atoms%ntype),hybrid%lcutwf(atoms%ntype))
  ALLOCATE (hybrid%ddist(dimension%jspd))
  hybrid%ddist = 1.0

  ! Generate lattice harmonics

  IF (.NOT.oneD%odd%d1) THEN
     CALL local_sym(atoms%lmaxd,atoms%lmax,sym%nop,sym%mrot,sym%tau,&
          atoms%nat,atoms%ntype,atoms%neq,cell%amat,cell%bmat,atoms%taual,&
          sphhar%nlhd,sphhar%memd,sphhar%ntypsd,.FALSE.,&
          atoms%nlhtyp,atoms%ntypsy,sphhar%nlh,sphhar%llh,sphhar%nmem,sphhar%mlh,sphhar%clnu)
     sym%nsymt = sphhar%ntypsd
     oneD%mrot1(:,:,:) = sym%mrot(:,:,:)
     oneD%tau1(:,:) = sym%tau(:,:)
  ELSE IF (oneD%odd%d1) THEN
     WRITE(*,*) 'Note: I would be surprised if lattice harmonics generation works'
     WRITE(*,*) 'Dimensioning of local arrays seems to be inconsistent with routine local_sym'
     CALL od_chisym(oneD%odd,oneD%mrot1,oneD%tau1,sym%zrfs,sym%invs,sym%invs2,cell%amat)
     ALLOCATE (nq1(atoms%nat),lmx1(atoms%nat),nlhtp1(atoms%nat))
     ii = 1
     DO i = 1,atoms%ntype
        DO j = 1,atoms%neq(i)
           nq1(ii) = 1
           lmx1(ii) = atoms%lmax(i)
           ii = ii + 1
        END DO
     END DO
     CALL local_sym(atoms%lmaxd,lmx1,sym%nop,sym%mrot,sym%tau,&
          atoms%nat,atoms%nat,nq1,cell%amat,cell%bmat,atoms%taual,&
          sphhar%nlhd,sphhar%memd,sphhar%ntypsd,.FALSE.,&
          nlhtp1,atoms%ntypsy,sphhar%nlh,sphhar%llh,sphhar%nmem,sphhar%mlh,sphhar%clnu)
     sym%nsymt = sphhar%ntypsd
     ii = 1
     DO i = 1,atoms%ntype
        atoms%nlhtyp(i) = nlhtp1(ii)
        ii = ii + atoms%neq(i)
     END DO
     DEALLOCATE (lmx1,nlhtp1)
  END IF

  ! Calculate additional symmetry information

  IF (atoms%n_u.GT.0) THEN
     CALL d_wigner(sym%nop,sym%mrot,cell%bmat,3,sym%d_wgn)
  END IF
  IF (.NOT.oneD%odd%d1) THEN
     CALL mapatom(sym,atoms,cell,input,noco)
     oneD%ngopr1(1:atoms%nat) = atoms%ngopr(1:atoms%nat)
     !     DEALLOCATE ( nq1 )
  ELSE
     CALL juDFT_error("The oneD version is broken here. Compare call to mapatom with old version")
     CALL mapatom(sym,atoms,cell,input,noco)
     CALL od_mapatom(oneD,atoms,sym,cell)
  END IF

  ! Missing xc functionals initializations
  IF (xcpot%igrd.NE.0) THEN
     ALLOCATE (stars%ft2_gfx(0:dimension%nn2d-1),stars%ft2_gfy(0:dimension%nn2d-1))
     ALLOCATE (oneD%pgft1x(0:oneD%odd%nn2d-1),oneD%pgft1xx(0:oneD%odd%nn2d-1),&
          oneD%pgft1xy(0:oneD%odd%nn2d-1),&
          oneD%pgft1y(0:oneD%odd%nn2d-1),oneD%pgft1yy(0:oneD%odd%nn2d-1))
  ELSE
     ALLOCATE (stars%ft2_gfx(0:1),stars%ft2_gfy(0:1))
     ALLOCATE (oneD%pgft1x(0:1),oneD%pgft1xx(0:1),oneD%pgft1xy(0:1),&
          oneD%pgft1y(0:1),oneD%pgft1yy(0:1))
  END IF
  oneD%odd%nq2 = oneD%odd%n2d
  oneD%odi%nq2 = oneD%odd%nq2

  ! Store structure data

  CALL storeStructureIfNew(input, atoms, cell, vacuum, oneD, sym)

  ! Generate stars

  IF (input%film.OR.(sym%namgrp.NE.'any ')) THEN
     CALL strgn1(stars,sym,atoms,vacuum,sphhar,input,cell,xcpot)
     IF (oneD%odd%d1) THEN
        CALL od_strgn1(xcpot,cell,sym,oneD)
     END IF
  ELSE
     CALL strgn2(stars,sym,atoms,vacuum,sphhar,input,cell,xcpot)
  END IF

  ! Other small stuff

  input%strho = .FALSE.

  INQUIRE(file="cdn1",exist=l_opti)
  if (noco%l_noco) INQUIRE(file="rhomat_inp",exist=l_opti)
  l_opti=.not.l_opti
  IF ((sliceplot%iplot).OR.(input%strho).OR.(input%swsp).OR.&
       (input%lflip).OR.(obsolete%l_f2u).OR.(obsolete%l_u2f).OR.(input%l_bmt)) l_opti = .TRUE.

  IF (.NOT.l_opti) THEN
     !      The following call to inpeig should not be required.
     !      CALL inpeig(atoms,cell,input,oneD%odd%d1,kpts,enpara)
  END IF

  CALL prp_qfft(stars,cell,noco,input)

  IF (input%gw.GE.1) THEN
     CALL write_gw(atoms%ntype,sym%nop,1,input%jspins,atoms%nat,&
          atoms%ncst,atoms%neq,atoms%lmax,sym%mrot,cell%amat,cell%bmat,input%rkmax,&
          atoms%taual,atoms%zatom,cell%vol,1.0,DIMENSION%neigd,atoms%lmaxd,&
          atoms%nlod,atoms%llod,atoms%nlo,atoms%llo,noco%l_soc)
  END IF

  CALL  prp_xcfft(stars,input,cell,xcpot)

  IF (.NOT.sliceplot%iplot) THEN
     CALL stepf(sym,stars,atoms,oneD,input,cell,vacuum)
     CALL convn(DIMENSION,atoms,stars)
     CALL efield(atoms,DIMENSION,stars,sym,vacuum,cell,input)
  END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End of input postprocessing (calculate missing parameters)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  CALL xmlFreeResources()

  WRITE(*,*) 'Reading of inp.xml file finished'

  DEALLOCATE(speciesNames, speciesNLO)

END SUBROUTINE r_inpXML

SUBROUTINE getIntegerSequenceFromString(string, sequence, count)

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
               STOP 'String has wrong syntax (in getIntegerSequenceFromString)'
            END IF
            singleNumber = .FALSE.
            comma = .TRUE.
            READ(string(start:i-1),*) lastNumber
            count = count + 1
            start = i+1
         CASE ('-')
            IF ((start.EQ.i).OR.(dash).OR.(comma)) THEN
               STOP 'String has wrong syntax (in getIntegerSequenceFromString)'
            END IF
            singleNumber = .FALSE.
            dash = .TRUE.
            READ(string(start:i-1),*) lastNumber
            start = i+1
         CASE DEFAULT
            STOP 'String has wrong syntax (in getIntegerSequenceFromString)'
      END SELECT
   END DO
   IF(start.GT.length) THEN
      STOP 'String has wrong syntax (in getIntegerSequenceFromString)'
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
               sequence(index) = lastNumber
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

END MODULE m_rinpXML
