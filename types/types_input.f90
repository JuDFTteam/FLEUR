!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_input
  use m_judfT
  use m_types_fleurinput
  IMPLICIT NONE
  private
     TYPE,extends(t_fleurinput):: t_input
      LOGICAL :: film=.false.
      INTEGER:: jspins=1
      INTEGER:: neig=0
      LOGICAL:: total
      REAL :: rkmax
      REAL :: zelec
      LOGICAL :: strho =.FALSE.
      LOGICAL :: cdinf =.false.
      LOGICAL :: vchk =.false.
      LOGICAL :: l_f =.false.
      LOGICAL :: eonly =.false.
      LOGICAL :: ctail =.true.
      INTEGER :: coretail_lmax =0 
      INTEGER :: itmax =9 
      REAL    :: minDistance=1.0e-5
      INTEGER :: maxiter=99
      INTEGER :: imix=7
      INTEGER :: gw=0
      INTEGER :: gw_neigd=0
      INTEGER :: qfix=0
      REAL    :: forcealpha =1.0 !< mixing parameter for geometry optimzer
      REAL    :: epsdisp =0.00001!< minimal displacement. If all displacements are < epsdisp stop
      REAL    :: epsforce =0.00001!< minimal force. If all forces <epsforce stop
      REAL    :: force_converged=0.00001
      INTEGER :: forcemix=3
      REAL    :: delgau =0.001  !TODO = tkb?
      REAL    :: alpha=0.05
      REAL    :: preconditioning_param=0.0
      REAL    :: spinf=2.0
      REAL    :: tkb=0.001
      LOGICAL :: gauss=.false.
      LOGICAL :: l_bmt=.false.
      !INTEGER:: scale
      INTEGER:: kcrel =0
      LOGICAL:: frcor =.false. !frozen core
      LOGICAL:: lflip=.false.
      LOGICAL:: score=.false.
      LOGICAL:: swsp=.false.
      LOGICAL:: tria=.false.
      LOGICAL:: integ=.false.
      LOGICAL:: pallst=.false.
      LOGICAL:: l_coreSpec=.false.
      LOGICAL:: l_wann=.false.
      LOGICAL:: secvar=.false.
      LOGICAL:: evonly=.false.
      LOGICAL:: l_inpXML=.true.
      REAL :: scaleCell=1.0
      REAL :: scaleA1=1.0
      REAL :: scaleA2=1.0
      REAL :: scaleC=1.0
      REAL :: ellow=-1.8
      REAL :: elup=1.0
      REAL :: fixed_moment = 0.0
      CHARACTER(LEN=100) :: comment="FLEUR calculation without a title"
      REAL, POINTER :: sigma !this is the difference in charge due to the electric field it points to the value stored in t_efield
      LOGICAL :: l_core_confpot=.true. !Former CPP_CORE
      LOGICAL :: l_useapw=.false.
      LOGICAL :: ldauLinMix=.false.
      REAL    :: ldauMixParam=0.1
      REAL    :: ldauSpinf=2.0
      LOGICAL :: l_rdmft=.false.
      REAL    :: rdmftOccEps=0.0
      INTEGER :: rdmftStatesBelow=0
      INTEGER :: rdmftStatesAbove=0
      INTEGER :: rdmftFunctional=0
    CONTAINS
      PROCEDURE :: read_xml=>read_xml_input
      procedure :: init
   END TYPE t_input
   public t_input
 CONTAINS
   SUBROUTINE read_xml_input(input,xml)
     USE m_types_xml
     use m_constants
     CLASS(t_input),INTENT(out):: input
     TYPE(t_xml),intent(in)    :: xml

     CHARACTER(len=100):: valueString,xpathA,xpathB
     INTEGER:: numberNodes,nodeSum
     
     !TODO! these switches should be in the inp-file
     input%l_core_confpot=.TRUE. !former CPP_CORE
     input%l_useapw=.FALSE.   !former CPP_APW
     input%comment =  xml%GetAttributeValue('/fleurInput/comment')
     input%rkmax = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/cutoffs/@Kmax'))

     xPathA = '/fleurInput/calculationSetup/cutoffs/@numbands'
     numberNodes = xml%GetNumberOfNodes(xPathA)
     input%neig = 0
      IF(numberNodes.EQ.1) THEN
         valueString = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA)))))
         IF(TRIM(ADJUSTL(valueString)).EQ.'all') THEN
            input%neig = -1
         ELSE
            READ(valueString,*) input%neig
         END IF
      END IF

     ! Read SCF loop parametrization

     input%itmax = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/scfLoop/@itmax'))
     input%minDistance = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/scfLoop/@minDistance'))
     input%maxiter = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/scfLoop/@maxIterBroyd'))
     valueString = TRIM(ADJUSTL(xml%GetAttributeValue('/fleurInput/calculationSetup/scfLoop/@imix')))
     SELECT CASE (valueString)
     CASE ('straight')
        input%imix = 0
     CASE ('Broyden1')
        input%imix = 3
     CASE ('Broyden2')
        input%imix = 5
     CASE ('Anderson')
        input%imix = 7
     CASE ("Pulay")
        input%imix = 9
     CASE ("pPulay")
        input%imix = 11
     CASE ("rPulay")
        input%imix = 13
     CASE ("aPulay")
        input%imix = 15
     CASE DEFAULT
        CALL juDFT_error('Error: unknown mixing scheme selected!')
     END SELECT

     input%alpha = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/scfLoop/@alpha'))
     input%preconditioning_param = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/scfLoop/@preconditioning_param'))
     input%spinf = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/scfLoop/@spinf'))
     ! Get parameters for core electrons
     input%ctail = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/coreElectrons/@ctail'))
     input%coretail_lmax = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/coreElectrons/@coretail_lmax'))
     input%frcor = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/coreElectrons/@frcor'))
     input%kcrel = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/coreElectrons/@kcrel'))
     ! Read in magnetism parameters
     input%jspins = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/magnetism/@jspins'))
     input%swsp = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/magnetism/@swsp'))
     input%lflip = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/magnetism/@lflip'))
     input%fixed_moment=evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/magnetism/@fixed_moment'))
     ! Read in optional expert modes switches
     xPathA = '/fleurInput/calculationSetup/expertModes'
     IF (xml%GetNumberOfNodes(xPathA)==1) THEN
        input%gw = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@gw'))
        input%secvar = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@secvar'))
     END IF
     ! Read in Brillouin zone integration parameters
     valueString = TRIM(ADJUSTL(xml%GetAttributeValue('/fleurInput/calculationSetup/bzIntegration/@mode')))
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
     numberNodes = xml%GetNumberOfNodes(xPathA)
     nodeSum = nodeSum + numberNodes
     IF (numberNodes.EQ.1) THEN
        input%tkb = evaluateFirstOnly(xml%GetAttributeValue(xPathA))
     END IF
     xPathA = '/fleurInput/calculationSetup/bzIntegration/@fermiSmearingTemp'
     numberNodes = xml%GetNumberOfNodes(xPathA)
     nodeSum = nodeSum + numberNodes
     IF (numberNodes.EQ.1) THEN
        input%tkb = evaluateFirstOnly(xml%GetAttributeValue(xPathA))
        input%tkb = boltzmann_Const * input%tkb
     END IF
     IF(nodeSum>1) THEN
        CALL juDFT_error('Error: Multiple fermi Smearing parameters provided in input file!')
     END IF
     xPathA = '/fleurInput/calculationSetup/bzIntegration/@valenceElectrons'
     numberNodes = xml%GetNumberOfNodes(xPathA)
     IF (numberNodes.EQ.1) THEN
        input%zelec = evaluateFirstOnly(xml%GetAttributeValue(xPathA))
     ELSE
        CALL juDFT_error('Error: Optionality of valence electrons in input file not yet implemented!')
     END IF
     input%film =  xml%GetNumberOfNodes('/fleurInput/cell/filmLattice')==1
     ! Read in optional geometry optimization parameters
     xPathA = '/fleurInput/calculationSetup/geometryOptimization'
     numberNodes = xml%GetNumberOfNodes(xPathA)
     IF (numberNodes.EQ.1) THEN
        input%l_f = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_f'))
        input%forcealpha = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@forcealpha'))
        input%epsdisp = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@epsdisp'))
        input%epsforce = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@epsforce'))
        input%forcemix = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@forcemix'))
        input%force_converged = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@force_converged'))
        input%qfix = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@qfix'))
     END IF
     ! Read in optional general LDA+U parameters
     xPathA = '/fleurInput/calculationSetup/ldaU'
     numberNodes = xml%GetNumberOfNodes(xPathA)
     IF (numberNodes.EQ.1) THEN
        input%ldauLinMix = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_linMix'))
        input%ldauMixParam = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@mixParam'))
        input%ldauSpinf = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@spinf'))
     END IF
     ! Read in RDMFT parameters
     xPathA = '/fleurInput/calculationSetup/rdmft'
     numberNodes = xml%GetNumberOfNodes(xPathA)
     IF (numberNodes.EQ.1) THEN
        input%l_rdmft = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_rdmft'))
        input%rdmftOccEps = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@occEps'))
        input%rdmftStatesBelow = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@statesBelow'))
        input%rdmftStatesAbove = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@statesAbove'))
        valueString = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@functional')))
        SELECT CASE (valueString)
        CASE ('Muller')
           input%rdmftFunctional = 1
        CASE DEFAULT
           STOP 'Error: unknown RDMFT functional selected!'
        END SELECT
     END IF
     ! Read in optional energy parameter limits
     xPathA = '/fleurInput/calculationSetup/energyParameterLimits'
     numberNodes = xml%GetNumberOfNodes(xPathA)
     IF (numberNodes.EQ.1) THEN
        input%ellow = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@ellow'))
        input%elup = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@elup'))
     END IF
     ! !! Start of output section
     xPathA = '/fleurInput/output'
     numberNodes = xml%GetNumberOfNodes(xPathA)
     IF (numberNodes.EQ.1) THEN
        ! Read in general output switches
         input%l_coreSpec = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@coreSpec'))
         input%l_wann = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@wannier'))
         ! Read in optional switches for checks
         xPathA = '/fleurInput/output/checks'
         numberNodes = xml%GetNumberOfNodes(xPathA)
         IF (numberNodes.EQ.1) THEN
            input%vchk = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@vchk'))
            input%cdinf = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@cdinf'))
         END IF
         ! Read in optional plotting parameters
         xPathA = '/fleurInput/output/plotting'
         numberNodes = xml%GetNumberOfNodes(xPathA)
         IF (numberNodes.EQ.1) THEN
            input%score = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@score'))
         END IF
         ! Read in optional specialOutput switches
         xPathA = '/fleurInput/output/specialOutput'
         numberNodes = xml%GetNumberOfNodes(xPathA)
         IF (numberNodes.EQ.1) THEN
            input%eonly = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@eonly'))
            input%l_bmt = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@bmt'))
         END IF
         ! Read in optional vacuumDOS parameters
         xPathA = '/fleurInput/output/vacuumDOS'
         numberNodes = xml%GetNumberOfNodes(xPathA)
         IF (numberNodes.EQ.1) THEN
            input%integ = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@integ'))
         END IF
         ! Read in optional chargeDensitySlicing parameters
         xPathA = '/fleurInput/output/chargeDensitySlicing'
         numberNodes = xml%GetNumberOfNodes(xPathA)
         IF (numberNodes.EQ.1) THEN
            input%pallst = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@pallst'))
         END IF
      END IF
    END SUBROUTINE read_xml_input

    subroutine init(input)
      class(t_input),intent(input)::input
    end subroutine init
END MODULE m_types_input
  
