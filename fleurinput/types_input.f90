!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_input
  USE m_judfT
  USE m_constants
  USE m_types_fleurinput_base
  IMPLICIT NONE
  PRIVATE
  PUBLIC:: t_input
  TYPE,EXTENDS(t_fleurinput_base):: t_input
  LOGICAL :: film=.FALSE.
  LOGICAL :: l_real
  INTEGER :: jspins=1
  INTEGER :: neig=0
  REAL    :: rkmax=0.0
  REAL    :: gmax
  REAL    :: zelec
  LOGICAL :: eig66(2)=.FALSE.
  LOGICAL :: strho =.FALSE.
  LOGICAL :: cdinf =.FALSE.
  LOGICAL :: vchk =.FALSE.
  LOGICAL :: l_f =.FALSE.
  INTEGER :: f_level = -1
  !     f_level ==-1: Original force calculation
  !     f_level == 0: Original force calculation with FORCES and POSCAR printout
  !     f_level == 1: Forces from coretails calculated over whole unit cell
  !     f_level == 2: Kinetic energy surface term evaluated with IR functions
  !     f_level == 3: Surface term for density and potential discontinuity at the MT boundaries
  !     level 3 needs a large gmax cutoff, which backfires at the rest of the calculation
  LOGICAL :: eonly =.FALSE.
  LOGICAL :: ctail =.TRUE.
  INTEGER :: coretail_lmax =0
  INTEGER :: itmax =15
  REAL    :: minDistance=1.0e-5
  INTEGER :: maxiter=15
  INTEGER :: imix=7
  INTEGER :: gw=0
  INTEGER :: gw_neigd=0
  INTEGER :: qfix=0
  REAL    :: forcealpha =1.0 !< mixing parameter for geometry optimzer
  REAL    :: epsdisp =0.00001!< minimal displacement. If all displacements are < epsdisp stop
  REAL    :: epsforce =0.00001!< minimal force. If all forces <epsforce stop
  REAL    :: force_converged=0.00001
  INTEGER :: forcemix=2
  REAL    :: alpha=0.05
  REAL    :: preconditioning_param=0.0
  REAL    :: spinf=2.0
  REAL    :: tkb=0.001
  INTEGER :: bz_integration=BZINT_METHOD_HIST
  LOGICAL :: l_bloechl=.FALSE. !Are the bloechl corrections used for bz_integration=BZINT_METHOD_TETRA
  LOGICAL :: l_bmt=.FALSE.
  !INTEGER:: scale
  INTEGER:: kcrel =0
  LOGICAL:: frcor =.FALSE. !frozen core
  LOGICAL:: lflip=.FALSE.
  LOGICAL:: score=.FALSE.
  LOGICAL:: swsp=.FALSE.
  LOGICAL:: integ=.FALSE.
  LOGICAL:: pallst=.FALSE.
  LOGICAL:: l_coreSpec=.FALSE.
  LOGICAL:: l_wann=.FALSE.
  LOGICAL:: l_kpts_fullbz=.FALSE.
  LOGICAL:: secvar=.FALSE.
  LOGICAL:: evonly=.FALSE.
  !     LOGICAL:: l_inpXML=.TRUE.
  REAL :: fixed_moment = 0.0
  LOGICAL :: l_onlyMtStDen=.FALSE.
  CHARACTER(LEN=100) :: comment="FLEUR calculation without a title"
  LOGICAL :: l_core_confpot=.TRUE. !Former CPP_CORE
  LOGICAL :: l_useapw=.FALSE.
  LOGICAL :: ldauLinMix=.FALSE.
  REAL    :: ldauMixParam=0.05
  REAL    :: ldauSpinf=1.0
  LOGICAL :: ldauAdjEnpara=.FALSE.
  LOGICAL :: ldauSpinoffd=.FALSE.
  LOGICAL :: l_rdmft=.FALSE.
  REAL    :: rdmftOccEps=0.0
  INTEGER :: rdmftStatesBelow=0
  INTEGER :: rdmftStatesAbove=0
  INTEGER :: rdmftFunctional=0
  INTEGER :: lResMax = -1
CONTAINS
  PROCEDURE :: read_xml=>read_xml_input
  PROCEDURE :: init => init_input
  PROCEDURE ::mpi_bc =>mpi_bc_input
END TYPE t_input

CONTAINS
SUBROUTINE mpi_bc_input(this,mpi_comm,irank)
   USE m_mpi_bc_tool
   CLASS(t_input),INTENT(INOUT)::this
   INTEGER,INTENT(IN):: mpi_comm
   INTEGER,INTENT(IN),OPTIONAL::irank
   INTEGER ::rank
   IF (PRESENT(irank)) THEN
      rank=irank
   ELSE
      rank=0
   END IF
   CALL mpi_bc(this%eig66(1),rank,mpi_comm)
   CALL mpi_bc(this%film,rank,mpi_comm)
   CALL mpi_bc(this%l_real,rank,mpi_comm)
   CALL mpi_bc(this%jspins,rank,mpi_comm)
   CALL mpi_bc(this%neig,rank,mpi_comm)
   CALL mpi_bc(this%rkmax,rank,mpi_comm)
   CALL mpi_bc(this%gmax,rank,mpi_comm)
   CALL mpi_bc(this%zelec,rank,mpi_comm)
   CALL mpi_bc(this%strho,rank,mpi_comm)
   CALL mpi_bc(this%cdinf,rank,mpi_comm)
   CALL mpi_bc(this%vchk,rank,mpi_comm)
   CALL mpi_bc(this%l_f,rank,mpi_comm)
   CALL mpi_bc(this%f_level,rank,mpi_comm)
   CALL mpi_bc(this%eonly,rank,mpi_comm)
   CALL mpi_bc(this%ctail,rank,mpi_comm)
   CALL mpi_bc(this%coretail_lmax,rank,mpi_comm)
   CALL mpi_bc(this%itmax,rank,mpi_comm)
   CALL mpi_bc(this%minDistance,rank,mpi_comm)
   CALL mpi_bc(this%maxiter,rank,mpi_comm)
   CALL mpi_bc(this%imix,rank,mpi_comm)
   CALL mpi_bc(this%gw,rank,mpi_comm)
   CALL mpi_bc(this%gw_neigd,rank,mpi_comm)
   CALL mpi_bc(this%qfix,rank,mpi_comm)
   CALL mpi_bc(this%forcealpha,rank,mpi_comm)
   CALL mpi_bc(this%epsdisp,rank,mpi_comm)
   CALL mpi_bc(this%epsforce,rank,mpi_comm)
   CALL mpi_bc(this%force_converged,rank,mpi_comm)
   CALL mpi_bc(this%forcemix,rank,mpi_comm)
   CALL mpi_bc(this%alpha,rank,mpi_comm)
   CALL mpi_bc(this%preconditioning_param,rank,mpi_comm)
   CALL mpi_bc(this%spinf,rank,mpi_comm)
   CALL mpi_bc(this%tkb,rank,mpi_comm)
   CALL mpi_bc(this%bz_integration,rank,mpi_comm)
   CALL mpi_bc(this%l_bloechl,rank,mpi_comm)
   CALL mpi_bc(this%l_bmt,rank,mpi_comm)
   CALL mpi_bc(this%kcrel,rank,mpi_comm)
   CALL mpi_bc(this%frcor,rank,mpi_comm)
   CALL mpi_bc(this%lflip,rank,mpi_comm)
   CALL mpi_bc(this%score,rank,mpi_comm)
   CALL mpi_bc(this%swsp,rank,mpi_comm)
   CALL mpi_bc(this%integ,rank,mpi_comm)
   CALL mpi_bc(this%pallst,rank,mpi_comm)
   CALL mpi_bc(this%l_coreSpec,rank,mpi_comm)
   CALL mpi_bc(this%l_wann,rank,mpi_comm)
   CALL mpi_bc(this%secvar,rank,mpi_comm)
   CALL mpi_bc(this%evonly,rank,mpi_comm)
   CALL mpi_bc(this%l_onlyMtStDen,rank,mpi_comm)
   !    call mpi_bc(this%l_inpXML,rank,mpi_comm)
   CALL mpi_bc(this%fixed_moment ,rank,mpi_comm)
   CALL mpi_bc(this%l_core_confpot,rank,mpi_comm)
   CALL mpi_bc(this%l_useapw,rank,mpi_comm)
   CALL mpi_bc(this%ldauLinMix,rank,mpi_comm)
   CALL mpi_bc(this%ldauMixParam,rank,mpi_comm)
   CALL mpi_bc(this%ldauSpinf,rank,mpi_comm)
   CALL mpi_bc(this%ldauAdjEnpara,rank,mpi_comm)
   CALL mpi_bc(this%ldauSpinoffd,rank,mpi_comm)
   CALL mpi_bc(this%l_rdmft,rank,mpi_comm)
   CALL mpi_bc(this%rdmftOccEps,rank,mpi_comm)
   CALL mpi_bc(this%rdmftStatesBelow,rank,mpi_comm)
   CALL mpi_bc(this%rdmftStatesAbove,rank,mpi_comm)
   CALL mpi_bc(this%rdmftFunctional,rank,mpi_comm)
END SUBROUTINE mpi_bc_input

SUBROUTINE read_xml_input(this,xml)
   USE m_types_xml
   USE m_constants
   CLASS(t_input),INTENT(inout):: this
   TYPE(t_xml),INTENT(INOUT)  ::xml

   CHARACTER(len=100):: valueString,xpathA,xpathB,xPathC
   INTEGER:: numberNodes,nodeSum, i, numberNodesB,numberNodesC

   !TODO! these switches should be in the inp-file
   !this%l_core_confpot=.TRUE. !former CPP_CORE !Done (A.N.).
   this%l_useapw=.FALSE.   !former CPP_APW
   this%comment =  xml%GetAttributeValue('/fleurInput/comment')
   DO i = 1, LEN(this%comment)
      IF(IACHAR(this%comment(i:i)).LT.32) this%comment(i:i) = ' '
   END DO
   this%comment = TRIM(ADJUSTL(this%comment))
   this%rkmax = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/cutoffs/@Kmax'))
   this%gmax = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/cutoffs/@Gmax'))

   xPathA = '/fleurInput/calculationSetup/cutoffs/@numbands'
   numberNodes = xml%GetNumberOfNodes(xPathA)
   this%neig = 0
   IF(numberNodes.EQ.1) THEN
      valueString = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA)))))
      IF(TRIM(ADJUSTL(valueString)).EQ.'all') THEN
         this%neig = -1
      ELSE
         READ(valueString,*) this%neig
      END IF
   END IF

 ! Read SCF loop parametrization

   this%itmax = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/scfLoop/@itmax'))
   this%minDistance = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/scfLoop/@minDistance'))
   this%maxiter = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/scfLoop/@maxIterBroyd'))
   valueString = TRIM(ADJUSTL(xml%GetAttributeValue('/fleurInput/calculationSetup/scfLoop/@imix')))
   SELECT CASE (valueString)
      CASE ('straight')
         this%imix = 0
      CASE ('Broyden1')
         this%imix = 3
      CASE ('Broyden2')
         this%imix = 5
      CASE ('Anderson')
         this%imix = 7
      CASE ("Pulay")
         this%imix = 9
      CASE ("pPulay")
         this%imix = 11
      CASE ("rPulay")
         this%imix = 13
      CASE ("aPulay")
         this%imix = 15
      CASE DEFAULT
         CALL juDFT_error('Error: unknown mixing scheme selected!')
   END SELECT

   this%alpha = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/scfLoop/@alpha'))
   this%preconditioning_param = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/scfLoop/@precondParam'))
   this%spinf = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/scfLoop/@spinf'))
   ! Get parameters for core electrons
   this%ctail = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/coreElectrons/@ctail'))
   this%coretail_lmax = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/coreElectrons/@coretail_lmax'))
   IF (xml%GetNumberOfNodes('/fleurInput/calculationSetup/coreElectrons/@l_core_confpot')==1) THEN
      this%l_core_confpot = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/coreElectrons/@l_core_confpot'))
   END IF
   this%frcor = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/coreElectrons/@frcor'))
   this%kcrel = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/coreElectrons/@kcrel'))
   ! Read in magnetism parameters
   this%jspins = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/magnetism/@jspins'))
   this%swsp = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/magnetism/@swsp'))
   this%lflip = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/magnetism/@lflip'))
   IF (xml%versionNumber>31) &
        this%l_onlyMtStDen=evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/magnetism/@l_onlyMtStDen'))
   this%fixed_moment=evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/magnetism/@fixed_moment'))
   ! Read in optional expert modes switches
   xPathA = '/fleurInput/calculationSetup/expertModes'
   IF (xml%GetNumberOfNodes(xPathA)==1) THEN
      xPathB = TRIM(ADJUSTL(xPathA))//'/@gw'
      xPathC = TRIM(ADJUSTL(xPathA))//'/@spex'
      numberNodesB = xml%GetNumberOfNodes(xPathB)
      numberNodesC = xml%GetNumberOfNodes(xPathC)
      IF((numberNodesB.EQ.1).AND.(numberNodesC.EQ.1)) THEN
         CALL juDFT_error("@gw and @spex specified. Choose only one!", calledby='types_input%read_xml_input')
      END IF
      IF (numberNodesB.EQ.1) this%gw = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))))
      IF (numberNodesC.EQ.1) this%gw = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathC))))
      this%secvar = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@secvar'))
   END IF
   ! Read in Brillouin zone integration parameters
   IF (xml%GetNumberOfNodes('/fleurInput/cell/bzIntegration/@mode')> 0) THEN
      valueString = TRIM(ADJUSTL(xml%GetAttributeValue('/fleurInput/cell/bzIntegration/@mode')))
   ELSE
      valueString = TRIM(ADJUSTL(xml%GetAttributeValue('/fleurInput/calculationSetup/bzIntegration/@mode')))
   END IF
   SELECT CASE (valueString)
      CASE ('hist')
         this%bz_integration = BZINT_METHOD_HIST
      CASE ('gauss')
         this%bz_integration = BZINT_METHOD_GAUSS
      CASE ('tria')
         this%bz_integration = BZINT_METHOD_TRIA
      CASE ('tetra')
         this%bz_integration = BZINT_METHOD_TETRA
      CASE DEFAULT
         CALL juDFT_error('Invalid bzIntegration mode selected!')
   END SELECT
   nodeSum = 0
   if (xml%GetNumberOfNodes('/fleurInput/cell/bzIntegration')>0) THEN
     xpathb='/fleurInput/cell/bzIntegration'
   else
     xpathb='/fleurInput/calculationSetup/bzIntegration'
   endif
   xpathA=trim(xpathb)//'/@fermiSmearingEnergy'
   numberNodes = xml%GetNumberOfNodes(xPathA)
   nodeSum = nodeSum + numberNodes
   IF (numberNodes.EQ.1) THEN
      this%tkb = evaluateFirstOnly(xml%GetAttributeValue(xPathA))
   END IF
   xpathA=trim(xpathb)//'/@fermiSmearingTemp'
   numberNodes = xml%GetNumberOfNodes(xPathA)
   nodeSum = nodeSum + numberNodes
   IF (numberNodes.EQ.1) THEN
      this%tkb = evaluateFirstOnly(xml%GetAttributeValue(xPathA))
      this%tkb = boltzmann_Const * this%tkb
   END IF
   IF(nodeSum>1) THEN
      CALL juDFT_error('Error: Multiple fermi Smearing parameters provided in input file!')
   END IF
   xpathA=trim(xpathb)//'/@valenceElectrons'
   numberNodes = xml%GetNumberOfNodes(xPathA)
   IF (numberNodes.EQ.1) THEN
      this%zelec = evaluateFirstOnly(xml%GetAttributeValue(xPathA))
   ELSE
      CALL juDFT_error('Error: Optionality of valence electrons in input file not yet implemented!')
   END IF
   xPathA = trim(xpathb)//'/@l_bloechl'
   IF (xml%versionNumber > 31) this%l_bloechl = evaluateFirstBoolOnly(xml%GetAttributeValue(xPathA))

   this%film =  xml%GetNumberOfNodes('/fleurInput/cell/filmLattice')==1
   ! Read in optional geometry optimization parameters
   xPathA = '/fleurInput/calculationSetup/geometryOptimization'
   numberNodes = xml%GetNumberOfNodes(xPathA)
   IF (numberNodes.EQ.1) THEN
      this%l_f = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_f'))
      if (xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@f_level')>0) this%f_level = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@f_level'))
      this%forcealpha = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@forcealpha'))
      this%epsdisp = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@epsdisp'))
      this%epsforce = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@epsforce'))
      this%force_converged = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@force_converged'))
      this%qfix = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@qfix'))
      valueString=xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@forcemix')
      SELECT CASE(TRIM(valuestring))
         CASE ("straight")
            this%forcemix = 0
         CASE ("CG")
            this%forcemix = 1
         CASE ("BFGS")
            this%forcemix = 2
         CASE default
            CALL juDFT_error("Illegal mixing scheme for force in inp.xml:"//TRIM(valuestring))
      END SELECT
   END IF
   ! Read in optional general LDA+U parameters
   xPathA = '/fleurInput/calculationSetup/ldaU'
   numberNodes = xml%GetNumberOfNodes(xPathA)
   IF (numberNodes.EQ.1) THEN
      this%ldauLinMix = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_linMix'))
      this%ldauMixParam = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@mixParam'))
      this%ldauSpinf = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@spinf'))
      this%ldauAdjEnpara = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_adjEnpara'))
      IF(xml%versionNumber>=35) THEN
        this%ldauSpinoffd = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_spinoffd'))
      ENDIF
   END IF
   ! Read in RDMFT parameters
   xPathA = '/fleurInput/calculationSetup/rdmft'
   numberNodes = xml%GetNumberOfNodes(xPathA)
   IF (numberNodes.EQ.1) THEN
      this%l_rdmft = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_rdmft'))
      this%rdmftOccEps = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@occEps'))
      this%rdmftStatesBelow = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@statesBelow'))
      this%rdmftStatesAbove = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@statesAbove'))
      valueString = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@functional')))
      SELECT CASE (valueString)
      CASE ('Muller')
         this%rdmftFunctional = 1
      CASE DEFAULT
         STOP 'Error: unknown RDMFT functional selected!'
      END SELECT
   END IF
   ! !! Start of output section
   xPathA = '/fleurInput/output'
   numberNodes = xml%GetNumberOfNodes(xPathA)
   IF (numberNodes.EQ.1) THEN
      ! Read in general output switches
      this%l_coreSpec = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@coreSpec'))
      this%l_wann = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@wannier'))
      IF (xml%versionNumber > 31)this%eig66(1) = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@eig66'))
      ! Read in optional switches for checks
      xPathA = '/fleurInput/output/checks'
      numberNodes = xml%GetNumberOfNodes(xPathA)
      IF (numberNodes.EQ.1) THEN
         this%vchk = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@vchk'))
         this%cdinf = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@cdinf'))
      END IF
      ! Read in optional plotting parameters
      xPathA = '/fleurInput/output/plotting'
      numberNodes = xml%GetNumberOfNodes(xPathA)
      IF (numberNodes.EQ.1) THEN
         !this%score = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@score'))
      END IF
      ! Read in optional specialOutput switches
      xPathA = '/fleurInput/output/specialOutput'
      numberNodes = xml%GetNumberOfNodes(xPathA)
      IF (numberNodes.EQ.1) THEN
         this%eonly = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@eonly'))
         this%l_bmt = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@bmt'))
      END IF
      ! Read in optional vacuumDOS parameters
      xPathA = '/fleurInput/output/vacuumDOS'
      numberNodes = xml%GetNumberOfNodes(xPathA)
      IF (numberNodes.EQ.1) THEN
         this%integ = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@integ'))
      END IF
      ! Read in optional chargeDensitySlicing parameters
      xPathA = '/fleurInput/output/chargeDensitySlicing'
      numberNodes = xml%GetNumberOfNodes(xPathA)
      IF (numberNodes.EQ.1) THEN
         this%pallst = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@pallst'))
      END IF
   END IF
END SUBROUTINE read_xml_input

SUBROUTINE init_input(input,noco,l_hybrid,invs,n_denmat,n_hia,nbasfcn)
   USE m_types_noco
   CLASS(t_input),INTENT(inout):: input
   TYPE(t_noco),INTENT(in)     :: noco
   LOGICAL, INTENT(IN)         :: invs
   INTEGER, INTENT(IN)         :: n_denmat, n_hia
   LOGICAL, INTENT(in)         :: l_hybrid
   INTEGER,INTENT(IN),OPTIONAL :: nbasfcn

   ! Generate missing general parameters
   INTEGER :: minNeigd
   minNeigd = MAX(5,NINT(0.75*input%zelec) + 1)
   IF (noco%l_soc.AND.(.NOT.noco%l_noco)) minNeigd = 2 * minNeigd
   IF (noco%l_soc.AND.noco%l_ss) minNeigd=(3*minNeigd)/2
   IF ((input%neig.NE.-1).AND.(input%neig.LT.minNeigd)) THEN
      IF (input%neig>0) THEN
         WRITE(*,*) 'numbands is too small. Setting parameter to default value.'
         WRITE(*,*) 'changed numbands (input%neig) to ',minNeigd
      ENDIF
      input%neig = minNeigd
   END IF
   IF(input%neig == -1. .AND. PRESENT(nbasfcn)) THEN
      input%neig = nbasfcn
   END IF
   IF (noco%l_noco) input%neig = 2*input%neig
   input%gw_neigd = MERGE(MAX(NINT(input%zelec)*10, 60),0, l_hybrid)

   IF(PRESENT(nbasfcn)) input%neig = MIN(input%neig, nbasfcn)

   input%l_real = invs.and..not.noco%l_noco.and..not.(noco%l_soc.and.n_denmat>0).and..not.n_hia>0
END SUBROUTINE init_input

END MODULE m_types_input
