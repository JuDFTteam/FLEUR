!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_input
  USE m_judft
  USE m_types_fleur_setup
  IMPLICIT NONE
  TYPE,EXTENDS(t_fleursetup):: t_input
     LOGICAL :: cdinf
     LOGICAL :: vchk
     LOGICAL :: l_f
     LOGICAL :: film
     LOGICAL :: ctail
     INTEGER :: coretail_lmax
     INTEGER :: qfix
     INTEGER :: jspins
     LOGICAL :: l_wann
     LOGICAL :: l_inpXML
     REAL    :: ellow
     REAL    :: elup
     REAL    :: rkmax
     REAL    :: zelec
     LOGICAL :: l_useapw  
     REAL,POINTER :: sigma !this is the difference in charge due to the electric field it points to the value stored in t_efield
     !Stuff for mixing
     INTEGER :: maxiter
     INTEGER :: imix
     REAL    :: preconditioning_param
     REAL    :: spinf
     LOGICAL :: ldauLinMix
     REAL    :: ldauMixParam
     REAL    :: ldauSpinf
     !Force stuff
     REAL    :: xa !< mixing parameter for geometry optimzer
     REAL    :: thetad !< Debey temperature for first step of geometry optimzer
     REAL    :: epsdisp !< minimal displacement. If all displacements are < epsdisp stop
     REAL    :: epsforce !< minimal force. If all forces <epsforce stop
     !Fermi
     REAL    :: delgau
     REAL    :: alpha
     REAL    :: tkb
     LOGICAL :: gauss
     LOGICAL :: tria
     REAL :: fixed_moment=0.0
     !Core 
     INTEGER:: kcrel
     LOGICAL:: frcor
     LOGICAL:: l_coreSpec
     LOGICAL:: l_core_confpot
     !Comment (not written to files, no broadcast)
     CHARACTER(len=:),ALLOCATABLE:: comment
   CONTAINS
     PROCEDURE,PASS :: broadcast=>broadcast_input
     PROCEDURE,PASS :: write=>WRITE_input
     PROCEDURE,PASS :: read=>READ_input
     PROCEDURE,PASS :: read_xml=>read_xml_input
     !PROCEDURE,PASS :: init=>init_input
  END TYPE t_input

CONTAINS
  SUBROUTINE broadcast_input(tt,mpi_comm,origin)
#ifdef CPP_MPI
    USE m_mpi_bc_tool
#endif    
    IMPLICIT NONE
    CLASS(t_input),INTENT(INOUT):: tt
    INTEGER,INTENT(IN)               :: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL      :: origin
    
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER :: pe,ierr,ntype,irank

    CALL MPI_COMM_RANK(mpi_comm,irank,ierr)
    
    IF (PRESENT(origin)) THEN
       pe=origin
    ELSE
       pe=0
    ENDIF

    CALL MPI_BCAST(tt%cdninf,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%vchk,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%lf,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%film,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%ctail,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%coretail_lmax,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%qfix,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%jspins,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%l_wann,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%l_inpXML,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%ellow,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%elup,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%rkmax,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%zelec,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%l_useapw,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    !Mixing
    CALL MPI_BCAST(tt%maxiter,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%imix,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%preconditioning_param,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%spinf,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%ldaulinmix,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%ldaumixparam,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%ldauspinf,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    !Force
    CALL MPI_BCAST(tt%xa,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%thetad,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%epsdisp,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%epsforce,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    !core
    CALL MPI_BCAST(tt%kcrel,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%frcor,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%l_corespec,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%l_core_confpot,1,MPI_LOGICAL,pe,mpi_comm,ierr)
#endif
      
  END SUBROUTINE broadcast_input

  SUBROUTINE write_input(tt, unit, iotype, v_list, iostat, iomsg)
    IMPLICIT NONE
    CLASS(t_input),INTENT(IN):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

   
    WRITE(unit,*,IOSTAT=iostat) '"input":{'
    CALL json_print(unit,"cdinf",tt%cdinf)
    CALL json_print(unit,"vchk",tt%vchk)
    CALL json_print(unit,"l_f",tt%l_f)
    CALL json_print(unit,"film",tt%film)
    CALL json_print(unit,"ctail",tt%ctail)
    CALL json_print(unit,"coretail_lmax",tt%coretail_lmax)
    CALL json_print(unit,"qfix",tt%qfix)
    CALL json_print(unit,"jspins",tt%jspins)
    CALL json_print(unit,"l_wann",tt%l_wann)
    CALL json_print(unit,"l_inpXML",tt%l_inpXML)
    CALL json_print(unit,"ellow",tt%ellow)
    CALL json_print(unit,"elup",tt%elup)
    CALL json_print(unit,"rkmax",tt%rkmax)
    CALL json_print(unit,"zelec",tt%zelec)
    CALL json_print(unit,"l_useapw  ",tt%l_useapw  )
    !Stuff for mixing
    CALL json_print(unit,"maxiter",tt%maxiter)
    CALL json_print(unit,"imix",tt%imix)
    CALL json_print(unit,"preconditioning_param",tt%preconditioning_param)
    CALL json_print(unit,"spinf",tt%spinf)
    CALL json_print(unit,"ldauLinMix",tt%ldauLinMix)
    CALL json_print(unit,"ldauMixParam",tt%ldauMixParam)
    CALL json_print(unit,"ldauSpinf",tt%ldauSpinf)
    !Force stuff
    CALL json_print(unit,"xa ",tt%xa )
    CALL json_print(unit,"thetad ",tt%thetad )
    CALL json_print(unit,"epsdisp ",tt%epsdisp )
    CALL json_print(unit,"epsforce",tt%epsforce)
    !Fermi
    CALL json_print(unit,"delgau",tt%delgau)
    CALL json_print(unit,"alpha",tt%alpha)
    CALL json_print(unit,"tkb",tt%tkb)
    CALL json_print(unit,"gauss",tt%gauss)
    CALL json_print(unit,"tria",tt%tria)
    CALL json_print(unit,"fixed_moment",tt%fixed_moment)
    !Core 
    CALL json_print(unit,"kcrel",tt%kcrel)
    CALL json_print(unit,"frcor",tt%frcor)
    CALL json_print(unit,"l_coreSpec",tt%l_coreSpec)
    CALL json_print(unit,"l_core_confpot",tt%l_core_confpot)
    WRITE(unit,*,IOSTAT=iostat) '}'
    
  END SUBROUTINE write_input
  SUBROUTINE read_input(tt, unit, iotype, v_list, iostat, iomsg)
    use m_json_tools
    IMPLICIT NONE
    CLASS(t_input),INTENT(INOUT):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    INTEGER :: ntype
    REAL,ALLOCATABLE::rt(:,:)
    CALL json_open_class("input",unit,iostat)
    IF (iostat.NE.0)   RETURN

    CALL json_read(unit,"cdinf",tt%cdinf)
    CALL json_read(unit,"vchk",tt%vchk)
    CALL json_read(unit,"l_f",tt%l_f)
    CALL json_read(unit,"film",tt%film)
    CALL json_read(unit,"ctail",tt%ctail)
    CALL json_read(unit,"coretail_lmax",tt%coretail_lmax)
    CALL json_read(unit,"qfix",tt%qfix)
    CALL json_read(unit,"jspins",tt%jspins)
    CALL json_read(unit,"l_wann",tt%l_wann)
    CALL json_read(unit,"l_inpXML",tt%l_inpXML)
    CALL json_read(unit,"ellow",tt%ellow)
    CALL json_read(unit,"elup",tt%elup)
    CALL json_read(unit,"rkmax",tt%rkmax)
    CALL json_read(unit,"zelec",tt%zelec)
    CALL json_read(unit,"l_useapw  ",tt%l_useapw  )
    !Stuff for mixing
    CALL json_read(unit,"maxiter",tt%maxiter)
    CALL json_read(unit,"imix",tt%imix)
    CALL json_read(unit,"preconditioning_param",tt%preconditioning_param)
    CALL json_read(unit,"spinf",tt%spinf)
    CALL json_read(unit,"ldauLinMix",tt%ldauLinMix)
    CALL json_read(unit,"ldauMixParam",tt%ldauMixParam)
    CALL json_read(unit,"ldauSpinf",tt%ldauSpinf)
    !Force stuff
    CALL json_read(unit,"xa ",tt%xa )
    CALL json_read(unit,"thetad ",tt%thetad )
    CALL json_read(unit,"epsdisp ",tt%epsdisp )
    CALL json_read(unit,"epsforce",tt%epsforce)
    !Fermi
    CALL json_read(unit,"delgau",tt%delgau)
    CALL json_read(unit,"alpha",tt%alpha)
    CALL json_read(unit,"tkb",tt%tkb)
    CALL json_read(unit,"gauss",tt%gauss)
    CALL json_read(unit,"tria",tt%tria)
    CALL json_read(unit,"fixed_moment",tt%fixed_moment)
    !Core 
    CALL json_read(unit,"kcrel",tt%kcrel)
    CALL json_read(unit,"frcor",tt%frcor)
    CALL json_read(unit,"l_coreSpec",tt%l_coreSpec)
    CALL json_read(unit,"l_core_confpot",tt%l_core_confpot)
    

    CALL json_close_class(unit,iostat)
    
  END SUBROUTINE read_input

 
  SUBROUTINE read_xml_input(tt)
    USE m_xmlIntWrapFort
    USE m_constants
    USE m_calculator
    USE m_inp_xml
    IMPLICIT NONE
    CLASS(t_input),INTENT(OUT):: tt


    INTEGER           :: i
    CHARACTER(len=200):: xpath,valueString,xpatha
    LOGICAL :: l_qfix
    
    INTEGER           :: numberNodes,nodesum

    !TODO! these switches should be in the inp-file
    tt%l_core_confpot=.TRUE. !former CPP_CORE
    tt%l_useapw=.FALSE.   !former CPP_APW
    ! Read general cutoff parameters
    tt%rkmax = evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/cutoffs/@Kmax'))
    tt%maxiter = evaluateFirstIntOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@maxIterBroyd'))
    valueString = TRIM(ADJUSTL(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@imix')))
    SELECT CASE (valueString)
    CASE ('straight')
       tt%imix = 0
    CASE ('Broyden1')
       tt%imix = 3
    CASE ('Broyden2')
       tt%imix = 5
    CASE ('Anderson')
       tt%imix = 7
    CASE DEFAULT
       STOP 'Error: unknown mixing scheme selected!'
    END SELECT
    tt%alpha = evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@alpha'))
    tt%preconditioning_param = evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@preconditioning_param'))
    tt%spinf = evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/scfLoop/@spinf'))

    ! Get parameters for core electrons
    tt%ctail = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/coreElectrons/@ctail'))
    tt%coretail_lmax = evaluateFirstIntOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/coreElectrons/@coretail_lmax'))
    tt%frcor = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/coreElectrons/@frcor'))
    tt%kcrel = evaluateFirstIntOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/coreElectrons/@kcrel'))

    
    tt%jspins = evaluateFirstIntOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/magnetism/@jspins'))
    tt%fixed_moment=evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/magnetism/@fixed_moment'))
    valueString = TRIM(ADJUSTL(xmlGetAttributeValue('/fleurInput/calculationSetup/bzIntegration/@mode')))
    SELECT CASE (valueString)
    CASE ('hist')
       tt%gauss = .FALSE.
       tt%tria = .FALSE.
    CASE ('gauss')
       tt%gauss = .TRUE.
       tt%tria = .FALSE.
    CASE ('tria')
       tt%gauss = .FALSE.
       tt%tria = .TRUE.
    CASE DEFAULT
       STOP 'Invalid bzIntegration mode selected!'
    END SELECT
    nodeSum = 0
    xPathA = '/fleurInput/calculationSetup/bzIntegration/@fermiSmearingEnergy'
    numberNodes = xmlGetNumberOfNodes(xPathA)
    nodeSum = nodeSum + numberNodes
    IF (numberNodes.EQ.1) THEN
       tt%tkb = evaluateFirstOnly(xmlGetAttributeValue(xPathA))
    END IF
    xPathA = '/fleurInput/calculationSetup/bzIntegration/@fermiSmearingTemp'
    numberNodes = xmlGetNumberOfNodes(xPathA)
    nodeSum = nodeSum + numberNodes
    IF (numberNodes.EQ.1) THEN
       tt%tkb = evaluateFirstOnly(xmlGetAttributeValue(xPathA))
       tt%tkb = boltzmannConst * tt%tkb
    END IF
    IF(nodeSum.GE.2) THEN
       STOP 'Error: Multiple fermi Smearing parameters provided in input file!'
    END IF

    xPathA = '/fleurInput/calculationSetup/bzIntegration/@valenceElectrons'
    numberNodes = xmlGetNumberOfNodes(xPathA)
    IF (numberNodes.EQ.1) THEN
       tt%zelec = evaluateFirstOnly(xmlGetAttributeValue(xPathA))
    ELSE
       STOP 'Error: Optionality of valence electrons in input file not yet implemented!'
    END IF
    

    xPathA = '/fleurInput/calculationSetup/geometryOptimization'
    numberNodes = xmlGetNumberOfNodes(xPathA)
    
    tt%l_f = .FALSE.
    tt%qfix = 0
    
    IF (numberNodes.EQ.1) THEN
       tt%l_f = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_f'))
       tt%xa = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@xa'))
       tt%thetad = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@thetad'))
       tt%epsdisp = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@epsdisp'))
       tt%epsforce = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@epsforce'))
       
       numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@qfix')
       IF (numberNodes.EQ.1) THEN
          tt%qfix = 1
          l_qfix = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@qfix'))
          IF (l_qfix) THEN
             tt%qfix = 2
          END IF
       END IF
    END IF

    xPathA = '/fleurInput/calculationSetup/ldaU'
    numberNodes = xmlGetNumberOfNodes(xPathA)
    IF (numberNodes.EQ.1) THEN
       tt%ldauLinMix = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_linMix'))
       tt%ldauMixParam = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@mixParam'))
       tt%ldauSpinf = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@spinf'))
    END IF
    


  ! Read in optional energy parameter limits

    xPathA = '/fleurInput/calculationSetup/energyParameterLimits'
    numberNodes = xmlGetNumberOfNodes(xPathA)
    
    IF (numberNodes.EQ.1) THEN
       tt%ellow = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@ellow'))
       tt%elup = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@elup'))
    END IF
    
    
    numberNodes = xmlGetNumberOfNodes('/fleurInput/cell/filmLattice')
    
    IF (numberNodes.EQ.1) THEN
       tt%film = .TRUE.
    END IF
    

    tt%l_coreSpec = .FALSE.
    tt%l_wann = .FALSE.
    
    tt%vchk = .FALSE.
    tt%cdinf = .FALSE.
    
    
    
    xPathA = '/fleurInput/output'
    numberNodes = xmlGetNumberOfNodes(xPathA)
    
    IF (numberNodes.EQ.1) THEN
       
       ! Read in general output switches
       tt%l_coreSpec = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@coreSpec'))
       tt%l_wann = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@wannier'))
       
       ! Read in optional switches for checks
       
       xPathA = '/fleurInput/output/checks'
       numberNodes = xmlGetNumberOfNodes(xPathA)
       
       IF (numberNodes.EQ.1) THEN
          tt%vchk = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@vchk'))
          tt%cdinf = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@cdinf'))
          
       END IF
       
    ENDIF
    xPathA = '/fleurInput/comment'
    tt%comment = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA)))))
    DO i = 1, LEN(tt%comment)
       IF (tt%comment(i:i).EQ.ACHAR(10)) tt%comment(i:i) = ' ' !remove line breaks
    END DO
  END SUBROUTINE read_xml_input
  
 

  
END MODULE m_types_input
