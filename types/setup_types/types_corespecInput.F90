!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_corespecInput
  USE m_judft
  USE m_types_fleur_setup
  USE m_json_tools

  IMPLICIT NONE
  TYPE,EXTENDS(t_fleursetup):: t_corespecInput
     INTEGER :: verb  ! output verbosity
     INTEGER :: atomType  ! atomic type used for calculation of core spectra
     CHARACTER(LEN=1) :: edge  ! edge character (K,L,M,N,O,P)
     INTEGER :: edgeidx(11)  ! l-j edges
     INTEGER :: lx  ! maximum lmax considered in spectra calculation
     REAL :: ek0  ! kinetic energy of incoming electrons
     REAL :: emn  ! energy spectrum lower bound
     REAL :: emx  ! energy spectrum upper bound
     REAL :: ein  ! energy spectrum increment
   CONTAINS
     PROCEDURE,PASS :: broadcast=>broadcast_corespecInput
     PROCEDURE,PASS :: WRITE=>WRITE_corespecInput
     PROCEDURE,PASS :: READ=>READ_corespecInput
     PROCEDURE,PASS :: read_xml=>read_xml_corespecInput
  END TYPE t_corespecInput

CONTAINS
  SUBROUTINE broadcast_corespecInput(tt,mpi_comm,origin)
#ifdef CPP_MPI
    USE m_mpi_bc_tool
#endif    
    IMPLICIT NONE
    CLASS(t_corespecInput),INTENT(INOUT):: tt
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
    CALL MPI_BCAST(tt%verb,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%atomType,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%edgeidx,11,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%lx,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%ek0,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%emn,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%emx,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%ein,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%edge,1,MPI_CHARACTER,pe,mpi_comm,ierr)
#endif

  END SUBROUTINE broadcast_corespecInput

  SUBROUTINE write_corespecInput(tt, unit, iotype, v_list, iostat, iomsg)
    IMPLICIT NONE
    CLASS(t_corespecInput),INTENT(IN):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    INTEGER:: ntype
    REAL,ALLOCATABLE::d_wgn(:,:)

    WRITE(unit,*,IOSTAT=iostat) '"corespecInput":{'

    CALL json_print(unit,"verb",tt%verb)
    CALL json_print(unit,"atomType",tt%atomType)
    CALL json_print(unit,"edgeidx",tt%edgeidx)
    CALL json_print(unit,"lx",tt%lx)
    CALL json_print(unit,"ek0",tt%ek0)
    CALL json_print(unit,"emn",tt%emn)
    CALL json_print(unit,"emx",tt%emx)
    CALL json_print(unit,"ein",tt%ein)
    CALL json_print(unit,"edge",tt%edge)

    WRITE(unit,*,IOSTAT=iostat) '}'

  END SUBROUTINE write_corespecInput
  SUBROUTINE read_corespecInput(tt, unit, iotype, v_list, iostat, iomsg)
    USE m_json_tools
    IMPLICIT NONE
    CLASS(t_corespecInput),INTENT(INOUT):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    INTEGER :: ntype
    REAL,ALLOCATABLE::rt(:,:)
    CALL json_open_class("corespecInput",unit,iostat)
    IF (iostat.NE.0)   RETURN

    CALL json_read(unit,"verb",tt%verb)
    CALL json_read(unit,"atomType",tt%atomType)
    CALL json_read(unit,"edgeidx",tt%edgeidx)
    CALL json_read(unit,"lx",tt%lx)
    CALL json_read(unit,"ek0",tt%ek0)
    CALL json_read(unit,"emn",tt%emn)
    CALL json_read(unit,"emx",tt%emx)
    CALL json_read(unit,"ein",tt%ein)
    CALL json_read(unit,"edge",tt%edge)


    CALL json_close_class(unit,iostat)

  END SUBROUTINE read_corespecInput


  SUBROUTINE read_xml_corespecInput(tt)
    USE m_xmlIntWrapFort
    USE m_constants
    USE m_calculator
    USE m_inp_xml
    IMPLICIT NONE
    CLASS(t_corespecInput),INTENT(OUT):: tt

    LOGICAL :: tempBool
    INTEGER:: numberNodes,numtokens
    CHARACTER(len=50):: xpathA,xpathb,valuestring
    
    xPathA = '/fleurInput/output/coreSpectrum'
    numberNodes = xmlGetNumberOfNodes(xPathA)

     IF (numberNodes.EQ.1) THEN
        tt%verb = 0
        tempBool = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@verbose'))
        IF(tempBool) tt%verb = 1
        tt%ek0 = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@eKin'))
        tt%atomType = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@atomType'))
        tt%lx = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@lmax'))
        tt%edge = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@edgeType')))
        tt%emn = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@eMin'))
        tt%emx = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@eMax'))
        tempInt = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@numPoints'))
        tt%ein = (tt%emx - tt%emn) / (tempInt - 1.0)
        xPathB = TRIM(ADJUSTL(xPathA))//'/edgeIndices'
        xPathB = TRIM(ADJUSTL(xPathB))//'/text()'
        valueString = xmlGetAttributeValue(TRIM(ADJUSTL(xPathB)))
        numTokens = countStringTokens(valueString)
        tt%edgeidx(:) = 0
        IF(numTokens.GT.SIZE(tt%edgeidx)) THEN
           CALL juDFT_error('More EELS edge indices provided than allowed.',calledby='r_inpXML')
        END IF
        DO i = 1, MAX(numTokens,SIZE(tt%edgeidx))
           tt%edgeidx(i) = evaluateFirstIntOnly(popFirstStringToken(valueString))
        END DO
     END IF


  END SUBROUTINE read_xml_corespecInput



END MODULE m_types_corespecInput
