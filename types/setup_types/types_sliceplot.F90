!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_sliceplot
  USE m_judft
  USE m_types_fleur_setup
  IMPLICIT NONE
  TYPE,EXTENDS(t_fleursetup):: t_sliceplot
     LOGICAL :: iplot
     LOGICAL :: slice
     LOGICAL :: plpot
     LOGICAL :: pallst
     INTEGER :: kk
     INTEGER :: nnne
     REAL    :: e1s
     REAL    :: e2s

   CONTAINS
     PROCEDURE,PASS :: broadcast=>broadcast_sliceplot
     PROCEDURE,PASS :: WRITE=>WRITE_sliceplot
     PROCEDURE,PASS :: READ=>READ_sliceplot
     PROCEDURE,PASS :: read_xml=>read_xml_sliceplot
  END TYPE t_sliceplot

CONTAINS
  SUBROUTINE broadcast_sliceplot(tt,mpi_comm,origin)
#ifdef CPP_MPI
    USE m_mpi_bc_tool
#endif    
    IMPLICIT NONE
    CLASS(t_sliceplot),INTENT(INOUT):: tt
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

    CALL MPI_BCAST(tt%iplot,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%slice,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%plpot,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%pallst,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%kk,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%nnne,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%e1s,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%e2s,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
#endif

  END SUBROUTINE broadcast_sliceplot

  SUBROUTINE write_sliceplot(tt, unit, iotype, v_list, iostat, iomsg)
    IMPLICIT NONE
    CLASS(t_sliceplot),INTENT(IN):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    INTEGER:: ntype
    REAL,ALLOCATABLE::d_wgn(:,:)

    WRITE(unit,*,IOSTAT=iostat) '"sliceplot":{'

    CALL json_print(unit,"iplot",tt%iplot)
    CALL json_print(unit,"slice",tt%slice)
    CALL json_print(unit,"plpot",tt%plpot)
    CALL json_print(unit,"pallst",tt%pallst)
    CALL json_print(unit,"kk",tt%kk)
    CALL json_print(unit,"nnne",tt%nnne)
    CALL json_print(unit,"e1s",tt%e1s)
    CALL json_print(unit,"e2s",tt%e2s)

    WRITE(unit,*,IOSTAT=iostat) '}'

  END SUBROUTINE write_sliceplot
  SUBROUTINE read_sliceplot(tt, unit, iotype, v_list, iostat, iomsg)
    USE m_json_tools
    IMPLICIT NONE
    CLASS(t_sliceplot),INTENT(INOUT):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    INTEGER :: ntype
    REAL,ALLOCATABLE::rt(:,:)
    CALL json_open_class("sliceplot",unit,iostat)
    IF (iostat.NE.0)   RETURN


    CALL json_read(unit,"iplot",tt%iplot)
    CALL json_read(unit,"slice",tt%slice)
    CALL json_read(unit,"plpot",tt%plpot)
    CALL json_read(unit,"pallst",tt%pallst)
    CALL json_read(unit,"kk",tt%kk)
    CALL json_read(unit,"nnne",tt%nnne)
    CALL json_read(unit,"e1s",tt%e1s)
    CALL json_read(unit,"e2s",tt%e2s)


    CALL json_close_class(unit,iostat)

  END SUBROUTINE read_sliceplot


  SUBROUTINE read_xml_sliceplot(tt)
    USE m_xmlIntWrapFort
    USE m_constants
    USE m_calculator
    USE m_inp_xml
    IMPLICIT NONE
    CLASS(t_sliceplot),INTENT(OUT):: tt

    sliceplot%slice = .FALSE.
    sliceplot%iplot = .FALSE.
    sliceplot%plpot = .FALSE.
    IF(xmlGetNumberOfNodes('/fleurInput/output/')==1)&
         sliceplot%slice = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/output/@slice'))
    IF(xmlGetNumberOfNodes('/fleurInput/output/plotting')==1) THEN
       sliceplot%iplot = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/output/plotting/@iplot'))
       sliceplot%plplot = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/output/plotting/@plplot'))
    END IF


    IF ((sliceplot%slice).AND.(xmlGetNumberOfNodes('/fleurInput/output/chargeDensitySlicing').EQ.0)) THEN
       CALL juDFT_error("slice is true but chargeDensitySlicing parameters are not set!", calledby = "r_inpXML")
    END IF

    IF (xmlGetNumberOfNodes('/fleurInput/output/chargeDensitySlicing').EQ.1) THEN
       sliceplot%kk = evaluateFirstIntOnly(xmlGetAttributeValue('/fleurInput/output/chargeDensitySlicing/@numkpt'))
       sliceplot%e1s = evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/output/chargeDensitySlicing/@minEigenval'))
       sliceplot%e2s = evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/output/chargeDensitySlicing/@maxEigenval'))
       sliceplot%nnne = evaluateFirstIntOnly(xmlGetAttributeValue('/fleurInput/output/chargeDensitySlicing/@nnne'))
       sliceplot%pallst = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/output/chargeDensitySlicing/@pallst'))
    END IF


  END SUBROUTINE read_xml_sliceplot



END MODULE m_types_sliceplot
