!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_vacuum
  USE m_judft
  USE m_types_fleur_setup
  use m_json_tools
  IMPLICIT NONE
  TYPE,EXTENDS(t_fleursetup):: t_vacuum
     INTEGER ::nmz
      INTEGER ::nmzxy
     INTEGER :: layers
     INTEGER :: nvac
     REAL :: delz
     REAL :: dvac
     INTEGER::nstars
     INTEGER:: nstm
     REAL :: tworkf
     REAL :: locx(2)
     REAL :: locy(2)
     LOGICAL ::starcoeff
     INTEGER, ALLOCATABLE :: izlay(:,:)
   CONTAINS
     PROCEDURE,PASS :: broadcast=>broadcast_vacuum
     PROCEDURE,PASS :: write=>WRITE_vacuum
     PROCEDURE,PASS :: read=>READ_vacuum
     PROCEDURE,PASS :: read_xml=>read_xml_vacuum
  END TYPE t_vacuum

CONTAINS
  SUBROUTINE broadcast_vacuum(tt,mpi_comm,origin)
    IMPLICIT NONE
    CLASS(t_vacuum),INTENT(INOUT):: tt
    INTEGER,INTENT(IN)               :: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL      :: origin

#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER :: pe,ierr,irank,n(2)
    CALL MPI_COMM_RANK(mpi_comm,irank,ierr)
    IF (PRESENT(origin)) THEN
       pe=origin
    ELSE
       pe=0
    ENDIF

    CALL MPI_BCAST(tt%nmz,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%nmzxy,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%layers,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%nvac,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%nstars,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%nstm,1,MPI_INTEGER,pe,mpi_comm,ierr)
    
    CALL MPI_BCAST(tt%delz,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%dvac,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%tworkf,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%locx,2,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%locy,2,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)

    CALL MPI_BCAST(tt%starcoeff,1,MPI_LOGICAL,pe,mpi_comm,ierr)

    IF (irank==pe) THEN
       n=SHAPE(tt%izlay)
       CALL MPI_BCAST(n,2,MPI_INTEGER,pe,mpi_comm,ierr)
    ELSE
       CALL MPI_BCAST(n,2,MPI_INTEGER,pe,mpi_comm,ierr)
       IF (ALLOCATED(tt%izlay)) DEALLOCATE(tt%izlay)
       ALLOCATE(tt%izlay(n(1),n(2)))
    ENDIF
    CALL MPI_BCAST(tt%izlay,SIZE(tt%izlay),MPI_INTEGER,pe,mpi_comm,ierr)
#endif
      
  END SUBROUTINE broadcast_vacuum

  SUBROUTINE write_vacuum(tt, unit, iotype, v_list, iostat, iomsg)
    IMPLICIT NONE
    CLASS(t_vacuum),INTENT(IN):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    WRITE(unit,*,IOSTAT=iostat) '"vacuum":{'


    CALL json_print(unit,"nmz",tt%nmz)
    CALL JSON_PRINT(unit,"nmzxy",tt%nmzxy)
    CALL JSON_PRINT(unit,"layers",tt%layers) 
    CALL JSON_PRINT(unit, "nvac",tt%nvac)
    CALL JSON_PRINT(unit, "delz",tt%delz)  
    CALL JSON_PRINT(unit, "dvac",tt%dvac)  
    CALL JSON_PRINT(unit, "nstars",tt%nstars)  
    CALL JSON_PRINT(unit, "nstm",tt%nstm)  
    CALL JSON_PRINT(unit, "tworkf",tt%tworkf)  
    CALL JSON_PRINT(unit, "locx",tt%locx)  
    CALL JSON_PRINT(unit, "locy",tt%locy)  
    CALL JSON_PRINT(unit, "starcoeff",tt%starcoeff)
    CALL JSON_PRINT(unit, "izlay_shape",SHAPE(tt%izlay))
    CALL JSON_PRINT(unit,"izlay",RESHAPE(tt%izlay,(/SIZE(tt%izlay)/)))
    
    WRITE(unit,*,IOSTAT=iostat) '}'
    
  END SUBROUTINE write_vacuum
  SUBROUTINE read_vacuum(tt, unit, iotype, v_list, iostat, iomsg)
    IMPLICIT NONE
    CLASS(t_vacuum),INTENT(INOUT):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    CHARACTER(len=40)::string
    INTEGER,ALLOCATABLE:: itemp(:),n(2)
    CALL json_open_class("vacuum",unit,iostat)
    IF (iostat.NE.0)   RETURN

    CALL json_read(unit,"nmz",tt%nmz)
    CALL JSON_READ(unit,"nmzxy",tt%nmzxy)
    CALL JSON_READ(unit,"layers",tt%layers) 
    CALL JSON_READ(unit, "nvac",tt%nvac)
    CALL JSON_READ(unit, "delz",tt%delz)  
    CALL JSON_READ(unit, "dvac",tt%dvac)  
    CALL JSON_READ(unit, "nstars",tt%nstars)  
    CALL JSON_READ(unit, "nstm",tt%nstm)  
    CALL JSON_READ(unit, "tworkf",tt%tworkf)  
    CALL JSON_READ(unit, "locx",tt%locx)  
    CALL JSON_READ(unit, "locy",tt%locy)  
    CALL JSON_READ(unit, "starcoeff",tt%starcoeff)
    CALL JSON_READ(unit, "izlay_shape",n)

    IF (ALLOCATED(tt%izlay)) DEALLOCATE(tt%izlay)
    ALLOCATE(tt%izlay(n(1),n(2)),itemp(n(1)*n(2)))
    CALL JSON_READ(unit,"izlay",itemp)
    tt%izlay=RESHAPE(itemp,n)

    CALL json_close_class(unit)
    
  END SUBROUTINE read_vacuum

 


  SUBROUTINE read_xml_vacuum(vacuum)
    USE m_xmlIntWrapFort
    USE m_constants
    USE m_calculator
    USE m_invs3
    IMPLICIT NONE
    CLASS(t_vacuum),INTENT(OUT):: vacuum


    LOGICAL::film,vacdos
    CHARACTER(len=200):: xpathA,valueString


    film=(xmlGetNumberOfNodes('/fleurInput/vacuum/filmLattice')==1)
    vacuum%nvac = 2

    vacuum%nmz = 250
    vacuum%delz = 25.0/vacuum%nmz
    vacuum%nmzxy = 100

    !DOS stuff for vacuum
    vacuum%layers = 1
    vacuum%starcoeff = .FALSE.
    vacuum%nstars = 0
    vacuum%locx = 0.0
    vacuum%locy = 0.0
    vacuum%nstm = 0
    vacuum%tworkf = 0.0
    !Check if vacdos=t
    vacdos=.FALSE.
    IF (xmlGetNumberOfNodes('/fleurInput/output')==1) &
         vacdos = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/output/@vacdos'))
    
    xPathA = '/fleurInput/output/vacuumDOS'
    IF (film.AND.vacdos.AND.xmlGetNumberOfNodes(xPathA).EQ.1) THEN
       vacuum%layers = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@layers'))
       vacuum%starcoeff = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@star'))
       vacuum%nstars = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nstars'))
       vacuum%locx(1) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locx1'))
       vacuum%locx(2) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locx2'))
       vacuum%locy(1) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locy1'))
       vacuum%locy(2) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locy2'))
       vacuum%nstm = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nstm'))
       vacuum%tworkf = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@tworkf'))
    END IF

    IF (vacuum%layers>1) THEN
       ALLOCATE(vacuum%izlay(vacuum%layers,2))
       CALL judft_error("Layer mode not supported in types_vacuum")
    ELSE
       ALLOCATE(vacuum%izlay(1,1))
       vacuum%izlay=1.0
    END IF
    
  END SUBROUTINE read_xml_vacuum


  SUBROUTINE init_vacuum(vacuum,sym)
    USE types_sym
    IMPLICIT NONE
    CLASS(t_vacuum),INTENT(OUT):: vacuum
    TYPE(t_sym),INTENT(IN)     :: sym


    IF (sym%zrfs.OR.sym%invs) vacuum%nvac = 1

  END SUBROUTINE init_vacuum

    
END MODULE m_types_vacuum
