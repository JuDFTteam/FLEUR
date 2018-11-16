!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_cell
  USE m_judft
  USE m_types_fleur_setup
  use m_json_tools
  IMPLICIT NONE
  TYPE,EXTENDS(t_fleursetup):: t_cell
     !vol of dtilde box
     REAL::omtil
     !2D area
     REAL::area
     !bravais matrix
     REAL::amat(3,3)
     !rez. bravais matrx
     REAL::bmat(3,3)
     !square of bbmat
     REAL::bbmat(3,3)
     !d-value
     REAL::z1
     !volume of cell
     REAL::vol
     !volume of interstitial
     REAL::volint
   CONTAINS
     PROCEDURE,PASS :: broadcast=>broadcast_cell
     PROCEDURE,PASS :: write=>WRITE_cell
     PROCEDURE,PASS :: read=>READ_cell
     PROCEDURE,PASS :: read_xml=>read_xml_cell
  END TYPE t_cell

CONTAINS
  SUBROUTINE broadcast_cell(tt,mpi_comm,origin)
    IMPLICIT NONE
    CLASS(t_cell),INTENT(INOUT):: tt
    INTEGER,INTENT(IN)               :: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL      :: origin

#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER :: pe,ierr
    
    IF (PRESENT(origin)) THEN
       pe=origin
    ELSE
       pe=0
    ENDIF

    CALL MPI_BCAST(tt%omtil,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%area,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%z1,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%vol,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%volint,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%c,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)


    CALL MPI_BCAST(tt%amat,9,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%bmat,9,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%bbmat,9,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
#endif
      
  END SUBROUTINE broadcast_cell

  SUBROUTINE write_cell(tt, unit, iotype, v_list, iostat, iomsg)
    IMPLICIT NONE
    CLASS(t_cell),INTENT(IN):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    WRITE(unit,*,IOSTAT=iostat) '"cell":{'


    call json_print(unit,"omtil",tt%omtil)
    CALL JSON_PRINT(unit,"area",tt%area)
    CALL JSON_PRINT(unit,"z1",tt%z1) 
    CALL JSON_PRINT(unit, "vol",tt%vol)
    CALL JSON_PRINT(unit, "volint",tt%volint)  
  
    CALL JSON_PRINT(unit,"amat",RESHAPE(tt%amat,(/9/)))
    CALL JSON_PRINT(unit,"bmat",reshape(tt%bmat,(/9/)))
    CALL JSON_PRINT(unit,"bbmat",reshape(tt%bbmat,(/9/)))
    
    WRITE(unit,*,IOSTAT=iostat) '}'
    
  END SUBROUTINE write_cell
  SUBROUTINE read_cell(tt, unit, iotype, v_list, iostat, iomsg)
    IMPLICIT NONE
    CLASS(t_cell),INTENT(INOUT):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    CHARACTER(len=40)::string
    real:: rtemp(9)
    CALL json_open_class("cell",unit,iostat)
    IF (iostat.NE.0)   RETURN

    call json_read(unit,"omtil",tt%omtil)
    call json_read(unit,"area",tt%area)
    call json_read(unit,"z1",tt%z1)
    call json_read(unit,"vol",tt%volint)
    call json_read(unit,"volint",tt%volint)
    call json_read(unit,"amat",rtemp)
    tt%amat=reshape(rtemp,(/3,3/))
    call json_read(unit,"bmat",rtemp)
    tt%bmat=reshape(rtemp,(/3,3/))
    call json_read(unit,"bbmat",rtemp)
    tt%bbmat=reshape(rtemp,(/3,3/))
    
    call json_close_class(unit)
    
  END SUBROUTINE read_cell

 


  SUBROUTINE read_xml_cell(tt)
    USE m_xmlIntWrapFort
    USE m_constants
    USE m_calculator
    USE m_invs3
    IMPLICIT NONE
    CLASS(t_cell),INTENT(OUT):: tt


    LOGICAL::film
    CHARACTER(len=200):: xpath,valueString
    CHARACTER(len=3)  :: latnam
    REAL              :: atemp, a1(3),a2(3),a3(3),latticeScale

    a1=0.0;a2=0.0;a3=0.0
    tt%z1 = 0.0

    film=(xmlGetNumberOfNodes('/fleurInput/cell/filmLattice')==1)
    IF (.NOT.(film.OR.(xmlGetNumberOfNodes('/fleurInput/cell/bulkLattice')==1))) CALL judft_error("Either filmLattice or bulkLattice must be specified in inp.xml",calledby="types_cell.F90")
    
    IF (film) THEN
       xPath = '/fleurInput/cell/filmLattice'
    ELSE
       xPath = '/fleurInput/cell/bulkLattice'
    END IF
    latticeScale = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@scale'))
    latnam = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@latnam')))

    IF (film) THEN
       tt%z1 = latticeScale*evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@dVac'))
       a3(3) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@dTilda'))
    END IF
    !read a1,a2,c if present
    IF (xmlGetNumberOfNodes(TRIM(ADJUSTL(xpath))//'/a1')==1) &
         a1(1) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xpath))//'/a1/@scale'))*&
         evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xpath))//'/a1'))
    
    IF (xmlGetNumberOfNodes(TRIM(ADJUSTL(xpath))//'/a2')==1) &
         a2(2) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xpath))//'/a2/@scale'))*&
         evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xpath))//'/a2'))
    
    IF (xmlGetNumberOfNodes(TRIM(ADJUSTL(xpath))//'/c')==1) &
         a3(3) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xpath))//'/c/@scale'))*&
         evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xpath))//'/c'))
    
    
    
    !Determine bravais matrix
    IF (latnam=='any') THEN
       valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/bravaisMatrix/row-1')))
       a1(1) = evaluateFirst(valueString)
       a1(2) = evaluateFirst(valueString)
       IF(.NOT.film) a1(3) = evaluateFirst(valueString)
       valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/bravaisMatrix/row-2')))
       a2(1) = evaluateFirst(valueString)
       a2(2) = evaluateFirst(valueString)
       IF(.NOT.film) THEN
          a2(3) = evaluateFirst(valueString)
          valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/bravaisMatrix/row-3')))
          a3(1) = evaluateFirst(valueString)
          a3(2) = evaluateFirst(valueString)
          a3(3) = evaluateFirst(valueString)
       ENDIF
    ELSEIF(latnam=='squ') THEN
       a2(2) = a1(1)
    ELSE IF (latnam.EQ.'c-b') THEN
       aTemp = a1(1)
       a1(1) = aTemp*0.5*SQRT(2.0)
       a1(2) = -aTemp*0.5
       a2(1) = aTemp*0.5*SQRT(2.0)
       a2(2) = aTemp*0.5
    ELSE IF (latnam.EQ.'hex') THEN
       aTemp = 0.5*a1(1)
       a1(1) = aTemp*SQRT(3.0)
       a1(2) = -aTemp
       a2(1) = a1(1)
       a2(2) = aTemp
    ELSE IF (latnam.EQ.'hx3') THEN
       aTemp = 0.5*a1(1)
       a1(1) = aTemp
       a1(2) = -aTemp*SQRT(3.0)
       a2(1) = a1(1)
       a2(2) = -a1(2)
    ELSE IF (latnam.EQ.'fcc') THEN
       aTemp = a1(1)
       a1(1) =       0.0 ; a1(2) = 0.5*aTemp ; a1(3) = 0.5*aTemp
       a2(1) = 0.5*aTemp ; a2(2) =       0.0 ; a2(3) = 0.5*aTemp
       a3(1) = 0.5*aTemp ; a3(2) = 0.5*aTemp ; a3(3) =       0.0
    ELSE IF (latnam.EQ.'bcc') THEN
       aTemp = a1(1)
       a1(1) =-0.5*aTemp ; a1(2) = 0.5*aTemp ; a1(3) = 0.5*aTemp
       a2(1) = 0.5*aTemp ; a2(2) =-0.5*aTemp ; a2(3) = 0.5*aTemp
       a3(1) = 0.5*aTemp ; a3(2) = 0.5*aTemp ; a3(3) =-0.5*aTemp
    ELSE IF ((latnam.EQ.'c-r').OR.(latnam.EQ.'p-r')) THEN
       IF (latnam.EQ.'c-r') THEN
          a1(2) = -a2(2)
          a2(1) =  a1(1)
       END IF
    ELSEIF (latnam.EQ.'obl') THEN
       valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xpath))//'/row-1')))
       a1(1) = evaluateFirst(valueString)
       a1(2) = evaluateFirst(valueString)
       valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xpath))//'/row-2')))
       a2(1) = evaluateFirst(valueString)
       a2(2) = evaluateFirst(valueString)
    ENDIF
    !Set Bravais matrix in type
    tt%amat(:,1) = a1(:) * latticeScale
    tt%amat(:,2) = a2(:) * latticeScale
    tt%amat(:,3) = a3(:) * latticeScale
    
    CALL inv3(tt%amat,tt%bmat,tt%omtil)
    tt%bmat(:,:) = tpi_const*tt%bmat(:,:)

    !Volumes
    tt%omtil = ABS(tt%omtil)
    IF (film) THEN
       tt%vol = (tt%omtil/tt%amat(3,3))*tt%z1
       tt%area = tt%omtil/tt%amat(3,3)
    ELSE
       tt%vol = tt%omtil
       tt%area = tt%amat(1,1)*tt%amat(2,2)-tt%amat(1,2)*tt%amat(2,1)
       IF (tt%area.LT.1.0e-7) THEN
          IF (latnam.EQ.'any') THEN
             tt%area = 1.
          ELSE
             CALL juDFT_error("area = 0",calledby ="r_inpXML")
          END IF
       END IF
    END IF
    tt%volint=-1 !could not be determined here
    
  END SUBROUTINE read_xml_cell
  
  
END MODULE m_types_cell
