!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_banddos
  USE m_judft
  USE m_types_fleur_setup
  USE m_json_tools
  IMPLICIT NONE
  TYPE,EXTENDS(t_fleursetup):: t_banddos
        LOGICAL :: dos
     LOGICAL :: band
     LOGICAL :: l_mcd
     LOGICAL :: l_orb
     LOGICAL :: vacdos
     INTEGER :: ndir
     INTEGER :: orbCompAtom
     REAL    :: e1_dos
     REAL    :: e2_dos
     REAL    :: sig_dos
     REAL    :: e_mcd_lo
     REAL    :: e_mcd_up
     LOGICAL :: unfoldband
     INTEGER :: s_cell_x
     INTEGER :: s_cell_y
     INTEGER :: s_cell_z     
   CONTAINS
     PROCEDURE,PASS :: broadcast=>broadcast_banddos
     PROCEDURE,PASS :: write=>WRITE_banddos
     PROCEDURE,PASS :: read=>READ_banddos
     PROCEDURE,PASS :: read_xml=>read_xml_banddos
     PROCEDURE,PASS :: init=>init_banddos
  END TYPE t_banddos

CONTAINS
  SUBROUTINE broadcast_banddos(tt,mpi_comm,origin)
#ifdef CPP_MPI
    USE m_mpi_bc_tool
#endif    
    IMPLICIT NONE
    CLASS(t_banddos),INTENT(INOUT):: tt
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


call MPI_BCAST(tt%dos,1,MPI_LOGICAL,pe,mpi_comm,ierr)
call MPI_BCAST(tt%band,1,MPI_LOGICAL,pe,mpi_comm,ierr)
call MPI_BCAST(tt%l_mcd,1,MPI_LOGICAL,pe,mpi_comm,ierr)
call MPI_BCAST(tt%l_orb,1,MPI_LOGICAL,pe,mpi_comm,ierr)
call MPI_BCAST(tt%vacdos,1,MPI_LOGICAL,pe,mpi_comm,ierr)
call MPI_BCAST(tt%ndir,1,MPI_INTEGER,pe,mpi_comm,ierr)
call MPI_BCAST(tt%orbCompAtom,1,MPI_INTEGER,pe,mpi_comm,ierr)
call MPI_BCAST(tt%e1_dos,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
call MPI_BCAST(tt%e2_dos,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
call MPI_BCAST(tt%sig_dos,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
call MPI_BCAST(tt%e_mcd_lo,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
call MPI_BCAST(tt%e_mcd_up,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
call MPI_BCAST(tt%unfoldband,1,MPI_LOGICAL,pe,mpi_comm,ierr)
call MPI_BCAST(tt%s_cell_x,1,MPI_INTEGER,pe,mpi_comm,ierr)
call MPI_BCAST(tt%s_cell_y,1,MPI_INTEGER,pe,mpi_comm,ierr)
call MPI_BCAST(tt%s_cell_z,1,MPI_INTEGER,pe,mpi_comm,ierr)
    
#endif
      
  END SUBROUTINE broadcast_banddos

  SUBROUTINE write_banddos(tt, unit, iotype, v_list, iostat, iomsg)
    IMPLICIT NONE
    CLASS(t_banddos),INTENT(IN):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    INTEGER:: ntype
    REAL,ALLOCATABLE::d_wgn(:,:)

    WRITE(unit,*,IOSTAT=iostat) '"banddos":{'
 
call json_print(unit,"dos",tt%dos)
call json_print(unit,"band",tt%band)
call json_print(unit,"l_mcd",tt%l_mcd)
call json_print(unit,"l_orb",tt%l_orb)
call json_print(unit,"vacdos",tt%vacdos)
call json_print(unit,"ndir",tt%ndir)
call json_print(unit,"orbCompAtom",tt%orbCompAtom)
call json_print(unit,"e1_dos",tt%e1_dos)
call json_print(unit,"e2_dos",tt%e2_dos)
call json_print(unit,"sig_dos",tt%sig_dos)
call json_print(unit,"e_mcd_lo",tt%e_mcd_lo)
call json_print(unit,"e_mcd_up",tt%e_mcd_up)
call json_print(unit,"unfoldband",tt%unfoldband)
call json_print(unit,"s_cell_x",tt%s_cell_x)
call json_print(unit,"s_cell_y",tt%s_cell_y)
call json_print(unit,"s_cell_z",tt%s_cell_z)
    WRITE(unit,*,IOSTAT=iostat) '}'
    
  END SUBROUTINE write_banddos
  SUBROUTINE read_banddos(tt, unit, iotype, v_list, iostat, iomsg)
    use m_json_tools
    IMPLICIT NONE
    CLASS(t_banddos),INTENT(INOUT):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    INTEGER :: ntype
    REAL,ALLOCATABLE::rt(:,:)
    CALL json_open_class("banddos",unit,iostat)
    IF (iostat.NE.0)   RETURN


call json_read(unit,"dos",tt%dos)
call json_read(unit,"band",tt%band)
call json_read(unit,"l_mcd",tt%l_mcd)
call json_read(unit,"l_orb",tt%l_orb)
call json_read(unit,"vacdos",tt%vacdos)
call json_read(unit,"ndir",tt%ndir)
call json_read(unit,"orbCompAtom",tt%orbCompAtom)
call json_read(unit,"e1_dos",tt%e1_dos)
call json_read(unit,"e2_dos",tt%e2_dos)
call json_read(unit,"sig_dos",tt%sig_dos)
call json_read(unit,"e_mcd_lo",tt%e_mcd_lo)
call json_read(unit,"e_mcd_up",tt%e_mcd_up)
call json_read(unit,"unfoldband",tt%unfoldband)
call json_read(unit,"s_cell_x",tt%s_cell_x)
call json_read(unit,"s_cell_y",tt%s_cell_y)
call json_read(unit,"s_cell_z",tt%s_cell_z)
    

    CALL json_close_class(unit,iostat)
    
  END SUBROUTINE read_banddos

 
  SUBROUTINE read_xml_banddos(tt)
    USE m_xmlIntWrapFort
    USE m_constants
    USE m_calculator
    USE m_inp_xml
    IMPLICIT NONE
    CLASS(t_banddos),INTENT(OUT):: tt

    CHARACTER(LEN=255) :: xPathA, xPathB
    INTEGER            :: n,i,na,numberNodes

    tt%l_mcd = .FALSE.
    tt%e_mcd_lo = -10.0
    tt%e_mcd_up = 0.0
    tt%unfoldband = .FALSE.
    
    tt%l_orb = .FALSE.
    tt%orbCompAtom = 0
    na=0
    DO n=1,xmlGetNumberOfNodes('/fleurInput/atomGroups/atomGroup')
       xpathA=inp_xml_xpath_for_group(n)
       if(xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/relPos')>0) xpatha=TRIM(ADJUSTL(xPathA))//'/relPos'
       if(xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/filmPos')>0) xpatha=TRIM(ADJUSTL(xPathA))//'/filmPos'

       DO i = 1,xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathA))) 
          na = na + 1
          WRITE(xPathB,*) TRIM(ADJUSTL(xPathA)),'[',i,']'
          IF(evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@orbcomp'))) THEN
             IF(tt%l_orb) THEN
                CALL juDFT_error("Multiple orbcomp flags set.", calledby = "r_inpXML")
             END IF
             tt%l_orb = .TRUE.
             tt%orbCompAtom = na
          END IF
       end do
    end DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Start of output section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    tt%dos = .FALSE.
    tt%band = .FALSE.
    tt%vacdos = .FALSE.
    xPathA = '/fleurInput/output'
    numberNodes = xmlGetNumberOfNodes(xPathA)

    IF (numberNodes.EQ.1) THEN

       ! Read in general output switches
       tt%dos = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@dos'))
       tt%band = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@band'))
       tt%vacdos = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@vacdos'))
       tt%l_mcd = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@mcd'))
       ! Read in optional densityOfStates output parameters

       xPathA = '/fleurInput/output/densityOfStates'
       numberNodes = xmlGetNumberOfNodes(xPathA)

       IF ((tt%dos).AND.(numberNodes.EQ.0)) THEN
          CALL juDFT_error("dos is true but densityOfStates parameters are not set!", calledby = "r_inpXML")
       END IF

       IF (numberNodes.EQ.1) THEN
          tt%ndir = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@ndir'))
          tt%e2_dos = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@minEnergy'))
          tt%e1_dos = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@maxEnergy'))
          tt%sig_dos = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@sigma'))
       END IF

       ! Read in optional vacuumDOS parameters

       xPathA = '/fleurInput/output/vacuumDOS'
       numberNodes = xmlGetNumberOfNodes(xPathA)

       IF ((tt%vacdos).AND.(numberNodes.EQ.0)) THEN
          CALL juDFT_error("vacdos is true but vacDOS parameters are not set!", calledby = "r_inpXML")
       END IF

       IF (tt%band) THEN
          tt%dos=.TRUE.
          tt%ndir = -4
          WRITE(*,*) 'band="T" --> Overriding "dos" and "ndir"!'
       ENDIF
       ! Read in optional magnetic circular dichroism parameters
       xPathA = '/fleurInput/output/magneticCircularDichroism'
       numberNodes = xmlGetNumberOfNodes(xPathA)

       IF ((tt%l_mcd).AND.(numberNodes.EQ.0)) THEN
          CALL juDFT_error("mcd is true but magneticCircularDichroism parameters are not set!", calledby = "r_inpXML")
       END IF

       IF (numberNodes.EQ.1) THEN
          tt%e_mcd_lo = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@energyLo'))
          tt%e_mcd_up = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@energyUp'))
       END IF

       ! Read in optional parameter for unfolding bandstructure of supercell
       xPathA = '/fleurInput/output/unfoldingBand'
       numberNodes = xmlGetNumberOfNodes(xPathA)

       IF (numberNodes.EQ.1) THEN
          tt%unfoldband = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@unfoldband'))
          tt%s_cell_x = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@supercellX'))
          tt%s_cell_y = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@supercellY'))
          tt%s_cell_z = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@supercellZ'))
       END IF
    ENDIF
  END SUBROUTINE read_xml_banddos
  
  SUBROUTINE init_banddos(banddos)
    IMPLICIT NONE
    CLASS(t_banddos),INTENT(INOUT):: banddos
    !TODO
  END SUBROUTINE init_banddos

  
END MODULE m_types_banddos
