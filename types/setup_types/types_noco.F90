!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_noco
  USE m_judft
  USE m_types_fleur_setup
  USE m_json_tools
  IMPLICIT NONE
  TYPE,EXTENDS(t_fleursetup):: t_noco
     LOGICAL:: l_noco
     LOGICAL:: l_ss
     LOGICAL:: l_mperp
     LOGICAL:: l_constr
     LOGICAL:: l_mtNocoPot
     REAL:: qss(3)
     REAL:: mix_b
     LOGICAL, ALLOCATABLE :: l_relax(:)
     REAL, ALLOCATABLE :: alphInit(:)
     REAL, ALLOCATABLE :: alph(:)
     REAL, ALLOCATABLE :: beta(:)
     REAL, ALLOCATABLE :: b_con(:,:)
     LOGICAL           :: l_soc
     LOGICAL           :: l_spav
     REAL              :: theta
     REAL              :: phi
     REAL,ALLOCATABLE  :: socscale(:)
   CONTAINS
     PROCEDURE,PASS :: broadcast=>broadcast_noco
     PROCEDURE,PASS :: write=>WRITE_noco
     PROCEDURE,PASS :: read=>READ_noco
     PROCEDURE,PASS :: read_xml=>read_xml_noco
  END TYPE t_noco

CONTAINS
  SUBROUTINE broadcast_noco(tt,mpi_comm,origin)
    IMPLICIT NONE
    CLASS(t_noco),INTENT(INOUT):: tt
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

    CALL MPI_BCAST(tt%l_noco,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%l_ss,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%l_mperp,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%l_constr,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%l_mtNocoPot,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%l_l_soc,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%l_l_spav,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%qss,3,MPI_DOUBLE_PRECSION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%mix_b,1,MPI_DOUBLE_PRECSION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%theta,1,MPI_DOUBLE_PRECSION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%phi,1,MPI_DOUBLE_PRECSION,pe,mpi_comm,ierr)

    CALL MPI_BC(tt%alphinit,pe,mpi_comm)
    CALL MPI_BC(tt%alph,,pe,mpi_comm)
    CALL MPI_BC(tt%beta,pe,mpi_comm)
    CALL MPI_BC(tt%b_con,pe,mpi_comm)
    CALL MPI_BC(tt%socscale,pe,mpi_comm)
    CALL MPI_BC(tt%l_relax,pe,mpi_comm)
  
#endif
      
  END SUBROUTINE broadcast_noco

  SUBROUTINE write_noco(tt, unit, iotype, v_list, iostat, iomsg)
    IMPLICIT NONE
    CLASS(t_noco),INTENT(IN):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    INTEGER:: ntype

    WRITE(unit,*,IOSTAT=iostat) '"noco":{'
    CALL json_print(unit,"l_noco",tt%l_noco)
    CALL json_print(unit,"l_ss",tt%l_noco)
    CALL json_print(unit,"l_mperp",tt%l_noco)
    CALL json_print(unit,"l_constr",tt%l_noco)
    CALL json_print(unit,"l_mtnocopot",tt%l_noco)
    CALL json_print(unit,"l_soc",tt%l_noco)
    CALL json_print(unit,"l_spav",tt%l_noco)
    CALL json_print(unit,"mix_b",tt%l_noco)
    CALL json_print(unit,"theta",tt%l_noco)
    CALL json_print(unit,"phi",tt%l_noco)
   
    CALL json_print(unit,"qss",tt%qss)
 
    CALL json_print(unit,"l_relax",tt%l_relax)
    CALL json_print(unit,"alphInit",tt%alphInit)
    CALL json_print(unit,"alph",tt%alph)
    CALL json_print(unit,"beta",tt%beta)
    CALL json_print(unit,"b_const",tt%b_con)
    CALL json_print(unit,"socscale",tt%socscale,.TRUE.)
     
    WRITE(unit,*,IOSTAT=iostat) '}'
    
  END SUBROUTINE write_noco
  SUBROUTINE read_noco(tt, unit, iotype, v_list, iostat, iomsg)
    use m_json_tools
    IMPLICIT NONE
    CLASS(t_noco),INTENT(INOUT):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    INTEGER :: ntype
    REAL,ALLOCATABLE::rt(:)
    CALL json_open_class("noco",unit,iostat)
    IF (iostat.NE.0)   RETURN
    
    call json_read(unit,"l_noco",tt%l_noco)
    call json_read(unit,"l_ss",tt%l_ss)
    call json_read(unit,"l_mperp",tt%l_mperp)
    CALL json_read(unit,"l_constr",tt%l_constr)
    call json_read(unit,"l_mtnocopot",tt%l_mtnocopot)
    call json_read(unit,"l_soc",tt%l_soc)
    call json_read(unit,"l_spav",tt%l_spav)
    
    call json_read(unit,"mix_b",tt%mix_b)
    call json_read(unit,"theta",tt%theta)
    call json_read(unit,"phi",tt%phi)
    
 
    call json_read(unit,"qss",rt)
    tt%qss=rt
    call json_read(unit,"l_relax",tt%l_relax)
    call json_read(unit,"alphInit",tt%alphInit)
    call json_read(unit,"alph",tt%alph)
    call json_read(unit,"beta",tt%beta)
    call json_read(unit,"b_const",tt%b_con)
    call json_read(unit,"socscale",tt%socscale)

    CALL json_close_class(unit,iostat)
    
  END SUBROUTINE read_noco

 
  SUBROUTINE read_xml_noco(tt)
    USE m_xmlIntWrapFort
    USE m_constants
    USE m_calculator
    USE m_inp_xml
    IMPLICIT NONE
    CLASS(t_noco),INTENT(OUT):: tt


    LOGICAL::film
    CHARACTER(len=200):: xpath,valueString
    INTEGER           :: ntype,numberNodes,n

    tt%l_mtNocoPot = .FALSE.
    ntype=xmlGetNumberOfNodes('/fleurInput/atomGroups/atomGroup')
    IF (ALLOCATED(tt%l_relax)) DEALLOCATE(tt%l_relax,tt%b_con,tt%alphInit,tt%alph,tt%beta,tt%socscale)
    ALLOCATE(tt%l_relax(ntype),tt%b_con(2,ntype))
    ALLOCATE(tt%alphInit(ntype),tt%alph(ntype),tt%beta(ntype))
    ALLOCATE(tt%socscale(ntype))

    tt%l_noco = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/magnetism/@l_noco'))

    xPath = '/fleurInput/calculationSetup/soc'
    numberNodes = xmlGetNumberOfNodes(xPath)
    tt%l_soc = .FALSE.
    tt%theta = 0.0
    tt%phi = 0.0

    IF (numberNodes.EQ.1) THEN
       tt%theta=evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@theta'))
       tt%phi=evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@phi'))
       tt%l_soc = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@l_soc'))
       tt%l_spav = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@spav'))
    END IF
    !read species depended stuff
    IF (tt%l_soc) THEN
       DO n=1,ntype
          xpath=inp_xml_speciesxpath_for_group(n)
          IF (xmlGetNumberOfNodes(TRIM(ADJUSTL(xPath))//'/special')==1) THEN
             tt%socscale(n)=evaluateFirstOnly(TRIM(ADJUSTL(xpath))//'/special/@socscale')
          ELSE
             tt%socscale(n)=1.0
          ENDIF
       END DO
    ENDIF
    ! Read in optional noco parameters if present

    xPath = '/fleurInput/calculationSetup/nocoParams'
    numberNodes = xmlGetNumberOfNodes(xPath)
    
    tt%l_ss = .FALSE.
    tt%l_mperp = .FALSE.
    tt%l_constr = .FALSE.
    tt%mix_b = 0.0
    tt%qss = 0.0
    
    tt%l_relax(:) = .FALSE.
    tt%alphInit(:) = 0.0
    tt%alph(:) = 0.0
    tt%beta(:) = 0.0
    tt%b_con(:,:) = 0.0

    IF ((tt%l_noco).AND.(numberNodes.NE.1)) THEN
       STOP 'Error: l_noco is true but no noco parameters set in xml input file!'
    END IF
    
    IF (tt%l_noco) THEN
       tt%l_ss = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@l_ss'))
       tt%l_mperp = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@l_mperp'))
       tt%l_constr = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@l_constr'))
       tt%l_mtNocoPot = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@l_mtNocoPot'))
       
       tt%mix_b = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@mix_b'))
       valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/qss')))
       READ(valueString,*) tt%qss(1), tt%qss(2), tt%qss(3)
    END IF
    !Now the group dependent stuff
    IF (tt%l_noco) THEN
       DO n=1,ntype
          xpath=inp_xml_xpath_for_group(n)
          IF (xmlGetNumberOfNodes(TRIM(ADJUSTL(xPath))//"/nocoParams")==1) THEN
             tt%l_relax(n) = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/nocoParams/@l_relax'))
            tt%alphInit(n) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/nocoParams/@alpha'))
            tt%alph(n) = tt%alphInit(n)
            tt%beta(n) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/nocoParams/@beta'))
            tt%b_con(1,n) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/nocoParams/@b_cons_x'))
            tt%b_con(2,n) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/nocoParams/@b_cons_y'))
         END IF
      END DO
   END IF
      
 END SUBROUTINE read_xml_noco
  
  
END MODULE m_types_noco
