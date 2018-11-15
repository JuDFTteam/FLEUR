!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_sym
  USE m_judft
  USE m_types_fleur_setup
  IMPLICIT NONE
  TYPE,EXTENDS(t_fleursetup):: t_sym
     !Symophic group
     LOGICAL ::symor
     INTEGER ::nsymt
     INTEGER :: nsym
     COMPLEX,ALLOCATABLE:: d_wgn(:,:,:,:)
     !2D-inv-sym
     LOGICAL ::invs2
     !Inversion-sym
     LOGICAL ::invs
     !Z-refls. sym
     LOGICAL ::zrfs
     !No of sym ops
     INTEGER ::nop
     !No of 2D-sym ops
     INTEGER ::nop2
     !Rot-matrices (3,3,nop)
     INTEGER,ALLOCATABLE::mrot(:,:,:)
     !inverse operation (nop)
     INTEGER,ALLOCATABLE::invtab(:)
     !translation vectors (3,nop)
     REAL,ALLOCATABLE::tau(:,:)
     INTEGER, ALLOCATABLE :: multab(:,:)
     INTEGER, ALLOCATABLE :: invsatnr(:)
     INTEGER, ALLOCATABLE :: invarop(:,:)
     INTEGER, ALLOCATABLE :: invarind(:)
   CONTAINS
     PROCEDURE,PASS :: broadcast=>broadcast_sym
     PROCEDURE,PASS :: write=>WRITE_sym
     PROCEDURE,PASS :: read=>READ_sym
     PROCEDURE,PASS :: read_xml=>read_xml_sym
  END TYPE t_sym

CONTAINS
  SUBROUTINE broadcast_sym(tt,mpi_comm,origin)
    IMPLICIT NONE
    CLASS(t_sym),INTENT(INOUT):: tt
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

    IF (irank==pe) THEN
       ntype=SIZE(tt%l_relax)
       CALL MPI_BCAST(ntype,1,MPI_INTEGER,pe,mpi_comm,ierr)
    ELSE
       CALL MPI_BCAST(ntype,1,MPI_INTEGER,pe,mpi_comm,ierr)
       IF (ALLOCATED(tt%l_relax)) &
            DEALLOCATE(tt%l_relax,tt%b_con,tt%alphInit,tt%alph,tt%beta,tt%socscale)
       ALLOCATE(tt%l_relax(ntype),tt%b_con(2,ntype))
       ALLOCATE(tt%alphInit(ntype),tt%alph(ntype),tt%beta(ntype))
       ALLOCATE(tt%socscale(ntype))
    ENDIF
    CALL MPI_BCAST(tt%alphinit,ntype,MPI_DOUBLE_PRECSION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%alph,ntype,MPI_DOUBLE_PRECSION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%beta,ntype,MPI_DOUBLE_PRECSION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%b_con,2*ntype,MPI_DOUBLE_PRECSION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%socscale,ntype,MPI_DOUBLE_PRECSION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%l_relax,ntype,MPI_LOGICAL,pe,mpi_comm,ierr)
  
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
    CALL json_print(unit,"ntype",size(tt%l_relax))
    
    CALL json_print(unit,"qss",tt%qss,',')
 
    CALL json_print(unit,"l_relax",tt%l_relax)
    CALL json_print(unit,"alphInit",tt%alphInit)
    CALL json_print(unit,"alph",tt%alph)
    CALL json_print(unit,"beta",tt%beta)
    CALL json_print(unit,"b_const",RESHAPE(tt%b_const,(/SIZE(tt%b_const)/))
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
    real,allocatable::rtemp(:)
    CALL json_open_class("noco",unit,iostat)
    IF (iostat.NE.0)   RETURN
    
    call json_read(unit,"l_noco",tt%l_noco)
    call json_read(unit,"l_ss",tt%l_ss)
    call json_read(unit,"l_mperp",tt%l_mperp)
    call json_read(unit,"l_constr",tt%l_constr)
    call json_read(unit,"l_mtnocopot",tt%l_mtnocopot)
    call json_read(unit,"l_soc",tt%l_soc)
    call json_read(unit,"l_spav",tt%l_spav)
    
    call json_read(unit,"mix_b",tt%mix_b)
    call json_read(unit,"theta",tt%theta)
    call json_read(unit,"phi",tt%phi)
    
    call json_read(unit,"ntype",ntype)

    call json_read(unit,"qss",tt%qss)

    IF (ALLOCATED(tt%l_relax)) DEALLOCATE(tt%l_relax,tt%b_con,tt%alphInit,tt%alph,tt%beta,tt%socscale)
    ALLOCATE(tt%l_relax(ntype),tt%b_con(2,ntype))
    ALLOCATE(tt%alphInit(ntype),tt%alph(ntype),tt%beta(ntype))
    ALLOCATE(tt%socscale(ntype))

    call json_read(unit,"l_relax",tt%l_relax)
    call json_read(unit,"alphInit",tt%l_alphInit)
    call json_read(unit,"alph",tt%l_alph)
    call json_read(unit,"beta",tt%l_beta)
    allocate(rtemp(2*ntype))
    call json_read_rarray(unit,"b_const",rtemp)
    tt%l_b_const=RESHAPE(rtemp,(/2,ntype/))
    call json_read(unit,"socscale",tt%l_socscale)

    CALL json_close_class(unit)
    
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
       tt%theta=evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@theta'))
       tt%phi=evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@phi'))
       tt%l_soc = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_soc'))
       tt%l_spav = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@spav'))
    END IF
    !read species depended stuff
    IF (tt%l_soc) THEN
       DO n=1,ntype
          xpath=inp_xml_speciesxpath_for_group(n)
          IF (xmlGetNumberOfNodes(TRIM(ADJUSTL(xPath))//'/special')==1) THEN
             tt%socscale(n)=evaluateFirstOnly(TRIM(ADJUSTL(xpath))//'/special/@socscale'))))
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
       tt%l_ss = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_ss'))
       tt%l_mperp = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_mperp'))
       tt%l_constr = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_constr'))
       tt%l_mtNocoPot = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_mtNocoPot'))
       
       tt%mix_b = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@mix_b'))
       valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/qss')))
       READ(valueString,*) tt%qss(1), tt%qss(2), tt%qss(3)
    END IF
    !Now the group dependent stuff
    IF (tt%l_noco) THEN
       DO n=1,ntype
          xpath=inp_xml_xpath_for_group(n)
          IF (xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathB))//"/nocoParams")==1) THEN
             tt%l_relax(n) = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/nocoParams/@l_relax'))
            tt%alphInit(n) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/nocoParams/@alpha'))
            tt%alph(n) = tt%alphInit(iType)
            tt%beta(n) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/nocoParams/@beta'))
            tt%b_con(1,n) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/nocoParams/@b_cons_x'))
            tt%b_con(2,n) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/nocoParams/@b_cons_y'))
         END IF
      END DO
   END IF
      
 END SUBROUTINE read_xml_noco
  
  
END MODULE m_types_cell
