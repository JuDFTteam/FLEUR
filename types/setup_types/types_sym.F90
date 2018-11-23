!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_sym
  USE m_judft
  USE m_types_fleur_setup
  USE m_json_tools
  IMPLICIT NONE
  TYPE,EXTENDS(t_fleursetup):: t_sym
     !Symophic group
     LOGICAL ::symor
     INTEGER ::nsymt
     INTEGER :: nsym
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
     COMPLEX,ALLOCATABLE  :: d_wgn(:,:,:,:)
   CONTAINS
     PROCEDURE,PASS :: broadcast=>broadcast_sym
     PROCEDURE,PASS :: write=>WRITE_sym
     PROCEDURE,PASS :: read=>READ_sym
     PROCEDURE,PASS :: read_xml=>read_xml_sym
     PROCEDURE,PASS :: init=>init_sym
  END TYPE t_sym

CONTAINS
  SUBROUTINE broadcast_sym(tt,mpi_comm,origin)
#ifdef CPP_MPI
    USE m_mpi_bc_tool
#endif    
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

    CALL MPI_BCAST(tt%l_symor,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%invs2,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%invs,1,MPI_LOGICAL,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%zrfs,1,MPI_LOGICAL,pe,mpi_comm,ierr)

    CALL MPI_BCAST(tt%nsymt,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%nsym,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%nop,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%nop2,1,MPI_INTEGER,pe,mpi_comm,ierr)
    
    
    CALL MPI_BC(tt%mrot,pe,mpi_comm)
    CALL MPI_BC(tt%invtab,pe,mpi_comm)
    CALL MPI_BC(tt%tau,pe,mpi_comm)
    CALL MPI_BC(tt%multab,pe,mpi_comm)

    CALL MPI_BC(tt%d_wgn,pe,mpi_comm)

    CALL MPI_BC(tt%invarop,pe,mpi_comm)
    CALL MPI_BC(tt%invarind,pe,mpi_comm)
    CALL MPI_BC(tt%invsatnr,pe,mpi_comm)
    
#endif
      
  END SUBROUTINE broadcast_sym

  SUBROUTINE write_sym(tt, unit, iotype, v_list, iostat, iomsg)
    IMPLICIT NONE
    CLASS(t_sym),INTENT(IN):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    INTEGER:: ntype
    REAL,ALLOCATABLE::d_wgn(:,:)

    WRITE(unit,*,IOSTAT=iostat) '"sym":{'
    CALL json_print(unit,"symor",tt%symor)
    CALL json_print(unit,"invs",tt%invs)
    CALL json_print(unit,"invs2",tt%invs2)
    CALL json_print(unit,"zrfs",tt%zrfs)
    CALL json_print(unit,"nsymt",tt%nsymt)
    CALL json_print(unit,"nsym",tt%nsym)
    CALL json_print(unit,"nop",tt%nop)
    CALL json_print(unit,"nop2",tt%nop2)

    
    CALL json_print(unit,"mrot",tt%mrot)
    CALL json_print(unit,"multab",tt%multab)
    CALL json_print(unit,"invtab",tt%invtab)
    CALL json_print(unit,"tau",tt%tau)
  
    ALLOCATE(d_wgn(2,SIZE(tt%d_wgn)))
    d_wgn(1,:)=RESHAPE(REAL(tt%d_wgn),(/SIZE(tt%d_wgn)/))
    d_wgn(2,:)=RESHAPE(aimag(tt%d_wgn),(/SIZE(tt%d_wgn)/))
    CALL json_print(unit,"d_wgn",d_wgn)

    
    CALL json_print(unit,"invarop",tt%invarop)
    CALL json_print(unit,"invarind",tt%invarind)
    CALL json_print(unit,"invsatnr",tt%invsatnr,.true.)

     
    WRITE(unit,*,IOSTAT=iostat) '}'
    
  END SUBROUTINE write_sym
  SUBROUTINE read_sym(tt, unit, iotype, v_list, iostat, iomsg)
    use m_json_tools
    IMPLICIT NONE
    CLASS(t_sym),INTENT(INOUT):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    INTEGER :: ntype
    REAL,ALLOCATABLE::rt(:,:)
    CALL json_open_class("sym",unit,iostat)
    IF (iostat.NE.0)   RETURN


    CALL json_read(unit,"symor",tt%symor)
    CALL json_read(unit,"invs",tt%invs)
    CALL json_read(unit,"invs2",tt%invs2)
    CALL json_read(unit,"zrfs",tt%zrfs)
    CALL json_read(unit,"nsymt",tt%nsymt)
    CALL json_read(unit,"nsym",tt%nsym)
    CALL json_read(unit,"nop",tt%nop)
    CALL json_read(unit,"nop2",tt%nop2)

    CALL json_read(unit,"mrot",tt%mrot)
    CALL json_read(unit,"multab",tt%multab)
    CALL json_read(unit,"invtab",tt%invtab)
    CALL json_read(unit,"tau",tt%tau)

    IF (ALLOCATED(tt%d_wgn)) DEALLOCATE(tt%d_wgn)
    CALL json_read(unit,"d_wgn",rt)
    ALLOCATE(tt%d_wgn(-3:3,-3:3,3,tt%nop))
    tt%d_wgn=CMPLX(RESHAPE(rt(1,:),SHAPE(tt%d_wgn)),RESHAPE(rt(2,:),SHAPE(tt%d_wgn)))
    
    CALL json_read(unit,"invarop",tt%invarop)
    CALL json_read(unit,"invarind",tt%invarind)
    CALL json_read(unit,"invsatnr",tt%invsatnr)
    

    CALL json_close_class(unit,iostat)
    
  END SUBROUTINE read_sym

 
  SUBROUTINE read_xml_sym(tt)
    USE m_xmlIntWrapFort
    USE m_constants
    USE m_calculator
    USE m_inp_xml
    USE m_symdata , ONLY : nammap, ord2, l_c2
    USE m_spg2set
    IMPLICIT NONE
    CLASS(t_sym),INTENT(OUT):: tt


    LOGICAL::film,invsym
    CHARACTER(len=200):: xpath,valueString,xpathB
    CHARACTER(len=4)  :: namgrp
    CHARACTER(len=3)  :: latnam
    INTEGER           :: ntype,numberNodes,n,i,n2spg
    INTEGER:: j,nop48
    INTEGER,ALLOCATABLE::invOps(:),multtab(:,:),optype(:)

    !Two different modes to specify symmetry in inp.xml
    !Either we give the sym-group (only in film case)
    !Or we list all operations (see case below)
    IF (xmlGetNumberOfNodes('/fleurInput/cell/symmetry')==1) THEN
       valueString = TRIM(ADJUSTL(xmlGetAttributeValue('/fleurInput/cell/symmetry/@spgrp')))
       READ(valueString,*) namgrp
       tt%invs = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/cell/symmetry/@invs'))
       tt%zrfs = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/cell/symmetry/@zrfs'))
       tt%invs2 = tt%invs.AND.tt%zrfs
       !Now determine the latnam from cell
       IF (xmlGetNumberOfNodes('/fleurInput/cell/bulkLattice')==1) THEN
          latnam=TRIM(ADJUSTL(xmlGetAttributeValue('/fleurInput/cell/bulkLattice/@latnam')))
          film=.FALSE.
       ELSE
          latnam=TRIM(ADJUSTL(xmlGetAttributeValue('/fleurInput/cell/filmLattice/@latnam')))
          film=.TRUE.
       END IF
       
       DO  n2spg=1, SIZE(nammap)
          IF (namgrp.EQ.nammap(i)) EXIT
       END DO
       IF (n2spg >size(nammap) ) THEN
          WRITE (*,*) 'Spacegroup ',namgrp,' not known! Choose one of:'
          WRITE (*,'(20(a4,1x))') (nammap(i),i=1,20)
          CALL juDFT_error("Could not determine spacegroup!", calledby = "r_inpXML")
       END IF
       IF ((n2spg.GE.13).AND.(n2spg.LE.17)) THEN
          IF (.NOT.((latnam.EQ.'hx3').OR.(latnam.EQ.'hex'))) THEN
             CALL juDFT_error("Use only hex or hx3 with p3, p3m1, p31m, p6 or p6m!", calledby ="r_inpXML")
          END IF
       END IF
       tt%nop = ord2(n2spg)
       IF (tt%invs) THEN
          tt%nop = 2*tt%nop
          IF (tt%zrfs.AND.(.NOT.l_c2(n2spg))) tt%nop = 2*tt%nop
       ELSE
          IF (tt%zrfs) tt%nop = 2*tt%nop
       END IF
       IF (ALLOCATED(tt%mrot)) DEALLOCATE(tt%mrot)
       ALLOCATE(tt%mrot(3,3,tt%nop))
       IF (ALLOCATED(tt%tau)) DEALLOCATE(tt%tau)
       ALLOCATE(tt%tau(3,tt%nop))
       CALL spg2set(tt%nop,tt%zrfs,tt%invs,namgrp,latnam,&
            tt%mrot,tt%tau,tt%nop2,tt%symor)
    ELSEIF(xmlGetNumberOfNodes('/fleurInput/cell/symmetryOperations')==1) THEN
       xpath='/fleurInput/cell/symmetryOperations'
       tt%nop = xmlGetNumberOfNodes(TRIM(ADJUSTL(xPath))//'/symOp')
       IF (ALLOCATED(tt%mrot)) DEALLOCATE(tt%mrot)
       ALLOCATE(tt%mrot(3,3,tt%nop))
       IF (ALLOCATED(tt%tau)) DEALLOCATE(tt%tau)
       ALLOCATE(tt%tau(3,tt%nop))
       tt%symor = .TRUE.
       DO i = 1, tt%nop
            WRITE(xPathB,*) TRIM(ADJUSTL(xPath))//'/symOp[',i,']/row-1'
            valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB)))))
            READ(valueString,*) tt%mrot(1,1,i), tt%mrot(1,2,i), tt%mrot(1,3,i), tt%tau(1,i)

            WRITE(xPathB,*) TRIM(ADJUSTL(xPath))//'/symOp[',i,']/row-2'
            valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB)))))
            READ(valueString,*) tt%mrot(2,1,i), tt%mrot(2,2,i), tt%mrot(2,3,i), tt%tau(2,i)

            WRITE(xPathB,*) TRIM(ADJUSTL(xPath))//'/symOp[',i,']/row-3'
            valueString = TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathB)))))
            READ(valueString,*) tt%mrot(3,1,i), tt%mrot(3,2,i), tt%mrot(3,3,i), tt%tau(3,i)

            IF ((tt%tau(1,i)**2 + tt%tau(2,i)**2 + tt%tau(3,i)**2).GT.1.e-8) THEN
               tt%symor = .FALSE.
            END IF
         END DO
         DO i=-3,3
            DO j=1,5
               WHERE(ABS(tt%tau-(1.*i)/j)<1E-5) tt%tau=(1.*i)/j
            ENDDO
         ENDDO
        
      ELSE
         CALL judft_error("No valid specification of symmetry in inp.xml")
      END IF
      
  END SUBROUTINE read_xml_sym
  
  SUBROUTINE init_sym(sym,cell,input)
    USE m_types_cell
    USE m_types_input
    USE m_closure
    USE m_symproperties
    IMPLICIT NONE
    CLASS(t_sym),INTENT(INOUT):: sym
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_input),INTENT(IN)  :: input

    LOGICAL :: invSym
    INTEGER :: nop48
    INTEGER,ALLOCATABLE::invOps(:),multtab(:,:),optype(:)
    
    nop48 = 48
    ALLOCATE (invOps(sym%nop),multtab(sym%nop,sym%nop),optype(nop48))
    CALL check_close(sym%nop,sym%mrot,sym%tau,&
         &                      multtab,invOps,optype)
    
    CALL symproperties(nop48,optype,input%film,sym%nop,multtab,cell%amat,&
         &                        sym%symor,sym%mrot,sym%tau,&
         &                        invSym,sym%invs,sym%zrfs,sym%invs2,sym%nop,sym%nop2)
    IF (.NOT.input%film) sym%nop2=sym%nop
    IF (input%film.AND.ANY(ABS(sym%tau(:,:sym%nop))>1E-5)) &
         CALL juDFT_error("nonsymmorphic symmetries not yet implemented for films!",calledby ="r_inpXML")
    sym%invs2 = sym%invs.AND.sym%zrfs
    !TODO
  END SUBROUTINE init_sym

  
END MODULE m_types_sym
