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
    
    !Arrays
    IF (irank.NE.pe.AND.ALLOCATED(tt%mrot)) THEN
       DEALLOCATE(tt%mrot,tt%invtab,tt%tau,tt%multab,tt%d_wgn)
    ENDIF
    
    CALL MPI_BCAST(tt%mrot,9*tt%nop,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%invtab,tt%nop,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%tau,3*tt%nop,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%multab,tt%nop**2,MPI_INTEGER,pe,mpi_comm,ierr)

    CALL MPI_BCAST(tt%d_wgn,147*tt%nop,MPI_DOUBLE_COMPLEX,pe,mpi_comm,ierr)

    !Atom dependent stuff
    IF (irank==pe) THEN
       nat=SIZE(tt%invsatnr)
       CALL MPI_BCAST(nat,1,MPI_INTEGER,pe,mpi_comm,ierr)
    ELSE
       CALL MPI_BCAST(nat,1,MPI_INTEGER,pe,mpi_comm,ierr)
       IF (ALLOCATED(tt%invarop)) DEALLOCATE(tt%invarop),tt%invarind,tt%invsatnr)
       ALLOCATE(tt%invarop(nat,tt%nop),tt%invarind(nat))
       ALLOCATE(tt%invsatnr(nat))
    ENDIF
    CALL MPI_BCAST(tt%invarop,nat*tt%nop,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%invarind,nat,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%invsatnr,nat,MPI_INTEGER,pe,mpi_comm,ierr)
    
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

    WRITE(unit,*,IOSTAT=iostat) '"sym":{'
    CALL json_print(unit,"l_symor",tt%l_symor)
    CALL json_print(unit,"invs",tt%invs)
    CALL json_print(unit,"invs2",tt%invs2)
    CALL json_print(unit,"zrfs",tt%zrfs)
    CALL json_print(unit,"nsymt",tt%nsymt)
    CALL json_print(unit,"nsym",tt%nsym)
    CALL json_print(unit,"nop",tt%nop)
    CALL json_print(unit,"nop2",tt%nop2)

    
    CALL json_print(unit,"mrot",RESHAPE(tt%mrot,(/SIZE(tt%mrot)/)))
    CALL json_print(unit,"multab",RESHAPE(tt%multab,(/SIZE(tt%multab)/)))
    CALL json_print(unit,"invtab",RESHAPE(tt%invtab,(/SIZE(tt%invtab)/)))
    CALL json_print(unit,"tau",RESHAPE(tt%tau,(/SIZE(tt%tau)/)))
  
    ALLOCATE(d_wgn(2,SIZE(tt%d_wgn))
    d_wgn(1,:)=REAL(tt%d_wgn)
    d_wgn(2,:)=AIMAG(tt%d_wgn)
    CALL json_print(unit,"d_wgn",RESHAPE(d_wgn,(/SHAPE(d_wgn)/))

    
    CALL json_print(unit,"nat",SIZE(tt%invsatnr))
    CALL json_print(unit,"invarop",RESHAPE(tt%invarop,(/SIZE(tt%invarop)/)))
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
    real,allocatable::rtemp(:)
    CALL json_open_class("sym",unit,iostat)
    IF (iostat.NE.0)   RETURN


    CALL json_read(unit,"l_symor",tt%l_symor)
    CALL json_read(unit,"invs",tt%invs)
    CALL json_read(unit,"invs2",tt%invs2)
    CALL json_read(unit,"zrfs",tt%zrfs)
    CALL json_read(unit,"nsymt",tt%nsymt)
    CALL json_read(unit,"nsym",tt%nsym)
    CALL json_read(unit,"nop",tt%nop)
    CALL json_read(unit,"nop2",tt%nop2)

    IF (ALLOCATED(tt%mrot)) DEALLOCATE(tt%mrot,tt%invtab,tt%tau,tt%multab)
    ALLOCATE(it(9*tt%nop))
    CALL json_read(unit,"mrot",it)
    ALLOCATE(tt%mrot(3,3,tt%nop))
    tt%mrot=RESHAPE(it,SHAPE(tt%mrot))
    DEALLOCATE(it)

    ALLOCATE(it(tt%nop**2))
    CALL json_read(unit,"multab",it)
    ALLOCATE(tt%multab(3,3,tt%nop))
    tt%multab=RESHAPE(it,SHAPE(tt%multab))
    DEALLOCATE(it)

    ALLOCATE (tt%invtab(tt%nop))
    CALL json_read(unit,"invtab",tt%invtab)

    ALLOCATE(rt(3*tt%nop))
    ALLOCATE(tt%tau(3,tt%nop))
    CALL json_read(unit,"tau",rt)
    tt%tau=RESHAPE(rt,SHAPE(tt%tau))
    DEALLOCATE(rt)

    IF (ALLOCATED(tt%d_wgn)) DEALLOCATE(tt%d_wgn)
    ALLOCATE(rt(2*147*tt%nop))
    CALL json_read(unit,"d_wgn",rt)
    ALLOCATE(tt%d_wgn(-3:3,-3:3,3,tt%nop))
    tt%d_wgn=CMPLX(RESHAPE(rt(::2),,SHAPE(tt%d_wgn)),RESHAPE(rt(2::2),,SHAPE(tt%d_wgn)))
    DEALLOCATE(rt)

    CALL json_read(unit,"nat",nat)
    IF (ALLOCATED(tt%invarop)) DEALLOCATE(tt%invarop,tt%invarind,tt%invsatnr)
    ALLOCATE(tt%invarind(nat))
    ALLOCATE(tt%invsatnr(nat))

    ALLOCATE(it(nat*tt%nop))
    ALLOCATE(tt%invarop(nat,tt%nop))
    CALL json_read(unit,"invarop",it)
    tt%invarop=RESHAPE(it,SHAPE(tt%invarop))
    DEALLOCATE(it)

    CALL json_read(unit,"invarind",tt%invarind)
    CALL json_read(unit,"invsatnr",tt%invsatnr)
    

    CALL json_close_class(unit)
    
  END SUBROUTINE read_sym

 
  SUBROUTINE read_xml_sym(tt)
    USE m_xmlIntWrapFort
    USE m_constants
    USE m_calculator
    USE m_inp_xml
    USE m_symdata , ONLY : nammap, ord2, l_c2
    IMPLICIT NONE
    CLASS(t_sym),INTENT(OUT):: tt


    LOGICAL::film
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
       ELSE
          latnam=TRIM(ADJUSTL(xmlGetAttributeValue('/fleurInput/cell/filmLattice/@latnam')))
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
         nop48 = 48
         ALLOCATE (invOps(tt%nop),multtab(tt%nop,tt%nop),optype(nop48))
         CALL check_close(tt%nop,tt%mrot,tt%tau,&
              &                      multtab,invOps,optype)

         CALL symproperties(nop48,optype,input%film,tt%nop,multtab,cell%amat,&
              &                        tt%symor,tt%mrot,tt%tau,&
              &                        invSym,tt%invs,tt%zrfs,tt%invs2,tt%nop,tt%nop2)
         IF (.not.input%film) tt%nop2=tt%nop
         IF (input%film.AND.ANY(ABS(tt%tau(:,:tt%nop))>1E-5)) &
              CALL juDFT_error("nonsymmorphic symmetries not yet implemented for films!",calledby ="r_inpXML")
         tt%invs2 = tt%invs.AND.tt%zrfs
      ELSE
         CALL judft_error("No valid specification of symmetry in inp.xml")
      END IF
      
  END SUBROUTINE read_xml_sym
  
  SUBROUTINE init_sym(sym)
    IMPLICIT NONE
    CLASS(t_sym),INTENT(INOUT):: sym
    !TODO
  END SUBROUTINE init_sym

  
END MODULE m_types_sym
