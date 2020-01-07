!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_vacuum
  USE m_juDFT
  USE m_types_fleurinput_base
  IMPLICIT NONE
  PRIVATE
  PUBLIC:: t_vacuum
  TYPE,EXTENDS(t_fleurinput_base):: t_vacuum
     !Stuff for the vacuum
  INTEGER ::nmz=250
  INTEGER ::nmzd=250
  INTEGER ::nmzxy=100
  INTEGER ::nmzxyd=100
  INTEGER :: layerd=1
  INTEGER :: layers=0
  INTEGER :: nvac=2
  INTEGER :: nvacd=2
  REAL :: delz=0.1
  REAL :: dvac=0.0
  INTEGER::nstars=0
  INTEGER:: nstm=0
  REAL :: tworkf=0.0
  REAL :: locx(2)=[0.,0.]
  REAL :: locy(2)=[0.,0.]
  LOGICAL ::starcoeff=.FALSE.
  INTEGER, ALLOCATABLE :: izlay(:, :)
CONTAINS
  PROCEDURE :: read_xml
  PROCEDURE :: mpi_bc => mpi_bc_vacuum
  PROCEDURE :: init =>vacuum_init
END TYPE t_vacuum
CONTAINS

SUBROUTINE mpi_bc_vacuum(this,mpi_comm,irank)
 USE m_mpi_bc_tool
 CLASS(t_vacuum),INTENT(INOUT)::this
 INTEGER,INTENT(IN):: mpi_comm
 INTEGER,INTENT(IN),OPTIONAL::irank
 INTEGER ::rank
 IF (PRESENT(irank)) THEN
    rank=irank
 ELSE
    rank=0
 END IF
 CALL mpi_bc(this%nmz,rank,mpi_comm)
 CALL mpi_bc(this%nmzd,rank,mpi_comm)
 CALL mpi_bc(this%nmzxy,rank,mpi_comm)
 CALL mpi_bc(this%nmzxyd,rank,mpi_comm)
 CALL mpi_bc(this%layerd,rank,mpi_comm)
 CALL mpi_bc(this%layers,rank,mpi_comm)
 CALL mpi_bc(this%nvac,rank,mpi_comm)
 CALL mpi_bc(this%nvacd,rank,mpi_comm)
 CALL mpi_bc(this%delz,rank,mpi_comm)
 CALL mpi_bc(this%dvac,rank,mpi_comm)
 CALL mpi_bc(this%nstars,rank,mpi_comm)
 CALL mpi_bc(this%nstm,rank,mpi_comm)
 CALL mpi_bc(this%tworkf,rank,mpi_comm)
 CALL mpi_bc(this%locx(1),rank,mpi_comm)
 CALL mpi_bc(this%locy(1),rank,mpi_comm)
 CALL mpi_bc(this%locx(2),rank,mpi_comm)
 CALL mpi_bc(this%locy(2),rank,mpi_comm)
 CALL mpi_bc(this%starcoeff,rank,mpi_comm)
 CALL mpi_bc(this%izlay,rank,mpi_comm)

END SUBROUTINE mpi_bc_vacuum
SUBROUTINE read_xml(this,xml)
 USE m_types_xml
 CLASS(t_vacuum),INTENT(INOUT)::this
 TYPE(t_xml),INTENT(INOUT)::xml
 CHARACTER(len=100)::xpatha


 IF (xml%GetNumberOfNodes('/fleurInput/cell/filmLattice')==1) THEN
    this%dvac = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/cell/filmLattice/@dVac'))
    !this%dtild = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/cell/filmLattice/@dTilda'))

    xPathA = '/fleurInput/output/vacuumDOS'
    IF (xml%GetNumberOfNodes(xpathA)==1) THEN
       this%layers = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@layers'))
       this%starcoeff = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@star'))
       this%nstars = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nstars'))
       this%locx(1) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locx1'))
       this%locx(2) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locx2'))
       this%locy(1) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locy1'))
       this%locy(2) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locy2'))
       this%nstm = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nstm'))
       this%tworkf = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@tworkf'))
       PRINT *,'Reading of layers not implemented '
    END IF
    this%layerd = this%layers
    ALLOCATE(this%izlay(this%layerd,2))
 ELSE
    this%layerd = this%layers
    ALLOCATE(this%izlay(this%layerd,2))
 ENDIF
END SUBROUTINE read_xml

SUBROUTINE vacuum_init(this,sym)
 USE m_types_sym
 CLASS(t_vacuum),INTENT(INOUT)::this
 TYPE(t_sym),INTENT(IN)::sym

 IF (sym%invs.OR.sym%zrfs) this%nvac=1
END SUBROUTINE vacuum_init

END MODULE m_types_vacuum
