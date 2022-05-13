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
  INTEGER :: nvac=2
  INTEGER :: nvacd=2
  INTEGER,allocatable :: mrot2(:,:)
  REAL,allocatable    :: tau2(:)
  REAL :: delz=0.1
  REAL :: dvac=0.0
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
 CALL mpi_bc(this%nvac,rank,mpi_comm)
 CALL mpi_bc(this%nvacd,rank,mpi_comm)
 CALL mpi_bc(this%delz,rank,mpi_comm)
 CALL mpi_bc(this%dvac,rank,mpi_comm)
 call mpi_bc(this%tau2,rank,mpi_comm)
 call mpi_bc(this%mrot2,rank,mpi_comm)

END SUBROUTINE mpi_bc_vacuum
SUBROUTINE read_xml(this,xml)
 USE m_types_xml
 CLASS(t_vacuum),INTENT(INOUT)::this
 TYPE(t_xml),INTENT(INOUT)::xml
 CHARACTER(len=100)::xpatha


 IF (xml%GetNumberOfNodes('/fleurInput/cell/filmLattice')==1) THEN
    this%dvac = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/cell/filmLattice/@dVac'))
    !this%dtild = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/cell/filmLattice/@dTilda'))
  endif
  END SUBROUTINE read_xml

SUBROUTINE vacuum_init(this,sym)
 USE m_types_sym
 CLASS(t_vacuum),INTENT(INOUT)::this
 TYPE(t_sym),INTENT(IN)::sym
 
 allocate(this%mrot2(2,2),this%tau2(2))  
 if (sym%nop>sym%nop2) THEN
   this%nvac=1
   this%mrot2(1:2,1:2) = sym%mrot(1:2,1:2,sym%nop2+1)
   this%tau2(1:2) = sym%tau(1:2,sym%invtab(sym%nop2+1))
 else
   this%mrot2(1:2,1:2) = sym%mrot(1:2,1:2,1)
   this%tau2(1:2) = sym%tau(1:2,1)
 endif   
   
END SUBROUTINE vacuum_init



END MODULE m_types_vacuum
