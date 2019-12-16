!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_field
  !*************************************************************
  !     This module contains definitions for electric and magnetic field -types
  !*************************************************************
  USE m_juDFT
  USE m_types_fleurinput_base
  IMPLICIT NONE
  PRIVATE
  TYPE:: t_efield
     REAL    :: zsigma  = 10.0  ! Distance to the charged plates
     REAL    :: sigma   ! charge at the plates
     REAL    :: sig_b(2)=  0.0  ! Extra charge for the top/bottom plate
     COMPLEX :: vslope  =  0.0  ! Dirichlet bnd. cond.: Slope
     REAL,    ALLOCATABLE :: sigEF(:,:,:) ! (nx, ny, nvac)
     COMPLEX, ALLOCATABLE :: rhoEF(:,:)   ! (g_||, nvac)
     COMPLEX, ALLOCATABLE :: C1(:), C2(:) ! Coeff. for Dirichlet bnd.cond.
     LOGICAL :: l_segmented = .FALSE.
     LOGICAL :: plot_charge = .FALSE. ! Plot charge as inputted
     LOGICAL :: plot_rho    = .FALSE. ! Plot Fourier-transformed charge
     LOGICAL :: autocomp    = .TRUE.  ! Auto-compensate film charge
     LOGICAL :: dirichlet = .FALSE. ! Dirichlet vs. Neumann boundary cond.
     LOGICAL :: l_dirichlet_coeff = .FALSE. ! For MPI, true if C1/C2 set
     LOGICAL :: l_eV =.FALSE. !Input in eV
     CHARACTER(len=50),ALLOCATABLE :: shapes(:)
  END TYPE t_efield

  TYPE,EXTENDS(t_fleurinput_base):: t_field
     TYPE(t_efield)   :: efield
     LOGICAL          :: l_b_field=.FALSE.
     REAL             :: b_field
     REAL,ALLOCATABLE :: b_field_mt(:)
   CONTAINS
     PROCEDURE :: init=>init_field
     PROCEDURE :: read_xml=>read_xml_field
     PROCEDURE :: mpi_bc=>mpi_bc_field
  END TYPE t_field

  PUBLIC t_field,t_efield
CONTAINS

  SUBROUTINE mpi_bc_field(this,mpi_comm,irank)
    USE m_mpi_bc_tool
    CLASS(t_field),INTENT(INOUT)::this
    INTEGER,INTENT(IN):: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL::irank
    INTEGER ::rank
    IF (PRESENT(irank)) THEN
       rank=irank
    ELSE
       rank=0
    END IF

    CALL mpi_bc(this%l_b_field,rank,mpi_comm)
    CALL mpi_bc(this%b_field,rank,mpi_comm)
    CALL mpi_bc(this%b_field_mt,rank,mpi_comm)

    CALL mpi_bc(this%efield%zsigma  ,rank,mpi_comm)
    CALL mpi_bc(this%efield%sigma   ,rank,mpi_comm)
    CALL mpi_bc(this%efield%sig_b(1),rank,mpi_comm)
    CALL mpi_bc(this%efield%sig_b(2),rank,mpi_comm)
    CALL mpi_bc(this%efield%vslope  ,rank,mpi_comm)
    CALL mpi_bc(this%efield%sigEF ,rank,mpi_comm)
    CALL mpi_bc(this%efield%rhoEF   ,rank,mpi_comm)
    CALL MPI_BC(THIS%EFIELD%C1,RANK,MPI_COMM)
    CALL mpi_bc(this%efield%C2 ,rank,mpi_comm)
    CALL mpi_bc(this%efield%l_segmented ,rank,mpi_comm)
    CALL mpi_bc(this%efield%plot_charge ,rank,mpi_comm)
    CALL mpi_bc(this%efield%plot_rho    ,rank,mpi_comm)
    CALL mpi_bc(this%efield%autocomp    ,rank,mpi_comm)
    CALL mpi_bc(this%efield%dirichlet ,rank,mpi_comm)
    CALL mpi_bc(this%efield%l_dirichlet_coeff ,rank,mpi_comm)
    CALL mpi_bc(this%efield%l_eV ,rank,mpi_comm)



  END SUBROUTINE mpi_bc_field

  SUBROUTINE init_field(this)
    !USE m_types_setup
    IMPLICIT NONE
    CLASS(t_field),INTENT(INOUT)::this
    !TYPE(t_input),INTENT(INOUT) ::input
    !input%sigma => sigma
    !this%efield%sigma=>sigma
  END SUBROUTINE init_field

  SUBROUTINE read_xml_field(this,xml)
    USE m_types_xml
    CLASS(t_field),INTENT(INOUT)::this
    TYPE(t_xml),INTENT(IN)::xml

    CHARACTER(len=100)::xpatha,xpathb
    INTEGER:: numberNodes,i
    
    xPathA = '/fleurInput/calculationSetup/fields'
    numberNodes = xml%GetNumberOfNodes(xPathA)
    IF (numberNodes.EQ.1) THEN
       IF (xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/@b_field')>0) THEN
          this%b_field=evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'//@b_field'))
          this%l_b_field=.TRUE.
       ENDIF
       this%efield%zsigma = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@zsigma'))
       this%efield%sig_b(1) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@sig_b_1'))
       this%efield%sig_b(2) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@sig_b_2'))
       this%efield%plot_charge = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@plot_charge'))
       this%efield%plot_rho = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@plot_rho'))
       this%efield%autocomp = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@autocomp'))
       this%efield%dirichlet = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@dirichlet'))
       this%efield%l_eV = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@eV'))
       
       numberNodes=xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/shape')
       ALLOCATE(this%efield%shapes(numberNodes))
       DO i=1,numberNodes
          WRITE(xPathB,"(a,a,i0,a)") TRIM(ADJUSTL(xpathA)),'/shape[',i,']'
          this%efield%shapes(i)=TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB)))))
       ENDDO
    ELSE
       ALLOCATE(this%efield%shapes(0))
    END IF
  END SUBROUTINE read_xml_field
END MODULE m_types_field
