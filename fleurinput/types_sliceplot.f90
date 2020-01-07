!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_sliceplot
  USE m_judft
  USE m_types_fleurinput_base
  IMPLICIT NONE
  PRIVATE
  INTEGER,PUBLIC,PARAMETER :: PLOT_XSF_FORMAT=1
  PUBLIC :: t_sliceplot,t_plot
  TYPE ,extends(t_fleurinput_base):: t_plot
    LOGICAL :: cartesian=.false.
    LOGICAL :: twodim=.true.
    CHARACTER(len=100):: filename="default"
    integer :: grid(3)=[30,30,30]
    real    :: zero(3)=[0.,0.,0.]
    real    :: vec1(3)=[1.,0.,0.]
    real    :: vec2(3)=[0.,1.,0.]
    real    :: vec3(3)=[0.,0.,1.]
  CONTAINS
     PROCEDURE :: read_xml=>read_xml_plot
     PROCEDURE :: mpi_bc=>mpi_bc_plot
  END TYPE t_plot


  TYPE,EXTENDS(t_fleurinput_base) ::t_sliceplot
     INTEGER :: iplot=0
     LOGICAL :: slice=.FALSE.
     LOGICAL :: plpot=.FALSE.
     INTEGER :: kk=0
     INTEGER :: nnne=0
     REAL    :: e1s=0.
     REAL    :: e2s=0.
     LOGICAL :: polar=.false.
     INTEGER :: format=PLOT_XSF_FORMAT
     TYPE(t_plot),ALLOCATABLE::plot(:)
   CONTAINS
     PROCEDURE :: read_xml=>read_xml_sliceplot
     PROCEDURE :: mpi_bc=>mpi_bc_sliceplot
  END TYPE t_sliceplot

CONTAINS
  SUBROUTINE mpi_bc_sliceplot(this,mpi_comm,irank)
    USE m_mpi_bc_tool
    CLASS(t_sliceplot),INTENT(INOUT)::this
    INTEGER,INTENT(IN):: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL::irank
    INTEGER ::rank,i
    IF (PRESENT(irank)) THEN
       rank=irank
    ELSE
       rank=0
    END IF

    CALL mpi_bc(this%iplot,rank,mpi_comm)
    CALL mpi_bc(this%slice,rank,mpi_comm)
    CALL mpi_bc(this%plpot,rank,mpi_comm)
    CALL mpi_bc(this%kk,rank,mpi_comm)
    CALL mpi_bc(this%nnne,rank,mpi_comm)
    CALL mpi_bc(this%e1s,rank,mpi_comm)
    CALL mpi_bc(this%e2s,rank,mpi_comm)
    CALL mpi_bc(this%polar,rank,mpi_comm)
    CALL mpi_bc(this%format,rank,mpi_comm)


    if (allocated(this%plot)) then
      DO i=1,size(this%plot)
        call this%plot(i)%mpi_bc(mpi_comm,irank)
      ENDDO
    ENDIF
  END SUBROUTINE mpi_bc_sliceplot
  SUBROUTINE mpi_bc_plot(this,mpi_comm,irank)
    USE m_mpi_bc_tool
    CLASS(t_plot),INTENT(INOUT)::this
    INTEGER,INTENT(IN):: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL::irank
    INTEGER ::rank,i
    IF (PRESENT(irank)) THEN
       rank=irank
    ELSE
       rank=0
    END IF

    CALL mpi_bc(this%twodim,rank,mpi_comm)
    CALL mpi_bc(this%cartesian,rank,mpi_comm)
    CALL mpi_bc(this%grid(1),rank,mpi_comm)
    CALL mpi_bc(this%grid(2),rank,mpi_comm)
    CALL mpi_bc(this%grid(3),rank,mpi_comm)
    CALL mpi_bc(rank,mpi_comm,this%vec1)
    CALL mpi_bc(rank,mpi_comm,this%vec2)
    CALL mpi_bc(rank,mpi_comm,this%vec3)

  END SUBROUTINE mpi_bc_plot

  SUBROUTINE read_xml_sliceplot(this,xml)
    USE m_types_xml
    CLASS(t_sliceplot),INTENT(inOUT)::this
    TYPE(t_xml),INTENT(INOUT)::xml

    CHARACTER(len=200)::xpatha
    INTEGER::numberNodes,i

    xPathA = '/fleurInput/output'
    numberNodes = xml%GetNumberOfNodes(xPathA)

    IF (numberNodes.EQ.1) THEN
       this%slice = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@slice'))
    ENDIF
    xPathA = '/fleurInput/output/plotting'
    numberNodes = xml%GetNumberOfNodes(xPathA)

    IF (numberNodes.EQ.1) THEN
       this%iplot = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@iplot'))
       this%polar = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@polar'))
       xPathA = '/fleurInput/output/plotting/plot'
       numberNodes = xml%GetNumberOfNodes(xPathA)
       allocate(this%plot(numberNodes))
       do i = 1, numberNodes
         write(xPathA,'(a,i0,a)') '/fleurInput/output/plotting/plot[',i,']'
         call xml%set_basepath(xPathA)
         call this%plot(i)%read_xml(xml)
         call xml%set_basepath('')
       end do
    END IF

    xPathA = '/fleurInput/output/chargeDensitySlicing'
    numberNodes = xml%GetNumberOfNodes(xPathA)

    IF ((this%slice).AND.(numberNodes.EQ.0)) THEN
       CALL juDFT_error("slice is true but chargeDensitySlicing parameters are not set!")
    END IF

    IF (numberNodes.EQ.1) THEN
       this%kk = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@numkpt'))
       this%e1s = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@minEigenval'))
       this%e2s = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@maxEigenval'))
       this%nnne = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nnne'))
    END IF
  END SUBROUTINE read_xml_sliceplot
  SUBROUTINE read_xml_plot(this,xml)
    USE m_types_xml
    CLASS(t_plot),INTENT(inOUT)::this
    TYPE(t_xml),INTENT(INOUT)::xml

    CHARACTER(len=200)::line
    real,allocatable :: x(:)
    INTEGER::numberNodes,i

    this%cartesian = evaluateFirstBoolOnly(xml%GetAttributeValue('/@cartesian'))
    this%twodim     = evaluateFirstBoolOnly(xml%GetAttributeValue('/@2D'))

    line=xml%GetAttributeValue('/@grid')
    call evaluateList(x,line)
    if (size(x)/=3) call judft_error("Wrong number of grid points")
    this%grid=x

    line=xml%GetAttributeValue('/@vec1')
    call evaluateList(x,line)
    if (size(x)/=3) call judft_error("Wrong number of coordinates for vec1")
    this%vec1=x
    line=xml%GetAttributeValue('/@vec2')
    call evaluateList(x,line)
    if (size(x)/=3) call judft_error("Wrong number of coordinates for vec2")
    this%vec2=x
    line=xml%GetAttributeValue('/@vec3')
    call evaluateList(x,line)
    if (size(x)/=3) call judft_error("Wrong number of coordinates for vec3")
    this%vec3=x

    this%filename=xml%GetAttributeValue('/@file')

  END SUBROUTINE read_xml_plot

END MODULE m_types_sliceplot
