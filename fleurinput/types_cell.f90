!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_cell
  USE m_judft
  USE m_types_fleurinput_base
  IMPLICIT NONE
  PRIVATE
  !> This type contains the basic information on the lattice-cell of the calculation
  !> To use it, you basically only have to provide cell%amat (and cell%z1 for films)
  !> and call its init routine.
  TYPE,EXTENDS(t_fleurinput_base):: t_cell
     !vol of dtilde box
     REAL::omtil
     !2D area
     REAL::area
     !bravais matrix
     REAL::amat(3, 3)
     !rez. bravais matrx
     REAL::bmat(3, 3)
     !square of bbmat
     REAL::bbmat(3, 3),aamat(3,3)
     !d-value
     REAL::z1=0.0
     !volume of cell
     REAL::vol
     !volume of interstitial
     REAL::volint
   CONTAINS
     PROCEDURE :: init
     PROCEDURE :: read_xml=>read_xml_cell
     PROCEDURE :: mpi_bc=>mpi_bc_cell
  END TYPE t_cell
  PUBLIC t_cell
CONTAINS
  subroutine mpi_bc_cell(this,mpi_comm,irank)
    use m_mpi_bc_tool
    class(t_cell),INTENT(INOUT)::this
    integer,INTENT(IN):: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL::irank
    INTEGER ::rank
    if (present(irank)) THEN
       rank=irank
    else
       rank=0
    end if

    call mpi_bc(this%omtil,rank,mpi_comm)
    call mpi_bc(this%area,rank,mpi_comm)
    call mpi_bc(rank,mpi_comm,this%amat)
    call mpi_bc(rank,mpi_comm,this%bmat)
    call mpi_bc(rank,mpi_comm,this%bbmat)
    call mpi_bc(rank,mpi_comm,this%aamat)
    call mpi_bc(this%z1,rank,mpi_comm)
    call mpi_bc(this%vol,rank,mpi_comm)
    call mpi_bc(this%volint,rank,mpi_comm)
  end subroutine mpi_bc_cell

  SUBROUTINE init(cell,volmts)
    !initialize cell, only input is cell%amat and cell%z1 in case of a film
    USE m_constants,ONLY:tpi_const
    CLASS (t_cell),INTENT(INOUT):: cell
    real,intent(in):: volmts !Volume of all MT-spheres

    CALL inv3(cell%amat,cell%bmat,cell%omtil)
    IF (cell%omtil<0) CALL judft_warn("Negative volume! You are using a left-handed coordinate system")
    cell%omtil=ABS(cell%omtil)

    cell%bmat=tpi_const*cell%bmat
    IF (cell%z1>0) THEN
       cell%vol = (cell%omtil/cell%amat(3,3))*cell%z1*2
       cell%area = cell%omtil/cell%amat(3,3)
    ELSE
       cell%vol = cell%omtil
       cell%area =ABS(cell%amat(1,1)*cell%amat(2,2)-cell%amat(1,2)*cell%amat(2,1))
       IF (cell%area < 1.0e-7) THEN
          cell%area = 1.
          CALL juDFT_warn("area = 0",calledby ="types_cell")
       END IF
     END IF

     cell%bbmat=matmul(cell%bmat,transpose(cell%bmat))
     cell%aamat=matmul(transpose(cell%amat),cell%amat)
     cell%volint = cell%vol
  cell%volint = cell%volint-volmts
   CONTAINS
     !This is a copy of the code in math/inv3
     !Put here to make library independent
     SUBROUTINE inv3(a,b,d)
       IMPLICIT NONE
       !     ..
       !     .. Arguments ..
       REAL, INTENT (IN)  :: a(3,3)
       REAL, INTENT (OUT) :: b(3,3)  ! inverse matrix
       REAL, INTENT (OUT) :: d       ! determinant
       !     ..
       d = a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1) + &
            a(2,1)*a(3,2)*a(1,3) - a(1,3)*a(2,2)*a(3,1) - &
            a(2,3)*a(3,2)*a(1,1) - a(2,1)*a(1,2)*a(3,3)
       b(1,1) = (a(2,2)*a(3,3)-a(2,3)*a(3,2))/d
       b(1,2) = (a(1,3)*a(3,2)-a(1,2)*a(3,3))/d
       b(1,3) = (a(1,2)*a(2,3)-a(2,2)*a(1,3))/d
       b(2,1) = (a(2,3)*a(3,1)-a(2,1)*a(3,3))/d
       b(2,2) = (a(1,1)*a(3,3)-a(3,1)*a(1,3))/d
       b(2,3) = (a(1,3)*a(2,1)-a(1,1)*a(2,3))/d
       b(3,1) = (a(2,1)*a(3,2)-a(2,2)*a(3,1))/d
       b(3,2) = (a(1,2)*a(3,1)-a(1,1)*a(3,2))/d
       b(3,3) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))/d

     END SUBROUTINE inv3

   END SUBROUTINE init

   SUBROUTINE  read_xml_cell(this,xml)
     use m_types_xml
     class(t_cell),intent(INout)::this
     type(t_xml),intent(in)   ::xml

     ! Read in lattice parameters
     character(len=200)::valueString,path
     REAL:: scale,dvac,dtild

     if (xml%GetNumberOfNodes('/fleurInput/cell/filmLattice')==1) then
        path= '/fleurInput/cell/filmLattice'
        this%z1=evaluateFirstOnly(xml%GetAttributeValue(trim(path)//'/@dVac'))/2
        dvac=this%z1*2
        dtild=evaluateFirstOnly(xml%GetAttributeValue(trim(path)//'/@dTilda'))
     else
        dvac=0.0
        path = '/fleurInput/cell/bulkLattice'
     endif

     scale=evaluateFirstOnly(xml%GetAttributeValue(trim(path)//'/@scale'))
     path=trim(path)//'/bravaisMatrix'
     valueString = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(path))//'/row-1')))
     this%amat(1,1) = evaluateFirst(valueString)
     this%amat(2,1) = evaluateFirst(valueString)
     this%amat(3,1) = evaluateFirst(valueString)
     valueString = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(path))//'/row-2')))
     this%amat(1,2) = evaluateFirst(valueString)
     this%amat(2,2) = evaluateFirst(valueString)
     this%amat(3,2) = evaluateFirst(valueString)
     valueString = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(path))//'/row-3')))
     this%amat(1,3) = evaluateFirst(valueString)
     this%amat(2,3) = evaluateFirst(valueString)
     this%amat(3,3) = evaluateFirst(valueString)

     IF (dvac>0) THEN
        this%amat(3,3)=dtild
     ENDIF
     this%amat=this%amat*scale

   END SUBROUTINE read_xml_cell

 END MODULE m_types_cell
