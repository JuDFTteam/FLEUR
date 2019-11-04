!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_mae
  USE m_types
  USE m_types_forcetheo
  USE m_judft
  IMPLICIT NONE
  PRIVATE
  TYPE,EXTENDS(t_forcetheo) :: t_forcetheo_mae
     INTEGER :: directions_done
     REAL,ALLOCATABLE:: theta(:)
     REAL,ALLOCATABLE:: phi(:)
     REAL,ALLOCATABLE:: evsum(:)
   CONTAINS
     PROCEDURE :: start   =>mae_start
     PROCEDURE :: next_job=>mae_next_job 
     PROCEDURE :: eval    =>mae_eval
     PROCEDURE :: postprocess => mae_postprocess
     PROCEDURE :: init   => mae_init !not overloaded
     PROCEDURE :: dist   => mae_dist !not overloaded
  END TYPE t_forcetheo_mae
  PUBLIC t_forcetheo_mae
CONTAINS
 

  SUBROUTINE mae_init(this,theta,phi,cell,sym)
    USE m_calculator
    USE m_socsym
    USE m_types
    IMPLICIT NONE
    CLASS(t_forcetheo_mae),INTENT(INOUT):: this
    TYPE(t_cell),INTENT(IN)             :: cell
    TYPE(t_sym),INTENT(IN)              :: sym
    REAL,INTENT(in)                     :: theta(:),phi(:)

    INTEGER::n
    LOGICAL::error(sym%nop)
    
    this%phi=phi
    this%theta=theta

    IF (SIZE(this%phi).NE.SIZE(this%theta)) CALL &
         judft_error("Lists for theta/phi must have the same length in MAE force theorem calculations")
    DO n=1,SIZE(this%phi)
       CALL soc_sym(sym%nop,sym%mrot,this%theta(n),this%phi(n),cell%amat,error)
       IF (ANY(error)) CALL judft_error("Force theory choice of SOC-SQA breaks symmetry")
    END DO
    ALLOCATE(this%evsum(SIZE(this%phi)))
    this%evsum=0
  END SUBROUTINE mae_init
    

  SUBROUTINE mae_start(this,potden,l_io)
    USE m_types_potden
    IMPLICIT NONE
    CLASS(t_forcetheo_mae),INTENT(INOUT):: this
    TYPE(t_potden) ,INTENT(INOUT)       :: potden
    LOGICAL,INTENT(IN)                  :: l_io
    this%directions_done=0
    CALL this%t_forcetheo%start(potden,l_io) !call routine of basis type
  END SUBROUTINE  mae_start


  LOGICAL FUNCTION mae_next_job(this,lastiter,atoms,noco)
    USE m_types_setup
    USE m_xmlOutput
    USE m_constants
    IMPLICIT NONE
    CLASS(t_forcetheo_mae),INTENT(INOUT):: this
    LOGICAL,INTENT(IN)                  :: lastiter
    TYPE(t_atoms),INTENT(IN)            :: atoms
    !Stuff that might be modified...
    TYPE(t_noco),INTENT(INOUT) :: noco
       IF (.NOT.lastiter) THEN
          mae_next_job=this%t_forcetheo%next_job(lastiter,atoms,noco)
          RETURN
       ENDIF
       !OK, now we start the MAE-loop
       this%directions_done=this%directions_done+1
       mae_next_job=(this%directions_done<=SIZE(this%phi)) !still angles to do
       IF (.NOT.mae_next_job) RETURN

       noco%theta=this%theta(this%directions_done)
       noco%phi=this%phi(this%directions_done)
       noco%l_soc=.true.
       IF (this%directions_done.NE.1.AND.this%l_io) CALL closeXMLElement('Forcetheorem_Loop_MAE')
       IF (this%l_io) CALL openXMLElementPoly('Forcetheorem_Loop_MAE',(/'No'/),(/this%directions_done/))
  END FUNCTION mae_next_job

  FUNCTION mae_eval(this,eig_id,atoms,kpts,sym,&
       cell,noco, input,mpi, oneD,enpara,v,results)RESULT(skip)
    USE m_types
    IMPLICIT NONE
    CLASS(t_forcetheo_mae),INTENT(INOUT):: this
    LOGICAL :: skip
    !Stuff that might be used...
    TYPE(t_mpi),INTENT(IN)         :: mpi
    
    TYPE(t_oneD),INTENT(IN)        :: oneD
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_noco),INTENT(IN)        :: noco
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_cell),INTENT(IN)        :: cell
    TYPE(t_kpts),INTENT(IN)        :: kpts
    TYPE(t_atoms),INTENT(IN)       :: atoms
    TYPE(t_enpara),INTENT(IN)      :: enpara
    TYPE(t_potden),INTENT(IN)      :: v
    TYPE(t_results),INTENT(IN)     :: results
    INTEGER,INTENT(IN)             :: eig_id
    IF (this%directions_done==0) THEN
       skip=.FALSE.
       RETURN
    ENDIF
    this%evsum(this%directions_done)=results%seigv/2.0 
    skip=.TRUE.
  END FUNCTION  mae_eval

  SUBROUTINE mae_postprocess(this)
    USE m_xmlOutput
    IMPLICIT NONE
    CLASS(t_forcetheo_mae),INTENT(INOUT):: this

    !Locals
    INTEGER:: n
    CHARACTER(LEN=12):: attributes(3)
    IF (this%directions_done==0) THEN
       RETURN
    ENDIF
    
    IF (this%l_io) THEN
       !Now output the results
       CALL closeXMLElement('Forcetheorem_Loop_MAE')
       CALL openXMLElementPoly('Forcetheorem_MAE',(/'Angles'/),(/SIZE(this%evsum)/))
       DO n=1,SIZE(this%evsum)
          WRITE(attributes(1),'(f12.7)') this%theta(n)
          WRITE(attributes(2),'(f12.7)') this%phi(n)
          WRITE(attributes(3),'(f12.7)') this%evsum(n)     
          CALL writeXMLElementForm('Angle',(/'theta ','phi   ','ev-sum'/),attributes,&
               RESHAPE((/5,3,6,12,12,12/),(/3,2/)))
       END DO
       CALL closeXMLElement('Forcetheorem_MAE')
    ENDIF
    CALL judft_end("Forcetheorem MAE")
  END SUBROUTINE mae_postprocess

  SUBROUTINE mae_dist(this,mpi)
    USE m_types_mpi
    IMPLICIT NONE
    CLASS(t_forcetheo_mae),INTENT(INOUT):: this
    TYPE(t_mpi),INTENT(in):: mpi

    INTEGER:: i,ierr
#ifdef CPP_MPI    
    INCLUDE 'mpif.h'
    IF (mpi%irank==0) i=SIZE(this%theta)
    call MPI_BCAST(i,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    IF (mpi%irank.NE.0) ALLOCATE(this%phi(i),this%theta(i),this%evsum(i));this%evsum=0.0
    CALL MPI_BCAST(this%phi,i,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(this%theta,i,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
#endif    
  END SUBROUTINE mae_dist
END MODULE m_types_mae
