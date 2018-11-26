!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_ssdisp

  USE m_types
  USE m_types_forcetheo
  USE m_judft
  TYPE,EXTENDS(t_forcetheo) :: t_forcetheo_ssdisp
     INTEGER :: q_done
     REAL,ALLOCATABLE:: qvec(:,:)
     REAL,ALLOCATABLE:: evsum(:)
   CONTAINS
     PROCEDURE :: start   =>ssdisp_start
     PROCEDURE :: next_job=>ssdisp_next_job 
     PROCEDURE :: eval    =>ssdisp_eval
     PROCEDURE :: postprocess => ssdisp_postprocess
     PROCEDURE :: init   => ssdisp_init !not overloaded
     PROCEDURE :: dist   => ssdisp_dist !not overloaded
  END TYPE t_forcetheo_ssdisp

CONTAINS

  SUBROUTINE ssdisp_init(this,q)
    USE m_calculator
    USE m_constants
    IMPLICIT NONE
    CLASS(t_forcetheo_ssdisp),INTENT(INOUT):: this
    REAL,INTENT(in)                     :: q(:,:)
    
    ALLOCATE(this%qvec(3,SIZE(q,2)))
    this%qvec=q
    
    ALLOCATE(this%evsum(SIZE(q,2)))
    this%evsum=0
  END SUBROUTINE ssdisp_init

  SUBROUTINE ssdisp_start(this,potden)
    USE m_types_potden
    IMPLICIT NONE
    CLASS(t_forcetheo_ssdisp),INTENT(INOUT):: this
    TYPE(t_potden) ,INTENT(INOUT)          :: potden
    this%q_done=0
    CALL this%t_forcetheo%start(potden) !call routine of basis type

    IF (SIZE(potden%pw,2)<2) RETURN
    !Average out magnetic part of potential/charge in INT+Vacuum
    potden%pw(:,1)=(potden%pw(:,1)+potden%pw(:,2))/2.0
    potden%pw(:,2)=potden%pw(:,1)
    IF (SIZE(potden%pw,2)==3) potden%pw(:,3)=0.0
    
    potden%vacz(:,:,1)=(potden%vacz(:,:,1)+potden%vacz(:,:,2))/2.0
    potden%vacxy(:,:,:,1)=(potden%vacxy(:,:,:,1)+potden%vacxy(:,:,:,2))/2.0
    potden%vacz(:,:,2)=potden%vacz(:,:,1)
    potden%vacxy(:,:,:,2)=potden%vacxy(:,:,:,1)
    
  END SUBROUTINE  ssdisp_start

  LOGICAL FUNCTION ssdisp_next_job(this,lastiter,noco)
    USE m_types_setup
    USE m_xmlOutput
    IMPLICIT NONE
    CLASS(t_forcetheo_ssdisp),INTENT(INOUT):: this
    LOGICAL,INTENT(IN)                  :: lastiter
    !Stuff that might be modified...
    TYPE(t_noco),INTENT(INOUT) :: noco
    IF (.NOT.lastiter) THEN
       ssdisp_next_job=this%t_forcetheo%next_job(lastiter,noco)
       RETURN
    ENDIF
    !OK, now we start the SSDISP-loop
    this%q_done=this%q_done+1
    ssdisp_next_job=(this%q_done<=SIZE(this%qvec,2)) !still q-vectors to do
    IF (.NOT.ssdisp_next_job) RETURN
    
    !Now modify the noco-file
    noco%qss=this%qvec(:,this%q_done)
    IF (this%q_done.NE.1) CALL closeXMLElement('Forcetheorem_Loop_SSDISP')
    CALL openXMLElementPoly('Forcetheorem_Loop_SSDISP',(/'Q-vec:'/),(/this%q_done/))
  END FUNCTION ssdisp_next_job

  SUBROUTINE ssdisp_postprocess(this)
    USE m_xmlOutput
    IMPLICIT NONE
    CLASS(t_forcetheo_ssdisp),INTENT(INOUT):: this

    !Locals
    INTEGER:: n,q
    CHARACTER(LEN=12):: attributes(4)
    IF (this%q_done==0) RETURN
    !Now output the results
    CALL closeXMLElement('Forcetheorem_Loop_SSDISP')
    CALL openXMLElementPoly('Forcetheorem_SSDISP',(/'qvectors'/),(/SIZE(this%evsum)/))
    DO q=1,SIZE(this%evsum)
       WRITE(attributes(1),'(i5)') q
       WRITE(attributes(2),'(f12.7)') this%evsum(q) 
       CALL writeXMLElementForm('Entry',(/'q     ','ev-sum'/),attributes(1:2),&
            RESHAPE((/1,6,5,12/),(/2,2/)))
    ENDDO
    CALL closeXMLElement('Forcetheorem_SSDISP')
    CALL judft_end("Forcetheorem:SpinSpiralDispersion")
  END SUBROUTINE ssdisp_postprocess

  SUBROUTINE ssdisp_dist(this,mpi)
    USE m_types_mpi
    IMPLICIT NONE
    CLASS(t_forcetheo_ssdisp),INTENT(INOUT):: this
    TYPE(t_mpi),INTENT(in):: mpi
    
    INTEGER:: q,ierr
#ifdef CPP_MPI    
    INCLUDE 'mpif.h'
    IF (mpi%irank==0) q=SIZE(this%qvec,2)
    CALL MPI_BCAST(q,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    IF (mpi%irank.NE.0) ALLOCATE(this%qvec(3,q),this%evsum(q));this%evsum=0.0
    CALL MPI_BCAST(this%qvec,3*q,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
#endif    
  END SUBROUTINE ssdisp_dist

  FUNCTION ssdisp_eval(this,eig_id,DIMENSION,atoms,kpts,sym,&
       cell,noco, input,mpi, oneD,enpara,v,results)RESULT(skip)
     USE m_types
     USE m_ssomat
    IMPLICIT NONE
    CLASS(t_forcetheo_ssdisp),INTENT(INOUT):: this
    LOGICAL :: skip
    !Stuff that might be used...
    TYPE(t_mpi),INTENT(IN)         :: mpi
    TYPE(t_dimension),INTENT(IN)   :: dimension
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
    skip=.FALSE.
    IF (this%q_done==0) RETURN
  
    this%evsum(this%q_done)=results%seigv
    skip=.TRUE.
  END FUNCTION  ssdisp_eval

  
END MODULE m_types_ssdisp
