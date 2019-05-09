!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_dmi

  USE m_types
  USE m_types_forcetheo
  USE m_judft
  TYPE,EXTENDS(t_forcetheo) :: t_forcetheo_dmi
     INTEGER :: q_done
     REAL,ALLOCATABLE:: qvec(:,:)
     REAL,ALLOCATABLE:: theta(:)
     REAL,ALLOCATABLE:: phi(:)
     REAL,ALLOCATABLE:: evsum(:,:)
   CONTAINS
     PROCEDURE :: start   =>dmi_start
     PROCEDURE :: next_job=>dmi_next_job 
     PROCEDURE :: eval    =>dmi_eval
     PROCEDURE :: postprocess => dmi_postprocess
     PROCEDURE :: init   => dmi_init !not overloaded
     PROCEDURE :: dist   => dmi_dist !not overloaded
  END TYPE t_forcetheo_dmi

CONTAINS

  SUBROUTINE dmi_init(this,q,theta_s,phi_s)
    USE m_calculator
    USE m_constants
    IMPLICIT NONE
    CLASS(t_forcetheo_dmi),INTENT(INOUT):: this
    REAL,INTENT(in)                     :: q(:,:)
    CHARACTER(len=*),INTENT(INOUT)      :: theta_s,phi_s

    CALL evaluateList(this%theta,theta_s)
    CALL evaluateList(this%phi,phi_s)

    IF (SIZE(this%phi).NE.SIZE(this%theta)) CALL &
         judft_error("Lists for theta/phi must have the same length in DMI force theorem calculations")

    ! use same definition of rotation angles as in noco-routines 
    this%theta=-this%theta
    this%phi=this%phi+pi_const
    
    
    ALLOCATE(this%qvec(3,SIZE(q,2)))
    this%qvec=q
    
    ALLOCATE(this%evsum(0:SIZE(this%phi),SIZE(q,2)))
    this%evsum=0
  END SUBROUTINE dmi_init

  SUBROUTINE dmi_start(this,potden,l_io)
    USE m_types_potden
    IMPLICIT NONE
    CLASS(t_forcetheo_dmi),INTENT(INOUT):: this
    TYPE(t_potden) ,INTENT(INOUT)       :: potden
    LOGICAL,INTENT(IN)                  :: l_io
    this%q_done=0
    CALL this%t_forcetheo%start(potden,l_io) !call routine of basis type
  END SUBROUTINE  dmi_start

  LOGICAL FUNCTION dmi_next_job(this,lastiter,noco)
    USE m_types_setup
    USE m_xmlOutput
    IMPLICIT NONE
    CLASS(t_forcetheo_dmi),INTENT(INOUT):: this
    LOGICAL,INTENT(IN)                  :: lastiter
    !Stuff that might be modified...
    TYPE(t_noco),INTENT(INOUT) :: noco
    INTEGER                    :: itype
    IF (.NOT.lastiter) THEN
       dmi_next_job=this%t_forcetheo%next_job(lastiter,noco)
       RETURN
    ENDIF
    !OK, now we start the DMI-loop
    this%q_done=this%q_done+1
    dmi_next_job=(this%q_done<=SIZE(this%qvec,2)) !still q-vectors to do
    IF (.NOT.dmi_next_job) RETURN
    
    !Now modify the noco-file
    noco%qss=this%qvec(:,this%q_done)
    !Modify the alpha-angles
    DO iType = 1,atoms%ntype
       noco%alph(iType) = noco%alphInit(iType) + tpi_const*dot_PRODUCT(noco%qss,atoms%taual(:,SUM(atoms%neq(:itype-1))+1))
    END DO
    IF (.NOT.this%l_io) RETURN
  
    IF (this%q_done.NE.1) CALL closeXMLElement('Forcetheorem_Loop_DMI')
    CALL openXMLElementPoly('Forcetheorem_Loop_DMI',(/'Q-vec'/),(/this%q_done/))
  END FUNCTION dmi_next_job

  SUBROUTINE dmi_postprocess(this)
    USE m_xmlOutput
    IMPLICIT NONE
    CLASS(t_forcetheo_dmi),INTENT(INOUT):: this

    !Locals
    INTEGER:: n,q
    CHARACTER(LEN=12):: attributes(4)
    IF (this%q_done==0) RETURN
    IF (this%l_io) THEN
       !Now output the results
       CALL closeXMLElement('Forcetheorem_Loop_DMI')
       CALL openXMLElementPoly('Forcetheorem_DMI',(/'qPoints','Angles '/),(/SIZE(this%evsum,2),SIZE(this%evsum,1)/))
       DO q=1,SIZE(this%evsum,2)
          WRITE(attributes(1),'(i5)') q
          WRITE(attributes(2),'(f12.7)') this%evsum(0,q) 
          CALL writeXMLElementForm('Entry',(/'q     ','ev-sum'/),attributes(1:2),&
               RESHAPE((/1,6,5,12/),(/2,2/)))
          DO n=1,SIZE(this%evsum,1)-1
             WRITE(attributes(2),'(f12.7)') this%theta(n)
             WRITE(attributes(3),'(f12.7)') this%phi(n)
             WRITE(attributes(4),'(f12.7)') this%evsum(n,q)     
             CALL writeXMLElementForm('Entry',(/'q     ','theta ','phi   ','ev-sum'/),attributes,&
                  RESHAPE((/1,5,3,6,5,12,12,12/),(/4,2/)))
          END DO
       ENDDO
       CALL closeXMLElement('Forcetheorem_DMI')
    ENDIF

    CALL judft_end("Forcetheorem DMI")
  END SUBROUTINE dmi_postprocess

  SUBROUTINE dmi_dist(this,mpi)
    USE m_types_mpi
    IMPLICIT NONE
    CLASS(t_forcetheo_dmi),INTENT(INOUT):: this
    TYPE(t_mpi),INTENT(in):: mpi
    
    INTEGER:: i,q,ierr
#ifdef CPP_MPI    
    INCLUDE 'mpif.h'
    IF (mpi%irank==0) i=SIZE(this%theta)
    call MPI_BCAST(i,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    IF (mpi%irank==0) q=SIZE(this%qvec,2)
    CALL MPI_BCAST(q,1,MPI_INTEGER,0,mpi%mpi_comm,ierr)
    IF (mpi%irank.NE.0) ALLOCATE(this%qvec(3,q),this%phi(i),this%theta(i),this%evsum(0:i,q));this%evsum=0.0
    CALL MPI_BCAST(this%phi,i,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(this%theta,i,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(this%qvec,3*q,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
#endif    
  END SUBROUTINE dmi_dist

  FUNCTION dmi_eval(this,eig_id,DIMENSION,atoms,kpts,sym,&
       cell,noco, input,mpi, oneD,enpara,v,results)RESULT(skip)
     USE m_types
     USE m_ssomat
    IMPLICIT NONE
    CLASS(t_forcetheo_dmi),INTENT(INOUT):: this
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
  
    this%evsum(0,this%q_done)=results%seigv
    CALL ssomat(this%evsum(1:,this%q_done),this%theta,this%phi,eig_id,DIMENSION,atoms,kpts,sym,&
       cell,noco, input,mpi, oneD,enpara,v,results) 
    skip=.TRUE.
  END FUNCTION  dmi_eval

  
END MODULE m_types_dmi
