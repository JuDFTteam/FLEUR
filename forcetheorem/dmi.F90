!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_dmi
  USE m_types
  USE m_types_forcetheo
  USE m_judft
  IMPLICIT NONE
  PRIVATE
  TYPE,EXTENDS(t_forcetheo) :: t_forcetheo_dmi
     INTEGER :: q_done
     REAL,ALLOCATABLE:: qvec(:,:)
     REAL,ALLOCATABLE:: theta(:)
     REAL,ALLOCATABLE:: phi(:)
     REAL,ALLOCATABLE:: evsum(:,:)
     REAL,ALLOCATABLE:: h_so(:,:,:)
   CONTAINS
     PROCEDURE :: start   =>dmi_start
     PROCEDURE :: next_job=>dmi_next_job
     PROCEDURE :: eval    =>dmi_eval
     PROCEDURE :: postprocess => dmi_postprocess
     PROCEDURE :: init   => dmi_init !not overloaded
     PROCEDURE :: dist   => dmi_dist !not overloaded
  END TYPE t_forcetheo_dmi
  PUBLIC t_forcetheo_dmi
CONTAINS


  SUBROUTINE dmi_init(this,q,theta,phi,ntype)
    USE m_calculator
    USE m_constants
    IMPLICIT NONE
    CLASS(t_forcetheo_dmi),INTENT(INOUT):: this
    REAL,INTENT(in)                     :: q(:,:)
    REAL,INTENT(IN)                     :: theta(:),phi(:)
    INTEGER,INTENT(IN)                  :: ntype

    this%theta=theta
    this%phi=phi

    IF (SIZE(this%phi).NE.SIZE(this%theta)) CALL &
         judft_error("Lists for theta/phi must have the same length in DMI force theorem calculations")

    ! use same definition of rotation angles as in noco-routines
    this%theta=-this%theta
    this%phi=this%phi+pi_const


    ALLOCATE(this%qvec(3,SIZE(q,2)))
    this%qvec=q

    ALLOCATE(this%evsum(0:SIZE(this%phi),SIZE(q,2)))
    ALLOCATE(this%h_so(0:ntype,0:SIZE(this%phi),SIZE(q,2)))

    this%evsum=0;this%h_so=0.0
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

  LOGICAL FUNCTION dmi_next_job(this,lastiter,atoms,noco,nococonv)
    USE m_types_setup
    USE m_xmlOutput
    USE m_constants
    USE m_types_nococonv
    IMPLICIT NONE
    CLASS(t_forcetheo_dmi),INTENT(INOUT):: this
    LOGICAL,INTENT(IN)                  :: lastiter
    TYPE(t_atoms),INTENT(IN)            :: atoms
    !Stuff that might be modified...
    TYPE(t_noco),INTENT(IN) :: noco
    TYPE(t_nococonv),INTENT(INOUT) :: nococonv
    INTEGER                 :: itype
    CHARACTER(LEN=12):: attributes(2)
    IF (.NOT.lastiter) THEN
       dmi_next_job=this%t_forcetheo%next_job(lastiter,atoms,noco,nococonv)
       RETURN
    ENDIF
    !OK, now we start the DMI-loop
    this%q_done=this%q_done+1
    dmi_next_job=(this%q_done<=SIZE(this%qvec,2)) !still q-vectors to do
    IF (.NOT.dmi_next_job) RETURN

    !Now modify the noco-file
    nococonv%qss=this%qvec(:,this%q_done)
    if (.not.noco%l_spav) call judft_warn("l_spav=T should be set in DMI force theorem mode")
    !Modify the alpha-angles
    DO iType = 1,atoms%ntype
       nococonv%alph(iType) = noco%alph_inp(iType) + tpi_const*dot_PRODUCT(nococonv%qss,atoms%taual(:,SUM(atoms%neq(:itype-1))+1))
    END DO
    IF (.NOT.this%l_io) RETURN

    IF (this%q_done.NE.1) CALL closeXMLElement('Forcetheorem_Loop')
    WRITE(attributes(1),'(a)') 'DMI'
    WRITE(attributes(2),'(i5)') this%q_done
    CALL openXMLElementPoly('Forcetheorem_Loop',(/'calculationType','No             '/),attributes)
  END FUNCTION dmi_next_job

  SUBROUTINE dmi_postprocess(this)
    USE m_xmlOutput
    IMPLICIT NONE
    CLASS(t_forcetheo_dmi),INTENT(INOUT):: this

    !Locals
    INTEGER:: n,q,i
    CHARACTER(LEN=12):: attributes(5)
    CHARACTER(LEN=16) :: atom_name
    IF (this%q_done==0) RETURN
    IF (this%l_io) THEN
       !Now output the results
       CALL closeXMLElement('Forcetheorem_Loop')
       attributes = ''
       WRITE(attributes(1),'(i5)') SIZE(this%evsum,2)
       WRITE(attributes(2),'(i5)') SIZE(this%evsum,1)
       WRITE(attributes(3),'(a)') 'Htr'
       CALL openXMLElement('Forcetheorem_DMI',(/'qPoints','Angles ','units  '/),attributes(:3))
       DO q=1,SIZE(this%evsum,2)
          WRITE(attributes(2),'(i5)') q
          WRITE(attributes(3),'(f12.7)') this%evsum(0,q)
          CALL writeXMLElementForm('Entry',(/'q     ','ev-sum'/),attributes(2:3),&
                                   RESHAPE((/1,6,5,12/),(/2,2/)))
          DO n=1,SIZE(this%evsum,1)-1
             WRITE(attributes(3),'(f12.7)') this%theta(n)
             WRITE(attributes(4),'(f12.7)') this%phi(n)
             WRITE(attributes(5),'(f12.7)') this%evsum(n,q)
             CALL writeXMLElementForm('Entry',(/'q     ','theta ','phi   ','ev-sum'/),attributes(2:),RESHAPE((/1,5,3,6,5,12,12,12/),(/4,2/)))
             write(attributes(5),'(f12.7)') this%h_so(0,n,q)
             CALL writeXMLElementForm('allAtoms',(/'q     ','theta ','phi   ','H_so  '/),attributes(2:),RESHAPE((/1,5,3,6,5,12,12,12/),(/4,2/)))
             DO i=1,size(this%h_so,1)-1
               write(attributes(5),'(f12.7)') this%h_so(i,n,q)
               write(attributes(1),'(i0)') i
               CALL writeXMLElementForm('singleAtom',(/'atomType','q       ','theta   ','phi     ','H_so    '/),attributes,RESHAPE((/8,1,5,3,6,5,5,12,12,12/),(/5,2/)))
             ENDDO
          END DO
       ENDDO
       CALL closeXMLElement('Forcetheorem_DMI')
    ENDIF

    CALL judft_end("Forcetheorem DMI")
  END SUBROUTINE dmi_postprocess

  SUBROUTINE dmi_dist(this,fmpi)
#ifdef CPP_MPI
    USE mpi
#endif    
    USE m_types_mpi
    IMPLICIT NONE
    CLASS(t_forcetheo_dmi),INTENT(INOUT):: this
    TYPE(t_mpi),INTENT(in):: fmpi

    INTEGER:: i,q,ierr,n
#ifdef CPP_MPI
    IF (fmpi%irank==0) i=SIZE(this%theta)
    call MPI_BCAST(i,1,MPI_INTEGER,0,fmpi%mpi_comm,ierr)
    IF (fmpi%irank==0) q=SIZE(this%qvec,2)
    CALL MPI_BCAST(q,1,MPI_INTEGER,0,fmpi%mpi_comm,ierr)
    IF (fmpi%irank.NE.0) ALLOCATE(this%qvec(3,q),this%phi(i),this%theta(i),this%evsum(0:i,q));this%evsum=0.0
    if (fmpi%irank==0) n=size(this%h_so,1)-1
    call MPI_BCAST(n,1,MPI_INTEGER,0,fmpi%mpi_comm,ierr)
    IF (fmpi%irank.NE.0) ALLOCATE(this%h_so(0:n,0:i,q));this%h_so=0.0
    CALL MPI_BCAST(this%phi,i,MPI_DOUBLE_PRECISION,0,fmpi%mpi_comm,ierr)
    CALL MPI_BCAST(this%theta,i,MPI_DOUBLE_PRECISION,0,fmpi%mpi_comm,ierr)
    CALL MPI_BCAST(this%qvec,3*q,MPI_DOUBLE_PRECISION,0,fmpi%mpi_comm,ierr)
#endif
  END SUBROUTINE dmi_dist

  FUNCTION dmi_eval(this,eig_id,atoms,kpts,sym,&
       cell,noco,nococonv, input,fmpi, oneD,enpara,v,results)RESULT(skip)
     USE m_types
     USE m_ssomat
    IMPLICIT NONE
    CLASS(t_forcetheo_dmi),INTENT(INOUT):: this
    LOGICAL :: skip
    !Stuff that might be used...
    TYPE(t_mpi),INTENT(IN)         :: fmpi

    TYPE(t_oneD),INTENT(IN)        :: oneD
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_noco),INTENT(IN)        :: noco
    TYPE(t_nococonv),INTENT(IN)    :: nococonv
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
    CALL ssomat(this%evsum(1:,this%q_done),this%h_so(:,:,this%q_done),this%theta,this%phi,eig_id,atoms,kpts,sym,&
       cell,noco,nococonv, input,fmpi, oneD,enpara,v,results)
    skip=.TRUE.
  END FUNCTION  dmi_eval


END MODULE m_types_dmi
