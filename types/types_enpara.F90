!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_enpara
  USE m_judft
  use m_types_enparaxml
  IMPLICIT NONE
  PRIVATE
  TYPE,extends(t_enparaxml):: t_enpara
     REAL, ALLOCATABLE CPP_MANAGED   :: el0(:,:,:)
     REAL                 :: evac(2,2)
     REAL, ALLOCATABLE    :: ello0(:,:,:)
     REAL, ALLOCATABLE    :: el1(:,:,:)
     REAL                 :: evac1(2,2)
     REAL, ALLOCATABLE    :: ello1(:,:,:)
     REAL, ALLOCATABLE    :: enmix(:)
     INTEGER, ALLOCATABLE :: skiplo(:,:)
     LOGICAL, ALLOCATABLE :: lchange(:,:,:)
     LOGICAL, ALLOCATABLE :: lchg_v(:,:)
     LOGICAL, ALLOCATABLE :: llochg(:,:,:)
     REAL                 :: epara_min
     LOGICAL              :: ready ! are the enpara's ok for calculation?
     LOGICAL              :: floating !floating energy parameters are relative to potential
   CONTAINS
     PROCEDURE :: init_enpara
     PROCEDURE :: update
     PROCEDURE :: read
     PROCEDURE :: write
     PROCEDURE :: mix
     PROCEDURE :: calcOutParams
  END TYPE t_enpara



  PUBLIC:: t_enpara

CONTAINS
  SUBROUTINE init_enpara(this,atoms,jspins,film,enparaXML)
    USE m_types_setup
    USE m_constants
    CLASS(t_enpara),INTENT(inout):: this
    TYPE(t_atoms),INTENT(IN)     :: atoms
    INTEGER,INTENT(IN)           :: jspins
    LOGICAL,INTENT(IN)           :: film
    TYPE(t_enparaxml),OPTIONAL   :: enparaxml
 
    LOGICAL :: l_enpara
    
    ALLOCATE(this%el0(0:atoms%lmaxd,atoms%ntype,jspins),this%el1(0:atoms%lmaxd,atoms%ntype,jspins))
    ALLOCATE(this%ello0(atoms%nlod,atoms%ntype,jspins),this%ello1(atoms%nlod,atoms%ntype,jspins))
    this%el0=-1E99
    this%ello0=-1E99
    this%evac0=-1E99


    ALLOCATE(this%llochg(atoms%nlod,atoms%ntype,jspins))
    ALLOCATE(this%lchg_v(2,jspins))
    ALLOCATE(this%skiplo(atoms%ntype,jspins))
    ALLOCATE(this%lchange(0:atoms%lmaxd,atoms%ntype,jspins))
    this%llochg=.FALSE.;this%lchg_v=.FALSE.;this%lchange=.FALSE.
    this%skiplo=0
    ALLOCATE(this%enmix(jspins))
    this%enmix=0.0

    this%ready=.FALSE.
    this%floating=.FALSE.

    ALLOCATE(this%qn_el(0:3,atoms%ntype,jspins))
    ALLOCATE(this%qn_ello(atoms%nlod,atoms%ntype,jspins))

    this%evac0=eVac0Default_const

    inquire(file="enpara",exist=l_enpara)
    if (l_enpara) then
       call this%read(atoms,jspins,film,.true.)
    else
       IF (PRESENT(enparaxml)) THEN
          this%qn_el=enparaxml%qn_el
          this%qn_ello=enparaxml%qn_ello
          this%evac0=enparaxml%evac0
       END IF
    endif
    
    
  END SUBROUTINE init_enpara

  !> This subroutine adjusts the energy parameters to the potential. In particular, it
  !! calculated them in case of qn_el>-1,qn_ello>-1
  !! Before this was done in lodpot.F
  SUBROUTINE update(enpara,mpi_comm,atoms,vacuum,input,v)
    USE m_types_atoms
    USE m_types_vacuum
    USE m_types_input
    USE m_xmlOutput
    USE m_types_potden
    USE m_find_enpara
    CLASS(t_enpara),INTENT(inout):: enpara
    INTEGER,INTENT(IN)          :: mpi_comm
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_vacuum),INTENT(IN)   :: vacuum
    TYPE(t_potden),INTENT(IN)   :: v


    LOGICAL ::  l_enpara
    LOGICAL ::  l_done(0:atoms%lmaxd,atoms%ntype,input%jspins)
    LOGICAL ::  lo_done(atoms%nlod,atoms%ntype,input%jspins)
    REAL    ::  vbar,vz0,rj
    INTEGER ::  n,jsp,l,ilo,j,ivac,irank
    CHARACTER(LEN=20)    :: attributes(5)
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER :: ierr
    CALL MPI_COMM_RANK(mpi_comm,irank,ierr)
#else
    irank=0
#endif
    IF (irank  == 0) CALL openXMLElement('energyParameters',(/'units'/),(/'Htr'/))

    l_done = .FALSE.;lo_done=.FALSE.
    DO jsp = 1,input%jspins
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP SHARED(atoms,enpara,jsp,l_done,irank,v,lo_done) &
       !$OMP PRIVATE(n,l,ilo)
       !! First calculate energy paramter from quantum numbers if these are given...
       !! l_done stores the index of those energy parameter updated
       DO n = 1, atoms%ntype
          DO l = 0,3
             IF( enpara%qn_el(l,n,jsp).ne.0)THEN 
                l_done(l,n,jsp) = .TRUE.
                enpara%el0(l,n,jsp)=find_enpara(.FALSE.,l,n,jsp,enpara%qn_el(l,n,jsp),atoms,irank,v%mt(:,0,n,jsp))
                IF( l .EQ. 3 ) THEN
                   enpara%el0(4:,n,jsp) = enpara%el0(3,n,jsp)
                   l_done(4:,n,jsp) = .TRUE.
                END IF
             ELSE 
                l_done(l,n,jsp) = .FALSE.
             END IF
          ENDDO ! l
          ! Now for the lo's
          DO ilo = 1, atoms%nlo(n)
             l = atoms%llo(ilo,n)
             IF( enpara%qn_ello(ilo,n,jsp).NE.0) THEN
                lo_done(ilo,n,jsp) = .TRUE.
                enpara%ello0(ilo,n,jsp)=find_enpara(.TRUE.,l,n,jsp,enpara%qn_ello(ilo,n,jsp),atoms,irank,v%mt(:,0,n,jsp))
             ELSE
                lo_done(ilo,n,jsp) = .FALSE.
             ENDIF
          ENDDO
       ENDDO ! n
       !$OMP END PARALLEL DO

       !!   Now check for floating energy parameters (not for those with l_done=T)
       IF (enpara%floating) THEN
          types_loop: DO n = 1,atoms%ntype 
             !
             !--->    determine energy parameters if lepr=1. the reference energy
             !--->    is the value of the l=0 potential at approximately rmt/4.
             !
             j = atoms%jri(n) - (LOG(4.0)/atoms%dx(n)+1.51)
             rj = atoms%rmt(n)*EXP(atoms%dx(n)* (j-atoms%jri(n)))
             vbar = v%mt(j,0,n,jsp)/rj
             IF (irank.EQ.0) THEN
                WRITE ( 6,'('' spin'',i2,'', atom type'',i3,'' ='',f12.6,''   r='',f8.5)') jsp,n,vbar,rj
             ENDIF
             DO l = 0,atoms%lmax(n)
                IF ( .NOT.l_done(l,n,jsp) ) THEN
                   enpara%el0(l,n,jsp) = vbar + enpara%el0(l,n,jsp)
                END IF
             ENDDO
             IF (atoms%nlo(n).GE.1) THEN
                DO ilo = 1,atoms%nlo(n)
                   IF ( .NOT. lo_done(ilo,n,jsp) ) THEN
                      enpara%ello0(ilo,n,jsp) = vbar + enpara%ello0(ilo,n,jsp)
                      !+apw+lo
                      IF (atoms%l_dulo(ilo,n)) THEN
                         enpara%ello0(ilo,n,jsp) = enpara%el0(atoms%llo(ilo,n),n,jsp)
                      ENDIF
                      !-apw+lo
                   END IF
                END DO
             ENDIF
          END DO types_loop
       ENDIF
       IF (input%film) THEN

          INQUIRE (file ='enpara',exist= l_enpara)
          IF(l_enpara) enpara%evac0(:,jsp) = enpara%evac(:,jsp)

          !--->    vacuum energy parameters: for floating: relative to potential
          !--->    at vacuum-interstitial interface (better for electric field)

          DO ivac = 1,vacuum%nvac
             vz0 = 0.0
             IF (enpara%floating) THEN
                vz0 = v%vacz(1,ivac,jsp)
                IF (irank.EQ.0) THEN
                   WRITE ( 6,'('' spin'',i2,'', vacuum   '',i3,'' ='',f12.6)') jsp,ivac,vz0 
                ENDIF
             ENDIF
             enpara%evac(ivac,jsp) = enpara%evac0(ivac,jsp) + vz0
             IF (.NOT.l_enpara) THEN
                enpara%evac(ivac,jsp) = v%vacz(vacuum%nmz,ivac,jsp) + enpara%evac0(ivac,jsp)
             END IF
             IF (irank.EQ.0) THEN
                attributes = ''
                WRITE(attributes(1),'(i0)') ivac
                WRITE(attributes(2),'(i0)') jsp
                WRITE(attributes(3),'(f16.10)') v%vacz(1,ivac,jsp)
                WRITE(attributes(4),'(f16.10)') v%vacz(vacuum%nmz,ivac,jsp)
                WRITE(attributes(5),'(f16.10)') enpara%evac(ivac,jsp)
                CALL writeXMLElementForm('vacuumEP',(/'vacuum','spin  ','vzIR  ','vzInf ','value '/),&
                     attributes(1:5),RESHAPE((/6+4,4,4,5,5+13,8,1,16,16,16/),(/5,2/)))
             END IF
          ENDDO
          IF (vacuum%nvac.EQ.1) THEN
             enpara%evac(2,jsp) = enpara%evac(1,jsp)
          END IF
       END IF
    END DO

!    enpara%ready=(ALL(enpara%el0>-1E99).AND.ALL(enpara%ello0>-1E99))
    enpara%epara_min=MIN(MINVAL(enpara%el0),MINVAL(enpara%ello0))
    
    IF (irank  == 0) CALL closeXMLElement('energyParameters')
  END SUBROUTINE update

  SUBROUTINE READ(enpara,atoms,jspins,film,l_required)
    USE m_types_setup
    IMPLICIT NONE
    CLASS(t_enpara),INTENT(INOUT):: enpara
    INTEGER, INTENT (IN)        :: jspins
    TYPE(t_atoms),INTENT(IN)    :: atoms
    LOGICAL,INTENT(IN)          :: film,l_required

    INTEGER :: n,l,lo,skip_t,io_err,jsp
    logical :: l_exist

    INQUIRE(file="enpara",exist=l_exist)
    IF (.NOT.l_exist.AND.l_required) CALL juDFT_error("No enpara file found")
    IF (.NOT.l_exist) RETURN

    OPEN(40,file="enpara",form="formatted",status="old")

    DO jsp=1,jspins

       !-->  first line contains mixing parameter!

       READ (40,FMT ='(48x,f10.6)',END=200) enpara%enmix(jsp)
       READ (40,*)                       ! skip next line
       IF (enpara%enmix(jsp).EQ.0.0) enpara%enmix(jsp) = 1.0
       WRITE (6,FMT=8001) jsp
       WRITE (6,FMT=8000)
       skip_t = 0
       DO n = 1,atoms%ntype
          READ (40,FMT=8040,END=200) (enpara%el0(l,n,jsp),l=0,3),&
               (enpara%lchange(l,n,jsp),l=0,3),enpara%skiplo(n,jsp)    
          WRITE (6,FMT=8140) n,(enpara%el0(l,n,jsp),l=0,3),&
               (enpara%lchange(l,n,jsp),l=0,3),enpara%skiplo(n,jsp)    
          !
          !--->    energy parameters for the local orbitals
          !
          IF (atoms%nlo(n).GE.1) THEN
             skip_t = skip_t + enpara%skiplo(n,jsp) * atoms%neq(n)
             READ (40,FMT=8039,END=200)  (enpara%ello0(lo,n,jsp),lo=1,atoms%nlo(n))
             READ (40,FMT=8038,END=200) (enpara%llochg(lo,n,jsp),lo=1,atoms%nlo(n))
             WRITE (6,FMT=8139)          (enpara%ello0(lo,n,jsp),lo=1,atoms%nlo(n))
             WRITE (6,FMT=8138)         (enpara%llochg(lo,n,jsp),lo=1,atoms%nlo(n))
          ELSEIF (enpara%skiplo(n,jsp).GT.0) THEN
             WRITE (6,*) "for atom",n," no LO's were specified"
             WRITE (6,*) 'but skiplo was set to',enpara%skiplo 
             CALL juDFT_error("No LO's but skiplo",calledby ="enpara",&
                  hint="If no LO's are set skiplo must be 0 in enpara")
          END IF
          !Integer values mean we have to set the qn-arrays
          enpara%qn_el(:,n,jsp)=0
          DO l=0,3
             IF (enpara%el0(l,n,jsp)==NINT(enpara%el0(l,n,jsp))) enpara%qn_el(l,n,jsp)=NINT(enpara%el0(l,n,jsp))
          ENDDO
          enpara%qn_ello(:,n,jsp)=0
          DO l=1,atoms%nlo(n)
             IF (enpara%ello0(l,n,jsp)==NINT(enpara%ello0(l,n,jsp))) enpara%qn_ello(l,n,jsp)=NINT(enpara%ello0(l,n,jsp))
          ENDDO
          !
          !--->    set the energy parameters with l>3 to the value of l=3
          !
          enpara%el0(4:,n,jsp) = enpara%el0(3,n,jsp)
       ENDDO   ! atoms%ntype

       IF (film) THEN
          enpara%lchg_v = .TRUE.
          READ (40,FMT=8050,END=200) enpara%evac0(1,jsp),enpara%lchg_v(1,jsp),enpara%evac0(2,jsp)
          WRITE (6,FMT=8150)         enpara%evac0(1,jsp),enpara%lchg_v(1,jsp),enpara%evac0(2,jsp)
       ENDIF
      ! IF (atoms%nlod.GE.1) THEN               
      !    WRITE (6,FMT=8090) jsp,skip_t
      !    WRITE (6,FMT=8091) 
      ! END IF
    END DO

    enpara%evac(:,:) = enpara%evac0(:,:)

    CLOSE(40)
    ! input formats

8038 FORMAT (14x,60(l1,8x))
8039 FORMAT (8x,60f9.5)
8040 FORMAT (8x,4f9.5,9x,4l1,9x,i3)
8050 FORMAT (19x,f9.5,9x,l1,15x,f9.5)

    ! output formats

8138 FORMAT (' --> change   ',60(l1,8x))
8139 FORMAT (' --> lo ',60f9.5)
8140 FORMAT (' -->',i3,1x,4f9.5,' change: ',4l1,' skiplo: ',i3)
8150 FORMAT ('  vacuum parameter=',f9.5,' change: ',l1, ' second vacuum=',f9.5)
8001 FORMAT ('READING enpara for spin: ',i1)
8000 FORMAT (/,' energy parameters:',/,t10,'s',t20, 'p',t30,'d',t37,'higher l - - -')
8090 FORMAT ('Spin: ',i1,' -- ',i3,'eigenvalues')
8091 FORMAT ('will be skipped for energyparameter computation')

    RETURN

200 WRITE (6,*) 'the end of the file enpara has been reached while'
    WRITE (6,*) 'reading the energy-parameters.'
    WRITE (6,*) 'possible reason: energy parameters have not been'
    WRITE (6,*) 'specified for all atom types.'
    WRITE (6,FMT='(a,i4)') 'the actual number of atom-types is: ntype=',atoms%ntype
    CALL juDFT_error ("unexpected end of file enpara reached while reading")
  END SUBROUTINE read


  SUBROUTINE WRITE(enpara, atoms,jspins,film)

    ! write enpara-file
    !
    USE m_types_setup
    IMPLICIT NONE
    CLASS(t_enpara),INTENT(IN) :: enpara
    INTEGER, INTENT (IN) :: jspins
    LOGICAL,INTENT(IN)   :: film
    TYPE(t_atoms),INTENT(IN) :: atoms

    INTEGER n,l,lo,jspin
    REAL el0Temp(0:3), ello0Temp(1:atoms%nlod)

    OPEN(unit=40,file="enpara",form="formatted",status="replace")

    DO jspin=1,jspins
       WRITE (40,FMT=8035) jspin,enpara%enmix(jspin)
       WRITE (40,FMT=8036)
8035   FORMAT (5x,'energy parameters          for spin ',i1,' mix=',f10.6)
8036   FORMAT (t6,'atom',t15,'s',t24,'p',t33,'d',t42,'f')
       DO n = 1,atoms%ntype
          el0Temp(0:3) = enpara%el0(0:3,n,jspin)
          DO l = 0, 3
             IF (enpara%qn_el(l,n,jspin).NE.0) el0Temp(l) =  REAL(enpara%qn_el(l,n,jspin))
          END DO
          WRITE (6,FMT=8040)  n, (el0Temp(l),l=0,3),&
               &                          (enpara%lchange(l,n,jspin),l=0,3),enpara%skiplo(n,jspin)
          WRITE (40,FMT=8040) n, (el0Temp(l),l=0,3),&
               &                          (enpara%lchange(l,n,jspin),l=0,3),enpara%skiplo(n,jspin)
          !--->    energy parameters for the local orbitals
          IF (atoms%nlo(n).GE.1) THEN
             ello0Temp(1:atoms%nlo(n)) = enpara%ello0(1:atoms%nlo(n),n,jspin)
             DO lo = 1, atoms%nlo(n)
                IF (enpara%qn_ello(lo,n,jspin).NE.0) ello0Temp(lo) = enpara%qn_ello(lo,n,jspin)
             END DO
             WRITE (6,FMT=8039) (ello0Temp(lo),lo=1,atoms%nlo(n))
             WRITE (6,FMT=8038) (enpara%llochg(lo,n,jspin),lo=1,atoms%nlo(n))
             WRITE (40,FMT=8039) (ello0Temp(lo),lo=1,atoms%nlo(n))
             WRITE (40,FMT=8038) (enpara%llochg(lo,n,jspin),lo=1,atoms%nlo(n))
          END IF

       ENDDO
8038   FORMAT (' --> change   ',60(l1,8x))
8039   FORMAT (' --> lo ',60f9.5)
8040   FORMAT (' -->',i3,1x,4f9.5,' change: ',4l1,' skiplo: ',i3)

       IF (film) THEN
          WRITE (40,FMT=8050) enpara%evac(1,jspin),enpara%lchg_v(1,jspin),enpara%evac(2,jspin)
          WRITE (6,FMT=8050)  enpara%evac(1,jspin),enpara%lchg_v(1,jspin),enpara%evac(2,jspin)
8050      FORMAT ('  vacuum parameter=',f9.5,' change: ',l1,&
               &           ' second vacuum=',f9.5)
       ENDIF
    ENDDO
    CLOSE(40)
    RETURN
  END SUBROUTINE WRITE





  SUBROUTINE mix(enpara,mpi_comm,atoms,vacuum,input,vr,vz)
    !------------------------------------------------------------------
    USE m_types_atoms
    USE m_types_input
    USE m_types_vacuum
    USE m_types_mpi
    IMPLICIT NONE
    CLASS(t_enpara),INTENT(INOUT)  :: enpara
    INTEGER,INTENT(IN)             :: mpi_comm
    TYPE(t_atoms),INTENT(IN)       :: atoms
    TYPE(t_vacuum),INTENT(IN)      :: vacuum
    TYPE(t_input),INTENT(IN)       :: input

    REAL,    INTENT(IN) :: vr(:,:,:)
    REAL,    INTENT(IN) :: vz(vacuum%nmzd,2)

    INTEGER ityp,j,l,lo,jsp,n
    REAL    vbar,maxdist,maxdist2
    INTEGER same(atoms%nlod),irank
    LOGICAL l_enpara
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER :: ierr
    CALL MPI_COMM_RANK(mpi_comm,irank,ierr)
#else
    irank=0
#endif

    IF (irank==0) THEN
       maxdist2=0.0
       DO jsp=1,SIZE(enpara%el0,3)
          maxdist=0.0
          DO ityp = 1,atoms%ntype
             !        look for LO's energy parameters equal to the LAPW (and previous LO) ones
             same = 0
             DO lo = 1,atoms%nlo(ityp)
                IF(enpara%el0(atoms%llo(lo,ityp),ityp,jsp).EQ.enpara%ello0(lo,ityp,jsp)) same(lo)=-1
                DO l = 1,lo-1
                   IF(atoms%llo(l,ityp).NE.atoms%llo(lo,ityp)) CYCLE
                   IF(enpara%ello0(l,ityp,jsp).EQ.enpara%ello0(lo,ityp,jsp).AND.same(lo).EQ.0) same(lo)=l
                ENDDO
             ENDDO
             !
             !--->   change energy parameters
             !
             DO l = 0,3
                WRITE(6,*) 'Type:',ityp,' l:',l
                WRITE(6,FMT=777) enpara%el0(l,ityp,jsp),enpara%el1(l,ityp,jsp),&
                     ABS(enpara%el0(l,ityp,jsp)-enpara%el1(l,ityp,jsp))
                maxdist=MAX(maxdist,ABS(enpara%el0(l,ityp,jsp)-enpara%el1(l,ityp,jsp)))
                IF ( enpara%lchange(l,ityp,jsp) ) THEN
                   maxdist2=MAX(maxdist2,ABS(enpara%el0(l,ityp,jsp)-enpara%el1(l,ityp,jsp)))
                   enpara%el0(l,ityp,jsp) =(1.0-enpara%enmix(jsp))*enpara%el0(l,ityp,jsp) + &
                        enpara%enmix(jsp)*enpara%el1(l,ityp,jsp)
                ENDIF
             ENDDO
             IF ( enpara%lchange(3,ityp,jsp) ) enpara%el0(4:,ityp,jsp) = enpara%el0(3,ityp,jsp)
             !
             !--->    determine and change local orbital energy parameters
             !
             DO lo = 1,atoms%nlo(ityp)
                IF (atoms%l_dulo(lo,ityp)) THEN
                   enpara%ello0(lo,ityp,jsp) =enpara%el0(atoms%llo(lo,ityp),ityp,jsp)
                ELSE
                   IF (enpara%llochg(lo,ityp,jsp) ) THEN
                      IF(same(lo).EQ.-1) THEN
                         enpara%ello0(lo,ityp,jsp) = enpara%el0(atoms%llo(lo,ityp),ityp,jsp)
                         CYCLE
                      ELSE IF(same(lo).GT.0) THEN
                         enpara%ello0(lo,ityp,jsp) = enpara%ello0(same(lo),ityp,jsp)
                         CYCLE
                      ENDIF
                   ENDIF
                   WRITE(6,*) 'Type:',ityp,' lo:',lo
                   WRITE(6,FMT=777) enpara%ello0(lo,ityp,jsp),enpara%ello1(lo,ityp,jsp),&
                        ABS(enpara%ello0(lo,ityp,jsp)-enpara%ello1(lo,ityp,jsp))
                   maxdist=MAX(maxdist,ABS(enpara%ello0(lo,ityp,jsp)-enpara%ello1(lo,ityp,jsp)))
                   IF (enpara%llochg(lo,ityp,jsp) ) THEN
                      maxdist2=MAX(maxdist2,ABS(enpara%ello0(lo,ityp,jsp)-enpara%ello1(lo,ityp,jsp)))
                      enpara%ello0(lo,ityp,jsp) =(1.0-enpara%enmix(jsp))*enpara%ello0(lo,ityp,jsp)+&
                           enpara%enmix(jsp)*enpara%ello1(lo,ityp,jsp)
                   ENDIF
                END IF
             END DO
             !Shift if floating energy parameters are used
             IF (enpara%floating) THEN
                j = atoms%jri(ityp) - (LOG(4.0)/atoms%dx(ityp)+1.51)
                vbar = vr(j,ityp,jsp)/( atoms%rmt(ityp)*EXP(atoms%dx(ityp)*(j-atoms%jri(ityp))) )
                enpara%el0(:,n,jsp)=enpara%el0(:,n,jsp)-vbar
             ENDIF
          END DO
          
          
          IF (input%film) THEN
             WRITE(6,*) 'Vacuum:'
             DO n=1,vacuum%nvac
                WRITE(6,FMT=777) enpara%evac(n,jsp),enpara%evac1(n,jsp),ABS(enpara%evac(n,jsp)-enpara%evac1(n,jsp))
                maxdist=MAX(maxdist,ABS(enpara%evac(n,jsp)-enpara%evac1(n,jsp)))
                IF (enpara%lchg_v(n,jsp) ) THEN
                   maxdist2=MAX(maxdist2,ABS(enpara%evac(n,jsp)-enpara%evac1(n,jsp)))
                   enpara%evac(n,jsp) =(1.0-enpara%enmix(jsp))*enpara%evac(n,jsp)+&
                        enpara%enmix(jsp)*enpara%evac1(n,jsp)
                END IF
             END DO
             IF (vacuum%nvac==1) enpara%evac(2,jsp) = enpara%evac(1,jsp)
             IF (enpara%floating) enpara%evac(:,jsp)=enpara%evac(:,jsp)-vz(:,jsp)
          ENDIF
          WRITE(6,'(a36,f12.6)') 'Max. mismatch of energy parameters:', maxdist
       END DO
       IF (maxdist2>1.0) CALL juDFT_warn&
            ("Energy parameter mismatch too large",hint&
            ="If any energy parameters calculated from the output "//&
            "differ from the input by more than 1Htr, chances are "//&
            "high that your initial setup was broken.")
    ENDIF
    INQUIRE(file='enpara',exist=l_enpara)
    IF (irank==0.AND.l_enpara) CALL enpara%WRITE(atoms,input%jspins,input%film)
#ifdef CPP_MPI
    CALL MPI_BCAST(enpara%el0,SIZE(enpara%el0),MPI_DOUBLE_PRECISION,0,mpi_comm,ierr)
    CALL MPI_BCAST(enpara%ello0,SIZE(enpara%ello0),MPI_DOUBLE_PRECISION,0,mpi_comm,ierr)
    CALL MPI_BCAST(enpara%evac,SIZE(enpara%evac),MPI_DOUBLE_PRECISION,0,mpi_comm,ierr)
#endif    
    RETURN
777 FORMAT('Old:',f8.5,' new:',f8.5,' diff:',f8.5)

  END SUBROUTINE mix

  SUBROUTINE calcOutParams(enpara,input,atoms,vacuum,regCharges)
    USE m_types_setup
    USE m_types_regionCharges
    IMPLICIT NONE
    CLASS(t_enpara),INTENT(INOUT)    :: enpara
    TYPE(t_input),INTENT(IN)         :: input
    TYPE(t_atoms),INTENT(IN)         :: atoms
    TYPE(t_vacuum),INTENT(IN)        :: vacuum
    TYPE(t_regionCharges),INTENT(IN) :: regCharges

    INTEGER :: ispin, n

    DO ispin = 1,input%jspins
       DO n=1,atoms%ntype
          enpara%el1(0:3,n,ispin)=regCharges%ener(0:3,n,ispin)/regCharges%sqal(0:3,n,ispin)
          IF (atoms%nlo(n)>0) enpara%ello1(:atoms%nlo(n),n,ispin)=regCharges%enerlo(:atoms%nlo(n),n,ispin)/regCharges%sqlo(:atoms%nlo(n),n,ispin)
       END DO
       IF (input%film) enpara%evac1(:vacuum%nvac,ispin)=regCharges%pvac(:vacuum%nvac,ispin)/regCharges%svac(:vacuum%nvac,ispin)
    END DO
  END SUBROUTINE calcOutParams

END MODULE m_types_enpara
