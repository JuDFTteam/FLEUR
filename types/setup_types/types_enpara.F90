!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_enpara
  USE m_judft
  USE m_types_fleur_setup
  use m_json_tools
  IMPLICIT NONE
  TYPE,EXTENDS(t_fleursetup):: t_enpara
     !energy parameters actually used
     REAL, ALLOCATABLE CPP_MANAGED   :: el0(:,:,:)
     REAL                 :: evac0(2,2)
     REAL, ALLOCATABLE    :: ello0(:,:,:)
     !energy parameters calculated from output
     REAL, ALLOCATABLE    :: el1(:,:,:)  
     REAL                 :: evac1(2,2)
     REAL, ALLOCATABLE    :: ello1(:,:,:)
     !energy parameters determined by these quantum numbers
     INTEGER,ALLOCATABLE  :: qn_el(:,:,:)    !>if these are .ne.0 they are understood as
     INTEGER,ALLOCATABLE  :: qn_ello(:,:,:)  !>quantum numbers
     
     REAL                 :: evac(2,2)
     REAL, ALLOCATABLE    :: enmix(:)
     INTEGER, ALLOCATABLE :: skiplo(:,:)
     LOGICAL, ALLOCATABLE :: lchange(:,:,:)
     LOGICAL, ALLOCATABLE :: lchg_v(:,:)
     LOGICAL, ALLOCATABLE :: llochg(:,:,:)
     REAL                 :: epara_min
     LOGICAL              :: ready ! are the enpara's ok for calculation?
     LOGICAL              :: floating !floating energy parameters are relative to potential
   CONTAINS
     PROCEDURE,PASS :: broadcast=>broadcast_enpara
     PROCEDURE,PASS :: write=>WRITE_enpara
     PROCEDURE,PASS :: read=>READ_enpara
     PROCEDURE,PASS :: read_xml=>read_xml_enpara

     PROCEDURE :: init
     PROCEDURE :: update
     !PROCEDURE :: read
     !PROCEDURE :: write
     PROCEDURE :: mix
     PROCEDURE :: calcOutParams
  END TYPE t_enpara
     
  END TYPE t_enpara

CONTAINS
  SUBROUTINE broadcast_enpara(tt,mpi_comm,origin)
    IMPLICIT NONE
    CLASS(t_enpara),INTENT(INOUT):: tt
    INTEGER,INTENT(IN)               :: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL      :: origin

#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER :: pe,ierr
    
    IF (PRESENT(origin)) THEN
       pe=origin
    ELSE
       pe=0
    ENDIF

    CALL mpi_bc(tt%el0,pe,mpi_comm)
    CALL mpi_bc(tt%evac0,pe,mpi_comm)
    CALL mpi_bc(tt%epara_min,pe,mpi_comm)
    CALL mpi_bc(tt%ready,pe,mpi_comm)
    CALL mpi_bc(tt%floating,pe,mpi_comm)
    CALL mpi_bc(tt%ello0,pe,mpi_comm)
    CALL mpi_bc(tt%qn_el,pe,mpi_comm)
    CALL mpi_bc(tt%qn_ello,pe,mpi_comm)
    CALL mpi_bc(tt%evac,pe,mpi_comm)
    CALL mpi_bc(tt%enmix,pe,mpi_comm)
    CALL mpi_bc(tt%skiplo,pe,mpi_comm)
    CALL mpi_bc(tt%lchange,pe,mpi_comm)
    CALL mpi_bc(tt%lchg_v,pe,mpi_comm)
    CALL mpi_bc(tt%llochg,pe,mpi_comm)
    
#endif
      
  END SUBROUTINE broadcast_enpara

  SUBROUTINE write_enpara(tt, unit, iotype, v_list, iostat, iomsg)
    IMPLICIT NONE
    CLASS(t_enpara),INTENT(IN):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    WRITE(unit,*,IOSTAT=iostat) '"enpara":{'

    
call json_print(unit,"el0",tt%el0)
call json_print(unit,"evac0",tt%evac0)
call json_print(unit,"epara_min",tt%epara_min)
call json_print(unit,"ready",tt%ready)
call json_print(unit,"floating",tt%floating)
call json_print(unit,"ello0",tt%ello0)
call json_print(unit,"qn_el",tt%qn_el)
call json_print(unit,"qn_ello",tt%qn_ello)
call json_print(unit,"evac",tt%evac)
call json_print(unit,"enmix",tt%enmix)
call json_print(unit,"skiplo",tt%skiplo)
call json_print(unit,"lchange",tt%lchange)
call json_print(unit,"lchg_v",tt%lchg_v)
call json_print(unit,"llochg",tt%llochg)

    WRITE(unit,*,IOSTAT=iostat) '}'
    
  END SUBROUTINE write_enpara
  SUBROUTINE read_enpara(tt, unit, iotype, v_list, iostat, iomsg)
    IMPLICIT NONE
    CLASS(t_enpara),INTENT(INOUT):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    CHARACTER(len=40)::string
    REAL,allocatable:: rtemp(:)
    CALL json_open_class("enpara",unit,iostat)
    IF (iostat.NE.0)   RETURN

  
call json_read(unit,"el0",tt%el0)
call json_read(unit,"evac0",tt%evac0)
call json_read(unit,"epara_min",tt%epara_min)
call json_read(unit,"ready",tt%ready)
call json_read(unit,"floating",tt%floating)
call json_read(unit,"ello0",tt%ello0)
call json_read(unit,"qn_el",tt%qn_el)
call json_read(unit,"qn_ello",tt%qn_ello)
call json_read(unit,"evac",tt%evac)
call json_read(unit,"enmix",tt%enmix)
call json_read(unit,"skiplo",tt%skiplo)
call json_read(unit,"lchange",tt%lchange)
call json_read(unit,"lchg_v",tt%lchg_v)
CALL json_read(unit,"llochg",tt%llochg)

    CALL json_close_class(unit,iostat)
    
  END SUBROUTINE read_enpara

 


  SUBROUTINE read_xml_enpara(tt)
    USE m_xmlIntWrapFort
    USE m_constants
    USE m_calculator
    USE m_inv3
    IMPLICIT NONE
    CLASS(t_enpara),INTENT(OUT):: tt


    ntype= xmlGetNumberOfNodes('/fleurInput/atomGroups/atomGroup')
    DO itype=1,ntype
       xpaths=inp_xml_speciesxpath_for_group(n)
       speciesEParams(0) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/energyParameters/@s'))
       speciesEParams(1) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/energyParameters/@p'))
       speciesEParams(2) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/energyParameters/@d'))
       speciesEParams(3) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/energyParameters/@f'))
       DO jsp = 1, input%jspins
          DO l = 0, 3
             enpara%el0(l,iType,jsp) = speciesEParams(l)
             IF (enpara%el0(l,iType,jsp)==NINT(enpara%el0(l,iType,jsp))) THEN
                enpara%qn_el(l,iType,jsp)=NINT(enpara%el0(l,iType,jsp))
                enpara%el0(l,iType,jsp)=0
             ELSE
                enpara%qn_el(l,iType,jsp)=0
             ENDIF
          END DO
          DO l = 4,atoms%lmax(iType)
             enpara%el0(l,iType,jsp) = enpara%el0(3,iType,jsp)
          END DO
       END DO
       !now the LOs
       DO ilo = 1,xmlGetNumberOfNodes(TRIM(ADJUSTL(xpaths))//'/lo')+xmlGetNumberOfNodes(TRIM(ADJUSTL(xpathg))//'/lo')
          IF (ilo>xmlGetNumberOfNodes(TRIM(ADJUSTL(xpaths))//'/lo')) THEN
             WRITE(xpath,*) TRIM(ADJUSTL(xpathg))//'/lo[',ilo-xmlGetNumberOfNodes(TRIM(ADJUSTL(xpaths))//'/lo'),']'
          ELSE
             WRITE(xpath,*) TRIM(ADJUSTL(xpaths))//'/lo[',ilo,']'
          END IF
          .....
          
    
  END SUBROUTINE read_xml_enpara


    SUBROUTINE init(this,atoms,jspins,l_defaults)
    USE m_types_setup
    USE m_constants
    CLASS(t_enpara),INTENT(inout):: this
    TYPE(t_atoms),INTENT(IN)     :: atoms
    INTEGER,INTENT(IN)           :: jspins
    LOGICAL,INTENT(IN),OPTIONAL  :: l_defaults

    INTEGER :: n,i,jsp,l

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

    IF (PRESENT(l_defaults)) THEN
       IF (.NOT.l_defaults) RETURN
    ENDIF
    !Set most simple defaults
    DO jsp=1,jspins
       DO n = 1,atoms%ntype
          IF ( atoms%nz(n) < 3 ) THEN
             this%qn_el(0:3,n,jsp) =  (/1,2,3,4/) 
          ELSEIF ( atoms%nz(n) < 11 ) THEN
             this%qn_el(0:3,n,jsp) =  (/2,2,3,4/) 
          ELSEIF ( atoms%nz(n) < 19 ) THEN
             this%qn_el(0:3,n,jsp) =  (/3,3,3,4/) 
          ELSEIF ( atoms%nz(n) < 31 ) THEN
             this%qn_el(0:3,n,jsp) =  (/4,4,3,4/) 
          ELSEIF ( atoms%nz(n) < 37 ) THEN
             this%qn_el(0:3,n,jsp) =  (/4,4,4,4/) 
          ELSEIF ( atoms%nz(n) < 49 ) THEN
             this%qn_el(0:3,n,jsp) =  (/5,5,4,4/) 
          ELSEIF ( atoms%nz(n) < 55 ) THEN
             this%qn_el(0:3,n,jsp) =  (/5,5,5,4/) 
          ELSEIF ( atoms%nz(n) < 72 ) THEN
             this%qn_el(0:3,n,jsp) =  (/6,6,5,4/) 
          ELSEIF ( atoms%nz(n) < 81 ) THEN
             this%qn_el(0:3,n,jsp) =  (/6,6,5,5/) 
          ELSEIF ( atoms%nz(n) < 87 ) THEN
             this%qn_el(0:3,n,jsp) =  (/6,6,6,5/) 
          ELSE
             this%qn_el(0:3,n,jsp) =  (/7,7,6,5/) 
          ENDIF
          
          DO i = 1, atoms%nlo(n)
             IF (atoms%llo(i,n)<0) THEN
                !llo might not be initialized
                !in this case defaults broken
                this%qn_ello(i,n,jsp) = 0
                this%skiplo(n,jsp) = 0
             ELSE
                this%qn_ello(i,n,jsp) = this%qn_el(atoms%llo(i,n),n,jsp) - 1
                this%skiplo(n,jsp) = this%skiplo(n,jsp) + (2*atoms%llo(i,n)+1)
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    this%evac0=eVac0Default_const

  END SUBROUTINE init
  !> This subroutine adjusts the energy parameters to the potential. In particular, it
  !! calculated them in case of qn_el>-1,qn_ello>-1
  !! Before this was done in lodpot.F
  SUBROUTINE update(enpara,mpi,atoms,vacuum,input,v)
    USE m_types_setup
    USE m_types_mpi
    USE m_xmlOutput
    USE m_types_potden
    USE m_find_enpara
    CLASS(t_enpara),INTENT(inout):: enpara
    TYPE(t_mpi),INTENT(IN)      :: mpi
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_vacuum),INTENT(IN)   :: vacuum
    TYPE(t_potden),INTENT(IN)   :: v


    LOGICAL ::  l_enpara
    LOGICAL ::  l_done(0:atoms%lmaxd,atoms%ntype,input%jspins)
    LOGICAL ::  lo_done(atoms%nlod,atoms%ntype,input%jspins)
    REAL    ::  vbar,vz0,rj
    INTEGER ::  n,jsp,l,ilo,j,ivac
    CHARACTER(LEN=20)    :: attributes(5)

    IF (mpi%irank  == 0) CALL openXMLElement('energyParameters',(/'units'/),(/'Htr'/))

    l_done = .FALSE.;lo_done=.FALSE.
    DO jsp = 1,input%jspins
       !$OMP PARALLEL DO DEFAULT(none) &
       !$OMP SHARED(atoms,enpara,jsp,l_done,mpi,v,lo_done) &
       !$OMP PRIVATE(n,l,ilo)
       !! First calculate energy paramter from quantum numbers if these are given...
       !! l_done stores the index of those energy parameter updated
       DO n = 1, atoms%ntype
          DO l = 0,3
             IF( enpara%qn_el(l,n,jsp).ne.0)THEN 
                l_done(l,n,jsp) = .TRUE.
                enpara%el0(l,n,jsp)=find_enpara(.FALSE.,l,n,jsp,enpara%qn_el(l,n,jsp),atoms,mpi,v%mt(:,0,n,jsp))
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
                enpara%ello0(ilo,n,jsp)=find_enpara(.TRUE.,l,n,jsp,enpara%qn_ello(ilo,n,jsp),atoms,mpi,v%mt(:,0,n,jsp))
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
             IF (mpi%irank.EQ.0) THEN
                WRITE ( 6,'('' spin'',i2,'', atom type'',i3,'' ='',f12.6,''   r='',f8.5)') jsp,n,vbar,rj
                WRITE (16,'('' spin'',i2,'', atom type'',i3,'' ='',f12.6,''   r='',f8.5)') jsp,n,vbar,rj
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
                IF (mpi%irank.EQ.0) THEN
                   WRITE ( 6,'('' spin'',i2,'', vacuum   '',i3,'' ='',f12.6)') jsp,ivac,vz0 
                   WRITE (16,'('' spin'',i2,'', vacuum   '',i3,'' ='',f12.6)') jsp,ivac,vz0
                ENDIF
             ENDIF
             enpara%evac(ivac,jsp) = enpara%evac0(ivac,jsp) + vz0
             IF (.NOT.l_enpara) THEN
                enpara%evac(ivac,jsp) = v%vacz(vacuum%nmz,ivac,jsp) + enpara%evac0(ivac,jsp)
             END IF
             IF (mpi%irank.EQ.0) THEN
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
    
    IF (mpi%irank  == 0) CALL closeXMLElement('energyParameters')
  END SUBROUTINE update
  SUBROUTINE mix(enpara,mpi,atoms,vacuum,input,vr,vz)
    !------------------------------------------------------------------
    USE m_types_setup
    USE m_types_mpi
    IMPLICIT NONE
    CLASS(t_enpara),INTENT(INOUT)  :: enpara
    TYPE(t_mpi),INTENT(IN)         :: mpi
    TYPE(t_atoms),INTENT(IN)       :: atoms
    TYPE(t_vacuum),INTENT(IN)      :: vacuum
    TYPE(t_input),INTENT(IN)       :: input

    REAL,    INTENT(IN) :: vr(:,:,:)
    REAL,    INTENT(IN) :: vz(vacuum%nmz,2)

    INTEGER ityp,j,l,lo,jsp,n
    REAL    vbar,maxdist,maxdist2
    INTEGER same(atoms%nlod)
    LOGICAL l_enpara
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER :: ierr
#endif    

    IF (mpi%irank==0) THEN
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
    IF (mpi%irank==0.AND.l_enpara) CALL enpara%WRITE(atoms,input%jspins,input%film)
#ifdef CPP_MPI
    CALL MPI_BCAST(enpara%el0,SIZE(enpara%el0),MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(enpara%ello0,SIZE(enpara%ello0),MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
    CALL MPI_BCAST(enpara%evac,SIZE(enpara%evac),MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
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
