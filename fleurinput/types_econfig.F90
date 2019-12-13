!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_econfig
  USE m_judft
  IMPLICIT NONE
  PRIVATE
  !This is used by t_atoms and does not extend t_fleurinput_base by itself
  TYPE:: t_econfig
     CHARACTER(len=200) :: coreconfig
     CHARACTER(len=100) :: valenceconfig
     INTEGER            :: num_core_states
     INTEGER            :: num_states
     INTEGER,ALLOCATABLE:: nprnc(:)
     INTEGER,ALLOCATABLE:: kappa(:)
     REAL,ALLOCATABLE   :: occupation(:,:)
     REAL               :: core_electrons
     REAL               :: valence_electrons
   CONTAINS
     PROCEDURE :: init_num,init_simple,init_nz,init_all
     GENERIC   :: init => init_num,init_simple,init_nz,init_all
     PROCEDURE :: set_occupation
     PROCEDURE :: set_initial_moment
     PROCEDURE :: broadcast
     PROCEDURE :: get_state_string
     PROCEDURE :: get_core
  END TYPE t_econfig
  PUBLIC :: t_econfig
CONTAINS

  SUBROUTINE get_core(econf,nst,nprnc,kappa,occupation,l_valence)
    CLASS(t_econfig),INTENT(IN)  :: econf
    INTEGER         ,INTENT(out) :: nst
    INTEGER         ,INTENT(out) :: nprnc(:),kappa(:)
    REAL            ,INTENT(out) :: occupation(:,:)
    LOGICAL,OPTIONAL,INTENT(IN)  :: l_valence

    nst=econf%num_core_states
    if (present(l_valence)) then
       if (l_valence) nst=econf%num_states
    endif
    nprnc(:nst)=econf%nprnc(:nst)
    kappa(:nst)=econf%kappa(:nst)
    occupation(:nst,:)=econf%occupation(:nst,:SIZE(occupation,2))
    if (size(occupation,2)==1) occupation=occupation*2
  END SUBROUTINE get_core



  FUNCTION get_state_string(econf,i)RESULT(str)
    CLASS(t_econfig),INTENT(IN):: econf
    INTEGER,INTENT(in)         :: i
    CHARACTER(len=7)::str

    CHARACTER(len=5)::s

    SELECT CASE (econf%kappa(i))
    CASE (-1)
       s="s1/2)"
    CASE (1)
       s="p1/2)"
    CASE (-2)
       s="p3/2)"
    CASE (2)
       s="d3/2)"
    CASE (-3)
       s="d5/2)"
    CASE (3)
       s="f5/2)"
    CASE (-4)
       s="f7/2)"
    CASE default
       call judft_error("Invalid reqest for string with kappa")
    END SELECT


    WRITE(str,"(a1,i1,a)") "(",econf%nprnc(i),s
  END FUNCTION get_state_string


  SUBROUTINE broadcast(econf,irank,mpi_comm)
    USE m_mpi_bc_tool
    CLASS(t_econfig),INTENT(INOUT):: econf
    INTEGER,INTENT(in)            :: irank,mpi_comm
#ifdef CPP_MPI
    CALL mpi_bc(econf%num_core_states,irank,mpi_comm)
    CALL mpi_bc(econf%num_states,irank,mpi_comm)
    CALL mpi_bc(econf%nprnc,irank,mpi_comm)
    CALL mpi_bc(econf%kappa,irank,mpi_comm)
    CALL mpi_bc(econf%occupation,Irank,mpi_comm)
    CALL mpi_bc(econf%core_electrons,irank,mpi_comm)
    CALL mpi_bc(econf%valence_electrons,irank,mpi_comm)
#endif
  END SUBROUTINE broadcast


  SUBROUTINE init_num(econf,nc,nz)
    USE m_constants
    CLASS(t_econfig),INTENT(OUT):: econf
    INTEGER,INTENT(in)          :: nc,nz

    CHARACTER(len=200)::core
    INTEGER :: n
    core=''
    IF (nc>54) CALL judft_error("Setting the core by number only until Xe")
    DO n=1,nc
       core=core//" "//coreStateList_const(n)
    ENDDO

    CALL econf%init(core,nz)
  END SUBROUTINE init_num


  SUBROUTINE init_simple(econf,str)
    CLASS(t_econfig),INTENT(OUT):: econf
    CHARACTER(len=*),INTENT(IN) :: str

    CHARACTER(len=200)::core,valence
    INTEGER :: n
    REAL,allocatable :: core_occ(:),valence_occ(:)

    n=INDEX(str,"|")
    IF (n==0) CALL judft_error(("Invalid econfig:"//TRIM(str)))

    IF (INDEX(str,"|")==1) THEN
       ! No core
       core=""
       ALLOCATE(core_occ(0))
    ELSE
       CALL convert_to_extended(str(1:n-1),core,core_occ)
    ENDIF
    CALL convert_to_extended(str(INDEX(str,"|")+1:),valence,valence_occ)

    CALL econf%init(core,valence)
    !Now set occupations
    IF (SIZE(core_occ)>0) core=econf%coreconfig(INDEX(econf%coreconfig,"("):) !remove noble gas
    DO n=1,SIZE(core_occ)
       CALL econf%set_occupation(core(:INDEX(core,")")),core_occ(n),-1.)
       core=core(INDEX(core,")")+1:)
    ENDDO
    valence=econf%valenceconfig
    DO n=1,SIZE(valence_occ)
       CALL econf%set_occupation(valence(:INDEX(valence,")")),valence_occ(n),-1.)
       valence=valence(INDEX(valence,")")+1:)
    ENDDO

  END SUBROUTINE init_simple

  SUBROUTINE init_all(econf,core,valence)
    CLASS(t_econfig),INTENT(OUT)   :: econf
    CHARACTER(len=*),INTENT(INOUT) :: core,valence


    INTEGER :: np(40),kap(40)

    econf%coreconfig=core
    econf%valenceconfig=valence

    CALL expand_noble_gas(core) !extend noble gas config

    IF (VERIFY(core,"(1234567spdf/) ")>0) call judft_error(("Invalid econfig:"//TRIM(core)))
    IF (VERIFY(valence,"(1234567spdf/) ")>0) CALL judft_error(("Invalid econfig:"//TRIM(valence)))

    econf%num_core_states=0
    DO WHILE (len_TRIM(core)>1)
       CALL extract_next(core,np(econf%num_core_states+1),kap(econf%num_core_states+1))
       econf%num_core_states=econf%num_core_states+1
    ENDDO
    econf%num_states=econf%num_core_states
    DO WHILE (len_TRIM(valence)>1)
       CALL extract_next(valence,np(econf%num_states+1),kap(econf%num_states+1))
       econf%num_states=econf%num_states+1
    ENDDO
    ALLOCATE(econf%nprnc(econf%num_states),econf%kappa(econf%num_states))
    econf%nprnc=np(:econf%num_states)
    econf%kappa=kap(:econf%num_states)
    ALLOCATE(econf%occupation(econf%num_states,2))

    CALL econf%set_occupation("(1s1/2)",-1.,-1.)
  END SUBROUTINE init_all


  SUBROUTINE init_nz(econf,core,nz)
    use m_constants
    CLASS(t_econfig),INTENT(OUT)   :: econf
    CHARACTER(len=*),INTENT(INOUT) :: core
    INTEGER,INTENT(IN)             :: nz

    INTEGER :: charge,np(40),kap(40)
    CHARACTER(len=7)::str

    IF (nz>54) CALL judft_warn("Specifying no explicit valence config for systems with f-states might lead to broken configs")

    econf%coreconfig=core

    CALL expand_noble_gas(core) !extend noble gas config

    IF (VERIFY(core,"(1234567spdf/) ")>0) call judft_error(("Invalid econfig:"//TRIM(core)))

    econf%num_core_states=0
    DO WHILE (len_TRIM(core)>1)
       CALL extract_next(core,np(econf%num_core_states+1),kap(econf%num_core_states+1))
       econf%num_core_states=econf%num_core_states+1
    ENDDO
    econf%num_states=econf%num_core_states
    !valence charge
    charge=nz-SUM(ABS(kap(:econf%num_core_states))**2)*2
    DO  WHILE(charge>0) !Add valence
       str=coreStateList_const(econf%num_states+1)
       PRINT*,econf%num_states,str,charge
       CALL extract_next(str,np(econf%num_states+1),kap(econf%num_states+1))
       econf%num_states=econf%num_states+1
       charge=charge-ABS(kap(econf%num_states)**2)*2
    ENDDO

    ALLOCATE(econf%nprnc(econf%num_states),econf%kappa(econf%num_states))
    econf%nprnc=np(:econf%num_states)
    econf%kappa=kap(:econf%num_states)
    ALLOCATE(econf%occupation(econf%num_states,2))

    CALL econf%set_occupation("(1s1/2)",-1.,-1.)
    !last level might be partially occupied
    IF (charge<0) THEN
       str=coreStateList_const(econf%num_states)
       CALL econf%set_occupation(str,ABS(charge)*1.,-1.)
    END IF
  END SUBROUTINE init_nz

  SUBROUTINE set_occupation(econf,state,u,d)
    CLASS(t_econfig),INTENT(INOUT):: econf
    CHARACTER(len=*),INTENT(IN) :: state
    REAL,INTENT(in)             :: u,d!up,down

    INTEGER :: n,np,kappa
    CHARACTER(len=7)::str
    REAL :: up,down
    up=u
    down=d

    str=ADJUSTL(state)

    IF (up==-1. .AND. down==-1) THEN
       !Set all defaults
       econf%occupation=0.0
       DO n=1,econf%num_states
          econf%occupation(n,:)=ABS(econf%kappa(n))
       END DO
    ELSE
       CALL extract_next(str,np,kappa)
       DO n=1,econf%num_states
          IF (np==econf%nprnc(n).AND.kappa==econf%kappa(n)) THEN
             IF (down==-1) THEN
                IF (up==-1) CALL judft_error("Invalid econfig:"//TRIM(str))
                econf%occupation(n,:)=up/2.0
                up=-1
             ELSE
                econf%occupation(n,1)=up
                econf%occupation(n,2)=down
                up=-1;down=-1
             END IF
          END IF
       END DO
    END IF
    econf%core_electrons=SUM(econf%occupation(:econf%num_core_states,:))
    econf%valence_electrons=SUM(econf%occupation(econf%num_core_states+1:econf%num_states,:))
  END SUBROUTINE set_occupation


  SUBROUTINE set_initial_moment(econf,bmu)
    CLASS(t_econfig),INTENT(INOUT):: econf
    real            ,INTENT(IN  ) :: bmu

    INTEGER:: n
    REAL   :: p,el

    DO n=1,econf%num_states
       IF (ABS(econf%occupation(n,1)-econf%occupation(n,2))>1E-5) CALL judft_error("Could not add a moment to an already polarized configuration")
    END DO

    p=bmu/2.0
    DO n=econf%num_states,1,-1
       el=econf%occupation(n,1)
       el=MAX(el,ABS(econf%kappa(n)-el)) !Maximal charge to redistribtue
       IF (ABS(p)>el) THEN
          p=p-SIGN(el,p)
          econf%occupation(n,1)=econf%occupation(n,1)+SIGN(el,p)
          econf%occupation(n,2)=econf%occupation(n,2)-SIGN(el,p)
       ELSE
          econf%occupation(n,1)=econf%occupation(n,1)+p
          econf%occupation(n,2)=econf%occupation(n,2)-p
          p=0.0
       ENDIF
       IF (ABS(p)<1E-8) EXIT
    END DO
  END SUBROUTINE set_initial_moment




  SUBROUTINE extract_next(str,n,kappa)
    CHARACTER(len=*),INTENT(INOUT) :: str
    INTEGER,INTENT(out)         :: n,kappa

    str=ADJUSTL(str)
    READ(str(2:2),"(i1)") n
    SELECT CASE (str(3:6))
    CASE ("s1/2")
       kappa=-1
    CASE ("p1/2")
       kappa=1
    CASE ("p3/2")
       kappa=-2
    CASE ("d3/2")
       kappa=2
    CASE ("d5/2")
       kappa=-3
    CASE ("f5/2")
       kappa=3
    CASE ("f7/2")
       kappa=-4
    CASE default
       call judft_error(("Invalid econfig:"//TRIM(str)))
    END SELECT
    str=ADJUSTL(str(8:))
  END SUBROUTINE extract_next

  SUBROUTINE expand_noble_gas(str)
    CHARACTER(len=*),INTENT(INOUT) :: str

    INTEGER:: n
    CHARACTER(len=2)::elem

    n=INDEX(str,"[")
    IF (n==0) RETURN
    elem=str(n+1:n+2)
    str=str(n+4:)
    SELECT CASE (elem)
    CASE ("He","he","HE")
       str="(1s1/2) "//str
    CASE ("Ne","ne","NE")
       str="(1s1/2) (2s1/2) (2p1/2) "//str
    CASE ("Ar","ar","AR")
       str="(1s1/2) (2s1/2) (2p1/2) (2p3/2) (3s1/2) (3p1/2) (3p3/2) "//str
    CASE ("Kr","kr","KR")
       str="(1s1/2) (2s1/2) (2p1/2) (2p3/2) (3s1/2) (3p1/2) (3p3/2) (4s1/2) (3d3/2) (3d5/2) (4p1/2) (4p3/2)"//str
    CASE ("Xe","xe","XE")
       str="(1s1/2) (2s1/2) (2p1/2) (2p3/2) (3s1/2) (3p1/2) (3p3/2) (4s1/2) (3d3/2) (3d5/2) (4p1/2) (4p3/2) (5s1/2) (4d3/2) (4d5/2) (5p1/2) (5p3/2) "//str
    CASE ("Rn","rn","RN")
       str="(1s1/2) (2s1/2) (2p1/2) (2p3/2) (3s1/2) (3p1/2) (3p3/2) (4s1/2) (3d3/2) (3d5/2) (4p1/2) (4p3/2) (5s1/2) (4d3/2) (4d5/2) (5p1/2) (5p3/2) (6s1/2) (4f5/2) (4f7/2) (5d3/2) (5d5/2) (6p1/2) (6p3/2) "//str
    CASE default
       call judft_error(("Invalid econfig:"//TRIM(str)))
    END SELECT
  END SUBROUTINE expand_noble_gas





    SUBROUTINE convert_to_extended(simple,extended,occupations)
      CHARACTER(len=*),INTENT(IN)  :: simple
      CHARACTER(len=*),INTENT(OUT) :: extended
      REAL,ALLOCATABLE,INTENT(OUT) :: occupations(:)

      CHARACTER(len=200)::conf
      REAL :: occupation(200) !this is the tmp local variable (no 's')
      INTEGER:: n,nn,occ
      CHARACTER:: n_ch,ch
      extended=""

      conf=ADJUSTL(simple)
      n=INDEX(conf,"[")
      IF (n>0) THEN
         !copy noble gs string
         nn=INDEX(conf,"]")
         IF (nn==0) CALL judft_error(("Invalid econfig:"//TRIM(simple)))
         extended=conf(n:nn)//" "
         conf=ADJUSTL(conf(nn+1:))
      ENDIF
      IF (VERIFY(conf,"01234567890spdf ")>0) CALL judft_error(("Invalid econfig:"//TRIM(simple)))
      n=1
      DO WHILE (len_TRIM(conf)>1)
         READ(conf,"(a1,a1,i2)") n_ch,ch,occ
         SELECT CASE (ch)
         CASE ('s')
            extended=trim(extended)//" ("//n_ch//"s1/2)"
            occupation(n)=occ
            n=n+1
         CASE ('p')
            extended=trim(extended)//" ("//n_ch//"p1/2) ("//n_ch//"p3/2)"
            occupation(n)=occ*1./3.
            occupation(n+1)=occ*2./3.
            n=n+2
         CASE('d')
            extended=trim(extended)//" ("//n_ch//"d3/2) ("//n_ch//"d5/2)"
            occupation(n)=occ*2./5.
            occupation(n+1)=occ*3./5.
            n=n+2
         CASE('f')
            extended=trim(extended)//" ("//n_ch//"f5/2) ("//n_ch//"f7/2)"
            occupation(n)=occ*3./7.
            occupation(n+1)=occ*4./7.
            n=n+2
         CASE default
            CALL judft_error(("Invalid econfig:"//TRIM(simple)))
         END SELECT
         IF (INDEX(conf," ")>0) THEN
            conf=conf(INDEX(conf," ")+1:)
         ELSE
            conf=''
         ENDIF
         conf=ADJUSTL(conf)
      END DO
      ALLOCATE(occupations(n-1))
      occupations=occupation(:n-1)
    END SUBROUTINE convert_to_extended






  END MODULE m_types_econfig
