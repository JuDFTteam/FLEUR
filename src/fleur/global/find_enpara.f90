!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_find_enpara
  USE m_judft
  IMPLICIT NONE
  PRIVATE
  PUBLIC:: find_enpara

CONTAINS

  !> Function to determine the energy parameter given the quantum number and the potential
  !! Different schemes are implemented. Nqn (main quantum number) is used as a switch.
  !! This code was previously in lodpot.f
  REAL FUNCTION find_enpara(lo,l,n,jsp,nqn,atoms,vr,e_lo,e_up,l_scalar_relativi)RESULT(e)
    USE m_types_setup
    USE m_radsra
    USE m_differ
    USE m_constants
    IMPLICIT NONE
    LOGICAL,INTENT(IN):: lo
    INTEGER,INTENT(IN):: l,n,nqn,jsp
    REAL,INTENT(OUT) :: e_lo,e_up
    TYPE(t_atoms),INTENT(IN)::atoms
    REAL,INTENT(IN):: vr(:)
    LOGICAL, OPTIONAL, INTENT(IN) :: l_scalar_relativi
    
    IF(PRESENT(l_scalar_relativi)) THEN
       e = priv_scalar_relativi(lo,l,n,jsp,nqn,atoms,vr,e_lo,e_up)
       RETURN
    END IF

    IF (nqn>0) e=priv_method1(lo,l,n,jsp,nqn,atoms,vr,e_lo,e_up)
    IF (nqn<0) e=priv_method2(lo,l,n,jsp,nqn,atoms,vr,e_lo,e_up)
  END FUNCTION find_enpara


  REAL FUNCTION priv_method1(lo,l,n,jsp,nqn,atoms,vr,e_lo,e_up)RESULT(e)
    USE m_types_setup
    USE m_radsra
    USE m_differ
    USE m_constants
    IMPLICIT NONE
    LOGICAL,INTENT(IN):: lo
    INTEGER,INTENT(IN):: l,n,nqn,jsp
    REAL,INTENT(OUT)  :: e_lo,e_up
    TYPE(t_atoms),INTENT(IN)::atoms
    REAL,INTENT(IN):: vr(:)


    INTEGER j,ilo,i
    INTEGER nodeu,node,ierr,msh
    REAL   lnd,e1, e_up_local, e_lo_local
    REAL   d,rn,fl,fn,fj,t2,rr,t1,ldmt,us,dus,c
    !     ..
    !     .. Local Arrays ..
    REAL, ALLOCATABLE :: f(:,:),vrd(:)
    CHARACTER(LEN=20)    :: attributes(6)
    c=c_light(1.0)
    !Core potential setup done for each n,l now
    d = EXP(atoms%dx(n))
    ! set up core-mesh
    rn = atoms%rmt(n)
    msh = atoms%jri(n)
    DO WHILE (rn < atoms%rmt(n) + 20.0)
       msh = msh + 1
       rn = rn*d
    ENDDO
    rn = atoms%rmsh(1,n)*( d**(msh-1) )
    ALLOCATE ( f(msh,2),vrd(msh) )
    f = 0.0
    vrd = 0.0
    ! extend core potential (linear with slope t1 / a.u.)
    vrd(:atoms%jri(n))=vr(:atoms%jri(n))
    t1=0.125
    t2 = vrd(atoms%jri(n))/atoms%rmt(n) - atoms%rmt(n)*t1
    rr = atoms%rmt(n)
    DO j = atoms%jri(n) + 1, msh
       rr = d*rr
       vrd(j) = rr*( t2 + rr*t1 )
    ENDDO

    node = nqn - (l+1)
    IF (node<0) CALL judft_error("Error in setup of energy-parameters",hint="This could e.g. happen if you try to use 1p-states")
    e = 0.0
    ! determine upper edge
    nodeu = -1
    
    DO WHILE ( nodeu <= node )
       CALL radsra(e,l,vr(:),atoms%rmsh(1,n),&
            atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
       IF ( nodeu.LE.node ) THEN
          e = e + 10.0
          IF (e>1E10) CALL judft_error("Determination of energy parameters did not converge",hint="Perhaps your potential is broken?")
       ENDIF
    ENDDO
    
    e_up_local = e
    
    DO WHILE (nodeu.GT.node)
       CALL radsra(e,l,vr(:),atoms%rmsh(1,n),&
            atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
       IF ( nodeu.GT.node ) THEN
          e = e - 10.0
          IF (e<-1E10) CALL judft_error("Determination of energy parameters did not converge",hint="Perhaps your potential is broken?")
       ENDIF
    ENDDO
    
    e_lo_local = e

    DO WHILE ( (e_up_local-e_lo_local).GT.1.0e-5)
       e  = (e_up_local+e_lo_local) / 2.0
       CALL radsra(e,l,vr(:),atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),c,us,dus,nodeu,f(:,1),f(:,2))
       IF (nodeu.GT.node) THEN
          e_up_local = e
       ELSE
          e_lo_local = e
       END IF
    END DO
    
    e_up = (e_up_local+e_lo_local) / 2.0

    IF (node /= 0) THEN
       ! determine lower edge
       
       e_up_local = e

       DO WHILE (nodeu.GE.node)
          CALL radsra(e,l,vr(:),atoms%rmsh(1,n),&
               atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
          e = e - 10.0       
       END DO
       
       e_lo_local = e

       DO WHILE ( (e_up_local-e_lo_local).GT.1.0e-5)
          e  = (e_up_local+e_lo_local) / 2.0
          CALL radsra(e,l,vr(:),atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),c,us,dus,nodeu,f(:,1),f(:,2))
          IF (nodeu.GE.node) THEN
             e_up_local = e
          ELSE
             e_lo_local = e
          END IF
       END DO
       
       e_lo = (e_up_local+e_lo_local) / 2.0       
    ELSE
       e_lo = -9.99
    ENDIF
    ! calculate core
    e  = (e_up+e_lo)/2
    
    fn = REAL(nqn) ; fl = REAL(l) ; fj = fl + 0.5
    CALL differ(fn,fl,fj,c,atoms%zatom(n),atoms%dx(n),atoms%rmsh(1,n),&
         rn,d,msh,vrd, e, f(:,1),f(:,2),ierr)
    IF (lo.AND.l>0) THEN
       e1  = (e_up+e_lo)/2
       fn = REAL(nqn) ; fl = REAL(l) ; fj = fl-0.5
       CALL differ(fn,fl,fj,c,atoms%zatom(n),atoms%dx(n),atoms%rmsh(1,n),&
            rn,d,msh,vrd, e1, f(:,1),f(:,2),ierr)
       e = (2.0*e + e1 ) / 3.0
    ENDIF
  END FUNCTION priv_method1

  REAL FUNCTION priv_method2(lo,l,n,jsp,nqn,atoms,vr,e_lo,e_up)RESULT(e)
    USE m_types_setup
    USE m_radsra
    USE m_differ
    USE m_constants
    IMPLICIT NONE
    LOGICAL,INTENT(IN):: lo
    INTEGER,INTENT(IN):: l,n,nqn,jsp
    REAL,INTENT(OUT)  :: e_lo,e_up
    TYPE(t_atoms),INTENT(IN)::atoms
    REAL,INTENT(IN):: vr(:)

    INTEGER j,ilo,i
    INTEGER nodeu,node,ierr,msh
    REAL   lnd,e_up_temp,e_lo_temp,large_e_step
    REAL   d,rn,fl,fn,fj,t2,rr,t1,ldmt,us,dus,c
    !     ..
    !     .. Local Arrays ..

    REAL, ALLOCATABLE :: f(:,:),vrd(:)
    CHARACTER(LEN=20)    :: attributes(6)

    c=c_light(1.0)

    !Core potential setup done for each n,l now
    d = EXP(atoms%dx(n))
    ! set up core-mesh
    rn = atoms%rmt(n)
    msh = atoms%jri(n)
    DO WHILE (rn < atoms%rmt(n) + 20.0)
       msh = msh + 1
       rn = rn*d
    ENDDO
    rn = atoms%rmsh(1,n)*( d**(msh-1) )
    ALLOCATE ( f(msh,2),vrd(msh) )
    f = 0.0
    vrd = 0.0
    ! extend core potential (linear with slope t1 / a.u.)
    vrd(:atoms%jri(n))=vr(:atoms%jri(n))
    t1=0.125
    t2 = vrd(atoms%jri(n))/atoms%rmt(n) - atoms%rmt(n)*t1
    rr = atoms%rmt(n)
    DO j = atoms%jri(n) + 1, msh
       rr = d*rr
       vrd(j) = rr*( t2 + rr*t1 )
    ENDDO
    ! search for branches
    node = ABS(nqn) - (l+1)
    IF (node<0) CALL judft_error("Error in setup of energy-parameters",hint="This could e.g. happen if you try to use 1p-states")
    e = 0.0 ! The initial value of e is arbitrary.
    large_e_step = 5.0 ! 5.0 Htr steps for coarse energy searches

    ! determine upper band edge
    ! Step 1: Coarse search for the band edge
    CALL radsra(e,l,vr(:),atoms%rmsh(1,n),&
         atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
    DO WHILE ( nodeu > node )
       e = e - large_e_step
       CALL radsra(e,l,vr(:),atoms%rmsh(1,n),&
            atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
    END DO
    DO WHILE ( nodeu <= node )
       e = e + large_e_step
       CALL radsra(e,l,vr(:),atoms%rmsh(1,n),&
            atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
    END DO
    e_up_temp = e
    e_lo_temp = e - large_e_step
    ! Step 2: Fine band edge determination by bisection search
    DO WHILE ((e_up_temp - e_lo_temp) > 1e-2)
       e = (e_up_temp + e_lo_temp) / 2.0
       CALL radsra(e,l,vr(:),atoms%rmsh(1,n),&
            atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
       IF (nodeu > node) THEN
          e_up_temp = e
       ELSE
          e_lo_temp = e
       END IF
    END DO
    e_up = (e_up_temp + e_lo_temp) / 2.0
    e    = e_up

    ! determine lower band edge
    IF (node == 0) THEN
       e_lo = -49.99
    ELSE
       ! Step 1: Coarse search for the band edge
       nodeu = node
       DO WHILE ( nodeu >= node )
          e = e - large_e_step
          CALL radsra(e,l,vr(:),atoms%rmsh(1,n),&
               atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
       ENDDO
       e_up_temp = e + large_e_step
       e_lo_temp = e
       ! Step 2: Fine band edge determination by bisection search
       DO WHILE ((e_up_temp - e_lo_temp) > 1e-2)
          e = (e_up_temp + e_lo_temp) / 2.0
          CALL radsra(e,l,vr(:),atoms%rmsh(1,n),&
               atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
          IF (nodeu < node) THEN
             e_lo_temp = e
          ELSE
             e_up_temp = e
          END IF
       END DO
       e_lo = (e_up_temp + e_lo_temp) / 2.0
    END IF


    ! determince notches by intersection
    ldmt= -99.0 !ldmt = logarithmic derivative @ MT boundary
    lnd = -l-1
    DO WHILE ( ABS(ldmt-lnd) .GE. 1E-07)
       e = (e_up+e_lo)/2
       CALL radsra(e,l,vr(:),atoms%rmsh(1,n),&
            atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
       ldmt = dus/us
       IF( ldmt .GT. lnd) THEN
          e_lo = e
       ELSE IF( ldmt .LT. lnd ) THEN
          e_up = e
          e_lo = e_lo
       END IF
    END DO

     END FUNCTION priv_method2

  REAL FUNCTION priv_scalar_relativi(lo,l,n,jsp,nqn,atoms,vr,e_lo,e_up)RESULT(e)
    USE m_types_setup
    USE m_radsra
    USE m_differ
    USE m_constants
    IMPLICIT NONE
    LOGICAL,INTENT(IN):: lo
    INTEGER,INTENT(IN):: l,n,nqn,jsp
    REAL,INTENT(OUT)  :: e_lo,e_up
    TYPE(t_atoms),INTENT(IN)::atoms
    REAL,INTENT(IN):: vr(:)


    INTEGER j,ilo,i
    INTEGER nodeu,node,ierr,msh
    REAL   lnd,e1
    REAL   d,rn,fl,fn,fj,t2,rr,t1,ldmt,us,dus,c
    !     ..
    !     .. Local Arrays ..
    REAL, ALLOCATABLE :: f(:,:),vrd(:)
    CHARACTER(LEN=20)    :: attributes(6)
    c=c_light(1.0)

    !Core potential setup done for each n,l now
    d = EXP(atoms%dx(n))
    ! set up core-mesh
    rn = atoms%rmt(n)
    msh = atoms%jri(n)
    DO WHILE (rn < atoms%rmt(n) + 8.0)
       msh = msh + 1
       rn = rn*d
    ENDDO
    rn = atoms%rmsh(1,n)*( d**(msh-1) )
    ALLOCATE ( f(msh,2),vrd(msh) )
    f = 0.0
    vrd = 0.0
    ! extend core potential (linear with slope t1 / a.u.)
    vrd(:atoms%jri(n))=vr(:atoms%jri(n))
    t1=0.125
    t2 = vrd(atoms%jri(n))/atoms%rmt(n) - atoms%rmt(n)*t1
    rr = atoms%rmt(n)
    DO j = atoms%jri(n) + 1, msh
       rr = d*rr
       vrd(j) = rr*( t2 + rr*t1 )
    ENDDO

    node = nqn - (l+1)
    IF (node<0) CALL judft_error("Error in setup of energy-parameters",hint="This could e.g. happen if you try to use 1p-states")
    e_up = 0.0
    ! determine upper edge
    nodeu = -1
    DO WHILE ( nodeu.LE.node )
       CALL radsra(e_up,l,vr(:),atoms%rmsh(1,n),&
            atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
       IF (nodeu.LE.node) THEN
          e_up = e_up + 5.0
       END IF
    ENDDO
    
    e_lo = e_up - 50.0

    DO WHILE ( nodeu.GT.node )
       CALL radsra(e_lo,l,vr(:),atoms%rmsh(1,n),&
            atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
       IF (nodeu.GT.node) THEN
          e_lo = e_lo - 50.0
       END IF
    ENDDO
    
    e_lo = e_lo - 0.1
    DO WHILE ( (e_up-e_lo).GT.1.0e-10)
       e  = (e_up+e_lo) / 2.0
       CALL radsra(e,l,vrd,atoms%rmsh(1,n),atoms%dx(n),msh,c,us,dus,nodeu,f(:,1),f(:,2))
       IF (nodeu.GT.node) THEN
          e_up = e
       ELSE
          e_lo = e
       END IF
    END DO
    e = (e_up+e_lo) / 2.0
    
  END FUNCTION priv_scalar_relativi


END MODULE m_find_enpara
