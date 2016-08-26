!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_lodpot
CONTAINS
  SUBROUTINE lodpot(mpi,atoms,sphhar,obsolete,vacuum,&
       input, vr,vz, enpara, enpara_new)
    !*********************************************************************
    !
    ! set el and evac from el0 & evac0 (depending on lepr)
    !
    !*********************************************************************
    USE m_constants, ONLY : c_light
    USE m_radsra
    USE m_differ
    USE m_types
    USE m_xmlOutput
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)      :: mpi
    TYPE(t_enpara),INTENT(IN)   :: enpara
    TYPE(t_obsolete),INTENT(IN) :: obsolete
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_vacuum),INTENT(IN)   :: vacuum
    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_atoms),INTENT(IN)    :: atoms
    Type(t_enpara),INTENT(OUT)  :: enpara_new
    !
    !     .. Scalar Arguments ..
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: vr(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd,input%jspins),vz(vacuum%nmzd,2,4)
    REAL :: el(0:atoms%lmaxd,atoms%ntypd,input%jspins) 
    REAL :: evac(2,input%jspins),ello(atoms%nlod,atoms%ntypd,input%jspins)

    !     ..
    !     .. Local Scalars ..
    INTEGER jsp,n,ivac,j,l,ilo,i
    INTEGER nodeu,node,ierr,msh
    REAL vbar,vz0,rj,e_up,e_lo,lnd,e_up_temp,e_lo_temp,large_e_step
    REAL   e,d,rn,fl,fn,fj,t2,rr,t1,ldmt,us,dus,c
    LOGICAL start
    !     ..
    !     .. Local Arrays .. 
    INTEGER nqn(0:atoms%lmaxd),nqn_lo(atoms%nlod)
    REAL, ALLOCATABLE :: f(:,:),vrd(:)
    LOGICAL l_done(0:atoms%lmaxd,atoms%ntypd,input%jspins)
    LOGICAL lo_done(atoms%nlod,atoms%ntypd,input%jspins)
    CHARACTER(len=1) :: ch(0:9)
    CHARACTER(LEN=20)    :: attributes(6)
    !     ..
  
    el(:,:,:) = 0.0 ; evac(:,:) = 0.0 ; ello(:,:,:) = 0.0
    l_done = .FALSE.
    c=c_light(1.0)
    IF (mpi%irank  == 0) CALL openXMLElement('energyParameters',(/'units'/),(/'Htr'/))
    IF ( obsolete%lepr == 0 ) THEN ! not for floating energy parameters
       e=1.0
       ch(0:9) = (/'s','p','d','f','g','h','i','j','k','l'/)
       DO jsp = 1,input%jspins
          IF( input%jspins .GT. 1 ) WRITE(6,'(A,i3)') ' Spin: ',jsp
!$OMP PARALLEL DO DEFAULT(none) &
!$OMP SHARED(atoms,enpara,jsp,l_done,mpi,vr,c,el,ch,lo_done,ello) &
!$OMP PRIVATE(n,nqn,nqn_lo,d,rn,msh,f,vrd,j,t1,t2,rr,l,node,nodeu,e,start,us,dus,e_lo,e_up) &
!$OMP PRIVATE(fl,fn,fj,ierr,attributes,large_e_step,e_up_temp,e_lo_temp,ldmt,lnd,ilo)
          DO n = 1, atoms%ntype
             ! check what to do ( for 'normal' energy parameters )
             nqn(0:3) = NINT( enpara%el0(0:3,n,jsp) )

             nqn_lo( 1:atoms%nlo(n)) = NINT( enpara%ello0(1:atoms%nlo(n),n,jsp) )

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

             ! extend core potential (linear with slope t1 / a.u.)

             DO j = 1, atoms%jri(n)
                vrd(j) = vr(j,0,n,jsp)
             ENDDO
             t1=0.125
             t2 = vrd(atoms%jri(n))/atoms%rmt(n) - atoms%rmt(n)*t1
             rr = atoms%rmt(n)
             DO j = atoms%jri(n) + 1, msh
                rr = d*rr
                vrd(j) = rr*( t2 + rr*t1 )
             ENDDO


             DO l = 0,3

                IF( ABS(enpara%el0(l,n,jsp) - nqn(l)) .LE. 1E-05&
                     &.AND. nqn(l) .GT. 0                        )THEN 

                   l_done(l,n,jsp) = .TRUE.
                   ! search for branches
                   node = nqn(l) - (l+1)
                                   e = 0.0 
                   ! determine upper edge
                   nodeu = -1 ; start = .TRUE.
                   DO WHILE ( nodeu <= node ) 
                      CALL radsra(e,l,vr(:,0,n,jsp),atoms%rmsh(1,n),&
                           atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
                      IF  ( ( nodeu > node ) .AND. start ) THEN
                         e = e - 1.0
                         nodeu = -1
                      ELSE
                         e = e + 0.01
                         start = .FALSE.
                      ENDIF
                   ENDDO

                   e_up = e
                   IF (node /= 0) THEN
                      ! determine lower edge
                      nodeu = node + 1
                      DO WHILE ( nodeu >= node ) 
                         CALL radsra(e,l,vr(:,0,n,jsp),atoms%rmsh(1,n),&
                              atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
                         e = e - 0.01
                      ENDDO
                      e_lo = e
                   ELSE
                      e_lo = -9.99 
                   ENDIF


                   ! calculate core

                   e  = (e_up+e_lo)/2
                   fn = REAL(nqn(l)) ; fl = REAL(l) ; fj = fl + 0.5
                   CALL differ(fn,fl,fj,c,atoms%zatom(n),atoms%dx(n),atoms%rmsh(1,n),&
                        rn,d,msh,vrd, e, f(:,1),f(:,2),ierr)
                   el(l,n,jsp) = e
                   IF (mpi%irank  == 0) THEN
                      attributes = ''
                      WRITE(attributes(1),'(i0)') n
                      WRITE(attributes(2),'(i0)') jsp
                      WRITE(attributes(3),'(i0,a1)') nqn(l), ch(l)
                      WRITE(attributes(4),'(f8.2)') e_lo
                      WRITE(attributes(5),'(f8.2)') e_up
                      WRITE(attributes(6),'(f16.10)') e
                      CALL writeXMLElementForm('atomicEP',(/'atomType     ','spin         ','branch       ',&
                                                            'branchLowest ','branchHighest','value        '/),&
                                               attributes,reshape((/12,4,6,12,13,5,6,1,3,8,8,16/),(/6,2/)))
                      WRITE(6,'(a6,i3,i2,a1,a12,f6.2,a3,f6.2,a13,f8.4)') '  Atom',n,nqn(l),ch(l),' branch from',&
                                                                         e_lo, ' to',e_up,' htr. ; e_l =',e
                   ENDIF
                   IF( l .EQ. 3 ) THEN
                      el(4:atoms%lmax(n),n,jsp) = el(3,n,jsp)
                      l_done(4:atoms%lmax(n),n,jsp) = .TRUE.
                   END IF

                ELSE IF( ABS(enpara%el0(l,n,jsp)- nqn(l)) .LE. 1E-05&
                     &.AND. nqn(l) .LT. 0 ) THEN

                   l_done(l,n,jsp) = .TRUE.
                   ! search for branches
                   node = ABS(nqn(l)) - (l+1)
                   e = 0.0 ! The initial value of e is arbitrary.
                   large_e_step = 5.0 ! 5.0 Htr steps for coarse energy searches

                   ! determine upper band edge
                   ! Step 1: Coarse search for the band edge
                   CALL radsra(e,l,vr(:,0,n,jsp),atoms%rmsh(1,n),&
                        atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
                   DO WHILE ( nodeu > node )
                      e = e - large_e_step
                      CALL radsra(e,l,vr(:,0,n,jsp),atoms%rmsh(1,n),&
                           atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
                   END DO
                   DO WHILE ( nodeu <= node )
                      e = e + large_e_step
                      CALL radsra(e,l,vr(:,0,n,jsp),atoms%rmsh(1,n),&
                           atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
                   END DO
                   e_up_temp = e
                   e_lo_temp = e - large_e_step
                   ! Step 2: Fine band edge determination by bisection search
                   DO WHILE ((e_up_temp - e_lo_temp) > 1e-2)
                      e = (e_up_temp + e_lo_temp) / 2.0
                      CALL radsra(e,l,vr(:,0,n,jsp),atoms%rmsh(1,n),&
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
                         CALL radsra(e,l,vr(:,0,n,jsp),atoms%rmsh(1,n),&
                              atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
                      ENDDO
                      e_up_temp = e + large_e_step
                      e_lo_temp = e
                      ! Step 2: Fine band edge determination by bisection search
                      DO WHILE ((e_up_temp - e_lo_temp) > 1e-2)
                         e = (e_up_temp + e_lo_temp) / 2.0
                         CALL radsra(e,l,vr(:,0,n,jsp),atoms%rmsh(1,n),&
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
                      CALL radsra(e,l,vr(:,0,n,jsp),atoms%rmsh(1,n),&
                           atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
                      ldmt = dus/us
                      IF( ldmt .GT. lnd) THEN
                         e_lo = e
                      ELSE IF( ldmt .LT. lnd ) THEN
                         e_up = e
                         e_lo = e_lo
                      END IF
                   END DO

                   IF (mpi%irank == 0) THEN
                      attributes = ''
                      WRITE(attributes(1),'(i0)') n
                      WRITE(attributes(2),'(i0)') jsp
                      WRITE(attributes(3),'(i0,a1)') ABS(nqn(l)), ch(l)
                      WRITE(attributes(4),'(f16.10)') ldmt
                      WRITE(attributes(5),'(f16.10)') e
                      CALL writeXMLElementForm('heAtomicEP',(/'atomType      ','spin          ','branch        ',&
                                                              'logDerivMT    ','value         '/),&
                                               attributes(1:5),reshape((/10,4,6,12,5+17,6,1,3,16,16/),(/5,2/)))
                      WRITE (6,'(a7,i3,i2,a1,a12,f7.2,a4,f7.2,a5)') "  Atom ",n,nqn(l),ch(l)," branch, D = ",&
                                                                    ldmt, " at ",e," htr."
                   ENDIF

                   el(l,n,jsp) = e

                   IF( l .EQ. 3 ) THEN
                      el(4:atoms%lmax(n),n,jsp) = el(3,n,jsp)
                      l_done(4:atoms%lmax(n),n,jsp) = .TRUE.
                   END IF
                ELSE 
                   l_done(l,n,jsp) = .FALSE.
                END IF
             ENDDO ! l

             ! Now for the lo's

             DO ilo = 1, atoms%nlo(n)
                l = atoms%llo(ilo,n)

                IF( ABS(enpara%ello0(ilo,n,jsp) - nqn_lo(ilo)) .LE. 1E-05&
                     &.AND. nqn_lo(ilo) .GT. 0) THEN

                   lo_done(ilo,n,jsp) = .TRUE.
                   ! search for branches
                   node = nqn_lo(ilo) - (l+1)

                   e = 0.0
                   ! determine upper edge
                   nodeu = -1 ; start = .TRUE.
                   DO WHILE ( nodeu <= node )
                      CALL radsra(e,l,vr(:,0,n,jsp),atoms%rmsh(1,n),&
                           atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
                      IF  ( ( nodeu > node ) .AND. start ) THEN
                         e = e - 1.0
                         nodeu = -1
                      ELSE
                         e = e + 0.01
                         start = .FALSE.
                      ENDIF
                   ENDDO
                   e_up = e

                   IF (node /= 0) THEN
                      ! determine lower edge
                      nodeu = node + 1
                      DO WHILE ( nodeu >= node )
                         CALL radsra(e,l,vr(:,0,n,jsp),atoms%rmsh(1,n),&
                              atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
                         e = e - 0.01
                      ENDDO
                      e_lo = e
                   ELSE
                      e_lo = -9.99
                   ENDIF

                   ! calculate core

                   e = (e_up+e_lo)/2
                   fn = REAL(nqn_lo(ilo)) ; fl = REAL(l) ; fj = fl + 0.5
                   CALL differ(fn,fl,fj,c,atoms%zatom(n),atoms%dx(n),atoms%rmsh(1,n),&
                        rn,d,msh,vrd, e, f(:,1),f(:,2),ierr)
                   ello(ilo,n,jsp) = e

                   IF (mpi%irank == 0) THEN
                      attributes = ''
                      WRITE(attributes(1),'(i0)') n
                      WRITE(attributes(2),'(i0)') jsp
                      WRITE(attributes(3),'(i0,a1)') nqn_lo(ilo), ch(l)
                      WRITE(attributes(4),'(f8.2)') e_lo
                      WRITE(attributes(5),'(f8.2)') e_up
                      WRITE(attributes(6),'(f16.10)') e
                      CALL writeXMLElementForm('loAtomicEP',(/'atomType     ','spin         ','branch       ',&
                                                              'branchLowest ','branchHighest','value        '/),&
                                               attributes,reshape((/10,4,6,12,13,5,6,1,3,8,8,16/),(/6,2/)))
                      WRITE(6,'(a6,i3,i2,a1,a12,f6.2,a3,f6.2,a13,f8.4)') '  Atom',n,nqn_lo(ilo),ch(l),' branch from',&
                                                                         e_lo,' to',e_up,' htr. ; e_l =',e
                   ENDIF


                ELSE IF(ABS(enpara%ello0(ilo,n,jsp) - nqn_lo(ilo)).LE. 1E-05 .AND. nqn_lo(ilo) .LT. 0  ) THEN

                   lo_done(ilo,n,jsp) = .TRUE.
                   ! search for branches
                   node = ABS(nqn_lo(ilo)) - (l+1)
                   e = 0.0 ! The initial value of e is arbitrary.
                   large_e_step = 5.0 ! 5.0 Htr steps for coarse energy searches

                   ! determine upper band edge
                   ! Step 1: Coarse search for the band edge
                   CALL radsra(e,l,vr(:,0,n,jsp),atoms%rmsh(1,n),&
                        atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
                   DO WHILE ( nodeu > node )
                      e = e - large_e_step
                      CALL radsra(e,l,vr(:,0,n,jsp),atoms%rmsh(1,n),&
                           atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
                   END DO
                   DO WHILE ( nodeu <= node )
                      e = e + large_e_step
                      CALL radsra(e,l,vr(:,0,n,jsp),atoms%rmsh(1,n),&
                           atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
                   END DO
                   e_up_temp = e
                   e_lo_temp = e - large_e_step
                   ! Step 2: Fine band edge determination by bisection search
                   DO WHILE ((e_up_temp - e_lo_temp) > 1e-2)
                      e = (e_up_temp + e_lo_temp) / 2.0
                      CALL radsra(e,l,vr(:,0,n,jsp),atoms%rmsh(1,n),&
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
                         CALL radsra(e,l,vr(:,0,n,jsp),atoms%rmsh(1,n),&
                              atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
                      ENDDO
                      e_up_temp = e + large_e_step
                      e_lo_temp = e
                      ! Step 2: Fine band edge determination by bisection search
                      DO WHILE ((e_up_temp - e_lo_temp) > 1e-2)
                         e = (e_up_temp + e_lo_temp) / 2.0
                         CALL radsra(e,l,vr(:,0,n,jsp),atoms%rmsh(1,n),&
                              atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
                         IF (nodeu < node) THEN
                            e_lo_temp = e
                         ELSE
                            e_up_temp = e
                         END IF
                      END DO
                      e_lo = (e_up_temp + e_lo_temp) / 2.0
                   END IF

                   ! determine energy with u'_l(e,R_MT)/u_l(e,R_MT)=-(l+1)
                   ldmt= -99.0 !ldmt = logarithmic derivative @ MT boundary
                   lnd = -l-1
                   DO WHILE (ABS(ldmt-lnd).GE.1e-7) 
                      e = (e_up + e_lo) / 2.0
                      CALL radsra(e,l,vr(:,0,n,jsp),atoms%rmsh(1,n),&
                           atoms%dx(n),atoms%jri(n),c, us,dus,nodeu,f(:,1),f(:,2))
                      ldmt = dus/us
                      IF(ldmt.GT.lnd) THEN
                         e_lo = e
                      ELSE
                         e_up = e
                      END IF
                   END DO


                   IF (mpi%irank == 0) THEN
                      attributes = ''
                      WRITE(attributes(1),'(i0)') n
                      WRITE(attributes(2),'(i0)') jsp
                      WRITE(attributes(3),'(i0,a1)') ABS(nqn_lo(ilo)), ch(l)
                      WRITE(attributes(4),'(f16.10)') ldmt
                      WRITE(attributes(5),'(f16.10)') e
                      CALL writeXMLElementForm('heloAtomicEP',(/'atomType      ','spin          ','branch        ',&
                                                              'logDerivMT    ','value         '/),&
                                               attributes(1:5),reshape((/8,4,6,12,5+17,6,1,3,16,16/),(/5,2/)))
                      WRITE (6,'(a6,i3,i2,a1,a12,f6.2,a4,f6.2,a5)') '  Atom',n,ABS(nqn_lo(ilo)),ch(l),&
                                                                    ' branch, D = ',ldmt,' at ',e,' htr.'
                   ENDIF

                   ello(ilo,n,jsp) = e                
                ELSE
                   lo_done(ilo,n,jsp) = .FALSE.
                ENDIF
             ENDDO

             DEALLOCATE ( f,vrd )
          ENDDO ! n
!$OMP END PARALLEL DO
       ENDDO   ! jsp
    ELSE
       l_done  = .FALSE.
       lo_done = .FALSE.
    ENDIF ! obsolete%lepr == 0
    !
    IF ((obsolete%lepr.EQ.1).AND.(mpi%irank.EQ.0)) THEN
       WRITE ( 6,'(//,'' Reference energies for energy parameters'')')
       WRITE ( 6,'('' ----------------------------------------'')')
       WRITE (16,'(//,'' Reference energies for energy parameters'')')
       WRITE (16,'('' ----------------------------------------'')')
    ENDIF
    !
    spins: DO jsp = 1,input%jspins
       types_loop: DO n = 1,atoms%ntype 
          
          !
          !--->    determine energy parameters if lepr=1. the reference energy
          !--->    is the value of the l=0 potential at approximately rmt/4.
          !
          IF (obsolete%lepr.EQ.1) THEN
             j = atoms%jri(n) - (LOG(4.0)/atoms%dx(n)+1.51)
             rj = atoms%rmt(n)*EXP(atoms%dx(n)* (j-atoms%jri(n)))
             vbar = vr(j,0,n,jsp)/rj
             IF (mpi%irank.EQ.0) THEN
                WRITE ( 6,'('' spin'',i2,'', atom type'',i3,'' ='',f12.6,''   r='',f8.5)') jsp,n,vbar,rj
                WRITE (16,'('' spin'',i2,'', atom type'',i3,'' ='',f12.6,''   r='',f8.5)') jsp,n,vbar,rj
             ENDIF
          ELSE
             vbar = 0.0
          END IF
          DO l = 0,atoms%lmax(n)
             IF ( .NOT.l_done(l,n,jsp) ) THEN
                el(l,n,jsp) = vbar + enpara%el0(l,n,jsp)
             END IF
          ENDDO

          IF (atoms%nlo(n).GE.1) THEN
             DO ilo = 1,atoms%nlo(n)
                IF ( .NOT. lo_done(ilo,n,jsp) ) THEN
                   ello(ilo,n,jsp) = vbar + enpara%ello0(ilo,n,jsp)
                   !+apw+lo
                   IF (atoms%l_dulo(ilo,n)) THEN
                      ello(ilo,n,jsp) = el(atoms%llo(ilo,n),n,jsp)
                   ENDIF
                   !-apw+lo
                END IF
             END DO
          ENDIF


       ENDDO types_loop

       IF (input%film) THEN
          !
          !--->    vacuum energy parameters: for lepr=1, relative to potential
          !--->    at vacuum-interstitial interface (better for electric field)
          !
          DO ivac = 1,vacuum%nvac
             vz0 = 0.0
             IF (obsolete%lepr.EQ.1) THEN
                vz0 = vz(1,ivac,jsp)
                IF (mpi%irank.EQ.0) THEN
                   WRITE ( 6,'('' spin'',i2,'', vacuum   '',i3,'' ='',f12.6)') jsp,ivac,vz0 
                   WRITE (16,'('' spin'',i2,'', vacuum   '',i3,'' ='',f12.6)') jsp,ivac,vz0
                ENDIF
             ENDIF
             evac(ivac,jsp) = enpara%evac0(ivac,jsp) + vz0
             IF (input%l_inpXML) THEN
                evac(ivac,jsp) = vz(vacuum%nmz,ivac,jsp) + enpara%evac0(ivac,jsp)
             END IF
             IF (mpi%irank.EQ.0) THEN
                attributes = ''
                WRITE(attributes(1),'(i0)') ivac
                WRITE(attributes(2),'(i0)') jsp
                WRITE(attributes(3),'(f16.10)') vz(1,ivac,jsp)
                WRITE(attributes(4),'(f16.10)') vz(vacuum%nmz,ivac,jsp)
                WRITE(attributes(5),'(f16.10)') evac(ivac,jsp)
                CALL writeXMLElementForm('vacuumEP',(/'vacuum','spin  ','vzIR  ','vzInf ','value '/),&
                                         attributes(1:5),reshape((/6+4,4,4,5,5+13,8,1,16,16,16/),(/5,2/)))
             END IF
          ENDDO
          IF (vacuum%nvac.EQ.1) THEN
             evac(2,jsp) = evac(1,jsp)
          END IF
       END IF
    ENDDO spins
    enpara_new=enpara
    enpara_new%el0=el
    enpara_new%ello0=ello
    enpara_new%evac0=evac
    IF (mpi%irank  == 0) CALL closeXMLElement('energyParameters')
    
  END SUBROUTINE lodpot
END MODULE m_lodpot
