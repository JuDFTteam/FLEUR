!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_jij

  USE m_types
  USE m_types_forcetheo
  USE m_judft
  TYPE,EXTENDS(t_forcetheo) :: t_forcetheo_jij
     INTEGER :: loopindex,no_loops
     INTEGER,ALLOCATABLE :: q_index(:),iatom(:),jatom(:)
     LOGICAL,ALLOCATABLE :: phase2(:)

     REAL,ALLOCATABLE:: qvec(:,:)
     REAL            :: thetaj
     REAL,ALLOCATABLE:: evsum(:)
   CONTAINS
     PROCEDURE :: start   =>jij_start
     PROCEDURE :: next_job=>jij_next_job 
     PROCEDURE :: eval    =>jij_eval
     PROCEDURE :: postprocess => jij_postprocess
     PROCEDURE :: init   => jij_init !not overloaded
     PROCEDURE :: dist   => jij_dist !not overloaded
  END TYPE t_forcetheo_jij

CONTAINS

  


  SUBROUTINE jij_init(this,qvec,thetaj,atoms)
    USE m_types_setup
    USE m_constants
    IMPLICIT NONE
    CLASS(t_forcetheo_jij),INTENT(INOUT):: this
    REAL,INTENT(in)                     :: qvec(:,:),thetaj
    TYPE(t_atoms),INTENT(IN)            :: atoms
  
    INTEGER:: n,na,ni,nj,j
    REAL,PARAMETER:: eps=1E-5

    this%qvec=qvec
    this%thetaj=thetaj

    !Max no of loops...
    n=atoms%nat**2*SIZE(this%qvec,2)+1
    ALLOCATE(this%q_index(n),this%iatom(n),this%jatom(n),this%phase2(n))
    

    !now construct the loops
    this%no_loops=0
    DO n=1,SIZE(this%qvec,2)
       DO ni=1,atoms%ntype
          IF (ABS(atoms%bmu(ni))<eps) CYCLE !no magnetic atom
          DO nj=ni,atoms%ntype
             IF (ABS(atoms%bmu(nj))<eps) CYCLE !no magnetic atom
             DO j=1,MERGE(1,2,ni==nj) !phase factor
                !new config found
                this%no_loops=this%no_loops+1
                this%q_index(this%no_loops)=n
                this%iatom(this%no_loops)=ni
                this%jatom(this%no_loops)=nj
                this%phase2(this%no_loops)=(j==2)
             ENDDO
          END DO
       END DO
    END DO

    ALLOCATE(this%evsum(this%no_loops))
    this%evsum=0
  END SUBROUTINE jij_init

  
  SUBROUTINE jij_dist(this,mpi)
    USE m_types_mpi
    IMPLICIT NONE
    CLASS(t_forcetheo_jij),INTENT(INOUT):: this
    TYPE(t_mpi),INTENT(in):: mpi
    
    INTEGER:: i,q,ierr
#ifdef CPP_MPI    
    INCLUDE 'mpif.h'
    call judft_error("jij has to be parallelized")
#endif    
  END SUBROUTINE jij_dist

  

  SUBROUTINE jij_start(this,potden,l_io)
    USE m_types_potden
    IMPLICIT NONE
    CLASS(t_forcetheo_jij),INTENT(INOUT):: this
    TYPE(t_potden) ,INTENT(INOUT)       :: potden
    LOGICAL,INTENT(IN)                  :: l_io
    this%loopindex=0
    CALL this%t_forcetheo%start(potden,l_io) !call routine of basis type
  END SUBROUTINE  jij_start

  LOGICAL FUNCTION jij_next_job(this,lastiter,atoms,noco)
    USE m_types_setup
    USE m_xmlOutput
    USE m_constants
    IMPLICIT NONE
    CLASS(t_forcetheo_jij),INTENT(INOUT):: this
    LOGICAL,INTENT(IN)                  :: lastiter
    TYPE(t_atoms),INTENT(IN)            :: atoms
    !Stuff that might be modified...
    TYPE(t_noco),INTENT(INOUT) :: noco

    !locals
    INTEGER:: n

    IF (.NOT.lastiter) THEN
       jij_next_job=this%t_forcetheo%next_job(lastiter,atoms,noco)
       RETURN
    ENDIF
    
    !OK, now we start the JIJ-loop

    this%loopindex=this%loopindex+1
    jij_next_job=(this%loopindex<=this%no_loops) !still loops to do...
    IF (.NOT.jij_next_job) RETURN

    ! Now set the noco-variable accordingly...

    noco%qss=this%qvec(:,this%q_index(this%loopindex))

!c     Determines the cone angles and phase shifts of the spin 
!c     vectors on magnetic atoms for the calculation of the
!c     interaction constants Jij from the Heisenberg model
!c                                   M. Lezaic 04

    noco%alph=0.0;noco%beta=0.0
    noco%beta(this%iatom(this%loopindex))=this%thetaj
    IF (this%phase2(this%loopindex)) noco%alph(this%iatom(this%loopindex))=pi_const*0.5
    noco%beta(this%jatom(this%loopindex))=this%thetaj

    !rotate according to q-vector
    DO n = 1,atoms%ntype
       noco%alph(n) = noco%alph(n) + tpi_const*DOT_PRODUCT(noco%qss,atoms%taual(:,SUM(atoms%neq(:n-1))+1))
    ENDDO
    
    IF (.NOT.this%l_io) RETURN
    IF (this%loopindex.NE.1) CALL closeXMLElement('Forcetheorem_Loop_JIJ')
    CALL openXMLElementPoly('Forcetheorem_Loop_JIJ',(/'Loop index:'/),(/this%loopindex/))
  END FUNCTION jij_next_job

  SUBROUTINE jij_postprocess(this)
    USE m_xmlOutput
    IMPLICIT NONE
    CLASS(t_forcetheo_jij),INTENT(INOUT):: this

    !Locals
    INTEGER:: n
    CHARACTER(LEN=18):: attributes(6)

    IF (this%loopindex==0) RETURN
    
    IF (.NOT.this%l_io) RETURN
  
    !Now output the results
    call closeXMLElement('Forcetheorem_Loop_JIJ')
    CALL openXMLElementPoly('Forcetheorem_JIJ',(/'Configs'/),(/this%no_loops/))
    DO n=1,this%no_loops
       WRITE(attributes(1),'(i5)') n
       WRITE(attributes(2),'(3(f5.3,1x))') this%qvec(:,this%q_index(n))
       WRITE(attributes(3),'(i0)') this%iatom(n)
       WRITE(attributes(4),'(i0)') this%jatom(n)
       WRITE(attributes(5),'(l1)') this%phase2(n)
       WRITE(attributes(6),'(f15.8)') this%evsum(n)
       CALL writeXMLElementForm('Config',(/'n     ','q     ','iatom ','jatom ','phase ','ev-sum'/),attributes,&
               RESHAPE((/1,1,5,5,5,6,4,18,4,4,2,15/),(/6,2/)))
    ENDDO
    CALL closeXMLElement('Forcetheorem_JIJ')
  END SUBROUTINE jij_postprocess


  FUNCTION jij_eval(this,eig_id,atoms,kpts,sym,&
       cell,noco, input,mpi, oneD,enpara,v,results)RESULT(skip)
     USE m_types
     USE m_ssomat
    IMPLICIT NONE
    CLASS(t_forcetheo_jij),INTENT(INOUT):: this
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
    skip=.FALSE.
    IF (this%loopindex==0) RETURN
  
    this%evsum(this%loopindex)=results%seigv
    skip=.TRUE.
  END FUNCTION  jij_eval

  

  SUBROUTINE priv_analyse_data()
!-------------------------------------------------------------------
!     calculates the
!     coupling constants J for all the couples of magnetic
!     atoms, for a given number of neighbouring shells
!                                   M. Lezaic 04
!-------------------------------------------------------------------

    USE m_constants, ONLY : pimach
    PRINT *,"jcoef2 has still to be reimplemented"
#ifdef CPP_NEVER
      USE m_nshell
#include"cpp_double.h"
      IMPLICIT NONE

c     .. Scalar arguments ..

      INTEGER, INTENT (IN)   :: ntypd,nmagn,nqpt,mtypes
      INTEGER, INTENT (IN)   :: natd,nop
      INTEGER, INTENT (INOUT)   :: nsh
      REAL,    INTENT (IN)   :: thetaJ
      LOGICAL, INTENT (IN)   :: invs,film

c     .. Array arguments ..

      INTEGER,  INTENT (IN) :: mrot(3,3,nop)
      INTEGER,  INTENT (IN) :: nmagtype(ntypd),magtype(ntypd),neq(ntypd)
      REAL,    INTENT (IN)  :: taual(3,natd)
      REAL,    INTENT (IN)  :: amat(3,3)

c     .. Temporary..
      INTEGER :: nshort !The number of non-zero interactions for a calculation 
c                        with the least square method
      REAL   Tc 
      REAL, ALLOCATABLE :: Cmat(:,:),DelE(:),work(:)
      INTEGER  lwork

c     .. Local scalars ..

      INTEGER n,nn,nnn,mu,nu,iatom,nqvect,atsh,phi,i,ii,wrJ
      INTEGER qcount,info,imt,limit,remt,sneq,nqcomp
      REAL    qn,sqsin,alph,beta,tpi,scp,J,rseig
      REAL    enmin,enmax,zcoord,deltaz
!     REAL    wei !(old version)
      INTEGER, PARAMETER   :: nmax=10,dims=(2*nmax+1)**3,shmax=192
      REAL, PARAMETER :: tol=0.000001 
 
c     ..Local Arrays..

      INTEGER w(nop,nqpt-1),itype(nmagn)
      INTEGER nat(dims),ierr(3)
      REAL seigv(nmagn,nmagn,nqpt-1,2),q(3,nop,nqpt-1),qss(3),Jw(shmax)
      REAL tauC1(3),tauC2(3),seigv0(nmagn,nmagn,2)
      REAL ReJq(nmagn,nmagn,nqpt-1),ImJq(nmagn,nmagn,nqpt-1)
      REAL M(nmagn),qssmin(3),qssmax(3),t(3),R(3,shmax,dims),lenR(dims)
      REAL Dabsq(3),divi(3)
      INTEGER IDabsq(3)

c     .. Intrinsic Functions ..

      INTRINSIC cos,sin

c     .. External Subroutines ..
      
      EXTERNAL CPP_LAPACK_dgels
#ifdef CPP_MPI
      INCLUDE 'mpif.h'
#endif

c-------------------------------------------------------------------
       OPEN(116,file='qptsinfo',status='unknown')
       OPEN(115,file='jconst',status='unknown')
        w=0
        nqvect=0
        nshort=0 
        ReJq=0.0
        ImJq=0.0
        sqsin=(sin(thetaJ))**2
        tpi = 2.0 * pimach()
        limit=nmagn-1
        IF (nmagn.gt.mtypes) limit=mtypes

            IF(nsh.LT.0)THEN
c...   Setup for a calculation using least squares method       
            WRITE(6,*) 'Jij calculated with the least squares method'
             nsh=-nsh
             nshort=nsh 
            ENDIF

      DO n=1,nqpt
         qcount=n-1
       mu=1
       DO imt=1,mtypes
        DO nu=mu,nmagn
         DO phi=1,2
         READ(114,*)
         READ(114,5000) itype(mu),alph,beta,M(mu) 
         READ(114,5000) itype(nu),alph,beta,M(nu) 
!        READ(114,5001) qss(1),qss(2),qss(3),wei, rseig !(old version)
         READ(114,5001) qss(1),qss(2),qss(3),rseig
 5000    FORMAT(3x,i4,4x,f14.10,1x,f14.10,4x,f8.5)
!5001    FORMAT(4(f14.10,1x),f20.10) !(old version)
 5001    FORMAT(3(f14.10,1x),f20.10)

c...     Aquire information on the maximal and minimal calculated energy
          IF((n*mu*nu*phi).EQ.1)THEN
             qssmin=qss
             enmin=rseig
             qssmax=qss
             enmax=rseig
          ELSE 
             IF(enmin.GE.rseig)THEN
                enmin=rseig
                qssmin=qss
                ENDIF
             IF(enmax.LE.rseig)THEN
                enmax=rseig
                qssmax=qss
                ENDIF
           ENDIF

               IF(n.EQ.1)THEN
               seigv0(mu,nu,phi)=rseig ! reference:energy of the collinear state 
               ELSE
               seigv(mu,nu,qcount,phi)=rseig ! energies of spin-spirals
               ENDIF 
         qn = ( qss(1)**2 + qss(2)**2 + qss(3)**2 )
         IF((mu.EQ.nu).OR.(invs.AND.(qn.GT.tol)))
     &   GOTO 33
         ENDDO !phi
   33   CONTINUE
        ENDDO !nu
        mu=mu+nmagtype(imt)
       ENDDO !imt
         IF(n.EQ.1)THEN
           IF(qn.GT.tol) THEN
           WRITE(*,*) 'jcoff2: The first energy in jenerg file should correspond 
     &to qss=0!'
#ifdef CPP_MPI
            CALL MPI_ABORT(MPI_COMM_WORLD,1,ierr)
#endif
           STOP
           ENDIF
         ELSE
         WRITE(116,*) qcount
c...
c...   Apply all the symmetry operations to the q-vectors from the irreducible wedge
c...
      DO nn=1,nop
          q(1,nn,qcount)=qss(1)*mrot(1,1,nn) + qss(2)*mrot(2,1,nn) +
     +                   qss(3)*mrot(3,1,nn)
          q(2,nn,qcount)=qss(1)*mrot(1,2,nn) + qss(2)*mrot(2,2,nn) +
     +                   qss(3)*mrot(3,2,nn)
          q(3,nn,qcount)=qss(1)*mrot(1,3,nn) + qss(2)*mrot(2,3,nn) +
     +                   qss(3)*mrot(3,3,nn)
             w(nn,qcount)=1
c...
c...    Eliminate the q-vectors which are repeating and for each q remove the -q
c... 
           DO nnn=1,nn-1
            DO ii=1,2
            Dabsq(:)=ABS(q(:,nn,qcount)+((-1)**ii)*q(:,nnn,qcount))
            IDabsq(:)=NINT(Dabsq(:))
            divi(:)=ABS(Dabsq(:)/FLOAT(IDabsq(:))-1.0) 
              IF(((Dabsq(1).LT.tol).OR.(divi(1).LT.tol)).AND.
     &           ((Dabsq(2).LT.tol).OR.(divi(2).LT.tol)).AND.
     &           ((Dabsq(3).LT.tol).OR.(divi(3).LT.tol)))THEN
                 w(nn,qcount)=0
                 GOTO 44
              ENDIF
             ENDDO ! ii 
           ENDDO !nnn
         nqvect=nqvect+1
        WRITE(116,5502) mrot(1,1,nn),mrot(1,2,nn),mrot(1,3,nn) 
        WRITE(116,5502) mrot(2,1,nn),mrot(2,2,nn),mrot(2,3,nn) 
        WRITE(116,5502) mrot(3,1,nn),mrot(3,2,nn),mrot(3,3,nn) 
        WRITE(116,5002) nqvect,q(1,nn,qcount),q(2,nn,qcount),
     &                         q(3,nn,qcount)
5002    FORMAT(i6,3(f14.10,1x))
5502    FORMAT(3(i3))
  44       CONTINUE 
      ENDDO !nn
c...      Now calculate Jq=Re(Jq)+i*Im(Jq)
              mu=1
           DO imt=1,mtypes
              ReJq(mu,mu,qcount)=-2.0*(seigv(mu,mu,qcount,1)
     &                -seigv0(mu,mu,1))/(M(mu)*M(mu)*sqsin)
              ImJq(mu,mu,qcount)=0.0
               DO remt=mu+1,mu+nmagtype(imt)-1
               ReJq(remt,remt,qcount)=ReJq(mu,mu,qcount)
               ImJq(remt,remt,qcount)=0.0
               ENDDO!remt
           mu=mu+nmagtype(imt)
           ENDDO !imt
              mu=1
           DO imt=1,limit
            DO nu=mu+1,nmagn
              ReJq(mu,nu,qcount)=((seigv0(mu,nu,2)-
     &        seigv(mu,nu,qcount,1))/(M(mu)*M(nu)*sqsin))
     &        -(0.5*M(mu)*ReJq(mu,mu,qcount)/M(nu))-
     &        (0.5*M(nu)*ReJq(nu,nu,qcount)/M(mu))
              IF(invs)THEN
               ImJq(mu,nu,qcount)=0.0
              ELSE  
               ImJq(mu,nu,qcount)=((seigv(mu,nu,qcount,2)
     &         -seigv(mu,nu,qcount,1))/
     &         (M(mu)*M(nu)*sqsin))-ReJq(mu,nu,qcount)
              ENDIF !invs
            ENDDO !nu
             mu=mu+nmagtype(imt)
           ENDDO !mu

         ENDIF !if(n.eq.1)   
      ENDDO !n   
      WRITE(6,5006)enmax,qssmax
      WRITE(6,5007)enmin,qssmin
 5006 FORMAT('Maximal calculated energy = ',f20.10,'htr',
     & ' for the spin-spiral vector qss = ',3(f14.10,1x))
 5007 FORMAT('Minimal calculated energy = ',f20.10,'htr',
     & ' for the spin-spiral vector qss = ',3(f14.10,1x))
      
      CLOSE(116)
      OPEN(117,file='shells',status='unknown')
      OPEN(118,file='MCinp',status='unknown')
         mu=1

      DO imt=1,mtypes
        DO nu=mu,nmagn
          WRITE(115,5004) itype(mu),itype(nu)
          WRITE(117,5004) itype(mu),itype(nu)
 5004     FORMAT('#',i4,i4)

          sneq=0
          do ii=1,itype(mu)-1
          sneq=sneq+neq(ii)
          enddo
          if(itype(mu).le.sneq) then
          itype(mu)=sneq+1
          endif

          do ii=itype(mu),itype(nu)-1
          sneq=sneq+neq(ii)
          enddo
          if(itype(nu).le.sneq) then
          itype(nu)=sneq+1
          endif

          t(:)=taual(:,itype(mu))-taual(:,itype(nu))
 
        deltaz=taual(3,itype(mu))-taual(3,itype(nu))   !Added for film 07/10 S. Schroeder 
!        zcoord=taual(3,itype(nu))*amat(3,3)-taual(3,itype(1))*amat(3,3)
        zcoord=taual(3,itype(nu))*amat(3,3)
          tauC1(:)=taual(1,itype(mu))*amat(:,1)
     &            +taual(2,itype(mu))*amat(:,2)
     &            +taual(3,itype(mu))*amat(:,3)
c...  Generate the shells of neighbouring atoms
      CALL nshell(
     >                  amat,t,nsh,dims,nmax,shmax,film,zcoord,
     <                  nat,R,lenR,nop,mrot,deltaz)

c ...
c ... For the case of calculation with the least squares method
c ... for one magnetic atom per unit cell
       IF (nshort.GT.0) THEN
       qcount=nqpt-1
       lwork=2*nshort 
       ALLOCATE (Cmat(qcount,nshort),DelE(qcount),work(lwork))
          Cmat=0.0
       IF (nshort.GE.nqpt)THEN 
        WRITE(*,*) ' Please supply the data for', nshort,
     & 'q-points different from zero' 
        STOP   
        ENDIF

          DO n=1,qcount
           DO nn=1,nshort
            DO atsh=1,nat(nn)
            scp=(q(1,1,n)*R(1,atsh,nn)    
     &          +q(2,1,n)*R(2,atsh,nn)
     &          +q(3,1,n)*R(3,atsh,nn))*tpi
            Cmat(n,nn)=Cmat(n,nn)-1.0+cos(scp)
            ENDDO
           ENDDO
          DelE(n)=ReJq(1,1,n)*2000 ! multiply by 2000 to get [mRy/muB**2]
          ENDDO

      IF (nu.gt.mu) STOP 'For more than one magnetic atom in the unit 
     & cell please set nshort=0'
       CALL dgels('N',qcount,nshort,1,Cmat,qcount,DelE,qcount,
     & work,lwork,info)

c      The routine dgels returns the solution, J(n), in the array DelE(n)  
       Tc=0.0
      DO n=1,nshort
      Tc=Tc+nat(n)*DelE(n) !Mean-field Tc=1/3*(Sum_i(J_0,i))
      WRITE(115,5005) n,lenR(n),DelE(n) ! J in units [mRy/muB**2]

c Now prepare and write the input for the Monte Carlo calculation:
            DO atsh=1,nat(n)
          tauC2(:)=(taual(1,itype(mu))-R(1,atsh,n))*amat(:,1)
     &            +(taual(2,itype(mu))-R(2,atsh,n))*amat(:,2)
     &            +(taual(3,itype(mu))-R(3,atsh,n))*amat(:,3)
         
       WRITE(118,5008) itype(mu),itype(nu),tauC1(1),tauC1(2),tauC1(3),
     &                 tauC2(1),tauC2(2),tauC2(3),DelE(n)
           ENDDO ! atsh

      ENDDO ! n
      DEALLOCATE (Cmat,DelE,work)
      ELSE
c... Perform the back-Fourier transform
          DO nnn=1,nsh
             wrJ=0
          DO atsh=1,nat(nnn)
          IF(atsh.gt.shmax) STOP 'jcoff2:increase shmax!' 
          J=0.0  
          DO n=1,nqpt-1
           DO nn=1,nop
            IF(w(nn,n).EQ.1)THEN
      scp=(q(1,nn,n)*R(1,atsh,nnn)
     &    +q(2,nn,n)*R(2,atsh,nnn)
     &    +q(3,nn,n)*R(3,atsh,nnn))*tpi
      J=J+(((cos(scp))*ReJq(mu,nu,n)
     &     +(sin(scp))*ImJq(mu,nu,n)))
            ENDIF
           ENDDO !nn
          ENDDO !n (qpts)
      J=(J/float(nqvect))*2000.0 ! J in units [mRy/muB**2]
      DO i=1,wrJ !A check for non-equivalent sub-shells
       IF(ABS(J-Jw(i)).LE.(tol))GOTO 55
      ENDDO
       WRITE(115,5005) nnn,lenR(nnn),J
       wrJ=wrJ+1
       Jw(wrJ)=J
   55 CONTINUE
 5005     FORMAT(i4,1x,2(f14.10,1x))

c Now prepare and write the input for the Monte Carlo calculation:
          tauC2(:)=(taual(1,itype(mu))-R(1,atsh,nnn))*amat(:,1)
     &            +(taual(2,itype(mu))-R(2,atsh,nnn))*amat(:,2)
     &            +(taual(3,itype(mu))-R(3,atsh,nnn))*amat(:,3)
         
       WRITE(118,5008) itype(mu),itype(nu),tauC1(1),tauC1(2),tauC1(3),
     &                 tauC2(1),tauC2(2),tauC2(3),J
         ENDDO !atsh (atoms in each shell)
c...  In case of only one magnetic atom per unit cell, calculate the mean-field Tc 
          IF(nmagn.EQ.1) Tc=Tc+nat(nnn)*J 
         ENDDO !nnn (shells)
      ENDIF !nshort
        ENDDO !nu
       mu=mu+nmagtype(imt)
      ENDDO !imt
          IF(nmagn.EQ.1) THEN
          Tc=157.889*M(1)*M(1)*Tc/3.0
          WRITE(115,*) '# Tc(mean field)= ',Tc
          ENDIF 
 5008     FORMAT(i4,i4,7(1x,f14.10))
      CLOSE(117)
      CLOSE(118)
      CLOSE(115)
#endif
    END SUBROUTINE priv_analyse_data
    

      END MODULE m_types_jij
