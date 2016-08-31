!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_radsrdn
      use m_juDFT
      CONTAINS
      SUBROUTINE radsrdn(
     >                  e,l,vr,r0,h,jri,c,
     <                  udn,dudn,ddnn,nodedn,fn,fni,
     >                  f0,du0,n)
C*********************************************************************
C     Calculates the nth energy derivative of the scalar relativistic
C     wavefuction for energy e and angular momentum l.
C
C     The large and small components
C       of the nth derivative of f0 are returned in fn,
C       of the inhomogeneous part of the nth derivative of the Dirac
C       equation ( approx. n * (n-1)st derivative of f0 ) are returned
C       in fni.
C     For non-relativistic case, choose large c ( >> 137.0359895).
C
C     f0 (= primitive, 0th derivative) and
C     the radial derivative du0 of f0,
C     are required.
C
C     The nth derivative of the Dirac equation is solved in radsrdn1.
C
C                  C. Friedrich   Apr. 2005
C*********************************************************************
      USE m_intgr, ONLY : intgr0
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: l,n
      INTEGER, INTENT (IN) :: jri
      INTEGER, INTENT (OUT):: nodedn
      REAL,    INTENT (IN) :: c
      REAL,    INTENT (IN) :: du0,e,h,r0
      REAL,    INTENT (OUT):: ddnn,dudn,udn
C     ..
C     .. Array Arguments ..
      REAL,    INTENT (IN) :: f0(:,:),vr(:)
      REAL,    INTENT (OUT):: fn(:,:),fni(:,:)
C     ..
C     .. Local Scalars ..
      INTEGER i,j,k
      REAL scprod,ovlap(9,9),fac
C     ..
C     .. Local Arrays ..
      REAL f(size(vr),2,0:n),help(size(vr))
C
C
      IF(n<1)  CALL juDFT_error("n<1",calledby ="radsrdn")
      f(:,:,0) = f0(:,:)
C
      DO i=1,n
C       Calculate
C           scprod = < f0(:,:) | f(:,:,i) >
C                  = - 1/2  SUM  ( i over j )  < f(:,:,j) | f(:,:,i-j) >   
C                          j=1,i ----fac-----  -----------ovlap---------
C           (This follows from d^n/de^n <f0|f0> = 0.)
        scprod=0
        DO j=1,int(i/2)
          help(:) = f(:,1,j)*f(:,1,i-j) + f(:,2,j)*f(:,2,i-j)
          CALL intgr0(help,r0,h,jri,ovlap(j,i-j))
          ovlap(i-j,j)=ovlap(j,i-j)
          fac=1        !
          DO k=i-j+1,i !
            fac=fac*k  ! Calculate
          ENDDO        ! fac = ( i over j ) = i! / j! / (i-j)!
          DO k=2,j     !
            fac=fac/k  !
          ENDDO        !
          IF(j.eq.i-j) fac=fac/2
          scprod=scprod-fac*ovlap(j,i-j)
        ENDDO
C
C        scprod=0 ! Instead of the nth derivative a linear combination
                  ! of derivatives of order <=n is calculated.
                  ! It can be shown, that the udot-coefficient in the LO
                  ! construction in setabc1lo is always zero with this choice.
                  ! This reduces the error arising from the inconsistency
                  ! between an "exact" skalar-relativistic udot and a
                  ! non-relativistic Hamiltonian. Uncomment this, if you
                  ! have set RELATIVISTIC_CORRECTIONS in radsrd. But note,
                  ! that then relativistic corrections are not fully treated
                  ! in the higher derivatives calculated here.
        CALL getfni(e,l,vr,r0,h,jri,c,i,f,fni)
        CALL radsrdn1(
     >             e,l,vr,r0,h,jri,c,
     <             udn,dudn,ddnn,nodedn,f(:,1,i),f(:,2,i),
     >             fni(:,1),fni(:,2),f0(:,1),f0(:,2),du0,scprod)
      ENDDO
      fn(:,:) = f(:,:,n)
      help(:) = fni(:,1)
      fni(:,1)=-fni(:,2)
      fni(:,2)= help(:)*c
      END SUBROUTINE radsrdn

C---------------------------------------------------------------------

      SUBROUTINE radsrdn1(
     >                  e,l,vr,r0,h,jri,c,
     <                  ud,dud,ddn,nodes,pe,qe,
     >                  pi,qi,phom,qhom,dus,scprod)
C*********************************************************************
C     Solves the "inhomogeneous" Dirac equation for energy e and angular
C     momentum l with pi (qi) on one side instead of zero.
C
C     The large (small) component is returned in pe (qe).
C     For non-relativistic case, choose large c ( >> 137.0359895).
C
C     The functions
C       pi,qi,
C       phom,qhom (= homogeneous solution, 0th derivative,
C                    as calculated in radsra),
C     as well as
C       the radial derivative dus of phom,
C       and scprod = < p(q)e | p(q)hom >
C     are required.
C
C     Modified from radsrd.    C. Friedrich   Apr. 2005
C*********************************************************************
C
      USE m_intgr, ONLY : intgr0
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: l
      INTEGER, INTENT (IN) :: jri
      INTEGER, INTENT (OUT):: nodes
      REAL,    INTENT (IN) :: c
      REAL,    INTENT (IN) :: dus,e,h,r0,scprod
      REAL,    INTENT (OUT):: ddn,dud,ud
C     ..
C     .. Array Arguments ..
      REAL,    INTENT (IN) :: pi(:),qi(:),vr(:)
      REAL,    INTENT (OUT):: pe(:),qe(:)
      REAL,    INTENT (IN) :: phom(:),qhom(:)
C     ..
C     .. Local Scalars ..
      REAL dr,drh,eps,erp,erq,fl1,p0,p1,p1p,q0,q1,q1p,r,rh,rh2,rm,rve,
     +     sk1,sk2,sk3,sk4,sl1,sl2,sl3,sl4,t,t1,t2,yn,zn,cin,cin2
      INTEGER i,it
C     ..
C     .. Local Arrays ..
      REAL pp(size(pi)),qp(size(pi))
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,exp,sign,sqrt
C     ..
C     .. Data statements ..
      DATA eps/1.e-06/
C     ..
      cin = 1.0/c
      cin2 = cin*cin
c
      IF (jri>size(pi))  CALL juDFT_error("dimension too small",  
     +     calledby ="radsrdn")
c--->    set up initial conditions
      fl1 = l* (l+1)
      pe(1) = 0
      qe(1) = 0
      pp(1) = r0*pi(1)
      qp(1) = r0*qi(1)
c--->    use 4th order runge-kutta to get first few mesh points
      dr = exp(h)
      drh = sqrt(dr)
      r = r0
      DO i = 1,5
         rh2 = drh*r
         rh = dr*r
         sk1 = h*pp(i)
         sl1 = h*qp(i)
         rve = 0.5* (vr(i)+vr(i+1)) - rh2*e
         rm = 2.*rh2 - cin2*rve
         yn = pe(i) + 0.5*sk1
         zn = qe(i) + 0.5*sl1
         t1 = 0.5*rh2*(pi(i)+pi(i+1))
         t2 = 0.5*rh2*(qi(i)+qi(i+1))
         sk2 = h* (rm*zn + yn              + t1)
         sl2 = h* (-zn   + (fl1/rm+rve)*yn + t2)
         yn = pe(i) + 0.5*sk2
         zn = qe(i) + 0.5*sl2
         sk3 = h* (rm*zn + yn              + t1)
         sl3 = h* (-zn   + (fl1/rm+rve)*yn + t2)
         rve = vr(i+1) - rh*e
         rm = 2.*rh - cin2*rve
         yn = pe(i) + sk3
         zn = qe(i) + sl3
         t1 = rh*pi(i+1)
         t2 = rh*qi(i+1)
         sk4 = h* (rm*zn + yn              + t1 )
         sl4 = h* (-zn   + (fl1/rm+rve)*yn + t2 )
         pe(i+1) = pe(i) + (sk1+2.*sk2+2.*sk3+sk4)/6.
         qe(i+1) = qe(i) + (sl1+2.*sl2+2.*sl3+sl4)/6.
         pp(i+1) = rm*qe(i+1) + pe(i+1)              + t1
         qp(i+1) = -qe(i+1)   + (fl1/rm+rve)*pe(i+1) + t2
         r = rh
      ENDDO
      nodes = 0
c--->    adams-bashforth-moulton predictor-corrector
      predictor: DO i = 6,jri - 1
         r = r*dr
c--->    predictor
         p0 = pe(i) + h* (4277.*pp(i)-7923.*pp(i-1)+9982.*pp(i-2)-
     +        7298.*pp(i-3)+2877.*pp(i-4)-475.*pp(i-5))/1440.
         q0 = qe(i) + h* (4277.*qp(i)-7923.*qp(i-1)+9982.*qp(i-2)-
     +        7298.*qp(i-3)+2877.*qp(i-4)-475.*qp(i-5))/1440.
c--->    evaluate derivatives at next point
         rve = vr(i+1) - r*e
         rm = 2.*r - cin2*rve
         t1 = r*pi(i+1)
         t2 = r*qi(i+1)
         p1p = rm*q0 + p0              + t1
         q1p = -q0   + (fl1/rm+rve)*p0 + t2
c--->    corrector
         corrector: DO it = 1,5
            p1 = pe(i) + h* (475.*p1p+1427.*pp(i)-798.*pp(i-1)+
     +           482.*pp(i-2)-173.*pp(i-3)+27.*pp(i-4))/1440.
            q1 = qe(i) + h* (475.*q1p+1427.*qp(i)-798.*qp(i-1)+
     +           482.*qp(i-2)-173.*qp(i-3)+27.*qp(i-4))/1440.
c--->    final evaluation
            p1p = rm*q1 + p1              + t1
            q1p = -q1   + (fl1/rm+rve)*p1 + t2
c--->    test quality of corrector and iterate if necessary
            erp = abs(p1-p0)/ (abs(p1)+abs(h*p1p))
            erq = abs(q1-q0)/ (abs(q1)+abs(h*p1p))
            IF (erp.LT.eps .AND. erq.LT.eps) EXIT corrector
            p0 = p1
            q0 = q1
         ENDDO corrector
         IF (it > 5)  CALL juDFT_error("Not converged.",calledby
     +        ="radsrdn")
c--->    store values
         pe(i+1) = p1
         qe(i+1) = q1
         pp(i+1) = p1p
         qp(i+1) = q1p
         nodes = nodes + 0.501*abs(sign(1.0,pe(i+1))-sign(1.0,pe(i)))
      ENDDO predictor
c--->    ensure < p(q)hom | p(q)e > = scprod
      DO i = 1,jri
         qe(i) = cin*qe(i)
      ENDDO
      DO i = 1,jri
         qp(i) = phom(i)*pe(i) + qhom(i)*qe(i)
      ENDDO
      CALL intgr0(qp,r0,h,jri,t)
      dud = (pp(jri)-pe(jri))/ (r*r)
      DO i = 1,jri
         pe(i) = pe(i) + (scprod-t)*phom(i)
         qe(i) = qe(i) + (scprod-t)*qhom(i)
      ENDDO
      ud = pe(jri)/r
      dud = dud + (scprod-t)*dus
      DO i = 1,jri
         qp(i) = pe(i)*pe(i) + qe(i)*qe(i)
      ENDDO
      CALL intgr0(qp,r0,h,jri,ddn)
      RETURN
      END SUBROUTINE radsrdn1
      
C---------------------------------------------------------------------

      SUBROUTINE getfni(e,l,vr,r0,h,jri,c,n,f,fni)
C*********************************************************************
C     Calculates the inhomogeneous part of the nth derivative of the
C     Dirac equation with energy e and angular momentum l.
C
C     It is
C
C     fni(:,1) = n*f(:,2,n-1)/c (Note that f(:,2)=Q/c)
C
C                                     l(l+1)
C     fni(:,2) =-n*f(:,1,n-1)* ( 1 + -------- )
C                                    (2Mrc)^2
C
C                  n      i           n! l(l+1)
C               + SUM (-1)   -----------------------------
C                 i=2        (n-i)! (2*M)^(n+1) r^2 c^(2n)
C
C     with M = 1 + (E-V)/(2*c^2).
C
C                  C. Friedrich   Apr. 2005
C*********************************************************************
C

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: n,jri,l
      REAL,    INTENT (IN) :: c,e,h,r0
C     ..
C     .. Array Arguments ..
      REAL,    INTENT (IN) :: f(:,:,0:),vr(:)
      REAL,    INTENT (OUT):: fni(:,:)
C     ..
C     .. Local Scalars ..
      INTEGER i
      REAL r
C     ..
C     .. Local Arrays ..
      REAL hlp(jri),fac(0:n),rmsh(jri)
C
      IF(n.lt.1)  CALL juDFT_error("getfni: n<1!",calledby="radsrdn")
      r=exp(h)
      rmsh(1)=r0
      DO i=2,jri
        rmsh(i)=rmsh(i-1)*r
      ENDDO
      fac(0)=1
      DO i=1,n
        fac(i)=fac(i-1)*i
      ENDDO
      hlp(:)      = 2*rmsh(:) + (e*rmsh(:)-vr(:jri))/c**2
      fni(:jri,1) = n*f(:jri,2,n-1) /c
      fni(:jri,2) =-n*f(:jri,1,n-1)
     &              * ( 1 + l*(l+1)* (hlp(:)*c)**(-2) )
      hlp(:)      = hlp(:)/rmsh(:)
      DO i=2,n
        r = fac(n)/fac(n-i) * (-1)**i * l*(l+1) / c**(2*i)
        fni(:jri,2) = fni(:jri,2) +
     &                f(:jri,1,n-i) * r / hlp(:)**(i+1) / rmsh(:)**2
      ENDDO
      END SUBROUTINE getfni
      
      END MODULE m_radsrdn
      
