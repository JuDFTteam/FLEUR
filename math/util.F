      MODULE m_util
      USE m_juDFT
c     error and warning codes for intgrf function
      INTEGER, PARAMETER :: NO_ERROR                  = 0
      INTEGER, PARAMETER :: NEGATIVE_EXPONENT_WARNING = 1
      INTEGER, PARAMETER :: NEGATIVE_EXPONENT_ERROR   = 2
c     return type of the pure intgrf function
      TYPE :: intgrf_out
        REAL    :: value    ! value of the integration
        INTEGER :: ierror   ! error code
      END TYPE intgrf_out

      INTERFACE derivative
        MODULE PROCEDURE derivative_t,derivative_nt
      END INTERFACE
      
      CONTAINS


c     Calculates Gaunt coefficients, i.e. the integrals of three spherical harmonics
c     integral ( conjg(Y(l1,m1)) * Y(l2,m2) * conjg(Y(l3,m3)) )
c     They are also the coefficients C(l1,l2,l3,m1,m2,m3) in
c     conjg(Y(l1,m1)) * Y(l2,m2) = sum(l3,m3) C(l1,l2,l3,m1,m2,m3) Y(l3,m3)
c     fac contains factorial up to maxfac, i.e. fac(i)= i! 
C     sfac contains square root of fac, i.e. sfac(i)= sqrt(i!) 

      FUNCTION gaunt(l1,l2,l3,m1,m2,m3,maxfac,fac,sfac)
      
      USE m_constants ,ONLY: pimach
      
      IMPLICIT NONE
      
      REAL                :: gaunt
      INTEGER, INTENT(IN) :: l1,l2,l3,m1,m2,m3,maxfac
      REAL   , INTENT(IN) :: fac(0:maxfac)
      REAL   , INTENT(IN) :: sfac(0:maxfac)

      gaunt = 0
      IF(m3.ne.m2-m1)   RETURN
      IF(abs(m1).gt.l1) RETURN
      IF(abs(m2).gt.l2) RETURN
      IF(abs(m3).gt.l3) RETURN
      IF(l3.lt.abs(l1-l2).or.l3.gt.l1+l2) RETURN
      gaunt = (-1)**(m1+m3) *
     &        sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)/pimach()/4)*
     &        wigner3j(l1,l2,l3,-m1,m2,-m3,maxfac,fac,sfac)*
     &        wigner3j(l1,l2,l3, 0,  0,  0,maxfac,fac,sfac)
      END FUNCTION gaunt
 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

c     Calculates the Wigner 3j symbols using Racah's formula

      FUNCTION wigner3j(l1,l2,l3,m1,m2,m3,maxfac,fac,sfac)

      IMPLICIT NONE

      REAL                :: wigner3j
      INTEGER, INTENT(IN) :: l1,l2,l3,m1,m2,m3,maxfac
      REAL   , INTENT(IN) :: fac(0:maxfac)
      REAL   , INTENT(IN) :: sfac(0:maxfac)
c     - local - 
      INTEGER             :: tmin,tmax,t,f1,f2,f3,f4,f5

      wigner3j = 0

c     The following IF clauses should be in the calling routine and commented here.
c      if(-m3.ne.m1+m2)  return
c      if(abs(m1).gt.l1) return
c      if(abs(m2).gt.l2) return
c      if(abs(m3).gt.l3) return
c      if(l3.lt.abs(l1-l2).or.l3.gt.l1+l2) return

      f1   = l3-l2+m1
      f2   = l3-l1-m2
      f3   = l1+l2-l3
      f4   = l1-m1
      f5   = l2+m2
      tmin = max(0,-f1,-f2) ! The arguments to fac (see below)
      tmax = min(f3,f4,f5)  ! must not be negative.

      ! The following line is only for testing and should be removed at a later time.      
      IF(tmax-tmin .ne. min(l1+m1,l1-m1,l2+m2,l2-m2,l3+m3,l3-m3,
     &                      l1+l2-l3,l1-l2+l3,-l1+l2+l3))
     &  STOP 'wigner3j: Number of terms incorrect.'
      
      IF(tmin.le.tmax) THEN
        DO t = tmin,tmax
          wigner3j = wigner3j + (-1)**t /
     &                          (  fac(t)    * fac(f1+t) * fac(f2+t) 
     &                            *fac(f3-t) * fac(f4-t) * fac(f5-t) )
        END DO
        wigner3j = wigner3j * (-1)**(l1-l2-m3) * sfac(l1+l2-l3)
     &                      * sfac(l1-l2+l3)   * sfac(-l1+l2+l3) 
     &                      / sfac(l1+l2+l3+1) *
     &                        sfac(l1+m1)      * sfac(l1-m1) *
     &                        sfac(l2+m2)      * sfac(l2-m2) *
     &                        sfac(l3+m3)      * sfac(l3-m3)
      END IF
      END FUNCTION wigner3j

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

c     Integrates function f numerically (Lagrange and Simpson integration)
c     on grid(itype) and is much faster than intgr. 
c     (Only normal outward integration.)
c     Before first use of this function it has to be initialized with 
c     intgrf_init.

      FUNCTION intgrf(f,jri,jmtd,rmsh,dx,ntype,itype,gridf)


      IMPLICIT NONE

      REAL                 :: intgrf
      INTEGER, INTENT(IN)  :: itype,ntype,jmtd
      INTEGER, INTENT(IN)  :: jri(ntype)
      REAL   , INTENT(IN)  :: dx(ntype),rmsh(jmtd,ntype)
      REAL   , INTENT(IN)  :: gridf(jmtd,ntype)
      REAL   , INTENT(IN)  :: f(*)
c     - local -
      TYPE(intgrf_out)     :: fct_res


      fct_res = pure_intgrf(f,jri,jmtd,rmsh,dx,ntype,itype,gridf)
      IF (fct_res%ierror == NEGATIVE_EXPONENT_WARNING) THEN
        write(6,*) 'intgrf: Warning!'//
     +               'Negative exponent x in extrapolation a+c*r**x'
      ELSEIF (fct_res%ierror == NEGATIVE_EXPONENT_ERROR) THEN
        write(6,*)
     +       'intgrf: Negative exponent x in extrapolation a+c*r**x'
        CALL juDFT_error(
     +          'intgrf: Negative exponent x in extrapolation a+c*r**x')
      END IF
      intgrf = fct_res%value

      END FUNCTION intgrf

c     pure wrapper for intgrf with same functionality
c     can be used within forall loops

      PURE FUNCTION pure_intgrf(f,jri,jmtd,rmsh,dx,ntype,itype,gridf)

      IMPLICIT NONE

      TYPE(intgrf_out)     :: pure_intgrf

      INTEGER, INTENT(IN)  :: itype,ntype,jmtd
      INTEGER, INTENT(IN)  :: jri(ntype)
      REAL   , INTENT(IN)  :: dx(ntype),rmsh(jmtd,ntype)
      REAL   , INTENT(IN)  :: gridf(jmtd,ntype)
      REAL   , INTENT(IN)  :: f(*)
c     - local -
      INTEGER              :: n
      REAL                 :: r1,h,a,x

      n  = jri(itype)
      r1 = rmsh(1,itype)
      h  = dx(itype)

      pure_intgrf%ierror = NO_ERROR

      ! integral from 0 to r1 approximated by leading term in power series expansion
      IF (f(1)*f(2).gt.1e-10.and.h.gt.0) THEN
        IF(f(2).eq.f(1)) THEN
          pure_intgrf%value = r1*f(1)
        ELSE
          x      = (f(3)-f(2)) / (f(2)-f(1))
          a      = (f(2)-x*f(1)) / (1-x)
          x      = log(x)/h
          IF(x.lt.0) THEN
            IF(x.gt.-1) THEN
              pure_intgrf%ierror = NEGATIVE_EXPONENT_WARNING
            ELSE IF(x.le.-1) THEN
              pure_intgrf%ierror = NEGATIVE_EXPONENT_ERROR
              RETURN
            END IF
          END IF

          pure_intgrf%value = r1*(f(1)+x*a) / (x+1)

!           x      = f(2) / f(1)
!           x      = log(x)/h
!           IF(x.lt.0) THEN
!             IF(x.gt.-1) write(6,'(A,ES9.1)') 'intgrf: Warning!&
!      &                                        Negative exponent x in& 
!      &                                        extrapolation c*r**x:',x
!             IF(x.le.-1) write(6,'(A,ES9.1)') 'intgrf: Negative exponent&
!      &                                        x in extrapolation&
!      &                                        c*r**x:',x
!             IF(x.le.-1) STOP                 'intgrf: Negative exponent& 
!      &                                        x in extrapolation &
!      &                                        c*r**x'
!           END IF
!           intgrf = (r1*f(1))/(x+1)

        END IF
      ELSE
        pure_intgrf%value = 0
      END IF 
      
      ! integrate from r(1) to r(n) by multiplying with gridf
      pure_intgrf%value = pure_intgrf%value 
     +                  + dot_product(gridf(:n,itype),f(:n))
      
      END FUNCTION pure_intgrf


c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

c     Initializes fast numerical integration intgrf (see below).
      
      SUBROUTINE intgrf_init(ntype,jmtd,jri,dx,rmsh,gridf)

      IMPLICIT NONE

      INTEGER, INTENT(IN)   :: ntype,jmtd
      INTEGER, INTENT(IN)   :: jri(ntype)
      REAL,    INTENT(IN)   :: dx(ntype),rmsh(jmtd,ntype)
      REAL, INTENT(OUT), ALLOCATABLE  :: gridf(:,:)

c     - local -
      INTEGER              :: i,j,itype
      INTEGER              :: n,nstep,n0 = 6
      INTEGER, PARAMETER   :: simpson(7) = (/41,216,27,272,27,216,41/)
      REAL                 :: r1,h,dr
      REAL                 :: r(7)     
      REAL, PARAMETER      :: lagrange(7,6)= reshape(
     &          (/19087.,65112.,-46461., 37504.,-20211., 6312., -863.,
     &             -863.,25128., 46989.,-16256.,  7299.,-2088.,  271.,
     &              271.,-2760., 30819., 37504., -6771., 1608., -191.,
     &             -191., 1608., -6771., 37504., 30819.,-2760.,  271.,
     &              271.,-2088.,  7299.,-16256., 46989.,25128., -863.,
     &             -863., 6312.,-20211., 37504.,-46461.,65112.,19087./), ! The last row is actually never used.
     &                                      (/7,6/) )

      n = jmtd
      ALLOCATE ( gridf(n,ntype) )

      gridf = 0

      DO itype = 1,ntype

        n  = jri(itype)
        r1 = rmsh(1,itype)
        h  = dx(itype)

        nstep = (n-1)/6
        n0    = n-6*nstep
        dr    = exp(h)

        ! Calculate Lagrange-integration coefficients from r(1) to r(n0)
        r(1)=r1
        IF(n0.gt.1) THEN
          DO i=2,7
            r(i) = r(i-1)*dr
          END DO
          DO i=1,7
            gridf(i,itype) = h/60480 * r(i) * sum(lagrange(i,1:n0-1))
          END DO
          r(1)  = r(n0)
        END IF

        ! Calculate Simpson-integration coefficients from r(n0) to r(n)
        DO i = 1,nstep
          DO j = 2,7
            r(j) = r(j-1)*dr
          END DO
          DO j = n0,n0+6
            gridf(j,itype) = gridf(j,itype) + h/140 * r(j-n0+1) *
     &                       simpson(j-n0+1)
          END DO
          n0   = n0 + 6
          r(1) = r(7)
        END DO

      END DO
        
      END SUBROUTINE intgrf_init

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Calculates the primitive of f, on grid(itypein):
c
c                   r
c     primf(r) = integral f(r') dr'   ( r = grid point )
c                   0
c
c     If itypein is negative, the primitive
c     
c                   R
c     primf(r) = integral f(r') dr'   ( R = MT sphere radius )
c                   r
c
c     is calculated instead.
c

c     -----------------------------

c     Fast calculation of primitive.
c     Only Lagrange integration is used

      SUBROUTINE primitivef(primf,fin,rmsh,dx,jri,jmtd,itypein,ntype)

      IMPLICIT NONE

c     - scalars -
      INTEGER, INTENT(IN)   :: itypein,jmtd,ntype

c     - arrays -
      INTEGER, INTENT(IN)   :: jri(ntype)
      REAL,    INTENT(OUT)  :: primf( jri( abs(itypein) ) )
      REAL,    INTENT(IN)   :: fin( jri( abs(itypein) ) )
      REAL,    INTENT(IN)   :: rmsh(jmtd,ntype),dx(ntype)

c     - local scalars -
      INTEGER               :: itype,n,i,n0
      REAL                  :: h,x,h1
      REAL                  :: intgr,r1,a,dr

c     - local arrays -
      REAL                  :: fr(7)
      REAL                  :: f( jri( abs(itypein) ))
      REAL                  :: r( jri( abs(itypein) ) )
      REAL, PARAMETER       :: lagrange(7,6)= reshape(
     &          (/19087.,65112.,-46461., 37504.,-20211., 6312., -863.,
     &             -863.,25128., 46989.,-16256.,  7299.,-2088.,  271.,
     &              271.,-2760., 30819., 37504., -6771., 1608., -191.,
     &             -191., 1608., -6771., 37504., 30819.,-2760.,  271.,
     &              271.,-2088.,  7299.,-16256., 46989.,25128., -863.,
     &             -863., 6312.,-20211., 37504.,-46461.,65112.,19087./),
     &                                       (/7,6/) )
    

      itype = abs(itypein)

      primf = 0

      n = jri(itype)
      h = dx(itype)

      IF(itypein.gt.0) THEN
        r1 = rmsh(1,itype)      ! perform outward integration
        f  = fin                ! (from 0 to r)
      ELSE
        r1 = rmsh(jri(itype),itype)         ! perform inward integration
        h  = -h                             ! (from MT sphere radius to r)
        f  = fin(n:1:-1)                    !
      END IF
      
      ! integral from 0 to r1 approximated by leading term in power series expansion (only if h>0)
      IF(h.gt.0.and.f(1)*f(2).gt.1e-10) THEN
        IF(f(2).eq.f(1)) THEN
          intgr = r1*f(1)
        ELSE
          x     = (f(3)-f(2))/(f(2)-f(1))
          a     = (f(2)-x*f(1)) / (1-x)
          x     = log(x)/h
          IF(x.lt.0) THEN
             IF(x>-1) WRITE(6,'(A,ES9.1)')
     +            '+intgr: Warning! Negative &exponent x in'//
     +        'extrapolation a+c*r**x:',x
             IF(x<=-1) WRITE(6,'(A,ES9.1)')'intgr: Negative exponent,'//
     +             'x in extrapolation a+c*r**x:',x
             IF(x<=-1) STOP 'intgr:Negative exponent x in extrapolation'
          END IF
          intgr = r1*(f(1)+x*a) / (x+1)
        END IF
      ELSE
        intgr = 0
      END IF

      primf(1) = intgr
      dr       = exp(h)
      r(1)     = r1
      n0       = 1
      h1       = h/60480
      
      ! Lagrange integration from r(n0) to r(n0+5)
 1    DO i=2,7
        r(i) = r(i-1)*dr
      END DO
      fr = f(n0:n0+6) * r(:7)
      DO i=1,6
        intgr = intgr + h1 * dot_product(lagrange(:,i),fr)
        IF(primf(n0+i).eq.0) primf(n0+i) = intgr ! avoid double-definition
      END DO
      IF(n0+12.le.n) THEN
        r(1) = r(7)
        n0   = n0 + 6
        GOTO 1
      ELSE IF(n0+6.lt.n) THEN
        r(1)  = r(n-5-n0)
        n0    = n-6
        intgr = primf(n-6)
        GOTO 1
      END IF

      IF(itypein.lt.0) THEN    !
        primf = -primf(n:1:-1) ! Inward integration
      END IF                    !

      END SUBROUTINE primitivef

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      !function modulo1 maps kpoint into first BZ
      FUNCTION modulo1(kpoint,nkpt3)

      IMPLICIT NONE

      INTEGER,INTENT(IN)  :: nkpt3(3)
      REAL   ,INTENT(IN)  :: kpoint(3)
      REAL                :: modulo1(3)
      INTEGER             :: help(3)

      modulo1 = kpoint*nkpt3
      help    = nint(modulo1)
      IF(any(abs(help-modulo1).gt.1e-10)) THEN
        WRITE(6,'(A,F5.3,2('','',F5.3),A)')'modulo1: argument (',kpoint,
     +                         ') is not an element of the k-point set.'
        CALL juDFT_error(
     +     'modulo1: argument not an element of k-point set.', 
     +     calledby = 'util:modulo1')
      END IF
      modulo1 = modulo(help,nkpt3)*1.0/nkpt3

      END FUNCTION modulo1
      


c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

c     Returns derivative of f in df.

      SUBROUTINE derivative_t(df,f,atoms,itype)
      USE m_types_atoms
      IMPLICIT NONE
      REAL,   INTENT(IN)   ::   f(:)
      REAL,   INTENT(OUT)  ::   df(:)
      TYPE(t_atoms),INTENT(IN)::atoms
      INTEGER,INTENT(IN)   ::   itype

      call derivative_nt(df,f,atoms%jmtd,atoms%jri,atoms%dx,atoms%rmsh,
     +     atoms%ntype,itype)
      END SUBROUTINE
      
      SUBROUTINE derivative_nt(df,f,jmtd,jri,dx,rmsh,ntype,itype)

      IMPLICIT NONE

      INTEGER,INTENT(IN)   :: ntype,itype,jmtd
      INTEGER,INTENT(IN)   :: jri(ntype)
      REAL,   INTENT(IN)   :: dx(ntype),rmsh(jmtd,ntype)
      REAL,   INTENT(IN)   :: f(jri(itype))
      REAL,   INTENT(OUT)  :: df(jri(itype))
      REAL                 :: x,h,r,d21,d32,d43,d31,d42,d41,df1,df2
      INTEGER              :: i,n
      
      n = jri(itype)
      h = dx(itype)
      r = rmsh(1,itype)
      ! use power series expansion a+c**x for first point
      IF(f(2).eq.f(1)) THEN
        df(1) = 0.0
      ELSE
        x     = (f(3)-f(2))/(f(2)-f(1))
        df(1) = (f(2)-f(1)) / (x-1) * log(x)/h / r
      END IF
      ! use Lagrange interpolation of 3rd order for all other points (and averaging)
      d21 = r * (exp(h)-1) ; d32 = d21 * exp(h) ; d43 = d32 * exp(h)
      d31 = d21 + d32        ; d42 = d32 + d43
      d41 = d31 + d43
      df(2) = - d32*d42 / (d21*d31*d41)         * f(1) 
     &        + ( 1.0/d21 - 1.0/d32 - 1.0/d42)  * f(2)
     &        + d21*d42 / (d31*d32*d43)         * f(3)
     &        - d21*d32 / (d41*d42*d43)         * f(4)
      df1   =   d32*d43 / (d21*d31*d41)         * f(1) 
     &        - d31*d43 / (d21*d32*d42)         * f(2)
     &        + ( 1.0/d31 + 1.0/d32 - 1.0/d43 ) * f(3)
     &        + d31*d32 / (d41*d42*d43)         * f(4)
      DO i = 3,n-2
        d21 = d32 ; d32 = d43 ; d43 = d43 * exp(h)
        d31 = d42 ; d42 = d42 * exp(h)
        d41 = d41 * exp(h)
        df2   = - d32*d42 / (d21*d31*d41)         * f(i-1) 
     &          + ( 1.0/d21 - 1.0/d32 - 1.0/d42)  * f(i)
     &          + d21*d42 / (d31*d32*d43)         * f(i+1)
     &          - d21*d32 / (d41*d42*d43)         * f(i+2)
        df(i) = ( df1 + df2 ) / 2
        df1   =   d32*d43 / (d21*d31*d41)         * f(i-1) 
     &          - d31*d43 / (d21*d32*d42)         * f(i)
     &          + ( 1.0/d31 + 1.0/d32 - 1.0/d43 ) * f(i+1)
     &          + d31*d32 / (d41*d42*d43)         * f(i+2)
      END DO
      df(n-1) = df1
      df(n)   = - d42*d43 / (d21*d31*d41)         * f(n-3) 
     &          + d41*d43 / (d21*d32*d42)         * f(n-2) 
     &          - d41*d42 / (d31*d32*d43)         * f(n-1) 
     &          + ( 1.0/d41 + 1.0/d42 + 1.0/d43 ) * f(n)

      END SUBROUTINE derivative_nt


c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Orders iarr(1:n) according to size and returns a correspondingly defined pointer in pnt

      SUBROUTINE iorderp(pnt,iarr,n)

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: n
      INTEGER, INTENT(OUT) :: pnt(1:n)
      INTEGER, INTENT(IN)  :: iarr(1:n)
      INTEGER              :: i,j,k

      DO i=1,n
        pnt(i) = i
        DO j=1,i-1
          IF(iarr(pnt(j)).gt.iarr(i)) THEN
            DO k=i,j+1,-1
              pnt(k) = pnt(k-1)
            END DO
            pnt(j) = i
            EXIT
          END IF
        END DO
      END DO

      END SUBROUTINE iorderp

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Orders rarr(1:n) according to size and returns a correspondingly defined pointer in pnt

      SUBROUTINE rorderp(pnt,rarr,n)

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: n
      INTEGER, INTENT(OUT) :: pnt(1:n)
      REAL,    INTENT(IN)  :: rarr(1:n)
      INTEGER              :: i,j,k

      DO i=1,n
        pnt(i) = i
        DO j=1,i-1
          IF(rarr(pnt(j)).gt.rarr(i)) THEN
            DO k=i,j+1,-1
              pnt(k) = pnt(k-1)
            END DO
            pnt(j) = i
            EXIT
          END IF
        END DO
      END DO

      END SUBROUTINE rorderp


c     Same as rorderp but divides the problem in halves np times (leading to 2**np intervals) and is
c     much faster than rorderp (devide and conquer algorithm).
c     There is an optimal np, while for larger np the overhead (also memory-wise) outweights the speed-up.
c     np = max(0,int(log(n*0.001)/log(2.0))) should be a safe choice.
      RECURSIVE SUBROUTINE rorderpf(pnt,rarr,n,np)

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: n,np
      INTEGER, INTENT(OUT) :: pnt(n)
      REAL   , INTENT(IN)  :: rarr(n)
      REAL                 :: rarr1(n)
      REAL                 :: ravg
      INTEGER              :: pnt1(n),pnt2(n)
      INTEGER              :: n1,n2,i
      
      IF(np.eq.0) THEN
        CALL rorderp(pnt,rarr,n)
        RETURN
      ELSE IF(np.lt.0) THEN
        STOP 'rorderpf: fourth argument must be non-negative (bug?).'
      END IF
      ravg = sum(rarr)/n
      ! first half
      n1 = 0
      DO i = 1,n
        IF(rarr(i).le.ravg) THEN
          n1        = n1 + 1
          rarr1(n1) = rarr(i)
          pnt1(n1)  = i
        END IF
      END DO
      CALL rorderpf(pnt2,rarr1,n1,np-1)
      pnt(:n1) = pnt1(pnt2(:n1))
      ! second half
      n2 = 0
      DO i = 1,n
        IF(rarr(i).gt.ravg) THEN
          n2        = n2 + 1
          rarr1(n2) = rarr(i)
          pnt1(n2)  = i
        END IF
      END DO
      CALL rorderpf(pnt2,rarr1,n2,np-1)
      pnt(n1+1:) = pnt1(pnt2(:n2))
      END  SUBROUTINE rorderpf


c     Calculates the spherical Bessel functions of orders 0 to l at x
c     by backward recurrence using j_l(x) = (2l+3)/x j_l+1(x) - j_l+2(x) .
c     (Starting points are calculated according to Zhang, Min, 
c      "Computation of Special Functions".)

      pure SUBROUTINE sphbessel(sphbes,x,l)

      IMPLICIT NONE

      INTEGER, INTENT(IN)    :: l
      REAL   , INTENT(IN)    :: x
      REAL   , INTENT(INOUT) :: sphbes(0:l)
      REAL                   :: s0,s1,f,f0,f1,cs
      INTEGER                :: ll,lsta,lmax,msta2

!      IF( x.lt.0 ) THEN
!        STOP 'sphbes: negative argument (bug?).' 
!      ELSE 
      IF( x.eq.0 ) THEN
        sphbes(0)  = 1.0
        DO ll = 1,l
          sphbes(ll) = 0.0
        END DO
        RETURN
      ENDIF
      sphbes(0) = sin(x) / x
      IF( l .eq. 0 ) RETURN
      sphbes(1) = ( sphbes(0) - cos(x) ) / x
!       IF(l.le.1) RETURN
      s0   = sphbes(0)
      s1   = sphbes(1)
      lsta = lsta1(x,200)      !
      lmax = l                 !
      IF(lsta.lt.l) THEN       !
        lmax            = lsta ! determine starting point lsta
        sphbes(lmax+1:) = 0.0  ! for backward recurrence
      ELSE                     !
        lsta = lsta2(x,l,15)   ! 
      END IF                    !
      f0 = 0.0                                                      !
      f1 = 1e-100                                                   !
      DO ll = lsta,0,-1                                             ! backward recurrence
        f  = f1 / x * (2*ll+3) - f0 ; IF(ll.le.lmax) sphbes(ll) = f ! with arbitrary start values
        f0 = f1                                                     !
        f1 = f                                                      !
      END DO                                                        !
      IF(abs(s0).gt.abs(s1)) THEN ; cs = s0 / f   !
      ELSE                        ; cs = s1 / f0  ! scale to correct values
      END IF                                      ! 
      sphbes = cs * sphbes                        !

      CONTAINS

c     Test starting point
      PURE FUNCTION lsta0(x,mp)

      IMPLICIT NONE

      INTEGER             :: lsta0
      INTEGER, INTENT(IN) :: mp
      REAL   , INTENT(IN) :: x
      REAL                :: f,lgx

      lgx   = log10(x)
      lsta0 = 0
      f     = lgx
      DO WHILE(f.gt.-mp)
        lsta0 = lsta0 + 1
        f     = f + lgx - log10(2.0*lsta0+1)
      END DO
      END FUNCTION lsta0

c     Returns starting point lsta1 for backward recurrence such that sphbes(lsta1) approx. 10^(-mp).
      PURE FUNCTION lsta1(x,mp)

      IMPLICIT NONE

      INTEGER             :: lsta1
      INTEGER, INTENT(IN) :: mp
      REAL   , INTENT(IN) :: x
      REAL                :: f0,f1,f
      INTEGER             :: n0,n1,nn,it

      n0 = int(1.1*x) + 1
      f0 = envj(n0,x) - mp
      n1 = n0 + 5
      f1 = envj(n1,x) - mp
      DO it = 1,20
        nn = n1 - (n1-n0) / (1.0-f0/f1)
        f  = envj(nn,x) - mp
        IF(abs(nn-n1).lt.1) EXIT
        n0 = n1
        f0 = f1
        n1 = nn
        f1 = f
      END DO
      lsta1 = nn
      END FUNCTION lsta1

c     Returns the starting point lsta2 for backward recurrence such that all sphbes(l) have mp significant digits. 
      PURE FUNCTION lsta2(x,l,mp)

      IMPLICIT NONE

      INTEGER             :: lsta2
      INTEGER, INTENT(IN) :: l,mp
      REAL   , INTENT(IN) :: x
      REAL                :: f0,f1,f,hmp,ejn,obj
      INTEGER             :: n0,n1,nn,it

      hmp = 0.5 * mp
      ejn = envj(l,x)
      IF( ejn.le.hmp ) THEN
        obj = mp
        n0  = int(1.1*x) + 1
      ELSE
        obj = hmp + ejn
        n0  = l
      END IF
      f0 = envj(n0,x) - obj
      n1 = n0 + 5
      f1 = envj(n1,x) - obj
      DO it = 1,20
        nn = n1 - (n1-n0) / (10-f0/f1)
        f  = envj(nn,x) - obj
        IF(abs(nn-n1).lt.1) EXIT
        n0 = n1
        f0 = f1
        n1 = nn
        f1 = f
      END DO
      lsta2 = nn + 10
      END FUNCTION lsta2


      PURE FUNCTION envj(l,x)

      IMPLICIT NONE

      REAL                :: envj
      REAL   , INTENT(IN) :: x
      INTEGER, INTENT(IN) :: l

      envj = 0.5 * log10(6.28*l) - l*log10(1.36*x/l)

      END FUNCTION envj

      END SUBROUTINE sphbessel



c     Returns the spherical harmonics Y_lm(^rvec)
c     for l = 0,...,ll in Y(1,...,(ll+1)**2).

      PURE SUBROUTINE harmonicsr(Y,rvec,ll)

      IMPLICIT NONE

      integer , intent(in)  :: ll
      real    , intent(in)  :: rvec(3)
      complex , intent(out) :: Y((ll+1)**2)
      complex               :: c
      real                  :: stheta,ctheta,sphi,cphi,r,rvec1(3)
      integer               :: l,m,lm
      complex , parameter   :: img = (0.0,1.0)
 
      Y(1) = 0.282094791773878
      IF(ll.eq.0) RETURN

      stheta = 0
      ctheta = 0
      sphi   = 0
      cphi   = 0
      r      = sqrt(sum(rvec**2))
      IF(r.gt.1e-16) THEN
        rvec1  = rvec / r
        ctheta = rvec1(3)
        stheta = sqrt(rvec1(1)**2+rvec1(2)**2)
        IF(stheta.gt.1e-16) THEN
          cphi = rvec1(1) / stheta
          sphi = rvec1(2) / stheta
        END IF
      ELSE
        Y(2:) = 0.0
        RETURN
      END IF

c     define Y,l,-l and Y,l,l
      r = Y(1)
      c = 1
      DO l=1,ll
        r           = r*stheta*sqrt(1.0+1.0/(2*l))
        c           = c * (cphi + img*sphi)
        Y(l**2+1)   = r*conjg(c)  ! l,-l
        Y((l+1)**2) = r*c*(-1)**l ! l,l
      END DO

c     define Y,l,-l+1 and Y,l,l-1
      Y(3) = 0.48860251190292*ctheta
      DO l=2,ll
        r          = sqrt(2.0*l+1) * ctheta
        Y(l**2+2)  = r*Y((l-1)**2+1) ! l,-l+1
        Y(l*(l+2)) = r*Y(l**2)       ! l,l-1
      END DO

c     define Y,l,m, |m|<l-1
      DO l=2,ll
        lm = l**2 + 2
        DO m=-l+2,l-2
          lm = lm + 1
          Y(lm) = sqrt((2.0*l+1)/(l+m)/(l-m)) * (
     &            sqrt(2.0*l-1)*ctheta           *Y(lm-2*l)-
     &            sqrt((l+m-1.0)*(l-m-1)/(2*l-3))*Y(lm-4*l+2) )
        END DO
      END DO

      END SUBROUTINE harmonicsr

c     Returns the complex error function.
      FUNCTION cerf(z)
      
      USE m_constants , ONLY: pimach
      
      IMPLICIT NONE
      COMPLEX, INTENT(IN) ::  z
      COMPLEX             ::  cerf
      COMPLEX             ::  z1,z2,c,d,delta
      REAL                ::  pi
      INTEGER             ::  i
      
      
      pi = pimach()
      z1 = z
      IF(real(z).lt.0) z1 = -z1
      IF(real(z1).lt.2.0) THEN ! McLaurin series
        z2   = z1**2
        i    = 0
        c    = z1
        cerf = z1
        DO
          i    = i + 1
          c    = -c * z2 / i
          cerf = cerf + c/(2*i+1)
          IF(abs(c/(2*i+1)).lt.1e-20) EXIT
        END DO
        cerf = cerf * 2/sqrt(pi)
      ELSE ! continued fraction using Lentz's method
        d    = 0.0
        c    = z1
        cerf = z1
        i    = 0
        DO
          i     = i + 1
          c     =   2*z1 + i/c
          d     = ( 2*z1 + i*d )**(-1)
          delta = c*d
          cerf  = cerf * delta
          IF(abs(1-delta).lt.1e-15) EXIT
          i     = i + 1
          c     =   z1 + i/c
          d     = ( z1 + i*d )**(-1)
          delta = c*d
          cerf  = cerf * delta
          IF(abs(1-delta).lt.1e-15) EXIT
          IF(i.eq.10000) 
     &        STOP 'cerf: Lentz method not converged after 10000 steps.'
        END DO
        cerf = 1 - exp(-z1**2) / cerf / sqrt(pi)        
      END IF
      IF(real(z).lt.0) cerf = -cerf
      END FUNCTION cerf

      FUNCTION chr(int)
      
      IMPLICIT NONE
      
      CHARACTER(5) :: chr
      INTEGER      :: int
      
      WRITE(chr,'(I5)') int
      END FUNCTION chr

      END MODULE m_util
