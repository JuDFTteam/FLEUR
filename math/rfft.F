      MODULE m_rfft
      use m_juDFT
      PRIVATE
      !Module also contains routines vrffti,vrfftf&vrfftb as private routines below
      PUBLIC rfft
      CONTAINS
      SUBROUTINE rfft(
     >                isn,n1d,n2d,n3d,n1,n2,n3,
     >                nw1,nw2,nw3,wsave,b,
     X                a)
c **********************************************************************
c
c isn = +1 : FFT a real, 3dim array "a" of dimensions n1*n2*n3 to a complex
c            array of size n1*n2*(n3+1)/2 [if n3 is odd] or n1*n2*(n3/2+1)
c            [n3 is even]. 
c isn = -1 : back-FFT of a complex array "a" of size n1*n2*(n3+1)/2 [odd] 
c            or n1*n2*(n3/2+1) [even] to a real, 3dim array 
c
c the actual array is assumed to be located between 
c                                  1 ... nw1+1    &     n1-nw1+1 ... n1
c                                  1 ... nw2+1    &     n2-nw2+1 ... n2
c                                  1 ... nw3+1    &     n3-nw3+1 ... n3
c and padded with zeros in between.
c if nw1 >= (n1-1)/2, no padding is assumed (-> aliasing errors!)
c
c                                                  G.Bihlmayer (UniWien)
c         
c **********************************************************************
      USE m_cfft
      IMPLICIT NONE

      INTEGER n1d,n2d,n3d,n1,n2,n3,nw1,nw2,nw3,isn
      REAL a(n1d,n2d,0:n3d),b(n1d,n2d,n3d),wsave(n3d+15)

      INTEGER i1,i2,i3,nup
      REAL factor
      LOGICAL l_nopad
c
c a ... array for FFT
c b ... work array
c wsave ... stores tables for r-2-c FFT
c n1,n2,n3 ... dimensions to be transformed
c nw1,nw2,nw3 ... actual dimensions of a before FFT
c n1d,n2d,n3d ... dimensions of a,b 
c
c
c check for input errors
c
      IF ((isn/=-1).AND.(isn.NE.1))  CALL juDFT_error("choose isn=+/- 1"
     +     ,calledby ="rfft")
      IF ((n1d.lt.n1).OR.(n2d.lt.n2).OR.(n3d.lt.n3)) THEN
        WRITE (6,*) 'n1d,n2d,n3d =',n1d,n2d,n3d
        WRITE (6,*) 'n1 ,n2 ,n3  =',n1 ,n2 ,n3 
        CALL juDFT_error("n(i) > n(i)d",calledby ="rfft")
      ENDIF
      IF ((n1.le.2*nw1+1).OR.
     +    (n2.le.2*nw2+1).OR.
     +    (n3.le.2*nw3+1)) THEN
c        WRITE (6,*) 'n1 ,n2 ,n3  =',n1 ,n2 ,n3 
c        WRITE (6,*) 'nw1,nw2,nw3 =',nw1,nw2,nw3
        l_nopad=.true.
      ELSE
        l_nopad=.false.
      ENDIF
c
c ******************** f o r w a r d - t r a n s f o r m *******************     
c      
      IF (isn.eq.1) THEN
c
c array a is assumed to be zero from (1,1,0) to (n1d,n2d,0) and the array 
c to be FFT'd starts at (1,1,1) as n1*n2*n3 real numbers. 
c first transform n1*n2 real sequences of lenghth n3 to n3/2 complex values
c
        CALL vrffti(n3,wsave) 
        IF (l_nopad) THEN
          CALL vrfftf(n1*n2,n3,a(1,1,1),b,n1d*n2d,wsave)
        ELSE
          DO i2=1,nw2+1
            CALL vrfftf(nw1+1,n3,a(1,i2,1),b,n1d*n2d,wsave)
            CALL vrfftf(nw1  ,n3,a(n1-nw1+1,i2,1),b,n1d*n2d,wsave)
          ENDDO
          DO i2=n2-nw2+1,n2
            CALL vrfftf(nw1+1,n3,a(1,i2,1),b,n1d*n2d,wsave)
            CALL vrfftf(nw1  ,n3,a(n1-nw1+1,i2,1),b,n1d*n2d,wsave)
          ENDDO
        ENDIF
c
c now we have the FFT'd array stored as described in vrfftf 
c (mixed real & compex data)
c remove the norm 1/sqrt(n3) (to be compatible with cfft)
c
        factor = sqrt(1.0*n3)
        !CALL CPP_BLAS_sscal(n1d*n2d*n3,factor,a(1,1,1),1)
        a=a*factor
c
c now, the real part of f(0) has to be moved to a(n1,n2,0) to get a purely
c complex array starting at a(1,1,0)
c
        DO i1=1,n1
          DO i2=1,n2
            a(i1,i2,0)=a(i1,i2,1) 
            a(i1,i2,1)=0.0 
          ENDDO
        ENDDO
c
      ENDIF
c
c ******************** g e n e r a l  p a r t *******************
c
c now perform n2*n3/2 and n1*n3/2 complex FFT's; a is assumed to be
c complex and starting at (1,1,0)
c
      IF (ABS((n3/2.)-NINT(n3/2.)).gt.0.1) THEN
        nup = n3
      ELSE
        nup = n3+1
        IF (n3+1>n3d)  CALL juDFT_error("n3 even & n3+1 > n3d" ,calledby
     +       ="rfft")
        a(:,:,n3+1)=0.0
      ENDIF

      IF (l_nopad) THEN
        DO i3=1,nup,2
          CALL cfft(a(1,1,i3-1),a(1,1,i3),n1*n2,n1,n1,isn)
          CALL cfft(a(1,1,i3-1),a(1,1,i3),n1*n2,n2,n1*n2,isn)
        ENDDO
      ELSE
        DO i3=1,nup,2
          CALL cfft(a(1,1,i3-1),a(1,1,i3),(nw2+1)*n1,n1,n1,isn)
          CALL cfft(a(1,n2-nw2+1,i3-1),a(1,n2-nw2+1,i3),
     +                                        nw2*n1,n1,n1,isn)
          CALL cfft(a(1,1,i3-1),a(1,1,i3),n1*n2,n2,n1*n2,isn)
        ENDDO
      ENDIF

c
c ******************** b a c k w a r d - t r a n s f o r m *******************     
c      
      IF (isn.eq.-1) THEN
c
c the real part of f(0) has to be moved to a(n1,n2,1) for a correct
c setup for vrfftb (see comments therein)
c
        DO i1=1,n1
          DO i2=1,n2
            a(i1,i2,1)=a(i1,i2,0)
            a(i1,i2,0)=0.0
          ENDDO
        ENDDO
c
c transform n1*n2 mixed real and complex sequences of lenghth n3/2
c to n3 real values
c
        CALL vrffti(n3,wsave)
        IF (l_nopad) THEN
          CALL vrfftb(n1*n2,n3,a(1,1,1),b,n1d*n2d,wsave)
        ELSE
          DO i2=1,nw2+1
            CALL vrfftb(nw1+1,n3,a(1,i2,1),b,n1d*n2d,wsave)
            CALL vrfftb(nw1  ,n3,a(n1-nw1+1,i2,1),b,n1d*n2d,wsave)
          ENDDO
          DO i2=n2-nw2+1,n2
            CALL vrfftb(nw1+1,n3,a(1,i2,1),b,n1d*n2d,wsave)
            CALL vrfftb(nw1  ,n3,a(n1-nw1+1,i2,1),b,n1d*n2d,wsave)
          ENDDO
        ENDIF
c
c remove the norm 1/sqrt(n3) (compatibility with cfft)
c
        factor = sqrt(1.0*n3)
        !CALL CPP_BLAS_sscal(n1d*n2d*n3,factor,a(1,1,1),1)
        a=a*factor
c
      ENDIF

      END SUBROUTINE rfft

      SUBROUTINE vrffti(n,wsave)
C***BEGIN PROLOGUE  VRFFTI
C***DATE WRITTEN   860701   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM,
C             MULTIPLE SEQUENCES
C***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST)
C***PURPOSE  Initialization for VRFFTF and VRFFTB.
C***DESCRIPTION
C
C  Subroutine VRFFTI initializes the array WSAVE which is used in
C  both VRFFTF and VRFFTB.  The prime factorization of N together with
C  a tabulation of certain trigonometric functions are computed and
C  stored in the array WSAVE.
C
C  Input Parameter
C
C  N       the length of the sequence to be transformed.  There is no
C          restriction on N.
C
C  Output Parameter
C
C  WSAVE   a work array which must be dimensioned at least N+15.
C          The same work array can be used for both VRFFTF and VRFFTB
C          as long as N remains unchanged.  Different WSAVE arrays
C          are required for different values of N.  The contents of
C          WSAVE must not be changed between calls of VRFFTF or VRFFTB.
C
C
C              * * * * * * * * * * * * * * * * * * * * *
C              *                                       *
C              *         PROGRAM SPECIFICATIONS        *
C              *                                       *
C              * * * * * * * * * * * * * * * * * * * * *
C
C
C     Dimension of    R(MDIMR,N), RT(MDIMR,N), WSAVE(N+15)
C     Arguments
C
C     Latest          AUGUST 1, 1985
C     Revision
C
C     Subprograms     VRFFTI, VRFTI1, VRFFTF, VRFTF1, VRADF2, VRADF3,
C     Required        VRADF4, VRADF5, VRADFG, VRFFTB, VRFTB1, VRADB2,
C                     VRADB3, VRADB4, VRADB5, VRADBG, PIMACH
C
C     Special         NONE
C     Conditions
C
C     Common          NONE
C     blocks
C
C     I/O             NONE
C
C     Precision       SINGLE
C
C     Specialist      ROLAND SWEET
C
C     Language        FORTRAN
C
C     History         WRITTEN BY LINDA LINDGREN AND ROLAND SWEET AT THE
C                     NATIONAL BUREAU OF STANDARDS (BOULDER).
C
C     Algorithm       A REAL VARIANT OF THE STOCKHAM AUTOSORT VERSION
C                     OF THE COOLEY-TUKEY FAST FOURIER TRANSFORM.
C
C     Portability     AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN 77.
C                     THE ONLY MACHINE DEPENDENT CONSTANT IS LOCATED IN
C                     THE FUNCTION PIMACH.
C
C     Required        COS,SIN
C     resident
C     Routines
C
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
C               pp. 51-83.
C***ROUTINES CALLED  VRFTI1
C***END PROLOGUE  VRFFTI
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION wsave(n+15)
C***FIRST EXECUTABLE STATEMENT  VRFFTI
      IF (n.EQ.1) RETURN
      CALL vrfti1(n,wsave(1),wsave(n+1))
      RETURN
      END subroutine
      SUBROUTINE vrfti1(n,wa,fac)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      USE m_constants, ONLY : pimach
      DIMENSION wa(n),fac(15),ntryh(4)
      DATA ntryh(1),ntryh(2),ntryh(3),ntryh(4)/4,2,3,5/

      nl = n
      nf = 0
      j = 0
   10 j = j + 1
      IF ( j <=4 ) THEN
        ntry = ntryh(j)
      ELSE
        ntry = ntry + 2
      ENDIF
   40 nq = nl/ntry
      nr = nl - ntry*nq
      IF ( nr /= 0 ) THEN
        GOTO 10
      ENDIF
      nf = nf + 1
      fac(nf+2) = ntry
      nl = nq
      IF (ntry.NE.2) GO TO 70
      IF (nf.EQ.1) GO TO 70
      DO 60 i = 2,nf
         ib = nf - i + 2
         fac(ib+2) = fac(ib+1)
   60 CONTINUE
      fac(3) = 2
   70 IF (nl.NE.1) GO TO 40
      fac(1) = n
      fac(2) = nf
      tpi = 2.*pimach()
      argh = tpi/real(n)
      is = 0
      nfm1 = nf - 1
      l1 = 1
      IF (nfm1.EQ.0) RETURN
      DO 100 k1 = 1,nfm1
         ip = fac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         ipm = ip - 1
         DO 90 j = 1,ipm
            ld = ld + l1
            i = is
            argld = real(ld)*argh
            fi = 0.
            DO 80 ii = 3,ido,2
               i = i + 2
               fi = fi + 1.
               arg = fi*argld
               wa(i-1) = cos(arg)
               wa(i) = sin(arg)
   80       CONTINUE
            is = is + ido
   90    CONTINUE
         l1 = l2
  100 CONTINUE
      RETURN
      END subroutine
      SUBROUTINE vrfftf(m,n,r,rt,mdimr,wsave)
C
C***BEGIN PROLOGUE  VRFFTF
C***DATE WRITTEN   850801   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM,
C             FOURIER ANALYSIS, FORWARD TRANSFORM, MULTIPLE SEQUENCES
C***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST)
C***PURPOSE  Forward real periodic transform, M sequences.
C***DESCRIPTION
C
C  Subroutine VRFFTF computes the Fourier coefficients (forward
C  transform) of a number of real periodic sequences.  Specifically,
C  for each sequence the subroutine claculates the independent
C  Fourier coefficients described below at output parameter R.
C
C  The array WSAVE which is used by subroutine VRFFTF must be
C  initialized by calling subroutine VRFFTI(N,WSAVE).
C
C
C  Input Parameters
C
C  M       the number of sequences to be transformed.
C
C  N       the length of the sequences to be transformed.  The method
C          is most efficient when N is a product of small primes,
C          however n may be any positive integer.
C
C  R       areal two-dimensional array of size MDIMX x N containing the
C          the sequences to be transformed.  The sequences are stored
C          in the ROWS of R.  Thus, the I-th sequence to be transformed,
C          X(I,J), J=0,1,...,N-1, is stored as
C
C               R(I,J) = X(I,J-1) , J=1, 2, . . . , N.
C
C  RT      a real two-dimensional work array of size MDIMX x N.
C
C  MDIMR   the row (or first) dimension of the arrays R and RT exactly
C          as they appear in the calling program.  This parameter is
C          used to specify the variable dimension of these arrays.
C
C  WSAVE   a real one-dimensional work array which must be dimensioned
C          at least N+15.  The WSAVE array must be initialized by
C          calling subroutine VRFFTI.  A different WSAVE array must be
C          used for each different value of N.  This initialization does
C          not have to be repeated so long as N remains unchanged.  The
C          same WSAVE array may be used by VRFFTF and VRFFTB.
C
C  Output Parameters
C
C  R       contains the Fourier coefficients F(K) for each of the M
C          input sequences.  Specifically, row I of R, R(I,J),
C          J=1,2,..,N, contains the independent Fourier coefficients
C          F(I,K), for the I-th input sequence stored as
C
C             R(I,1) = REAL( F(I,0) ),
C                    = SQRT(1/N)*SUM(J=0,N-1)[ X(I,J) ],
C
C             R(I,2*K) = REAL( F(I,K) )
C                      = SQRT(1/N)*SUM(J=0,N-1)[X(I,J)*COS(2J*K*PI/N)]
C
C             R(I,2*K+1) = IMAG( F(I,K) )
C                        =-SQRT(1/N)*SUM(J=0,N-1)[X(I,J)*SIN(2J*K*PI/N)]
C
C                   for K = 1, 2, . . . , M-1,
C
C              and, when N is even,
C
C              R(I,N) = REAL( F(I,N/2) ).
C                     = SQRT(1/N)*SUM(J=0,N-1)[ (-1)**J*X(I,J) ].
C
C  WSAVE   contains results which must not be destroyed between calls
C          to VRFFTF or VRFFTB.
C
C  -----------------------------------------------------------------
C
C  NOTE  -  A call of VRFFTF followed immediately by a call of
C           of VRFFTB will return the original sequences R.  Thus,
C           VRFFTB is the correctly normalized inverse of VRFFTF.
C
C  -----------------------------------------------------------------
C
C  VRFFTF is a straightforward extension of the subprogram RFFTF to
C  handle M simultaneous sequences.  RFFTF was originally developed
C  by P. N. Swarztrauber of NCAR.
C
C
C              * * * * * * * * * * * * * * * * * * * * *
C              *                                       *
C              *         PROGRAM SPECIFICATIONS        *
C              *                                       *
C              * * * * * * * * * * * * * * * * * * * * *
C
C
C     Dimension of    R(MDIMR,N), RT(MDIMR,N), WSAVE(N+15)
C     Arguments
C
C     Latest          AUGUST 1, 1985
C     Revision
C
C     Subprograms     VRFFTI, VRFTI1, VRFFTF, VRFTF1, VRADF2, VRADF3,
C     Required        VRADF4, VRADF5, VRADFG, VRFFTB, VRFTB1, VRADB2,
C                     VRADB3, VRADB4, VRADB5, VRADBG, PIMACH
C
C     Special         NONE
C     Conditions
C
C     Common          NONE
C     blocks
C
C     I/O             NONE
C
C     Precision       SINGLE
C
C     Specialist      ROLAND SWEET
C
C     Language        FORTRAN
C
C     History         WRITTEN BY LINDA LINDGREN AND ROLAND SWEET AT THE
C                     NATIONAL BUREAU OF STANDARDS (BOULDER).
C
C     Algorithm       A REAL VARIANT OF THE STOCKHAM AUTOSORT VERSION
C                     OF THE COOLEY-TUKEY FAST FOURIER TRANSFORM.
C
C     Portability     AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN 77.
C                     THE ONLY MACHINE DEPENDENT CONSTANT IS LOCATED IN
C                     THE FUNCTION PIMACH.
C
C     Required        COS,SIN
C     resident
C     Routines
C
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
C               pp. 51-83.
C***ROUTINES CALLED  VRFTF1
C***END PROLOGUE  VRFFTF
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION r(mdimr,n),rt(mdimr,n),wsave(n+15)
C***FIRST EXECUTABLE STATEMENT  VRFFTF
      IF (n.EQ.1) RETURN
      CALL vrftf1(m,n,r,rt,mdimr,wsave(1),wsave(n+1))
      RETURN
      END subroutine
      SUBROUTINE vradf2(mp,ido,l1,cc,ch,mdimc,wa1)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION ch(mdimc,ido,2,l1),cc(mdimc,ido,l1,2),wa1(ido)

      DO 20 k = 1,l1
         DO 10 m = 1,mp
            ch(m,1,1,k) = cc(m,1,k,1) + cc(m,1,k,2)
            ch(m,ido,2,k) = cc(m,1,k,1) - cc(m,1,k,2)
   10    CONTINUE
   20 CONTINUE
      IF ( ido ==  2 ) THEN
         GOTO 70
      ELSEIF ( ido < 2 ) THEN
         GOTO 100
      ENDIF
      idp2 = ido + 2
      DO 60 k = 1,l1
         DO 50 i = 3,ido,2
            ic = idp2 - i
            DO 40 m = 1,mp
               ch(m,i,1,k) = cc(m,i,k,1) +
     +                       (wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*
     +                       cc(m,i-1,k,2))
               ch(m,ic,2,k) = (wa1(i-2)*cc(m,i,k,2)-
     +                        wa1(i-1)*cc(m,i-1,k,2)) - cc(m,i,k,1)
               ch(m,i-1,1,k) = cc(m,i-1,k,1) +
     +                         (wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*
     +                         cc(m,i,k,2))
               ch(m,ic-1,2,k) = cc(m,i-1,k,1) -
     +                          (wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*
     +                          cc(m,i,k,2))
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
      IF (mod(ido,2).EQ.1) RETURN
   70 DO 90 k = 1,l1
         DO 80 m = 1,mp
            ch(m,1,2,k) = -cc(m,ido,k,2)
            ch(m,ido,1,k) = cc(m,ido,k,1)
   80    CONTINUE
   90 CONTINUE
  100 RETURN
      END subroutine
      SUBROUTINE vradf3(mp,ido,l1,cc,ch,mdimc,wa1,wa2)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      USE m_constants, ONLY : pimach
      DIMENSION ch(mdimc,ido,3,l1),cc(mdimc,ido,l1,3),wa1(ido),wa2(ido)

      arg = 2.*pimach()/3.
      taur = cos(arg)
      taui = sin(arg)
      DO 20 k = 1,l1
         DO 10 m = 1,mp
            ch(m,1,1,k) = cc(m,1,k,1) + (cc(m,1,k,2)+cc(m,1,k,3))
            ch(m,1,3,k) = taui* (cc(m,1,k,3)-cc(m,1,k,2))
            ch(m,ido,2,k) = cc(m,1,k,1) +
     +                      taur* (cc(m,1,k,2)+cc(m,1,k,3))
   10    CONTINUE
   20 CONTINUE
      IF (ido.EQ.1) RETURN
      idp2 = ido + 2
      DO 50 k = 1,l1
         DO 40 i = 3,ido,2
            ic = idp2 - i
            DO 30 m = 1,mp
               ch(m,i-1,1,k) = cc(m,i-1,k,1) +
     +                         ((wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i,
     +                         k,2))+ (wa2(i-2)*cc(m,i-1,k,
     +                         3)+wa2(i-1)*cc(m,i,k,3)))
               ch(m,i,1,k) = cc(m,i,k,1) +
     +                       ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,
     +                       2))+ (wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*cc(m,
     +                       i-1,k,3)))
               ch(m,i-1,3,k) = (cc(m,i-1,k,1)+
     +                         taur* ((wa1(i-2)*cc(m,i-1,k,
     +                         2)+wa1(i-1)*cc(m,i,k,2))+ (wa2(i-2)*cc(m,
     +                         i-1,k,3)+wa2(i-1)*cc(m,i,k,3)))) +
     +                         (taui* ((wa1(i-2)*cc(m,i,k,
     +                         2)-wa1(i-1)*cc(m,i-1,k,
     +                         2))- (wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*cc(m,
     +                         i-1,k,3))))
               ch(m,ic-1,2,k) = (cc(m,i-1,k,1)+
     +                          taur* ((wa1(i-2)*cc(m,i-1,k,
     +                          2)+wa1(i-1)*cc(m,i,k,
     +                          2))+ (wa2(i-2)*cc(m,i-1,k,
     +                          3)+wa2(i-1)*cc(m,i,k,3)))) -
     +                          (taui* ((wa1(i-2)*cc(m,i,k,
     +                          2)-wa1(i-1)*cc(m,i-1,k,
     +                          2))- (wa2(i-2)*cc(m,i,k,
     +                          3)-wa2(i-1)*cc(m,i-1,k,3))))
               ch(m,i,3,k) = (cc(m,i,k,1)+taur*
     +                       ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,
     +                       2))+ (wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*cc(m,
     +                       i-1,k,3)))) + (taui* ((wa2(i-2)*cc(m,i-1,k,
     +                       3)+wa2(i-1)*cc(m,i,k,3))- (wa1(i-2)*cc(m,
     +                       i-1,k,2)+wa1(i-1)*cc(m,i,k,2))))
               ch(m,ic,2,k) = (taui* ((wa2(i-2)*cc(m,i-1,k,
     +                        3)+wa2(i-1)*cc(m,i,k,3))- (wa1(i-2)*cc(m,
     +                        i-1,k,2)+wa1(i-1)*cc(m,i,k,2)))) -
     +                        (cc(m,i,k,1)+taur* ((wa1(i-2)*cc(m,i,k,
     +                        2)-wa1(i-1)*cc(m,i-1,k,
     +                        2))+ (wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*cc(m,
     +                        i-1,k,3))))
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      RETURN
      END subroutine
      SUBROUTINE vradf4(mp,ido,l1,cc,ch,mdimc,wa1,wa2,wa3)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION cc(mdimc,ido,l1,4),ch(mdimc,ido,4,l1),wa1(ido),wa2(ido),
     +          wa3(ido)

      hsqt2 = sqrt(2.)/2.
      DO 20 k = 1,l1
         DO 10 m = 1,mp
            ch(m,1,1,k) = (cc(m,1,k,2)+cc(m,1,k,4)) +
     +                    (cc(m,1,k,1)+cc(m,1,k,3))
            ch(m,ido,4,k) = (cc(m,1,k,1)+cc(m,1,k,3)) -
     +                      (cc(m,1,k,2)+cc(m,1,k,4))
            ch(m,ido,2,k) = cc(m,1,k,1) - cc(m,1,k,3)
            ch(m,1,3,k) = cc(m,1,k,4) - cc(m,1,k,2)
   10    CONTINUE
   20 CONTINUE
      IF ( ido ==  2 ) THEN
         GOTO 70
      ELSEIF ( ido < 2 ) THEN
         GOTO 100
      ENDIF
      idp2 = ido + 2
      DO 60 k = 1,l1
         DO 50 i = 3,ido,2
            ic = idp2 - i
            DO 40 m = 1,mp
               ch(m,i-1,1,k) = ((wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i,
     +                         k,2))+ (wa3(i-2)*cc(m,i-1,k,
     +                         4)+wa3(i-1)*cc(m,i,k,4))) +
     +                         (cc(m,i-1,k,1)+ (wa2(i-2)*cc(m,i-1,k,
     +                         3)+wa2(i-1)*cc(m,i,k,3)))
               ch(m,ic-1,4,k) = (cc(m,i-1,k,1)+
     +                          (wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i,
     +                          k,3))) - ((wa1(i-2)*cc(m,i-1,k,
     +                          2)+wa1(i-1)*cc(m,i,k,2))+
     +                          (wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i,
     +                          k,4)))
               ch(m,i,1,k) = ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,
     +                       2))+ (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m,
     +                       i-1,k,4))) + (cc(m,i,k,1)+
     +                       (wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*cc(m,i-1,k,
     +                       3)))
               ch(m,ic,4,k) = ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,
     +                        k,2))+ (wa3(i-2)*cc(m,i,k,
     +                        4)-wa3(i-1)*cc(m,i-1,k,4))) -
     +                        (cc(m,i,k,1)+ (wa2(i-2)*cc(m,i,k,
     +                        3)-wa2(i-1)*cc(m,i-1,k,3)))
               ch(m,i-1,3,k) = ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,
     +                         k,2))- (wa3(i-2)*cc(m,i,k,
     +                         4)-wa3(i-1)*cc(m,i-1,k,4))) +
     +                         (cc(m,i-1,k,1)- (wa2(i-2)*cc(m,i-1,k,
     +                         3)+wa2(i-1)*cc(m,i,k,3)))
               ch(m,ic-1,2,k) = (cc(m,i-1,k,1)-
     +                          (wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i,
     +                          k,3))) - ((wa1(i-2)*cc(m,i,k,
     +                          2)-wa1(i-1)*cc(m,i-1,k,2))-
     +                          (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m,i-1,
     +                          k,4)))
               ch(m,i,3,k) = ((wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i,k,
     +                       4))- (wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,
     +                       i,k,2))) + (cc(m,i,k,1)-
     +                       (wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*cc(m,i-1,k,
     +                       3)))
               ch(m,ic,2,k) = ((wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i,
     +                        k,4))- (wa1(i-2)*cc(m,i-1,k,
     +                        2)+wa1(i-1)*cc(m,i,k,2))) -
     +                        (cc(m,i,k,1)- (wa2(i-2)*cc(m,i,k,
     +                        3)-wa2(i-1)*cc(m,i-1,k,3)))
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
      IF (mod(ido,2).EQ.1) RETURN
   70 CONTINUE
      DO 90 k = 1,l1
         DO 80 m = 1,mp
            ch(m,ido,1,k) = (hsqt2* (cc(m,ido,k,2)-cc(m,ido,k,4))) +
     +                      cc(m,ido,k,1)
            ch(m,ido,3,k) = cc(m,ido,k,1) -
     +                      (hsqt2* (cc(m,ido,k,2)-cc(m,ido,k,4)))
            ch(m,1,2,k) = (-hsqt2* (cc(m,ido,k,2)+cc(m,ido,k,4))) -
     +                    cc(m,ido,k,3)
            ch(m,1,4,k) = (-hsqt2* (cc(m,ido,k,2)+cc(m,ido,k,4))) +
     +                    cc(m,ido,k,3)
   80    CONTINUE
   90 CONTINUE
  100 RETURN
      END subroutine
      SUBROUTINE vradf5(mp,ido,l1,cc,ch,mdimc,wa1,wa2,wa3,wa4)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      USE m_constants, ONLY : pimach
      DIMENSION cc(mdimc,ido,l1,5),ch(mdimc,ido,5,l1),wa1(ido),wa2(ido),
     +          wa3(ido),wa4(ido)

      arg = 2.*pimach()/5.
      tr11 = cos(arg)
      ti11 = sin(arg)
      tr12 = cos(2.*arg)
      ti12 = sin(2.*arg)
      DO 20 k = 1,l1
         DO 10 m = 1,mp
            ch(m,1,1,k) = cc(m,1,k,1) + (cc(m,1,k,5)+cc(m,1,k,2)) +
     +                    (cc(m,1,k,4)+cc(m,1,k,3))
            ch(m,ido,2,k) = cc(m,1,k,1) +
     +                      tr11* (cc(m,1,k,5)+cc(m,1,k,2)) +
     +                      tr12* (cc(m,1,k,4)+cc(m,1,k,3))
            ch(m,1,3,k) = ti11* (cc(m,1,k,5)-cc(m,1,k,2)) +
     +                    ti12* (cc(m,1,k,4)-cc(m,1,k,3))
            ch(m,ido,4,k) = cc(m,1,k,1) +
     +                      tr12* (cc(m,1,k,5)+cc(m,1,k,2)) +
     +                      tr11* (cc(m,1,k,4)+cc(m,1,k,3))
            ch(m,1,5,k) = ti12* (cc(m,1,k,5)-cc(m,1,k,2)) -
     +                    ti11* (cc(m,1,k,4)-cc(m,1,k,3))
   10    CONTINUE
   20 CONTINUE
      IF (ido.EQ.1) RETURN
      idp2 = ido + 2
      DO 50 k = 1,l1
         DO 40 i = 3,ido,2
            ic = idp2 - i
            DO 30 m = 1,mp
               ch(m,i-1,1,k) = cc(m,i-1,k,1) +
     +                         ((wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i,
     +                         k,2))+ (wa4(i-2)*cc(m,i-1,k,
     +                         5)+wa4(i-1)*cc(m,i,k,5))) +
     +                         ((wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i,
     +                         k,3))+ (wa3(i-2)*cc(m,i-1,k,
     +                         4)+wa3(i-1)*cc(m,i,k,4)))
               ch(m,i,1,k) = cc(m,i,k,1) +
     +                       ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,
     +                       2))+ (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m,
     +                       i-1,k,5))) + ((wa2(i-2)*cc(m,i,k,
     +                       3)-wa2(i-1)*cc(m,i-1,k,3))+
     +                       (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m,i-1,k,
     +                       4)))
               ch(m,i-1,3,k) = cc(m,i-1,k,1) +
     +                         tr11* (wa1(i-2)*cc(m,i-1,k,2)+
     +                         wa1(i-1)*cc(m,i,k,2)+
     +                         wa4(i-2)*cc(m,i-1,k,5)+
     +                         wa4(i-1)*cc(m,i,k,5)) +
     +                         tr12* (wa2(i-2)*cc(m,i-1,k,3)+
     +                         wa2(i-1)*cc(m,i,k,3)+
     +                         wa3(i-2)*cc(m,i-1,k,4)+
     +                         wa3(i-1)*cc(m,i,k,4)) +
     +                         ti11* (wa1(i-2)*cc(m,i,k,2)-
     +                         wa1(i-1)*cc(m,i-1,k,2)-
     +                         (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m,i-1,
     +                         k,5))) + ti12* (wa2(i-2)*cc(m,i,k,3)-
     +                         wa2(i-1)*cc(m,i-1,k,3)-
     +                         (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m,i-1,
     +                         k,4)))
               ch(m,ic-1,2,k) = cc(m,i-1,k,1) +
     +                          tr11* (wa1(i-2)*cc(m,i-1,k,2)+
     +                          wa1(i-1)*cc(m,i,k,2)+
     +                          wa4(i-2)*cc(m,i-1,k,5)+
     +                          wa4(i-1)*cc(m,i,k,5)) +
     +                          tr12* (wa2(i-2)*cc(m,i-1,k,3)+
     +                          wa2(i-1)*cc(m,i,k,3)+
     +                          wa3(i-2)*cc(m,i-1,k,4)+
     +                          wa3(i-1)*cc(m,i,k,4)) -
     +                          (ti11* (wa1(i-2)*cc(m,i,k,
     +                          2)-wa1(i-1)*cc(m,i-1,k,
     +                          2)- (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m,
     +                          i-1,k,5)))+ti12* (wa2(i-2)*cc(m,i,k,
     +                          3)-wa2(i-1)*cc(m,i-1,k,
     +                          3)- (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m,
     +                          i-1,k,4))))
               ch(m,i,3,k) = (cc(m,i,k,1)+tr11*
     +                       ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,
     +                       2))+ (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m,
     +                       i-1,k,5)))+tr12* ((wa2(i-2)*cc(m,i,k,
     +                       3)-wa2(i-1)*cc(m,i-1,k,3))+ (wa3(i-2)*cc(m,
     +                       i,k,4)-wa3(i-1)*cc(m,i-1,k,4)))) +
     +                       (ti11* ((wa4(i-2)*cc(m,i-1,k,
     +                       5)+wa4(i-1)*cc(m,i,k,5))- (wa1(i-2)*cc(m,
     +                       i-1,k,2)+wa1(i-1)*cc(m,i,k,2)))+
     +                       ti12* ((wa3(i-2)*cc(m,i-1,k,
     +                       4)+wa3(i-1)*cc(m,i,k,4))- (wa2(i-2)*cc(m,
     +                       i-1,k,3)+wa2(i-1)*cc(m,i,k,3))))
               ch(m,ic,2,k) = (ti11* ((wa4(i-2)*cc(m,i-1,k,
     +                        5)+wa4(i-1)*cc(m,i,k,5))- (wa1(i-2)*cc(m,
     +                        i-1,k,2)+wa1(i-1)*cc(m,i,k,2)))+
     +                        ti12* ((wa3(i-2)*cc(m,i-1,k,
     +                        4)+wa3(i-1)*cc(m,i,k,4))- (wa2(i-2)*cc(m,
     +                        i-1,k,3)+wa2(i-1)*cc(m,i,k,3)))) -
     +                        (cc(m,i,k,1)+tr11* ((wa1(i-2)*cc(m,i,k,
     +                        2)-wa1(i-1)*cc(m,i-1,k,
     +                        2))+ (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m,
     +                        i-1,k,5)))+tr12* ((wa2(i-2)*cc(m,i,k,
     +                        3)-wa2(i-1)*cc(m,i-1,k,
     +                        3))+ (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m,
     +                        i-1,k,4))))
               ch(m,i-1,5,k) = (cc(m,i-1,k,1)+
     +                         tr12* ((wa1(i-2)*cc(m,i-1,k,
     +                         2)+wa1(i-1)*cc(m,i,k,2))+ (wa4(i-2)*cc(m,
     +                         i-1,k,5)+wa4(i-1)*cc(m,i,k,5)))+
     +                         tr11* ((wa2(i-2)*cc(m,i-1,k,
     +                         3)+wa2(i-1)*cc(m,i,k,3))+ (wa3(i-2)*cc(m,
     +                         i-1,k,4)+wa3(i-1)*cc(m,i,k,4)))) +
     +                         (ti12* ((wa1(i-2)*cc(m,i,k,
     +                         2)-wa1(i-1)*cc(m,i-1,k,
     +                         2))- (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m,
     +                         i-1,k,5)))-ti11* ((wa2(i-2)*cc(m,i,k,
     +                         3)-wa2(i-1)*cc(m,i-1,k,
     +                         3))- (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m,
     +                         i-1,k,4))))
               ch(m,ic-1,4,k) = (cc(m,i-1,k,1)+
     +                          tr12* ((wa1(i-2)*cc(m,i-1,k,
     +                          2)+wa1(i-1)*cc(m,i,k,
     +                          2))+ (wa4(i-2)*cc(m,i-1,k,
     +                          5)+wa4(i-1)*cc(m,i,k,5)))+
     +                          tr11* ((wa2(i-2)*cc(m,i-1,k,
     +                          3)+wa2(i-1)*cc(m,i,k,
     +                          3))+ (wa3(i-2)*cc(m,i-1,k,
     +                          4)+wa3(i-1)*cc(m,i,k,4)))) -
     +                          (ti12* ((wa1(i-2)*cc(m,i,k,
     +                          2)-wa1(i-1)*cc(m,i-1,k,
     +                          2))- (wa4(i-2)*cc(m,i,k,
     +                          5)-wa4(i-1)*cc(m,i-1,k,5)))-
     +                          ti11* ((wa2(i-2)*cc(m,i,k,
     +                          3)-wa2(i-1)*cc(m,i-1,k,
     +                          3))- (wa3(i-2)*cc(m,i,k,
     +                          4)-wa3(i-1)*cc(m,i-1,k,4))))
               ch(m,i,5,k) = (cc(m,i,k,1)+tr12*
     +                       ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,
     +                       2))+ (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m,
     +                       i-1,k,5)))+tr11* ((wa2(i-2)*cc(m,i,k,
     +                       3)-wa2(i-1)*cc(m,i-1,k,3))+ (wa3(i-2)*cc(m,
     +                       i,k,4)-wa3(i-1)*cc(m,i-1,k,4)))) +
     +                       (ti12* ((wa4(i-2)*cc(m,i-1,k,
     +                       5)+wa4(i-1)*cc(m,i,k,5))- (wa1(i-2)*cc(m,
     +                       i-1,k,2)+wa1(i-1)*cc(m,i,k,2)))-
     +                       ti11* ((wa3(i-2)*cc(m,i-1,k,
     +                       4)+wa3(i-1)*cc(m,i,k,4))- (wa2(i-2)*cc(m,
     +                       i-1,k,3)+wa2(i-1)*cc(m,i,k,3))))
               ch(m,ic,4,k) = (ti12* ((wa4(i-2)*cc(m,i-1,k,
     +                        5)+wa4(i-1)*cc(m,i,k,5))- (wa1(i-2)*cc(m,
     +                        i-1,k,2)+wa1(i-1)*cc(m,i,k,2)))-
     +                        ti11* ((wa3(i-2)*cc(m,i-1,k,
     +                        4)+wa3(i-1)*cc(m,i,k,4))- (wa2(i-2)*cc(m,
     +                        i-1,k,3)+wa2(i-1)*cc(m,i,k,3)))) -
     +                        (cc(m,i,k,1)+tr12* ((wa1(i-2)*cc(m,i,k,
     +                        2)-wa1(i-1)*cc(m,i-1,k,
     +                        2))+ (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m,
     +                        i-1,k,5)))+tr11* ((wa2(i-2)*cc(m,i,k,
     +                        3)-wa2(i-1)*cc(m,i-1,k,
     +                        3))+ (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m,
     +                        i-1,k,4))))
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      RETURN
      END subroutine
      SUBROUTINE vradfg(mp,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,mdimc,wa)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      USE m_constants, ONLY : pimach
      DIMENSION ch(mdimc,ido,l1,ip),cc(mdimc,ido,ip,l1),
     +          c1(mdimc,ido,l1,ip),c2(mdimc,idl1,ip),
     +          ch2(mdimc,idl1,ip),wa(ido)

      tpi = 2.*pimach()
      arg = tpi/real(ip)
      dcp = cos(arg)
      dsp = sin(arg)
      ipph = (ip+1)/2
      ipp2 = ip + 2
      idp2 = ido + 2
      nbd = (ido-1)/2
      IF (ido.EQ.1) GO TO 250
      DO 20 ik = 1,idl1
         DO 10 m = 1,mp
            ch2(m,ik,1) = c2(m,ik,1)
   10    CONTINUE
   20 CONTINUE
      DO 50 j = 2,ip
         DO 40 k = 1,l1
            DO 30 m = 1,mp
               ch(m,1,k,j) = c1(m,1,k,j)
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      IF (nbd.GT.l1) GO TO 100
      is = -ido
      DO 90 j = 2,ip
         is = is + ido
         idij = is
         DO 80 i = 3,ido,2
            idij = idij + 2
            DO 70 k = 1,l1
               DO 60 m = 1,mp
                  ch(m,i-1,k,j) = wa(idij-1)*c1(m,i-1,k,j) +
     +                            wa(idij)*c1(m,i,k,j)
                  ch(m,i,k,j) = wa(idij-1)*c1(m,i,k,j) -
     +                          wa(idij)*c1(m,i-1,k,j)
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
      GO TO 150
  100 is = -ido
      DO 140 j = 2,ip
         is = is + ido
         DO 130 k = 1,l1
            idij = is
            DO 120 i = 3,ido,2
               idij = idij + 2
               DO 110 m = 1,mp
                  ch(m,i-1,k,j) = wa(idij-1)*c1(m,i-1,k,j) +
     +                            wa(idij)*c1(m,i,k,j)
                  ch(m,i,k,j) = wa(idij-1)*c1(m,i,k,j) -
     +                          wa(idij)*c1(m,i-1,k,j)
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
  150 IF (nbd.LT.l1) GO TO 200
      DO 190 j = 2,ipph
         jc = ipp2 - j
         DO 180 k = 1,l1
            DO 170 i = 3,ido,2
               DO 160 m = 1,mp
                  c1(m,i-1,k,j) = ch(m,i-1,k,j) + ch(m,i-1,k,jc)
                  c1(m,i-1,k,jc) = ch(m,i,k,j) - ch(m,i,k,jc)
                  c1(m,i,k,j) = ch(m,i,k,j) + ch(m,i,k,jc)
                  c1(m,i,k,jc) = ch(m,i-1,k,jc) - ch(m,i-1,k,j)
  160          CONTINUE
  170       CONTINUE
  180    CONTINUE
  190 CONTINUE
      GO TO 280
  200 DO 240 j = 2,ipph
         jc = ipp2 - j
         DO 230 i = 3,ido,2
            DO 220 k = 1,l1
               DO 210 m = 1,mp
                  c1(m,i-1,k,j) = ch(m,i-1,k,j) + ch(m,i-1,k,jc)
                  c1(m,i-1,k,jc) = ch(m,i,k,j) - ch(m,i,k,jc)
                  c1(m,i,k,j) = ch(m,i,k,j) + ch(m,i,k,jc)
                  c1(m,i,k,jc) = ch(m,i-1,k,jc) - ch(m,i-1,k,j)
  210          CONTINUE
  220       CONTINUE
  230    CONTINUE
  240 CONTINUE
      GO TO 280
  250 DO 270 ik = 1,idl1
         DO 260 m = 1,mp
            c2(m,ik,1) = ch2(m,ik,1)
  260    CONTINUE
  270 CONTINUE
  280 DO 310 j = 2,ipph
         jc = ipp2 - j
         DO 300 k = 1,l1
            DO 290 m = 1,mp
               c1(m,1,k,j) = ch(m,1,k,j) + ch(m,1,k,jc)
               c1(m,1,k,jc) = ch(m,1,k,jc) - ch(m,1,k,j)
  290       CONTINUE
  300    CONTINUE
  310 CONTINUE
C
      ar1 = 1.
      ai1 = 0.
      DO 370 l = 2,ipph
         lc = ipp2 - l
         ar1h = dcp*ar1 - dsp*ai1
         ai1 = dcp*ai1 + dsp*ar1
         ar1 = ar1h
         DO 330 ik = 1,idl1
            DO 320 m = 1,mp
               ch2(m,ik,l) = c2(m,ik,1) + ar1*c2(m,ik,2)
               ch2(m,ik,lc) = ai1*c2(m,ik,ip)
  320       CONTINUE
  330    CONTINUE
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         DO 360 j = 3,ipph
            jc = ipp2 - j
            ar2h = dc2*ar2 - ds2*ai2
            ai2 = dc2*ai2 + ds2*ar2
            ar2 = ar2h
            DO 350 ik = 1,idl1
               DO 340 m = 1,mp
                  ch2(m,ik,l) = ch2(m,ik,l) + ar2*c2(m,ik,j)
                  ch2(m,ik,lc) = ch2(m,ik,lc) + ai2*c2(m,ik,jc)
  340          CONTINUE
  350       CONTINUE
  360    CONTINUE
  370 CONTINUE
      DO 400 j = 2,ipph
         DO 390 ik = 1,idl1
            DO 380 m = 1,mp
               ch2(m,ik,1) = ch2(m,ik,1) + c2(m,ik,j)
  380       CONTINUE
  390    CONTINUE
  400 CONTINUE
C
      IF (ido.LT.l1) GO TO 440
      DO 430 k = 1,l1
         DO 420 i = 1,ido
            DO 410 m = 1,mp
               cc(m,i,1,k) = ch(m,i,k,1)
  410       CONTINUE
  420    CONTINUE
  430 CONTINUE
      GO TO 480
  440 DO 470 i = 1,ido
         DO 460 k = 1,l1
            DO 450 m = 1,mp
               cc(m,i,1,k) = ch(m,i,k,1)
  450       CONTINUE
  460    CONTINUE
  470 CONTINUE
  480 DO 510 j = 2,ipph
         jc = ipp2 - j
         j2 = j + j
         DO 500 k = 1,l1
            DO 490 m = 1,mp
               cc(m,ido,j2-2,k) = ch(m,1,k,j)
               cc(m,1,j2-1,k) = ch(m,1,k,jc)
  490       CONTINUE
  500    CONTINUE
  510 CONTINUE
      IF (ido.EQ.1) RETURN
      IF (nbd.LT.l1) GO TO 560
      DO 550 j = 2,ipph
         jc = ipp2 - j
         j2 = j + j
         DO 540 k = 1,l1
            DO 530 i = 3,ido,2
               ic = idp2 - i
               DO 520 m = 1,mp
                  cc(m,i-1,j2-1,k) = ch(m,i-1,k,j) + ch(m,i-1,k,jc)
                  cc(m,ic-1,j2-2,k) = ch(m,i-1,k,j) - ch(m,i-1,k,jc)
                  cc(m,i,j2-1,k) = ch(m,i,k,j) + ch(m,i,k,jc)
                  cc(m,ic,j2-2,k) = ch(m,i,k,jc) - ch(m,i,k,j)
  520          CONTINUE
  530       CONTINUE
  540    CONTINUE
  550 CONTINUE
      RETURN
  560 DO 600 j = 2,ipph
         jc = ipp2 - j
         j2 = j + j
         DO 590 i = 3,ido,2
            ic = idp2 - i
            DO 580 k = 1,l1
               DO 570 m = 1,mp
                  cc(m,i-1,j2-1,k) = ch(m,i-1,k,j) + ch(m,i-1,k,jc)
                  cc(m,ic-1,j2-2,k) = ch(m,i-1,k,j) - ch(m,i-1,k,jc)
                  cc(m,i,j2-1,k) = ch(m,i,k,j) + ch(m,i,k,jc)
                  cc(m,ic,j2-2,k) = ch(m,i,k,jc) - ch(m,i,k,j)
  570          CONTINUE
  580       CONTINUE
  590    CONTINUE
  600 CONTINUE
      RETURN
      END subroutine
      SUBROUTINE vrftf1(m,n,c,ch,mdimc,wa,fac)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION ch(mdimc,n),c(mdimc,n),wa(n),fac(15)

      nf = fac(2)
      na = 1
      l2 = n
      iw = n
      DO 110 k1 = 1,nf
         kh = nf - k1
         ip = fac(kh+3)
         l1 = l2/ip
         ido = n/l2
         idl1 = ido*l1
         iw = iw - (ip-1)*ido
         na = 1 - na
         IF (ip.NE.4) GO TO 20
         ix2 = iw + ido
         ix3 = ix2 + ido
         IF (na.NE.0) GO TO 10
         CALL vradf4(m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2),wa(ix3))
         GO TO 100
   10    CALL vradf4(m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2),wa(ix3))
         GO TO 100
   20    IF (ip.NE.2) GO TO 40
         IF (na.NE.0) GO TO 30
         CALL vradf2(m,ido,l1,c,ch,mdimc,wa(iw))
         GO TO 100
   30    CALL vradf2(m,ido,l1,ch,c,mdimc,wa(iw))
         GO TO 100
   40    IF (ip.NE.3) GO TO 60
         ix2 = iw + ido
         IF (na.NE.0) GO TO 50
         CALL vradf3(m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2))
         GO TO 100
   50    CALL vradf3(m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2))
         GO TO 100
   60    IF (ip.NE.5) GO TO 80
         ix2 = iw + ido
         ix3 = ix2 + ido
         ix4 = ix3 + ido
         IF (na.NE.0) GO TO 70
         CALL vradf5(m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         GO TO 100
   70    CALL vradf5(m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         GO TO 100
   80    IF (ido.EQ.1) na = 1 - na
         IF (na.NE.0) GO TO 90
         CALL vradfg(m,ido,ip,l1,idl1,c,c,c,ch,ch,mdimc,wa(iw))
         na = 1
         GO TO 100
   90    CALL vradfg(m,ido,ip,l1,idl1,ch,ch,ch,c,c,mdimc,wa(iw))
         na = 0
  100    l2 = l1
  110 CONTINUE
      scale = sqrt(1./n)
      IF (na.EQ.1) GO TO 140
      DO 130 j = 1,n
         DO 120 i = 1,m
            c(i,j) = scale*ch(i,j)
  120    CONTINUE
  130 CONTINUE
      RETURN
  140 DO 160 j = 1,n
         DO 150 i = 1,m
            c(i,j) = scale*c(i,j)
  150    CONTINUE
  160 CONTINUE
      RETURN
      END subroutine

      SUBROUTINE vrfftb(m,n,r,rt,mdimr,wsave)
C***BEGIN PROLOGUE  VRFFTB
C***DATE WRITTEN   850801   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM,
C             FOURIER SYNTHESIS, BACKWARD TRANSFORM, MULTIPLE SEQUENCES
C***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST)
C***PURPOSE  Backward real periodic transform, M sequences.
C***DESCRIPTION
C
C  Subroutine VRFFTB computes the synthesis (backward transform) of a
C  number of real periodic sequences from their Fourier coefficients.
C  Specifically, for each set of independent Fourier coefficients
C  F(K), the corresponding real periodic sequence is computed.
C
C  The array WSAVE which is used by subroutine VRFFTB must be
C  initialized by calling subroutine VRFFTI(N,WSAVE).
C
C
C  Input Parameters
C
C  M       the number of sets of coefficients.
C
C  N       the length of the sequences of coefficients to be
C          transformed.  The method is most efficient when N is a
C          product of small primes, however n may be any positive
C          integer.
C
C  R       areal two-dimensional array of size MDIMX x N containing the
C          coefficients to be transformed.  Each set of coefficients
C          F(K), K\0,1,..,N-1, is stored as a ROW of R.  Specifically,
C          the I-th set of independent Fourier coefficients is stored
C
C                R(I,1) = REAL( F(I,0) ),
C
C                R(I,2*K) = REAL( F(I,K) )
C
C                R(I,2*K+1) = IMAG( F(I,K) )
C
C                   for K = 1, 2, . . . , M-1,
C
C                and, when N is even,
C
C                R(I,N) = REAL( F(I,N/2) ).
C
C  RT      a real two-dimensional work array of size MDIMX x N.
C
C  MDIMR   the row (or first) dimension of the arrays R and RT exactly
C          as they appear in the calling program.  This parameter is
C          used to specify the variable dimension of these arrays.
C
C  WSAVE   a real one-dimensional work array which must be dimensioned
C          at least N+15.  The WSAVE array must be initialized by
C          calling subroutine VRFFTI.  A different WSAVE array must be
C          used for each different value of N.  This initialization does
C          not have to be repeated so long as N remains unchanged.  The
C          same WSAVE array may be used by VRFFTB and VRFFTB.
C
C  Output Parameters
C
C  R       contains M real periodic sequences corresponding to the given
C          coefficients.  Specifically, the I-th row of R contains the
C          real periodic sequence corresponding to the I-th set of
C          independent Fourier coefficients F(I,K) stored as
C
C               R(I,J) = X(I,J-1) ,   J = 1, 2, . . . , N, where
C
C               X(I,J) = SQRT(1/N)* F(I,0) + (-1)**J*F(I,N/2)
C                        + 2*SUM(K=1,M)[ REAL(F(I,2K))*COS(2K*J*PI/N)
C                        - IMAG(F(I,2K+1))*SIN(2K*J*PI/N) ]  ,
C
C                 when N is even, and
C
C               X(I,J) = SQRT(1/N)* F(I,0) +
C                        2*SUM(K=1,M)[ REAL(F(I,2K))*COS(2K*J*PI/N)
C                        - IMAG(F(I,2K+1))*SIN(2K*J*PI/N) ]  ,
C
C                 when N is odd.
C
C  WSAVE   contains results which must not be destroyed between calls
C          to VRFFTF or VRFFTB.
C
C  -----------------------------------------------------------------
C
C  NOTE  -  A call of VRFFTF followed immediately by a call of
C           of VRFFTB will return the original sequences R.  Thus,
C           VRFFTB is the correctly normalized inverse of VRFFTF.
C
C  -----------------------------------------------------------------
C
C  VRFFTB is a straightforward extension of the subprogram RFFTB to
C  handle M simultaneous sequences.  RFFTB was originally developed
C  by P. N. Swarztrauber of NCAR.
C
C
C              * * * * * * * * * * * * * * * * * * * * *
C              *                                       *
C              *         PROGRAM SPECIFICATIONS        *
C              *                                       *
C              * * * * * * * * * * * * * * * * * * * * *
C
C
C     Dimension of    R(MDIMR,N), RT(MDIMR,N), WSAVE(N+15)
C     Arguments
C
C     Latest          AUGUST 1, 1985
C     Revision
C
C     Subprograms     VRFFTI, VRFTI1, VRFFTF, VRFTF1, VRADF2, VRADF3,
C     required        VRADF4, VRADF5, VRADFG, VRFFTB, VRFTB1, VRADB2,
C                     VRADB3, VRADB4, VRADB5, VRADBG, PIMACH
C
C     Special         NONE
C     Conditions
C
C     Common          NONE
C     blocks
C
C     I/O             NONE
C
C     Precision       SINGLE
C
C     Specialist      ROLAND SWEET
C
C     Language        FORTRAN
C
C     History         WRITTEN BY LINDA LINDGREN AND ROLAND SWEET AT THE
C                     NATIONAL BUREAU OF STANDARDS (BOULDER).
C
C     Algorithm       A REAL VARIANT OF THE STOCKHAM AUTOSORT VERSION
C                     OF THE COOLEY-TUKEY FAST FOURIER TRANSFORM.
C
C     Portability     AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN 77.
C                     THE ONLY MACHINE DEPENDENT CONSTANT IS LOCATED IN
C                     THE FUNCTION PIMACH.
C
C     Required        COS,SIN
C     resident
C     Routines
C
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
C               pp. 51-83.
C***ROUTINES CALLED  VRFTB1
C***END PROLOGUE  VRFFTB
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n, m, mdimr
      REAL :: r(mdimr,n),rt(mdimr,n),wsave(n+15)

      IF (n.EQ.1) RETURN
      CALL vrftb1(m,n,r,rt,mdimr,wsave(1),wsave(n+1))
      RETURN

      CONTAINS

      SUBROUTINE vradb2(mp,ido,l1,cc,ch,mdimc,wa1)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: mp, ido, l1, mdimc
      REAL, INTENT(IN)    :: cc(mdimc,ido,2,l1), wa1(ido)
      REAL, INTENT(INOUT) :: ch(mdimc,ido,l1,2)

      INTEGER :: i, m, k, ic, idp2

      DO 20 k = 1,l1
         DO 10 m = 1,mp
            ch(m,1,k,1) = cc(m,1,1,k) + cc(m,ido,2,k)
            ch(m,1,k,2) = cc(m,1,1,k) - cc(m,ido,2,k)
   10    CONTINUE
   20 CONTINUE
      IF ( ido ==  2 ) THEN
         GOTO 70
      ELSEIF ( ido < 2 ) THEN
         GOTO 100
      ENDIF
      idp2 = ido + 2
      DO 60 k = 1,l1
         DO 50 i = 3,ido,2
            ic = idp2 - i
            DO 40 m = 1,mp
               ch(m,i-1,k,1) = cc(m,i-1,1,k) + cc(m,ic-1,2,k)
               ch(m,i,k,1) = cc(m,i,1,k) - cc(m,ic,2,k)
               ch(m,i-1,k,2) = wa1(i-2)* (cc(m,i-1,1,k)-
     +                         cc(m,ic-1,2,k)) - wa1(i-1)*
     +                         (cc(m,i,1,k)+cc(m,ic,2,k))
               ch(m,i,k,2) = wa1(i-2)* (cc(m,i,1,k)+cc(m,ic,2,k)) +
     +                       wa1(i-1)* (cc(m,i-1,1,k)-cc(m,ic-1,2,k))
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
      IF (mod(ido,2).EQ.1) RETURN
   70 DO 90 k = 1,l1
         DO 80 m = 1,mp
            ch(m,ido,k,1) = cc(m,ido,1,k) + cc(m,ido,1,k)
            ch(m,ido,k,2) = - (cc(m,1,2,k)+cc(m,1,2,k))
   80    CONTINUE
   90 CONTINUE
  100 RETURN
      END SUBROUTINE vradb2

      SUBROUTINE vradb3(mp,ido,l1,cc,ch,mdimc,wa1,wa2)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      USE m_constants, ONLY : pimach
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: mp, ido, l1, mdimc
      REAL, INTENT(IN)    :: cc(mdimc,ido,3,l1), wa1(ido), wa2(ido)
      REAL, INTENT(INOUT) :: ch(mdimc,ido,l1,3)

      INTEGER :: i, k, idp2, ic, m
      REAL :: arg, taur, taui

      arg = 2.*pimach()/3.
      taur = cos(arg)
      taui = sin(arg)
      DO 20 k = 1,l1
         DO 10 m = 1,mp
            ch(m,1,k,1) = cc(m,1,1,k) + 2.*cc(m,ido,2,k)
            ch(m,1,k,2) = cc(m,1,1,k) + (2.*taur)*cc(m,ido,2,k) -
     +                    (2.*taui)*cc(m,1,3,k)
            ch(m,1,k,3) = cc(m,1,1,k) + (2.*taur)*cc(m,ido,2,k) +
     +                    2.*taui*cc(m,1,3,k)
   10    CONTINUE
   20 CONTINUE
      IF (ido.EQ.1) RETURN
      idp2 = ido + 2
      DO 50 k = 1,l1
         DO 40 i = 3,ido,2
            ic = idp2 - i
            DO 30 m = 1,mp
               ch(m,i-1,k,1) = cc(m,i-1,1,k) +
     +                         (cc(m,i-1,3,k)+cc(m,ic-1,2,k))
               ch(m,i,k,1) = cc(m,i,1,k) + (cc(m,i,3,k)-cc(m,ic,2,k))
               ch(m,i-1,k,2) = wa1(i-2)* ((cc(m,i-1,1,k)+taur* (cc(m,
     +                         i-1,3,k)+cc(m,ic-1,2,k)))-
     +                         (taui* (cc(m,i,3,k)+cc(m,ic,2,k)))) -
     +                         wa1(i-1)* ((cc(m,i,1,k)+taur* (cc(m,i,3,
     +                         k)-cc(m,ic,2,k)))+ (taui* (cc(m,i-1,3,
     +                         k)-cc(m,ic-1,2,k))))
               ch(m,i,k,2) = wa1(i-2)* ((cc(m,i,1,k)+taur* (cc(m,i,3,
     +                       k)-cc(m,ic,2,k)))+ (taui* (cc(m,i-1,3,
     +                       k)-cc(m,ic-1,2,k)))) +
     +                       wa1(i-1)* ((cc(m,i-1,1,k)+taur* (cc(m,i-1,
     +                       3,k)+cc(m,ic-1,2,k)))-
     +                       (taui* (cc(m,i,3,k)+cc(m,ic,2,k))))
               ch(m,i-1,k,3) = wa2(i-2)* ((cc(m,i-1,1,k)+taur* (cc(m,
     +                         i-1,3,k)+cc(m,ic-1,2,k)))+
     +                         (taui* (cc(m,i,3,k)+cc(m,ic,2,k)))) -
     +                         wa2(i-1)* ((cc(m,i,1,k)+taur* (cc(m,i,3,
     +                         k)-cc(m,ic,2,k)))- (taui* (cc(m,i-1,3,
     +                         k)-cc(m,ic-1,2,k))))
               ch(m,i,k,3) = wa2(i-2)* ((cc(m,i,1,k)+taur* (cc(m,i,3,
     +                       k)-cc(m,ic,2,k)))- (taui* (cc(m,i-1,3,
     +                       k)-cc(m,ic-1,2,k)))) +
     +                       wa2(i-1)* ((cc(m,i-1,1,k)+taur* (cc(m,i-1,
     +                       3,k)+cc(m,ic-1,2,k)))+
     +                       (taui* (cc(m,i,3,k)+cc(m,ic,2,k))))
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      RETURN
      END SUBROUTINE vradb3

      SUBROUTINE vradb4(mp,ido,l1,cc,ch,mdimc,wa1,wa2,wa3)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: mp, ido, l1, mdimc
      REAL, INTENT(IN)    :: cc(mdimc,ido,4,l1), wa1(ido), wa2(ido),
     +                       wa3(ido)
      REAL, INTENT(INOUT) :: ch(mdimc,ido,l1,4)

      INTEGER :: k, m, i, idp2, ic
      REAL :: sqrt2

      sqrt2 = sqrt(2.)
      DO 20 k = 1,l1
         DO 10 m = 1,mp
            ch(m,1,k,3) = (cc(m,1,1,k)+cc(m,ido,4,k)) -
     +                    (cc(m,ido,2,k)+cc(m,ido,2,k))
            ch(m,1,k,1) = (cc(m,1,1,k)+cc(m,ido,4,k)) +
     +                    (cc(m,ido,2,k)+cc(m,ido,2,k))
            ch(m,1,k,4) = (cc(m,1,1,k)-cc(m,ido,4,k)) +
     +                    (cc(m,1,3,k)+cc(m,1,3,k))
            ch(m,1,k,2) = (cc(m,1,1,k)-cc(m,ido,4,k)) -
     +                    (cc(m,1,3,k)+cc(m,1,3,k))
   10    CONTINUE
   20 CONTINUE
      IF ( ido ==  2 ) THEN
         GOTO 70
      ELSEIF ( ido < 2 ) THEN
         GOTO 100
      ENDIF
      idp2 = ido + 2
      DO 60 k = 1,l1
         DO 50 i = 3,ido,2
            ic = idp2 - i
            DO 40 m = 1,mp
               ch(m,i-1,k,1) = (cc(m,i-1,1,k)+cc(m,ic-1,4,k)) +
     +                         (cc(m,i-1,3,k)+cc(m,ic-1,2,k))
               ch(m,i,k,1) = (cc(m,i,1,k)-cc(m,ic,4,k)) +
     +                       (cc(m,i,3,k)-cc(m,ic,2,k))
               ch(m,i-1,k,2) = wa1(i-2)* ((cc(m,i-1,1,k)-cc(m,ic-1,4,
     +                         k))- (cc(m,i,3,k)+cc(m,ic,2,k))) -
     +                         wa1(i-1)* ((cc(m,i,1,k)+cc(m,ic,4,k))+
     +                         (cc(m,i-1,3,k)-cc(m,ic-1,2,k)))
               ch(m,i,k,2) = wa1(i-2)* ((cc(m,i,1,k)+cc(m,ic,4,k))+
     +                       (cc(m,i-1,3,k)-cc(m,ic-1,2,k))) +
     +                       wa1(i-1)* ((cc(m,i-1,1,k)-cc(m,ic-1,4,k))-
     +                       (cc(m,i,3,k)+cc(m,ic,2,k)))
               ch(m,i-1,k,3) = wa2(i-2)* ((cc(m,i-1,1,k)+cc(m,ic-1,4,
     +                         k))- (cc(m,i-1,3,k)+cc(m,ic-1,2,k))) -
     +                         wa2(i-1)* ((cc(m,i,1,k)-cc(m,ic,4,k))-
     +                         (cc(m,i,3,k)-cc(m,ic,2,k)))
               ch(m,i,k,3) = wa2(i-2)* ((cc(m,i,1,k)-cc(m,ic,4,k))-
     +                       (cc(m,i,3,k)-cc(m,ic,2,k))) +
     +                       wa2(i-1)* ((cc(m,i-1,1,k)+cc(m,ic-1,4,k))-
     +                       (cc(m,i-1,3,k)+cc(m,ic-1,2,k)))
               ch(m,i-1,k,4) = wa3(i-2)* ((cc(m,i-1,1,k)-cc(m,ic-1,4,
     +                         k))+ (cc(m,i,3,k)+cc(m,ic,2,k))) -
     +                         wa3(i-1)* ((cc(m,i,1,k)+cc(m,ic,4,k))-
     +                         (cc(m,i-1,3,k)-cc(m,ic-1,2,k)))
               ch(m,i,k,4) = wa3(i-2)* ((cc(m,i,1,k)+cc(m,ic,4,k))-
     +                       (cc(m,i-1,3,k)-cc(m,ic-1,2,k))) +
     +                       wa3(i-1)* ((cc(m,i-1,1,k)-cc(m,ic-1,4,k))+
     +                       (cc(m,i,3,k)+cc(m,ic,2,k)))
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
      IF (mod(ido,2).EQ.1) RETURN
   70 CONTINUE
      DO 90 k = 1,l1
         DO 80 m = 1,mp
            ch(m,ido,k,1) = (cc(m,ido,1,k)+cc(m,ido,3,k)) +
     +                      (cc(m,ido,1,k)+cc(m,ido,3,k))
            ch(m,ido,k,2) = sqrt2* ((cc(m,ido,1,k)-cc(m,ido,3,k))-
     +                      (cc(m,1,2,k)+cc(m,1,4,k)))
            ch(m,ido,k,3) = (cc(m,1,4,k)-cc(m,1,2,k)) +
     +                      (cc(m,1,4,k)-cc(m,1,2,k))
            ch(m,ido,k,4) = -sqrt2* ((cc(m,ido,1,k)-cc(m,ido,3,k))+
     +                      (cc(m,1,2,k)+cc(m,1,4,k)))
   80    CONTINUE
   90 CONTINUE
  100 RETURN
      END SUBROUTINE vradb4

      SUBROUTINE vradb5(mp,ido,l1,cc,ch,mdimc,wa1,wa2,wa3,wa4)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      USE m_constants, ONLY : pimach
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: mp, ido, l1, mdimc
      REAL, INTENT(IN)    :: cc(mdimc,ido,5,l1), wa1(ido), wa2(ido),
     +                       wa3(ido),wa4(ido)
      REAL, INTENT(INOUT) :: ch(mdimc,ido,l1,5)

      INTEGER :: i, k, m, ic, idp2
      REAL :: arg, tr11, ti11, tr12, ti12

      arg = 2.*pimach()/5.
      tr11 = cos(arg)
      ti11 = sin(arg)
      tr12 = cos(2.*arg)
      ti12 = sin(2.*arg)
      DO 20 k = 1,l1
         DO 10 m = 1,mp
            ch(m,1,k,1) = cc(m,1,1,k) + 2.*cc(m,ido,2,k) +
     +                    2.*cc(m,ido,4,k)
            ch(m,1,k,2) = (cc(m,1,1,k)+tr11*2.*cc(m,ido,2,k)+
     +                    tr12*2.*cc(m,ido,4,k)) -
     +                    (ti11*2.*cc(m,1,3,k)+ti12*2.*cc(m,1,5,k))
            ch(m,1,k,3) = (cc(m,1,1,k)+tr12*2.*cc(m,ido,2,k)+
     +                    tr11*2.*cc(m,ido,4,k)) -
     +                    (ti12*2.*cc(m,1,3,k)-ti11*2.*cc(m,1,5,k))
            ch(m,1,k,4) = (cc(m,1,1,k)+tr12*2.*cc(m,ido,2,k)+
     +                    tr11*2.*cc(m,ido,4,k)) +
     +                    (ti12*2.*cc(m,1,3,k)-ti11*2.*cc(m,1,5,k))
            ch(m,1,k,5) = (cc(m,1,1,k)+tr11*2.*cc(m,ido,2,k)+
     +                    tr12*2.*cc(m,ido,4,k)) +
     +                    (ti11*2.*cc(m,1,3,k)+ti12*2.*cc(m,1,5,k))
   10    CONTINUE
   20 CONTINUE
      IF (ido.EQ.1) RETURN
      idp2 = ido + 2
      DO 50 k = 1,l1
         DO 40 i = 3,ido,2
            ic = idp2 - i
            DO 30 m = 1,mp
               ch(m,i-1,k,1) = cc(m,i-1,1,k) +
     +                         (cc(m,i-1,3,k)+cc(m,ic-1,2,k)) +
     +                         (cc(m,i-1,5,k)+cc(m,ic-1,4,k))
               ch(m,i,k,1) = cc(m,i,1,k) + (cc(m,i,3,k)-cc(m,ic,2,k)) +
     +                       (cc(m,i,5,k)-cc(m,ic,4,k))
               ch(m,i-1,k,2) = wa1(i-2)* ((cc(m,i-1,1,k)+tr11* (cc(m,
     +                         i-1,3,k)+cc(m,ic-1,2,k))+tr12* (cc(m,i-1,
     +                         5,k)+cc(m,ic-1,4,k)))-
     +                         (ti11* (cc(m,i,3,k)+cc(m,ic,2,
     +                         k))+ti12* (cc(m,i,5,k)+cc(m,ic,4,k)))) -
     +                         wa1(i-1)* ((cc(m,i,1,k)+tr11* (cc(m,i,3,
     +                         k)-cc(m,ic,2,k))+tr12* (cc(m,i,5,k)-cc(m,
     +                         ic,4,k)))+ (ti11* (cc(m,i-1,3,k)-cc(m,
     +                         ic-1,2,k))+ti12* (cc(m,i-1,5,k)-cc(m,
     +                         ic-1,4,k))))
               ch(m,i,k,2) = wa1(i-2)* ((cc(m,i,1,k)+tr11* (cc(m,i,3,
     +                       k)-cc(m,ic,2,k))+tr12* (cc(m,i,5,k)-cc(m,
     +                       ic,4,k)))+ (ti11* (cc(m,i-1,3,k)-cc(m,ic-1,
     +                       2,k))+ti12* (cc(m,i-1,5,k)-cc(m,ic-1,4,
     +                       k)))) + wa1(i-1)* ((cc(m,i-1,1,
     +                       k)+tr11* (cc(m,i-1,3,k)+cc(m,ic-1,2,
     +                       k))+tr12* (cc(m,i-1,5,k)+cc(m,ic-1,4,k)))-
     +                       (ti11* (cc(m,i,3,k)+cc(m,ic,2,
     +                       k))+ti12* (cc(m,i,5,k)+cc(m,ic,4,k))))
               ch(m,i-1,k,3) = wa2(i-2)* ((cc(m,i-1,1,k)+tr12* (cc(m,
     +                         i-1,3,k)+cc(m,ic-1,2,k))+tr11* (cc(m,i-1,
     +                         5,k)+cc(m,ic-1,4,k)))-
     +                         (ti12* (cc(m,i,3,k)+cc(m,ic,2,
     +                         k))-ti11* (cc(m,i,5,k)+cc(m,ic,4,k)))) -
     +                         wa2(i-1)* ((cc(m,i,1,k)+tr12* (cc(m,i,3,
     +                         k)-cc(m,ic,2,k))+tr11* (cc(m,i,5,k)-cc(m,
     +                         ic,4,k)))+ (ti12* (cc(m,i-1,3,k)-cc(m,
     +                         ic-1,2,k))-ti11* (cc(m,i-1,5,k)-cc(m,
     +                         ic-1,4,k))))
               ch(m,i,k,3) = wa2(i-2)* ((cc(m,i,1,k)+tr12* (cc(m,i,3,
     +                       k)-cc(m,ic,2,k))+tr11* (cc(m,i,5,k)-cc(m,
     +                       ic,4,k)))+ (ti12* (cc(m,i-1,3,k)-cc(m,ic-1,
     +                       2,k))-ti11* (cc(m,i-1,5,k)-cc(m,ic-1,4,
     +                       k)))) + wa2(i-1)* ((cc(m,i-1,1,
     +                       k)+tr12* (cc(m,i-1,3,k)+cc(m,ic-1,2,
     +                       k))+tr11* (cc(m,i-1,5,k)+cc(m,ic-1,4,k)))-
     +                       (ti12* (cc(m,i,3,k)+cc(m,ic,2,
     +                       k))-ti11* (cc(m,i,5,k)+cc(m,ic,4,k))))
               ch(m,i-1,k,4) = wa3(i-2)* ((cc(m,i-1,1,k)+tr12* (cc(m,
     +                         i-1,3,k)+cc(m,ic-1,2,k))+tr11* (cc(m,i-1,
     +                         5,k)+cc(m,ic-1,4,k)))+
     +                         (ti12* (cc(m,i,3,k)+cc(m,ic,2,
     +                         k))-ti11* (cc(m,i,5,k)+cc(m,ic,4,k)))) -
     +                         wa3(i-1)* ((cc(m,i,1,k)+tr12* (cc(m,i,3,
     +                         k)-cc(m,ic,2,k))+tr11* (cc(m,i,5,k)-cc(m,
     +                         ic,4,k)))- (ti12* (cc(m,i-1,3,k)-cc(m,
     +                         ic-1,2,k))-ti11* (cc(m,i-1,5,k)-cc(m,
     +                         ic-1,4,k))))
               ch(m,i,k,4) = wa3(i-2)* ((cc(m,i,1,k)+tr12* (cc(m,i,3,
     +                       k)-cc(m,ic,2,k))+tr11* (cc(m,i,5,k)-cc(m,
     +                       ic,4,k)))- (ti12* (cc(m,i-1,3,k)-cc(m,ic-1,
     +                       2,k))-ti11* (cc(m,i-1,5,k)-cc(m,ic-1,4,
     +                       k)))) + wa3(i-1)* ((cc(m,i-1,1,
     +                       k)+tr12* (cc(m,i-1,3,k)+cc(m,ic-1,2,
     +                       k))+tr11* (cc(m,i-1,5,k)+cc(m,ic-1,4,k)))+
     +                       (ti12* (cc(m,i,3,k)+cc(m,ic,2,
     +                       k))-ti11* (cc(m,i,5,k)+cc(m,ic,4,k))))
               ch(m,i-1,k,5) = wa4(i-2)* ((cc(m,i-1,1,k)+tr11* (cc(m,
     +                         i-1,3,k)+cc(m,ic-1,2,k))+tr12* (cc(m,i-1,
     +                         5,k)+cc(m,ic-1,4,k)))+
     +                         (ti11* (cc(m,i,3,k)+cc(m,ic,2,
     +                         k))+ti12* (cc(m,i,5,k)+cc(m,ic,4,k)))) -
     +                         wa4(i-1)* ((cc(m,i,1,k)+tr11* (cc(m,i,3,
     +                         k)-cc(m,ic,2,k))+tr12* (cc(m,i,5,k)-cc(m,
     +                         ic,4,k)))- (ti11* (cc(m,i-1,3,k)-cc(m,
     +                         ic-1,2,k))+ti12* (cc(m,i-1,5,k)-cc(m,
     +                         ic-1,4,k))))
               ch(m,i,k,5) = wa4(i-2)* ((cc(m,i,1,k)+tr11* (cc(m,i,3,
     +                       k)-cc(m,ic,2,k))+tr12* (cc(m,i,5,k)-cc(m,
     +                       ic,4,k)))- (ti11* (cc(m,i-1,3,k)-cc(m,ic-1,
     +                       2,k))+ti12* (cc(m,i-1,5,k)-cc(m,ic-1,4,
     +                       k)))) + wa4(i-1)* ((cc(m,i-1,1,
     +                       k)+tr11* (cc(m,i-1,3,k)+cc(m,ic-1,2,
     +                       k))+tr12* (cc(m,i-1,5,k)+cc(m,ic-1,4,k)))+
     +                       (ti11* (cc(m,i,3,k)+cc(m,ic,2,
     +                       k))+ti12* (cc(m,i,5,k)+cc(m,ic,4,k))))
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      RETURN
      END SUBROUTINE vradb5

      SUBROUTINE vradbg(mp,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,
     +                  mdimc,wa)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      USE m_constants, ONLY : pimach
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: mp, ido, ip, l1, idl1, mdimc
      REAL, INTENT(IN)    :: cc(mdimc,ido,ip,l1), wa(ido)
      REAL, INTENT(INOUT) :: ch(mdimc,ido,l1,ip), ch2(mdimc,idl1,ip),
     +                       c1(mdimc,ido,l1,ip), c2(mdimc,idl1,ip)

      INTEGER :: i, j, j2, k, l, m, ik, is, idp2, nbd, ipp2, ipph,
     +           ic, jc, lc, idij
      REAL :: ai1, ar1, ar2, ai2, ar1h, ar2h, tpi, arg, dcp, dsp,
     +        dc2, ds2

      tpi = 2.*pimach()
      arg = tpi/real(ip)
      dcp = cos(arg)
      dsp = sin(arg)
      idp2 = ido + 2
      nbd = (ido-1)/2
      ipp2 = ip + 2
      ipph = (ip+1)/2
      IF (ido.LT.l1) GO TO 40
      DO 30 k = 1,l1
         DO 20 i = 1,ido
            DO 10 m = 1,mp
               ch(m,i,k,1) = cc(m,i,1,k)
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
      GO TO 80
   40 DO 70 i = 1,ido
         DO 60 k = 1,l1
            DO 50 m = 1,mp
               ch(m,i,k,1) = cc(m,i,1,k)
   50       CONTINUE
   60    CONTINUE
   70 CONTINUE
   80 DO 110 j = 2,ipph
         jc = ipp2 - j
         j2 = j + j
         DO 100 k = 1,l1
            DO 90 m = 1,mp
               ch(m,1,k,j) = cc(m,ido,j2-2,k) + cc(m,ido,j2-2,k)
               ch(m,1,k,jc) = cc(m,1,j2-1,k) + cc(m,1,j2-1,k)
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE
      IF (ido.EQ.1) GO TO 210
      IF (nbd.LT.l1) GO TO 160
      DO 150 j = 2,ipph
         jc = ipp2 - j
         DO 140 k = 1,l1
            DO 130 i = 3,ido,2
               ic = idp2 - i
               DO 120 m = 1,mp
                  ch(m,i-1,k,j) = cc(m,i-1,2*j-1,k) + cc(m,ic-1,2*j-2,k)
                  ch(m,i-1,k,jc) = cc(m,i-1,2*j-1,k) -
     +                             cc(m,ic-1,2*j-2,k)
                  ch(m,i,k,j) = cc(m,i,2*j-1,k) - cc(m,ic,2*j-2,k)
                  ch(m,i,k,jc) = cc(m,i,2*j-1,k) + cc(m,ic,2*j-2,k)
  120          CONTINUE
  130       CONTINUE
  140    CONTINUE
  150 CONTINUE
      GO TO 210
  160 DO 200 j = 2,ipph
         jc = ipp2 - j
         DO 190 i = 3,ido,2
            ic = idp2 - i
            DO 180 k = 1,l1
               DO 170 m = 1,mp
                  ch(m,i-1,k,j) = cc(m,i-1,2*j-1,k) + cc(m,ic-1,2*j-2,k)
                  ch(m,i-1,k,jc) = cc(m,i-1,2*j-1,k) -
     +                             cc(m,ic-1,2*j-2,k)
                  ch(m,i,k,j) = cc(m,i,2*j-1,k) - cc(m,ic,2*j-2,k)
                  ch(m,i,k,jc) = cc(m,i,2*j-1,k) + cc(m,ic,2*j-2,k)
  170          CONTINUE
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE
  210 ar1 = 1.
      ai1 = 0.
      DO 270 l = 2,ipph
         lc = ipp2 - l
         ar1h = dcp*ar1 - dsp*ai1
         ai1 = dcp*ai1 + dsp*ar1
         ar1 = ar1h
         DO 230 ik = 1,idl1
            DO 220 m = 1,mp
               c2(m,ik,l) = ch2(m,ik,1) + ar1*ch2(m,ik,2)
               c2(m,ik,lc) = ai1*ch2(m,ik,ip)
  220       CONTINUE
  230    CONTINUE
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         DO 260 j = 3,ipph
            jc = ipp2 - j
            ar2h = dc2*ar2 - ds2*ai2
            ai2 = dc2*ai2 + ds2*ar2
            ar2 = ar2h
            DO 250 ik = 1,idl1
               DO 240 m = 1,mp
                  c2(m,ik,l) = c2(m,ik,l) + ar2*ch2(m,ik,j)
                  c2(m,ik,lc) = c2(m,ik,lc) + ai2*ch2(m,ik,jc)
  240          CONTINUE
  250       CONTINUE
  260    CONTINUE
  270 CONTINUE
      DO 300 j = 2,ipph
         DO 290 ik = 1,idl1
            DO 280 m = 1,mp
               ch2(m,ik,1) = ch2(m,ik,1) + ch2(m,ik,j)
  280       CONTINUE
  290    CONTINUE
  300 CONTINUE
      DO 330 j = 2,ipph
         jc = ipp2 - j
         DO 320 k = 1,l1
            DO 310 m = 1,mp
               ch(m,1,k,j) = c1(m,1,k,j) - c1(m,1,k,jc)
               ch(m,1,k,jc) = c1(m,1,k,j) + c1(m,1,k,jc)
  310       CONTINUE
  320    CONTINUE
  330 CONTINUE
      IF (ido.EQ.1) GO TO 430
      IF (nbd.LT.l1) GO TO 380
      DO 370 j = 2,ipph
         jc = ipp2 - j
         DO 360 k = 1,l1
            DO 350 i = 3,ido,2
               DO 340 m = 1,mp
                  ch(m,i-1,k,j) = c1(m,i-1,k,j) - c1(m,i,k,jc)
                  ch(m,i-1,k,jc) = c1(m,i-1,k,j) + c1(m,i,k,jc)
                  ch(m,i,k,j) = c1(m,i,k,j) + c1(m,i-1,k,jc)
                  ch(m,i,k,jc) = c1(m,i,k,j) - c1(m,i-1,k,jc)
  340          CONTINUE
  350       CONTINUE
  360    CONTINUE
  370 CONTINUE
      GO TO 430
  380 DO 420 j = 2,ipph
         jc = ipp2 - j
         DO 410 i = 3,ido,2
            DO 400 k = 1,l1
               DO 390 m = 1,mp
                  ch(m,i-1,k,j) = c1(m,i-1,k,j) - c1(m,i,k,jc)
                  ch(m,i-1,k,jc) = c1(m,i-1,k,j) + c1(m,i,k,jc)
                  ch(m,i,k,j) = c1(m,i,k,j) + c1(m,i-1,k,jc)
                  ch(m,i,k,jc) = c1(m,i,k,j) - c1(m,i-1,k,jc)
  390          CONTINUE
  400       CONTINUE
  410    CONTINUE
  420 CONTINUE
  430 CONTINUE
      IF (ido.EQ.1) RETURN
      DO 450 ik = 1,idl1
         DO 440 m = 1,mp
            c2(m,ik,1) = ch2(m,ik,1)
  440    CONTINUE
  450 CONTINUE
      DO 480 j = 2,ip
         DO 470 k = 1,l1
            DO 460 m = 1,mp
               c1(m,1,k,j) = ch(m,1,k,j)
  460       CONTINUE
  470    CONTINUE
  480 CONTINUE
      IF (nbd.GT.l1) GO TO 530
      is = -ido
      DO 520 j = 2,ip
         is = is + ido
         idij = is
         DO 510 i = 3,ido,2
            idij = idij + 2
            DO 500 k = 1,l1
               DO 490 m = 1,mp
                  c1(m,i-1,k,j) = wa(idij-1)*ch(m,i-1,k,j) -
     +                            wa(idij)*ch(m,i,k,j)
                  c1(m,i,k,j) = wa(idij-1)*ch(m,i,k,j) +
     +                          wa(idij)*ch(m,i-1,k,j)
  490          CONTINUE
  500       CONTINUE
  510    CONTINUE
  520 CONTINUE
      GO TO 580
  530 is = -ido
      DO 570 j = 2,ip
         is = is + ido
         DO 560 k = 1,l1
            idij = is
            DO 550 i = 3,ido,2
               idij = idij + 2
               DO 540 m = 1,mp
                  c1(m,i-1,k,j) = wa(idij-1)*ch(m,i-1,k,j) -
     +                            wa(idij)*ch(m,i,k,j)
                  c1(m,i,k,j) = wa(idij-1)*ch(m,i,k,j) +
     +                          wa(idij)*ch(m,i-1,k,j)
  540          CONTINUE
  550       CONTINUE
  560    CONTINUE
  570 CONTINUE
  580 RETURN
      END SUBROUTINE vradbg

      SUBROUTINE vrftb1(m,n,c,ch,mdimc,wa,fac)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: m, n, mdimc
      REAL, INTENT(IN)    :: wa(n), fac(15)
      REAL, INTENT(INOUT) :: c(mdimc,n), ch(mdimc,n)

      INTEGER :: i, j, nf, na, l1, iw, k1, ip, ipo, idl1, ix2, ix3, ix4
      INTEGER :: ido, l2
      REAL    :: scale

      nf = fac(2)
      na = 0
      l1 = 1
      iw = 1
      DO 160 k1 = 1,nf
         ip = fac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idl1 = ido*l1
         IF (ip.NE.4) GO TO 30
         ix2 = iw + ido
         ix3 = ix2 + ido
         IF (na.NE.0) GO TO 10
         CALL vradb4(m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2),wa(ix3))
         GO TO 20
   10    CALL vradb4(m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2),wa(ix3))
   20    na = 1 - na
         GO TO 150
   30    IF (ip.NE.2) GO TO 60
         IF (na.NE.0) GO TO 40
         CALL vradb2(m,ido,l1,c,ch,mdimc,wa(iw))
         GO TO 50
   40    CALL vradb2(m,ido,l1,ch,c,mdimc,wa(iw))
   50    na = 1 - na
         GO TO 150
   60    IF (ip.NE.3) GO TO 90
         ix2 = iw + ido
         IF (na.NE.0) GO TO 70
         CALL vradb3(m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2))
         GO TO 80
   70    CALL vradb3(m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2))
   80    na = 1 - na
         GO TO 150
   90    IF (ip.NE.5) GO TO 120
         ix2 = iw + ido
         ix3 = ix2 + ido
         ix4 = ix3 + ido
         IF (na.NE.0) GO TO 100
         CALL vradb5(m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         GO TO 110
  100    CALL vradb5(m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  110    na = 1 - na
         GO TO 150
  120    IF (na.NE.0) GO TO 130
         CALL vradbg(m,ido,ip,l1,idl1,c,c,c,ch,ch,mdimc,wa(iw))
         GO TO 140
  130    CALL vradbg(m,ido,ip,l1,idl1,ch,ch,ch,c,c,mdimc,wa(iw))
  140    IF (ido.EQ.1) na = 1 - na
  150    l1 = l2
         iw = iw + (ip-1)*ido
  160 CONTINUE
      scale = sqrt(1./n)
      IF (na.EQ.0) GO TO 190
      DO 180 j = 1,n
         DO 170 i = 1,m
            c(i,j) = scale*ch(i,j)
  170    CONTINUE
  180 CONTINUE
      RETURN
  190 DO 210 j = 1,n
         DO 200 i = 1,m
            c(i,j) = scale*c(i,j)
  200    CONTINUE
  210 CONTINUE
      RETURN
      END SUBROUTINE vrftb1

      END SUBROUTINE vrfftb


      END MODULE m_rfft
