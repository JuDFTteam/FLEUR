!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rfft
   use m_juDFT
   PRIVATE
! odule also contains routines vrffti,vrfftf&vrfftb as private routines below
   PUBLIC rfft
CONTAINS
   SUBROUTINE rfft( &
      isn,n1d,n2d,n3d,n1,n2,n3, &
      nw1,nw2,nw3,wsave,b, &
      a)
! **********************************************************************

! isn = +1 : FFT a real, 3dim array "a" of dimensions n1*n2*n3 to a complex
!            array of size n1*n2*(n3+1)/2 [if n3 is odd] or n1*n2*(n3/2+1)
!            [n3 is even].
! isn = -1 : back-FFT of a complex array "a" of size n1*n2*(n3+1)/2 [odd]
!            or n1*n2*(n3/2+1) [even] to a real, 3dim array

! the actual array is assumed to be located between
!                                  1 ... nw1+1    &     n1-nw1+1 ... n1
!                                  1 ... nw2+1    &     n2-nw2+1 ... n2
!                                  1 ... nw3+1    &     n3-nw3+1 ... n3
! and padded with zeros in between.
! if nw1 >= (n1-1)/2, no padding is assumed (-> aliasing errors!)

!                                                  G.Bihlmayer (UniWien)

! **********************************************************************
      USE m_cfft
      IMPLICIT NONE

      INTEGER :: n1d,n2d,n3d,n1,n2,n3,nw1,nw2,nw3,isn
      REAL :: a(n1d,n2d,0:n3d),b(n1d,n2d,n3d),wsave(n3d+15)

      INTEGER :: i1,i2,i3,nup
      REAL :: factor
      LOGICAL :: l_nopad

! a ... array for FFT
! b ... work array
! wsave ... stores tables for r-2-c FFT
! n1,n2,n3 ... dimensions to be transformed
! nw1,nw2,nw3 ... actual dimensions of a before FFT
! n1d,n2d,n3d ... dimensions of a,b

! check for input errors

      IF ((isn/=-1) .AND. (isn /= 1))  CALL juDFT_error("choose isn=+/- 1" &
                                                        ,calledby ="rfft")
      IF ((n1d < n1) .OR. (n2d < n2) .OR. (n3d < n3)) THEN
         WRITE (6,*) 'n1d,n2d,n3d =',n1d,n2d,n3d
         WRITE (6,*) 'n1 ,n2 ,n3  =',n1 ,n2 ,n3
         CALL juDFT_error("n(i) > n(i)d",calledby ="rfft")
      ENDIF
      IF ((n1 <= 2*nw1+1) .OR. &
          (n2 <= 2*nw2+1) .OR. &
          (n3 <= 2*nw3+1)) THEN
         !        WRITE (6,*) 'n1 ,n2 ,n3  =',n1 ,n2 ,n3
         !        WRITE (6,*) 'nw1,nw2,nw3 =',nw1,nw2,nw3
         l_nopad= .TRUE.
      ELSE
         l_nopad= .FALSE.
      ENDIF

! ******************** f o r w a r d - t r a n s f o r m *******************

      IF (isn == 1) THEN

         ! array a is assumed to be zero from (1,1,0) to (n1d,n2d,0) and the array
         ! to be FFT'd starts at (1,1,1) as n1*n2*n3 real numbers.
         ! first transform n1*n2 real sequences of lenghth n3 to n3/2 complex values

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

         ! now we have the FFT'd array stored as described in vrfftf
         ! (mixed real & compex data)
         ! remove the norm 1/sqrt(n3) (to be compatible with cfft)

         factor = sqrt(1.0*n3)
         ! ALL CPP_BLAS_sscal(n1d*n2d*n3,factor,a(1,1,1),1)
         a=a*factor

         ! now, the real part of f(0) has to be moved to a(n1,n2,0) to get a purely
         ! complex array starting at a(1,1,0)

         DO i1=1,n1
            DO i2=1,n2
               a(i1,i2,0)=a(i1,i2,1)
               a(i1,i2,1)=0.0
            ENDDO
         ENDDO

      ENDIF

! ******************** g e n e r a l  p a r t *******************

! now perform n2*n3/2 and n1*n3/2 complex FFT's; a is assumed to be
! complex and starting at (1,1,0)

      IF (ABS((n3/2.)-NINT(n3/2.)) > 0.1) THEN
         nup = n3
      ELSE
         nup = n3+1
         IF (n3+1>n3d)  CALL juDFT_error("n3 even & n3+1 > n3d" ,calledby &
                                         ="rfft")
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
            CALL cfft(a(1,n2-nw2+1,i3-1),a(1,n2-nw2+1,i3), &
                      nw2*n1,n1,n1,isn)
            CALL cfft(a(1,1,i3-1),a(1,1,i3),n1*n2,n2,n1*n2,isn)
         ENDDO
      ENDIF

! ******************** b a c k w a r d - t r a n s f o r m *******************

      IF (isn == -1) THEN

         ! the real part of f(0) has to be moved to a(n1,n2,1) for a correct
         ! setup for vrfftb (see comments therein)

         DO i1=1,n1
            DO i2=1,n2
               a(i1,i2,1)=a(i1,i2,0)
               a(i1,i2,0)=0.0
            ENDDO
         ENDDO

         ! transform n1*n2 mixed real and complex sequences of lenghth n3/2
         ! to n3 real values

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

         ! remove the norm 1/sqrt(n3) (compatibility with cfft)

         factor = sqrt(1.0*n3)
         ! ALL CPP_BLAS_sscal(n1d*n2d*n3,factor,a(1,1,1),1)
         a=a*factor

      ENDIF

   END SUBROUTINE rfft

   SUBROUTINE vrffti(n,wsave)
!***BEGIN PROLOGUE  VRFFTI
!***DATE WRITTEN   860701   (YYMMDD)
!***REVISION DATE  900509   (YYMMDD)
!***CATEGORY NO.  J1A1
!***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM,
!             MULTIPLE SEQUENCES
!***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST)
!***PURPOSE  Initialization for VRFFTF and VRFFTB.
!***DESCRIPTION

!  Subroutine VRFFTI initializes the array WSAVE which is used in
!  both VRFFTF and VRFFTB.  The prime factorization of N together with
!  a tabulation of certain trigonometric functions are computed and
!  stored in the array WSAVE.

!  Input Parameter

!  N       the length of the sequence to be transformed.  There is no
!          restriction on N.

!  Output Parameter

!  WSAVE   a work array which must be dimensioned at least N+15.
!          The same work array can be used for both VRFFTF and VRFFTB
!          as long as N remains unchanged.  Different WSAVE arrays
!          are required for different values of N.  The contents of
!          WSAVE must not be changed between calls of VRFFTF or VRFFTB.

!              * * * * * * * * * * * * * * * * * * * * *
!              *                                       *
!              *         PROGRAM SPECIFICATIONS        *
!              *                                       *
!              * * * * * * * * * * * * * * * * * * * * *

!     Dimension of    R(MDIMR,N), RT(MDIMR,N), WSAVE(N+15)
!     Arguments

!     Latest          AUGUST 1, 1985
!     Revision

!     Subprograms     VRFFTI, VRFTI1, VRFFTF, VRFTF1, VRADF2, VRADF3,
!     Required        VRADF4, VRADF5, VRADFG, VRFFTB, VRFTB1, VRADB2,
!                     VRADB3, VRADB4, VRADB5, VRADBG, PIMACH

!     Special         NONE
!     Conditions

!     Common          NONE
!     blocks

!     I/O             NONE

!     Precision       SINGLE

!     Specialist      ROLAND SWEET

!     Language        FORTRAN

!     History         WRITTEN BY LINDA LINDGREN AND ROLAND SWEET AT THE
!                     NATIONAL BUREAU OF STANDARDS (BOULDER).

!     Algorithm       A REAL VARIANT OF THE STOCKHAM AUTOSORT VERSION
!                     OF THE COOLEY-TUKEY FAST FOURIER TRANSFORM.

!     Portability     AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN 77.
!                     THE ONLY MACHINE DEPENDENT CONSTANT IS LOCATED IN
!                     THE FUNCTION PIMACH.

!     Required        COS,SIN
!     resident
!     Routines

!***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
!               pp. 51-83.
!***ROUTINES CALLED  VRFTI1
!***END PROLOGUE  VRFFTI

!     VRFFTPK, VERSION 1, AUGUST 1985

      DIMENSION wsave(n+15)
!***FIRST EXECUTABLE STATEMENT  VRFFTI
      IF (n == 1) RETURN
      CALL vrfti1(n,wsave(1),wsave(n+1))
      RETURN
   END subroutine
   SUBROUTINE vrfti1(n,wa,fac)

!     VRFFTPK, VERSION 1, AUGUST 1985

      USE m_constants, ONLY : pimach
      DIMENSION wa(n),fac(15),ntryh(4)
      DATA ntryh(1),ntryh(2),ntryh(3),ntryh(4)/4,2,3,5/

      nl = n
      nf = 0
      j = 0
10    j = j + 1
      IF ( j <=4 ) THEN
         ntry = ntryh(j)
      ELSE
         ntry = ntry + 2
      ENDIF
40    nq = nl/ntry
      nr = nl - ntry*nq
      IF ( nr /= 0 ) THEN
         GOTO 10
      ENDIF
      nf = nf + 1
      fac(nf+2) = ntry
      nl = nq
      IF (ntry /= 2) GO TO 70
      IF (nf == 1) GO TO 70
      DO 60 i = 2,nf
         ib = nf - i + 2
         fac(ib+2) = fac(ib+1)
60    END DO
      fac(3) = 2
70    IF (nl /= 1) GO TO 40
      fac(1) = n
      fac(2) = nf
      tpi = 2.*pimach()
      argh = tpi/real(n)
      is = 0
      nfm1 = nf - 1
      l1 = 1
      IF (nfm1 == 0) RETURN
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
80          END DO
            is = is + ido
90       END DO
         l1 = l2
100   END DO
      RETURN
   END subroutine
   SUBROUTINE vrfftf(m,n,r,rt,mdimr,wsave)

!***BEGIN PROLOGUE  VRFFTF
!***DATE WRITTEN   850801   (YYMMDD)
!***REVISION DATE  900509   (YYMMDD)
!***CATEGORY NO.  J1A1
!***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM,
!             FOURIER ANALYSIS, FORWARD TRANSFORM, MULTIPLE SEQUENCES
!***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST)
!***PURPOSE  Forward real periodic transform, M sequences.
!***DESCRIPTION

!  Subroutine VRFFTF computes the Fourier coefficients (forward
!  transform) of a number of real periodic sequences.  Specifically,
!  for each sequence the subroutine claculates the independent
!  Fourier coefficients described below at output parameter R.

!  The array WSAVE which is used by subroutine VRFFTF must be
!  initialized by calling subroutine VRFFTI(N,WSAVE).

!  Input Parameters

!  M       the number of sequences to be transformed.

!  N       the length of the sequences to be transformed.  The method
!          is most efficient when N is a product of small primes,
!          however n may be any positive integer.

!  R       areal two-dimensional array of size MDIMX x N containing the
!          the sequences to be transformed.  The sequences are stored
!          in the ROWS of R.  Thus, the I-th sequence to be transformed,
!          X(I,J), J=0,1,...,N-1, is stored as

!               R(I,J) = X(I,J-1) , J=1, 2, . . . , N.

!  RT      a real two-dimensional work array of size MDIMX x N.

!  MDIMR   the row (or first) dimension of the arrays R and RT exactly
!          as they appear in the calling program.  This parameter is
!          used to specify the variable dimension of these arrays.

!  WSAVE   a real one-dimensional work array which must be dimensioned
!          at least N+15.  The WSAVE array must be initialized by
!          calling subroutine VRFFTI.  A different WSAVE array must be
!          used for each different value of N.  This initialization does
!          not have to be repeated so long as N remains unchanged.  The
!          same WSAVE array may be used by VRFFTF and VRFFTB.

!  Output Parameters

!  R       contains the Fourier coefficients F(K) for each of the M
!          input sequences.  Specifically, row I of R, R(I,J),
!          J=1,2,..,N, contains the independent Fourier coefficients
!          F(I,K), for the I-th input sequence stored as

!             R(I,1) = REAL( F(I,0) ),
!                    = SQRT(1/N)*SUM(J=0,N-1)[ X(I,J) ],

!             R(I,2*K) = REAL( F(I,K) )
!                      = SQRT(1/N)*SUM(J=0,N-1)[X(I,J)*COS(2J*K*PI/N)]

!             R(I,2*K+1) = IMAG( F(I,K) )
!                        =-SQRT(1/N)*SUM(J=0,N-1)[X(I,J)*SIN(2J*K*PI/N)]

!                   for K = 1, 2, . . . , M-1,

!              and, when N is even,

!              R(I,N) = REAL( F(I,N/2) ).
!                     = SQRT(1/N)*SUM(J=0,N-1)[ (-1)**J*X(I,J) ].

!  WSAVE   contains results which must not be destroyed between calls
!          to VRFFTF or VRFFTB.

!  -----------------------------------------------------------------

!  NOTE  -  A call of VRFFTF followed immediately by a call of
!           of VRFFTB will return the original sequences R.  Thus,
!           VRFFTB is the correctly normalized inverse of VRFFTF.

!  -----------------------------------------------------------------

!  VRFFTF is a straightforward extension of the subprogram RFFTF to
!  handle M simultaneous sequences.  RFFTF was originally developed
!  by P. N. Swarztrauber of NCAR.

!              * * * * * * * * * * * * * * * * * * * * *
!              *                                       *
!              *         PROGRAM SPECIFICATIONS        *
!              *                                       *
!              * * * * * * * * * * * * * * * * * * * * *

!     Dimension of    R(MDIMR,N), RT(MDIMR,N), WSAVE(N+15)
!     Arguments

!     Latest          AUGUST 1, 1985
!     Revision

!     Subprograms     VRFFTI, VRFTI1, VRFFTF, VRFTF1, VRADF2, VRADF3,
!     Required        VRADF4, VRADF5, VRADFG, VRFFTB, VRFTB1, VRADB2,
!                     VRADB3, VRADB4, VRADB5, VRADBG, PIMACH

!     Special         NONE
!     Conditions

!     Common          NONE
!     blocks

!     I/O             NONE

!     Precision       SINGLE

!     Specialist      ROLAND SWEET

!     Language        FORTRAN

!     History         WRITTEN BY LINDA LINDGREN AND ROLAND SWEET AT THE
!                     NATIONAL BUREAU OF STANDARDS (BOULDER).

!     Algorithm       A REAL VARIANT OF THE STOCKHAM AUTOSORT VERSION
!                     OF THE COOLEY-TUKEY FAST FOURIER TRANSFORM.

!     Portability     AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN 77.
!                     THE ONLY MACHINE DEPENDENT CONSTANT IS LOCATED IN
!                     THE FUNCTION PIMACH.

!     Required        COS,SIN
!     resident
!     Routines

!***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
!               pp. 51-83.
!***ROUTINES CALLED  VRFTF1
!***END PROLOGUE  VRFFTF

!     VRFFTPK, VERSION 1, AUGUST 1985

      DIMENSION r(mdimr,n),rt(mdimr,n),wsave(n+15)
!***FIRST EXECUTABLE STATEMENT  VRFFTF
      IF (n == 1) RETURN
      CALL vrftf1(m,n,r,rt,mdimr,wsave(1),wsave(n+1))
      RETURN
   END subroutine
   SUBROUTINE vradf2(mp,ido,l1,cc,ch,mdimc,wa1)

!     VRFFTPK, VERSION 1, AUGUST 1985

      DIMENSION ch(mdimc,ido,2,l1),cc(mdimc,ido,l1,2),wa1(ido)

      DO 20 k = 1,l1
         DO 10 m = 1,mp
            ch(m,1,1,k) = cc(m,1,k,1) + cc(m,1,k,2)
            ch(m,ido,2,k) = cc(m,1,k,1) - cc(m,1,k,2)
10       END DO
20    END DO
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
               ch(m,i,1,k) = cc(m,i,k,1) + &
                             (wa1(i-2)*cc(m,i,k,2)-wa1(i-1)* &
                              cc(m,i-1,k,2))
               ch(m,ic,2,k) = (wa1(i-2)*cc(m,i,k,2)- &
                               wa1(i-1)*cc(m,i-1,k,2)) - cc(m,i,k,1)
               ch(m,i-1,1,k) = cc(m,i-1,k,1) + &
                               (wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)* &
                                cc(m,i,k,2))
               ch(m,ic-1,2,k) = cc(m,i-1,k,1) - &
                                (wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)* &
                                 cc(m,i,k,2))
40          END DO
50       END DO
60    END DO
      IF (mod(ido,2) == 1) RETURN
70    DO 90 k = 1,l1
         DO 80 m = 1,mp
            ch(m,1,2,k) = -cc(m,ido,k,2)
            ch(m,ido,1,k) = cc(m,ido,k,1)
80       END DO
90    END DO
100   RETURN
   END subroutine
   SUBROUTINE vradf3(mp,ido,l1,cc,ch,mdimc,wa1,wa2)

!     VRFFTPK, VERSION 1, AUGUST 1985

      USE m_constants, ONLY : pimach
      DIMENSION ch(mdimc,ido,3,l1),cc(mdimc,ido,l1,3),wa1(ido),wa2(ido)

      arg = 2.*pimach()/3.
      taur = cos(arg)
      taui = sin(arg)
      DO 20 k = 1,l1
         DO 10 m = 1,mp
            ch(m,1,1,k) = cc(m,1,k,1) + (cc(m,1,k,2)+cc(m,1,k,3))
            ch(m,1,3,k) = taui* (cc(m,1,k,3)-cc(m,1,k,2))
            ch(m,ido,2,k) = cc(m,1,k,1) + &
                            taur* (cc(m,1,k,2)+cc(m,1,k,3))
10       END DO
20    END DO
      IF (ido == 1) RETURN
      idp2 = ido + 2
      DO 50 k = 1,l1
         DO 40 i = 3,ido,2
            ic = idp2 - i
            DO 30 m = 1,mp
               ch(m,i-1,1,k) = cc(m,i-1,k,1) + &
               ((wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i, &
               k,2))+ (wa2(i-2)*cc(m,i-1,k, &
               &                          3)+wa2(i-1)*cc(m,i,k,3)))
               ch(m,i,1,k) = cc(m,i,k,1) + &
               ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k, &
               &                        2))+ (wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*cc(m, &
               i-1,k,3)))
               ch(m,i-1,3,k) = (cc(m,i-1,k,1)+ &
               taur* ((wa1(i-2)*cc(m,i-1,k, &
               &                          2)+wa1(i-1)*cc(m,i,k,2))+ (wa2(i-2)*cc(m, &
               i-1,k,3)+wa2(i-1)*cc(m,i,k,3)))) + &
               (taui* ((wa1(i-2)*cc(m,i,k, &
               &                          2)-wa1(i-1)*cc(m,i-1,k, &
               &                          2))- (wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*cc(m, &
               i-1,k,3))))
               ch(m,ic-1,2,k) = (cc(m,i-1,k,1)+ &
               taur* ((wa1(i-2)*cc(m,i-1,k, &
               &                           2)+wa1(i-1)*cc(m,i,k, &
               &                           2))+ (wa2(i-2)*cc(m,i-1,k, &
               &                           3)+wa2(i-1)*cc(m,i,k,3)))) - &
               (taui* ((wa1(i-2)*cc(m,i,k, &
               &                           2)-wa1(i-1)*cc(m,i-1,k, &
               &                           2))- (wa2(i-2)*cc(m,i,k, &
               &                           3)-wa2(i-1)*cc(m,i-1,k,3))))
               ch(m,i,3,k) = (cc(m,i,k,1)+taur* &
               ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k, &
               &                        2))+ (wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*cc(m, &
               i-1,k,3)))) + (taui* ((wa2(i-2)*cc(m,i-1,k, &
               &                        3)+wa2(i-1)*cc(m,i,k,3))- (wa1(i-2)*cc(m, &
               i-1,k,2)+wa1(i-1)*cc(m,i,k,2))))
               ch(m,ic,2,k) = (taui* ((wa2(i-2)*cc(m,i-1,k, &
               &                         3)+wa2(i-1)*cc(m,i,k,3))- (wa1(i-2)*cc(m, &
               i-1,k,2)+wa1(i-1)*cc(m,i,k,2)))) - &
               (cc(m,i,k,1)+taur* ((wa1(i-2)*cc(m,i,k, &
               &                         2)-wa1(i-1)*cc(m,i-1,k, &
               &                         2))+ (wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*cc(m, &
               i-1,k,3))))
30          END DO
40       END DO
50    END DO
      RETURN
   END subroutine
   SUBROUTINE vradf4(mp,ido,l1,cc,ch,mdimc,wa1,wa2,wa3)

!     VRFFTPK, VERSION 1, AUGUST 1985

      DIMENSION cc(mdimc,ido,l1,4),ch(mdimc,ido,4,l1),wa1(ido),wa2(ido), &
         wa3(ido)

      hsqt2 = sqrt(2.)/2.
      DO 20 k = 1,l1
         DO 10 m = 1,mp
            ch(m,1,1,k) = (cc(m,1,k,2)+cc(m,1,k,4)) + &
                          (cc(m,1,k,1)+cc(m,1,k,3))
            ch(m,ido,4,k) = (cc(m,1,k,1)+cc(m,1,k,3)) - &
                            (cc(m,1,k,2)+cc(m,1,k,4))
            ch(m,ido,2,k) = cc(m,1,k,1) - cc(m,1,k,3)
            ch(m,1,3,k) = cc(m,1,k,4) - cc(m,1,k,2)
10       END DO
20    END DO
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
               ch(m,i-1,1,k) = ((wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i, &
               k,2))+ (wa3(i-2)*cc(m,i-1,k, &
               &                          4)+wa3(i-1)*cc(m,i,k,4))) + &
               (cc(m,i-1,k,1)+ (wa2(i-2)*cc(m,i-1,k, &
               &                          3)+wa2(i-1)*cc(m,i,k,3)))
               ch(m,ic-1,4,k) = (cc(m,i-1,k,1)+ &
               (wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i, &
               k,3))) - ((wa1(i-2)*cc(m,i-1,k, &
               &                           2)+wa1(i-1)*cc(m,i,k,2))+ &
               (wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i, &
               k,4)))
               ch(m,i,1,k) = ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k, &
               &                        2))+ (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m, &
               i-1,k,4))) + (cc(m,i,k,1)+ &
               (wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*cc(m,i-1,k, &
               &                        3)))
               ch(m,ic,4,k) = ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1, &
               k,2))+ (wa3(i-2)*cc(m,i,k, &
               &                         4)-wa3(i-1)*cc(m,i-1,k,4))) - &
               (cc(m,i,k,1)+ (wa2(i-2)*cc(m,i,k, &
               &                         3)-wa2(i-1)*cc(m,i-1,k,3)))
               ch(m,i-1,3,k) = ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1, &
               k,2))- (wa3(i-2)*cc(m,i,k, &
               &                          4)-wa3(i-1)*cc(m,i-1,k,4))) + &
               (cc(m,i-1,k,1)- (wa2(i-2)*cc(m,i-1,k, &
               &                          3)+wa2(i-1)*cc(m,i,k,3)))
               ch(m,ic-1,2,k) = (cc(m,i-1,k,1)- &
               (wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i, &
               k,3))) - ((wa1(i-2)*cc(m,i,k, &
               &                           2)-wa1(i-1)*cc(m,i-1,k,2))- &
               (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m,i-1, &
               k,4)))
               ch(m,i,3,k) = ((wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i,k, &
               &                        4))- (wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m, &
               i,k,2))) + (cc(m,i,k,1)- &
               (wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*cc(m,i-1,k, &
               &                        3)))
               ch(m,ic,2,k) = ((wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i, &
               k,4))- (wa1(i-2)*cc(m,i-1,k, &
               &                         2)+wa1(i-1)*cc(m,i,k,2))) - &
               (cc(m,i,k,1)- (wa2(i-2)*cc(m,i,k, &
               &                         3)-wa2(i-1)*cc(m,i-1,k,3)))
40          END DO
50       END DO
60    END DO
      IF (mod(ido,2) == 1) RETURN
70    CONTINUE
      DO 90 k = 1,l1
         DO 80 m = 1,mp
            ch(m,ido,1,k) = (hsqt2* (cc(m,ido,k,2)-cc(m,ido,k,4))) + &
                            cc(m,ido,k,1)
            ch(m,ido,3,k) = cc(m,ido,k,1) - &
                            (hsqt2* (cc(m,ido,k,2)-cc(m,ido,k,4)))
            ch(m,1,2,k) = (-hsqt2* (cc(m,ido,k,2)+cc(m,ido,k,4))) - &
                          cc(m,ido,k,3)
            ch(m,1,4,k) = (-hsqt2* (cc(m,ido,k,2)+cc(m,ido,k,4))) + &
                          cc(m,ido,k,3)
80       END DO
90    END DO
100   RETURN
   END subroutine
   SUBROUTINE vradf5(mp,ido,l1,cc,ch,mdimc,wa1,wa2,wa3,wa4)

!     VRFFTPK, VERSION 1, AUGUST 1985

      USE m_constants, ONLY : pimach
      DIMENSION cc(mdimc,ido,l1,5),ch(mdimc,ido,5,l1),wa1(ido),wa2(ido), &
         wa3(ido),wa4(ido)

      arg = 2.*pimach()/5.
      tr11 = cos(arg)
      ti11 = sin(arg)
      tr12 = cos(2.*arg)
      ti12 = sin(2.*arg)
      DO 20 k = 1,l1
         DO 10 m = 1,mp
            ch(m,1,1,k) = cc(m,1,k,1) + (cc(m,1,k,5)+cc(m,1,k,2)) + &
                          (cc(m,1,k,4)+cc(m,1,k,3))
            ch(m,ido,2,k) = cc(m,1,k,1) + &
                            tr11* (cc(m,1,k,5)+cc(m,1,k,2)) + &
                            tr12* (cc(m,1,k,4)+cc(m,1,k,3))
            ch(m,1,3,k) = ti11* (cc(m,1,k,5)-cc(m,1,k,2)) + &
                          ti12* (cc(m,1,k,4)-cc(m,1,k,3))
            ch(m,ido,4,k) = cc(m,1,k,1) + &
                            tr12* (cc(m,1,k,5)+cc(m,1,k,2)) + &
                            tr11* (cc(m,1,k,4)+cc(m,1,k,3))
            ch(m,1,5,k) = ti12* (cc(m,1,k,5)-cc(m,1,k,2)) - &
                          ti11* (cc(m,1,k,4)-cc(m,1,k,3))
10       END DO
20    END DO
      IF (ido == 1) RETURN
      idp2 = ido + 2
      DO 50 k = 1,l1
         DO 40 i = 3,ido,2
            ic = idp2 - i
            DO 30 m = 1,mp
               ch(m,i-1,1,k) = cc(m,i-1,k,1) + &
               ((wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i, &
               k,2))+ (wa4(i-2)*cc(m,i-1,k, &
               &                          5)+wa4(i-1)*cc(m,i,k,5))) + &
               ((wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i, &
               k,3))+ (wa3(i-2)*cc(m,i-1,k, &
               &                          4)+wa3(i-1)*cc(m,i,k,4)))
               ch(m,i,1,k) = cc(m,i,k,1) + &
               ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k, &
               &                        2))+ (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m, &
               i-1,k,5))) + ((wa2(i-2)*cc(m,i,k, &
               &                        3)-wa2(i-1)*cc(m,i-1,k,3))+ &
               (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m,i-1,k, &
               &                        4)))
               ch(m,i-1,3,k) = cc(m,i-1,k,1) + &
                               tr11* (wa1(i-2)*cc(m,i-1,k,2)+ &
                                      wa1(i-1)*cc(m,i,k,2)+ &
                                      wa4(i-2)*cc(m,i-1,k,5)+ &
                                      wa4(i-1)*cc(m,i,k,5)) + &
                               tr12* (wa2(i-2)*cc(m,i-1,k,3)+ &
                                      wa2(i-1)*cc(m,i,k,3)+ &
                                      wa3(i-2)*cc(m,i-1,k,4)+ &
                                      wa3(i-1)*cc(m,i,k,4)) + &
                               ti11* (wa1(i-2)*cc(m,i,k,2)- &
                                      wa1(i-1)*cc(m,i-1,k,2)- &
                                      (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m,i-1, &
                                                                        k,5))) + ti12* (wa2(i-2)*cc(m,i,k,3)- &
                                                                                        wa2(i-1)*cc(m,i-1,k,3)- &
                                                                                        (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m,i-1, &
                                                                                                                          k,4)))
               ch(m,ic-1,2,k) = cc(m,i-1,k,1) + &
               tr11* (wa1(i-2)*cc(m,i-1,k,2)+ &
               wa1(i-1)*cc(m,i,k,2)+ &
               wa4(i-2)*cc(m,i-1,k,5)+ &
               wa4(i-1)*cc(m,i,k,5)) + &
               tr12* (wa2(i-2)*cc(m,i-1,k,3)+ &
               wa2(i-1)*cc(m,i,k,3)+ &
               wa3(i-2)*cc(m,i-1,k,4)+ &
               wa3(i-1)*cc(m,i,k,4)) - &
               (ti11* (wa1(i-2)*cc(m,i,k, &
               &                           2)-wa1(i-1)*cc(m,i-1,k, &
               &                           2)- (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m, &
               i-1,k,5)))+ti12* (wa2(i-2)*cc(m,i,k, &
               &                           3)-wa2(i-1)*cc(m,i-1,k, &
               &                           3)- (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m, &
               i-1,k,4))))
               ch(m,i,3,k) = (cc(m,i,k,1)+tr11* &
               ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k, &
               &                        2))+ (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m, &
               i-1,k,5)))+tr12* ((wa2(i-2)*cc(m,i,k, &
               &                        3)-wa2(i-1)*cc(m,i-1,k,3))+ (wa3(i-2)*cc(m, &
               i,k,4)-wa3(i-1)*cc(m,i-1,k,4)))) + &
               (ti11* ((wa4(i-2)*cc(m,i-1,k, &
               &                        5)+wa4(i-1)*cc(m,i,k,5))- (wa1(i-2)*cc(m, &
               i-1,k,2)+wa1(i-1)*cc(m,i,k,2)))+ &
               ti12* ((wa3(i-2)*cc(m,i-1,k, &
               &                        4)+wa3(i-1)*cc(m,i,k,4))- (wa2(i-2)*cc(m, &
               i-1,k,3)+wa2(i-1)*cc(m,i,k,3))))
               ch(m,ic,2,k) = (ti11* ((wa4(i-2)*cc(m,i-1,k, &
               &                         5)+wa4(i-1)*cc(m,i,k,5))- (wa1(i-2)*cc(m, &
               i-1,k,2)+wa1(i-1)*cc(m,i,k,2)))+ &
               ti12* ((wa3(i-2)*cc(m,i-1,k, &
               &                         4)+wa3(i-1)*cc(m,i,k,4))- (wa2(i-2)*cc(m, &
               i-1,k,3)+wa2(i-1)*cc(m,i,k,3)))) - &
               (cc(m,i,k,1)+tr11* ((wa1(i-2)*cc(m,i,k, &
               &                         2)-wa1(i-1)*cc(m,i-1,k, &
               &                         2))+ (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m, &
               i-1,k,5)))+tr12* ((wa2(i-2)*cc(m,i,k, &
               &                         3)-wa2(i-1)*cc(m,i-1,k, &
               &                         3))+ (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m, &
               i-1,k,4))))
               ch(m,i-1,5,k) = (cc(m,i-1,k,1)+ &
               tr12* ((wa1(i-2)*cc(m,i-1,k, &
               &                          2)+wa1(i-1)*cc(m,i,k,2))+ (wa4(i-2)*cc(m, &
               i-1,k,5)+wa4(i-1)*cc(m,i,k,5)))+ &
               tr11* ((wa2(i-2)*cc(m,i-1,k, &
               &                          3)+wa2(i-1)*cc(m,i,k,3))+ (wa3(i-2)*cc(m, &
               i-1,k,4)+wa3(i-1)*cc(m,i,k,4)))) + &
               (ti12* ((wa1(i-2)*cc(m,i,k, &
               &                          2)-wa1(i-1)*cc(m,i-1,k, &
               &                          2))- (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m, &
               i-1,k,5)))-ti11* ((wa2(i-2)*cc(m,i,k, &
               &                          3)-wa2(i-1)*cc(m,i-1,k, &
               &                          3))- (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m, &
               i-1,k,4))))
               ch(m,ic-1,4,k) = (cc(m,i-1,k,1)+ &
               tr12* ((wa1(i-2)*cc(m,i-1,k, &
               &                           2)+wa1(i-1)*cc(m,i,k, &
               &                           2))+ (wa4(i-2)*cc(m,i-1,k, &
               &                           5)+wa4(i-1)*cc(m,i,k,5)))+ &
               tr11* ((wa2(i-2)*cc(m,i-1,k, &
               &                           3)+wa2(i-1)*cc(m,i,k, &
               &                           3))+ (wa3(i-2)*cc(m,i-1,k, &
               &                           4)+wa3(i-1)*cc(m,i,k,4)))) - &
               (ti12* ((wa1(i-2)*cc(m,i,k, &
               &                           2)-wa1(i-1)*cc(m,i-1,k, &
               &                           2))- (wa4(i-2)*cc(m,i,k, &
               &                           5)-wa4(i-1)*cc(m,i-1,k,5)))- &
               ti11* ((wa2(i-2)*cc(m,i,k, &
               &                           3)-wa2(i-1)*cc(m,i-1,k, &
               &                           3))- (wa3(i-2)*cc(m,i,k, &
               &                           4)-wa3(i-1)*cc(m,i-1,k,4))))
               ch(m,i,5,k) = (cc(m,i,k,1)+tr12* &
               ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k, &
               &                        2))+ (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m, &
               i-1,k,5)))+tr11* ((wa2(i-2)*cc(m,i,k, &
               &                        3)-wa2(i-1)*cc(m,i-1,k,3))+ (wa3(i-2)*cc(m, &
               i,k,4)-wa3(i-1)*cc(m,i-1,k,4)))) + &
               (ti12* ((wa4(i-2)*cc(m,i-1,k, &
               &                        5)+wa4(i-1)*cc(m,i,k,5))- (wa1(i-2)*cc(m, &
               i-1,k,2)+wa1(i-1)*cc(m,i,k,2)))- &
               ti11* ((wa3(i-2)*cc(m,i-1,k, &
               &                        4)+wa3(i-1)*cc(m,i,k,4))- (wa2(i-2)*cc(m, &
               i-1,k,3)+wa2(i-1)*cc(m,i,k,3))))
               ch(m,ic,4,k) = (ti12* ((wa4(i-2)*cc(m,i-1,k, &
               &                         5)+wa4(i-1)*cc(m,i,k,5))- (wa1(i-2)*cc(m, &
               i-1,k,2)+wa1(i-1)*cc(m,i,k,2)))- &
               ti11* ((wa3(i-2)*cc(m,i-1,k, &
               &                         4)+wa3(i-1)*cc(m,i,k,4))- (wa2(i-2)*cc(m, &
               i-1,k,3)+wa2(i-1)*cc(m,i,k,3)))) - &
               (cc(m,i,k,1)+tr12* ((wa1(i-2)*cc(m,i,k, &
               &                         2)-wa1(i-1)*cc(m,i-1,k, &
               &                         2))+ (wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m, &
               i-1,k,5)))+tr11* ((wa2(i-2)*cc(m,i,k, &
               &                         3)-wa2(i-1)*cc(m,i-1,k, &
               &                         3))+ (wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m, &
               i-1,k,4))))
30          END DO
40       END DO
50    END DO
      RETURN
   END subroutine
   SUBROUTINE vradfg(mp,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,mdimc,wa)

!     VRFFTPK, VERSION 1, AUGUST 1985

      USE m_constants, ONLY : pimach
      DIMENSION ch(mdimc,ido,l1,ip),cc(mdimc,ido,ip,l1), &
         c1(mdimc,ido,l1,ip),c2(mdimc,idl1,ip), &
         ch2(mdimc,idl1,ip),wa(ido)

      tpi = 2.*pimach()
      arg = tpi/real(ip)
      dcp = cos(arg)
      dsp = sin(arg)
      ipph = (ip+1)/2
      ipp2 = ip + 2
      idp2 = ido + 2
      nbd = (ido-1)/2
      IF (ido == 1) GO TO 250
      DO 20 ik = 1,idl1
         DO 10 m = 1,mp
            ch2(m,ik,1) = c2(m,ik,1)
10       END DO
20    END DO
      DO 50 j = 2,ip
         DO 40 k = 1,l1
            DO 30 m = 1,mp
               ch(m,1,k,j) = c1(m,1,k,j)
30          END DO
40       END DO
50    END DO
      IF (nbd > l1) GO TO 100
      is = -ido
      DO 90 j = 2,ip
         is = is + ido
         idij = is
         DO 80 i = 3,ido,2
            idij = idij + 2
            DO 70 k = 1,l1
               DO 60 m = 1,mp
                  ch(m,i-1,k,j) = wa(idij-1)*c1(m,i-1,k,j) + &
                                  wa(idij)*c1(m,i,k,j)
                  ch(m,i,k,j) = wa(idij-1)*c1(m,i,k,j) - &
                                wa(idij)*c1(m,i-1,k,j)
60             END DO
70          END DO
80       END DO
90    END DO
      GO TO 150
100   is = -ido
      DO 140 j = 2,ip
         is = is + ido
         DO 130 k = 1,l1
            idij = is
            DO 120 i = 3,ido,2
               idij = idij + 2
               DO 110 m = 1,mp
                  ch(m,i-1,k,j) = wa(idij-1)*c1(m,i-1,k,j) + &
                                  wa(idij)*c1(m,i,k,j)
                  ch(m,i,k,j) = wa(idij-1)*c1(m,i,k,j) - &
                                wa(idij)*c1(m,i-1,k,j)
110            END DO
120         END DO
130      END DO
140   END DO
150   IF (nbd < l1) GO TO 200
      DO 190 j = 2,ipph
         jc = ipp2 - j
         DO 180 k = 1,l1
            DO 170 i = 3,ido,2
               DO 160 m = 1,mp
                  c1(m,i-1,k,j) = ch(m,i-1,k,j) + ch(m,i-1,k,jc)
                  c1(m,i-1,k,jc) = ch(m,i,k,j) - ch(m,i,k,jc)
                  c1(m,i,k,j) = ch(m,i,k,j) + ch(m,i,k,jc)
                  c1(m,i,k,jc) = ch(m,i-1,k,jc) - ch(m,i-1,k,j)
160            END DO
170         END DO
180      END DO
190   END DO
      GO TO 280
200   DO 240 j = 2,ipph
         jc = ipp2 - j
         DO 230 i = 3,ido,2
            DO 220 k = 1,l1
               DO 210 m = 1,mp
                  c1(m,i-1,k,j) = ch(m,i-1,k,j) + ch(m,i-1,k,jc)
                  c1(m,i-1,k,jc) = ch(m,i,k,j) - ch(m,i,k,jc)
                  c1(m,i,k,j) = ch(m,i,k,j) + ch(m,i,k,jc)
                  c1(m,i,k,jc) = ch(m,i-1,k,jc) - ch(m,i-1,k,j)
210            END DO
220         END DO
230      END DO
240   END DO
      GO TO 280
250   DO 270 ik = 1,idl1
         DO 260 m = 1,mp
            c2(m,ik,1) = ch2(m,ik,1)
260      END DO
270   END DO
280   DO 310 j = 2,ipph
         jc = ipp2 - j
         DO 300 k = 1,l1
            DO 290 m = 1,mp
               c1(m,1,k,j) = ch(m,1,k,j) + ch(m,1,k,jc)
               c1(m,1,k,jc) = ch(m,1,k,jc) - ch(m,1,k,j)
290         END DO
300      END DO
310   END DO

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
320         END DO
330      END DO
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
340            END DO
350         END DO
360      END DO
370   END DO
      DO 400 j = 2,ipph
         DO 390 ik = 1,idl1
            DO 380 m = 1,mp
               ch2(m,ik,1) = ch2(m,ik,1) + c2(m,ik,j)
380         END DO
390      END DO
400   END DO

      IF (ido < l1) GO TO 440
      DO 430 k = 1,l1
         DO 420 i = 1,ido
            DO 410 m = 1,mp
               cc(m,i,1,k) = ch(m,i,k,1)
410         END DO
420      END DO
430   END DO
      GO TO 480
440   DO 470 i = 1,ido
         DO 460 k = 1,l1
            DO 450 m = 1,mp
               cc(m,i,1,k) = ch(m,i,k,1)
450         END DO
460      END DO
470   END DO
480   DO 510 j = 2,ipph
         jc = ipp2 - j
         j2 = j + j
         DO 500 k = 1,l1
            DO 490 m = 1,mp
               cc(m,ido,j2-2,k) = ch(m,1,k,j)
               cc(m,1,j2-1,k) = ch(m,1,k,jc)
490         END DO
500      END DO
510   END DO
      IF (ido == 1) RETURN
      IF (nbd < l1) GO TO 560
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
520            END DO
530         END DO
540      END DO
550   END DO
      RETURN
560   DO 600 j = 2,ipph
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
570            END DO
580         END DO
590      END DO
600   END DO
      RETURN
   END subroutine
   SUBROUTINE vrftf1(m,n,c,ch,mdimc,wa,fac)

!     VRFFTPK, VERSION 1, AUGUST 1985

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
         IF (ip /= 4) GO TO 20
         ix2 = iw + ido
         ix3 = ix2 + ido
         IF (na /= 0) GO TO 10
         CALL vradf4(m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2),wa(ix3))
         GO TO 100
10       CALL vradf4(m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2),wa(ix3))
         GO TO 100
20       IF (ip /= 2) GO TO 40
         IF (na /= 0) GO TO 30
         CALL vradf2(m,ido,l1,c,ch,mdimc,wa(iw))
         GO TO 100
30       CALL vradf2(m,ido,l1,ch,c,mdimc,wa(iw))
         GO TO 100
40       IF (ip /= 3) GO TO 60
         ix2 = iw + ido
         IF (na /= 0) GO TO 50
         CALL vradf3(m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2))
         GO TO 100
50       CALL vradf3(m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2))
         GO TO 100
60       IF (ip /= 5) GO TO 80
         ix2 = iw + ido
         ix3 = ix2 + ido
         ix4 = ix3 + ido
         IF (na /= 0) GO TO 70
         CALL vradf5(m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         GO TO 100
70       CALL vradf5(m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         GO TO 100
80       IF (ido == 1) na = 1 - na
         IF (na /= 0) GO TO 90
         CALL vradfg(m,ido,ip,l1,idl1,c,c,c,ch,ch,mdimc,wa(iw))
         na = 1
         GO TO 100
90       CALL vradfg(m,ido,ip,l1,idl1,ch,ch,ch,c,c,mdimc,wa(iw))
         na = 0
100      l2 = l1
110   END DO
      scale = sqrt(1./n)
      IF (na == 1) GO TO 140
      DO 130 j = 1,n
         DO 120 i = 1,m
            c(i,j) = scale*ch(i,j)
120      END DO
130   END DO
      RETURN
140   DO 160 j = 1,n
         DO 150 i = 1,m
            c(i,j) = scale*c(i,j)
150      END DO
160   END DO
      RETURN
   END subroutine

   SUBROUTINE vrfftb(m,n,r,rt,mdimr,wsave)
!***BEGIN PROLOGUE  VRFFTB
!***DATE WRITTEN   850801   (YYMMDD)
!***REVISION DATE  900509   (YYMMDD)
!***CATEGORY NO.  J1A1
!***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM,
!             FOURIER SYNTHESIS, BACKWARD TRANSFORM, MULTIPLE SEQUENCES
!***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST)
!***PURPOSE  Backward real periodic transform, M sequences.
!***DESCRIPTION

!  Subroutine VRFFTB computes the synthesis (backward transform) of a
!  number of real periodic sequences from their Fourier coefficients.
!  Specifically, for each set of independent Fourier coefficients
!  F(K), the corresponding real periodic sequence is computed.

!  The array WSAVE which is used by subroutine VRFFTB must be
!  initialized by calling subroutine VRFFTI(N,WSAVE).

!  Input Parameters

!  M       the number of sets of coefficients.

!  N       the length of the sequences of coefficients to be
!          transformed.  The method is most efficient when N is a
!          product of small primes, however n may be any positive
!          integer.

!  R       areal two-dimensional array of size MDIMX x N containing the
!          coefficients to be transformed.  Each set of coefficients
!          F(K), K\0,1,..,N-1, is stored as a ROW of R.  Specifically,
!          the I-th set of independent Fourier coefficients is stored

!                R(I,1) = REAL( F(I,0) ),

!                R(I,2*K) = REAL( F(I,K) )

!                R(I,2*K+1) = IMAG( F(I,K) )

!                   for K = 1, 2, . . . , M-1,

!                and, when N is even,

!                R(I,N) = REAL( F(I,N/2) ).

!  RT      a real two-dimensional work array of size MDIMX x N.

!  MDIMR   the row (or first) dimension of the arrays R and RT exactly
!          as they appear in the calling program.  This parameter is
!          used to specify the variable dimension of these arrays.

!  WSAVE   a real one-dimensional work array which must be dimensioned
!          at least N+15.  The WSAVE array must be initialized by
!          calling subroutine VRFFTI.  A different WSAVE array must be
!          used for each different value of N.  This initialization does
!          not have to be repeated so long as N remains unchanged.  The
!          same WSAVE array may be used by VRFFTB and VRFFTB.

!  Output Parameters

!  R       contains M real periodic sequences corresponding to the given
!          coefficients.  Specifically, the I-th row of R contains the
!          real periodic sequence corresponding to the I-th set of
!          independent Fourier coefficients F(I,K) stored as

!               R(I,J) = X(I,J-1) ,   J = 1, 2, . . . , N, where

!               X(I,J) = SQRT(1/N)* F(I,0) + (-1)**J*F(I,N/2)
!                        + 2*SUM(K=1,M)[ REAL(F(I,2K))*COS(2K*J*PI/N)
!                        - IMAG(F(I,2K+1))*SIN(2K*J*PI/N) ]  ,

!                 when N is even, and

!               X(I,J) = SQRT(1/N)* F(I,0) +
!                        2*SUM(K=1,M)[ REAL(F(I,2K))*COS(2K*J*PI/N)
!                        - IMAG(F(I,2K+1))*SIN(2K*J*PI/N) ]  ,

!                 when N is odd.

!  WSAVE   contains results which must not be destroyed between calls
!          to VRFFTF or VRFFTB.

!  -----------------------------------------------------------------

!  NOTE  -  A call of VRFFTF followed immediately by a call of
!           of VRFFTB will return the original sequences R.  Thus,
!           VRFFTB is the correctly normalized inverse of VRFFTF.

!  -----------------------------------------------------------------

!  VRFFTB is a straightforward extension of the subprogram RFFTB to
!  handle M simultaneous sequences.  RFFTB was originally developed
!  by P. N. Swarztrauber of NCAR.

!              * * * * * * * * * * * * * * * * * * * * *
!              *                                       *
!              *         PROGRAM SPECIFICATIONS        *
!              *                                       *
!              * * * * * * * * * * * * * * * * * * * * *

!     Dimension of    R(MDIMR,N), RT(MDIMR,N), WSAVE(N+15)
!     Arguments

!     Latest          AUGUST 1, 1985
!     Revision

!     Subprograms     VRFFTI, VRFTI1, VRFFTF, VRFTF1, VRADF2, VRADF3,
!     required        VRADF4, VRADF5, VRADFG, VRFFTB, VRFTB1, VRADB2,
!                     VRADB3, VRADB4, VRADB5, VRADBG, PIMACH

!     Special         NONE
!     Conditions

!     Common          NONE
!     blocks

!     I/O             NONE

!     Precision       SINGLE

!     Specialist      ROLAND SWEET

!     Language        FORTRAN

!     History         WRITTEN BY LINDA LINDGREN AND ROLAND SWEET AT THE
!                     NATIONAL BUREAU OF STANDARDS (BOULDER).

!     Algorithm       A REAL VARIANT OF THE STOCKHAM AUTOSORT VERSION
!                     OF THE COOLEY-TUKEY FAST FOURIER TRANSFORM.

!     Portability     AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN 77.
!                     THE ONLY MACHINE DEPENDENT CONSTANT IS LOCATED IN
!                     THE FUNCTION PIMACH.

!     Required        COS,SIN
!     resident
!     Routines

!***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
!               pp. 51-83.
!***ROUTINES CALLED  VRFTB1
!***END PROLOGUE  VRFFTB

!     VRFFTPK, VERSION 1, AUGUST 1985

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n, m, mdimr
      REAL :: r(mdimr,n),rt(mdimr,n),wsave(n+15)

      IF (n == 1) RETURN
      CALL vrftb1(m,n,r,rt,mdimr,wsave(1),wsave(n+1))
      RETURN

   CONTAINS

      SUBROUTINE vradb2(mp,ido,l1,cc,ch,mdimc,wa1)

!     VRFFTPK, VERSION 1, AUGUST 1985

         IMPLICIT NONE
         INTEGER, INTENT(IN) :: mp, ido, l1, mdimc
         REAL, INTENT(IN)    :: cc(mdimc,ido,2,l1), wa1(ido)
         REAL, INTENT(INOUT) :: ch(mdimc,ido,l1,2)

         INTEGER :: i, m, k, ic, idp2

         DO 20 k = 1,l1
            DO 10 m = 1,mp
               ch(m,1,k,1) = cc(m,1,1,k) + cc(m,ido,2,k)
               ch(m,1,k,2) = cc(m,1,1,k) - cc(m,ido,2,k)
10          END DO
20       END DO
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
                  ch(m,i-1,k,2) = wa1(i-2)* (cc(m,i-1,1,k)- &
                                             cc(m,ic-1,2,k)) - wa1(i-1)* &
                                  (cc(m,i,1,k)+cc(m,ic,2,k))
                  ch(m,i,k,2) = wa1(i-2)* (cc(m,i,1,k)+cc(m,ic,2,k)) + &
                                wa1(i-1)* (cc(m,i-1,1,k)-cc(m,ic-1,2,k))
40             END DO
50          END DO
60       END DO
         IF (mod(ido,2) == 1) RETURN
70       DO 90 k = 1,l1
            DO 80 m = 1,mp
               ch(m,ido,k,1) = cc(m,ido,1,k) + cc(m,ido,1,k)
               ch(m,ido,k,2) = - (cc(m,1,2,k)+cc(m,1,2,k))
80          END DO
90       END DO
100      RETURN
      END SUBROUTINE vradb2

      SUBROUTINE vradb3(mp,ido,l1,cc,ch,mdimc,wa1,wa2)

!     VRFFTPK, VERSION 1, AUGUST 1985

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
               ch(m,1,k,2) = cc(m,1,1,k) + (2.*taur)*cc(m,ido,2,k) - &
                             (2.*taui)*cc(m,1,3,k)
               ch(m,1,k,3) = cc(m,1,1,k) + (2.*taur)*cc(m,ido,2,k) + &
               &                     2.*taui*cc(m,1,3,k)
10          END DO
20       END DO
         IF (ido == 1) RETURN
         idp2 = ido + 2
         DO 50 k = 1,l1
            DO 40 i = 3,ido,2
               ic = idp2 - i
               DO 30 m = 1,mp
                  ch(m,i-1,k,1) = cc(m,i-1,1,k) + &
                                  (cc(m,i-1,3,k)+cc(m,ic-1,2,k))
                  ch(m,i,k,1) = cc(m,i,1,k) + (cc(m,i,3,k)-cc(m,ic,2,k))
                  ch(m,i-1,k,2) = wa1(i-2)* ((cc(m,i-1,1,k)+taur* (cc(m, &
                                                                      i-1,3,k)+cc(m,ic-1,2,k)))- &
                                             (taui* (cc(m,i,3,k)+cc(m,ic,2,k)))) - &
                                  wa1(i-1)* ((cc(m,i,1,k)+taur* (cc(m,i,3, &
                                                                    k)-cc(m,ic,2,k)))+ (taui* (cc(m,i-1,3, &
                                                                                                  k)-cc(m,ic-1,2,k))))
                  ch(m,i,k,2) = wa1(i-2)* ((cc(m,i,1,k)+taur* (cc(m,i,3, &
                  k)-cc(m,ic,2,k)))+ (taui* (cc(m,i-1,3, &
                  k)-cc(m,ic-1,2,k)))) + &
                  wa1(i-1)* ((cc(m,i-1,1,k)+taur* (cc(m,i-1, &
                  &                        3,k)+cc(m,ic-1,2,k)))- &
                  (taui* (cc(m,i,3,k)+cc(m,ic,2,k))))
                  ch(m,i-1,k,3) = wa2(i-2)* ((cc(m,i-1,1,k)+taur* (cc(m, &
                                                                      i-1,3,k)+cc(m,ic-1,2,k)))+ &
                                             (taui* (cc(m,i,3,k)+cc(m,ic,2,k)))) - &
                                  wa2(i-1)* ((cc(m,i,1,k)+taur* (cc(m,i,3, &
                                                                    k)-cc(m,ic,2,k)))- (taui* (cc(m,i-1,3, &
                                                                                                  k)-cc(m,ic-1,2,k))))
                  ch(m,i,k,3) = wa2(i-2)* ((cc(m,i,1,k)+taur* (cc(m,i,3, &
                  k)-cc(m,ic,2,k)))- (taui* (cc(m,i-1,3, &
                  k)-cc(m,ic-1,2,k)))) + &
                  wa2(i-1)* ((cc(m,i-1,1,k)+taur* (cc(m,i-1, &
                  &                        3,k)+cc(m,ic-1,2,k)))+ &
                  (taui* (cc(m,i,3,k)+cc(m,ic,2,k))))
30             END DO
40          END DO
50       END DO
         RETURN
      END SUBROUTINE vradb3

      SUBROUTINE vradb4(mp,ido,l1,cc,ch,mdimc,wa1,wa2,wa3)

!     VRFFTPK, VERSION 1, AUGUST 1985

         IMPLICIT NONE
         INTEGER, INTENT(IN) :: mp, ido, l1, mdimc
         REAL, INTENT(IN)    :: cc(mdimc,ido,4,l1), wa1(ido), wa2(ido), &
                                wa3(ido)
         REAL, INTENT(INOUT) :: ch(mdimc,ido,l1,4)

         INTEGER :: k, m, i, idp2, ic
         REAL :: sqrt2

         sqrt2 = sqrt(2.)
         DO 20 k = 1,l1
            DO 10 m = 1,mp
               ch(m,1,k,3) = (cc(m,1,1,k)+cc(m,ido,4,k)) - &
                             (cc(m,ido,2,k)+cc(m,ido,2,k))
               ch(m,1,k,1) = (cc(m,1,1,k)+cc(m,ido,4,k)) + &
                             (cc(m,ido,2,k)+cc(m,ido,2,k))
               ch(m,1,k,4) = (cc(m,1,1,k)-cc(m,ido,4,k)) + &
                             (cc(m,1,3,k)+cc(m,1,3,k))
               ch(m,1,k,2) = (cc(m,1,1,k)-cc(m,ido,4,k)) - &
                             (cc(m,1,3,k)+cc(m,1,3,k))
10          END DO
20       END DO
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
                  ch(m,i-1,k,1) = (cc(m,i-1,1,k)+cc(m,ic-1,4,k)) + &
                                  (cc(m,i-1,3,k)+cc(m,ic-1,2,k))
                  ch(m,i,k,1) = (cc(m,i,1,k)-cc(m,ic,4,k)) + &
                                (cc(m,i,3,k)-cc(m,ic,2,k))
                  ch(m,i-1,k,2) = wa1(i-2)* ((cc(m,i-1,1,k)-cc(m,ic-1,4, &
                                                               k))- (cc(m,i,3,k)+cc(m,ic,2,k))) - &
                                  wa1(i-1)* ((cc(m,i,1,k)+cc(m,ic,4,k))+ &
                                             (cc(m,i-1,3,k)-cc(m,ic-1,2,k)))
                  ch(m,i,k,2) = wa1(i-2)* ((cc(m,i,1,k)+cc(m,ic,4,k))+ &
                                           (cc(m,i-1,3,k)-cc(m,ic-1,2,k))) + &
                                wa1(i-1)* ((cc(m,i-1,1,k)-cc(m,ic-1,4,k))- &
                                           (cc(m,i,3,k)+cc(m,ic,2,k)))
                  ch(m,i-1,k,3) = wa2(i-2)* ((cc(m,i-1,1,k)+cc(m,ic-1,4, &
                                                               k))- (cc(m,i-1,3,k)+cc(m,ic-1,2,k))) - &
                                  wa2(i-1)* ((cc(m,i,1,k)-cc(m,ic,4,k))- &
                                             (cc(m,i,3,k)-cc(m,ic,2,k)))
                  ch(m,i,k,3) = wa2(i-2)* ((cc(m,i,1,k)-cc(m,ic,4,k))- &
                                           (cc(m,i,3,k)-cc(m,ic,2,k))) + &
                                wa2(i-1)* ((cc(m,i-1,1,k)+cc(m,ic-1,4,k))- &
                                           (cc(m,i-1,3,k)+cc(m,ic-1,2,k)))
                  ch(m,i-1,k,4) = wa3(i-2)* ((cc(m,i-1,1,k)-cc(m,ic-1,4, &
                                                               k))+ (cc(m,i,3,k)+cc(m,ic,2,k))) - &
                                  wa3(i-1)* ((cc(m,i,1,k)+cc(m,ic,4,k))- &
                                             (cc(m,i-1,3,k)-cc(m,ic-1,2,k)))
                  ch(m,i,k,4) = wa3(i-2)* ((cc(m,i,1,k)+cc(m,ic,4,k))- &
                                           (cc(m,i-1,3,k)-cc(m,ic-1,2,k))) + &
                                wa3(i-1)* ((cc(m,i-1,1,k)-cc(m,ic-1,4,k))+ &
                                           (cc(m,i,3,k)+cc(m,ic,2,k)))
40             END DO
50          END DO
60       END DO
         IF (mod(ido,2) == 1) RETURN
70       CONTINUE
         DO 90 k = 1,l1
            DO 80 m = 1,mp
               ch(m,ido,k,1) = (cc(m,ido,1,k)+cc(m,ido,3,k)) + &
                               (cc(m,ido,1,k)+cc(m,ido,3,k))
               ch(m,ido,k,2) = sqrt2* ((cc(m,ido,1,k)-cc(m,ido,3,k))- &
                                       (cc(m,1,2,k)+cc(m,1,4,k)))
               ch(m,ido,k,3) = (cc(m,1,4,k)-cc(m,1,2,k)) + &
                               (cc(m,1,4,k)-cc(m,1,2,k))
               ch(m,ido,k,4) = -sqrt2* ((cc(m,ido,1,k)-cc(m,ido,3,k))+ &
                                        (cc(m,1,2,k)+cc(m,1,4,k)))
80          END DO
90       END DO
100      RETURN
      END SUBROUTINE vradb4

      SUBROUTINE vradb5(mp,ido,l1,cc,ch,mdimc,wa1,wa2,wa3,wa4)

!     VRFFTPK, VERSION 1, AUGUST 1985

         USE m_constants, ONLY : pimach
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: mp, ido, l1, mdimc
         REAL, INTENT(IN)    :: cc(mdimc,ido,5,l1), wa1(ido), wa2(ido), &
                                wa3(ido),wa4(ido)
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
               ch(m,1,k,1) = cc(m,1,1,k) + 2.*cc(m,ido,2,k) + &
               &                     2.*cc(m,ido,4,k)
               ch(m,1,k,2) = (cc(m,1,1,k)+tr11*2.*cc(m,ido,2,k)+ &
                              tr12*2.*cc(m,ido,4,k)) - &
                             (ti11*2.*cc(m,1,3,k)+ti12*2.*cc(m,1,5,k))
               ch(m,1,k,3) = (cc(m,1,1,k)+tr12*2.*cc(m,ido,2,k)+ &
                              tr11*2.*cc(m,ido,4,k)) - &
                             (ti12*2.*cc(m,1,3,k)-ti11*2.*cc(m,1,5,k))
               ch(m,1,k,4) = (cc(m,1,1,k)+tr12*2.*cc(m,ido,2,k)+ &
                              tr11*2.*cc(m,ido,4,k)) + &
                             (ti12*2.*cc(m,1,3,k)-ti11*2.*cc(m,1,5,k))
               ch(m,1,k,5) = (cc(m,1,1,k)+tr11*2.*cc(m,ido,2,k)+ &
                              tr12*2.*cc(m,ido,4,k)) + &
                             (ti11*2.*cc(m,1,3,k)+ti12*2.*cc(m,1,5,k))
10          END DO
20       END DO
         IF (ido == 1) RETURN
         idp2 = ido + 2
         DO 50 k = 1,l1
            DO 40 i = 3,ido,2
               ic = idp2 - i
               DO 30 m = 1,mp
                  ch(m,i-1,k,1) = cc(m,i-1,1,k) + &
                                  (cc(m,i-1,3,k)+cc(m,ic-1,2,k)) + &
                                  (cc(m,i-1,5,k)+cc(m,ic-1,4,k))
                  ch(m,i,k,1) = cc(m,i,1,k) + (cc(m,i,3,k)-cc(m,ic,2,k)) + &
                                (cc(m,i,5,k)-cc(m,ic,4,k))
                  ch(m,i-1,k,2) = wa1(i-2)* ((cc(m,i-1,1,k)+tr11* (cc(m, &
                  i-1,3,k)+cc(m,ic-1,2,k))+tr12* (cc(m,i-1, &
                  &                          5,k)+cc(m,ic-1,4,k)))- &
                  (ti11* (cc(m,i,3,k)+cc(m,ic,2, &
                  k))+ti12* (cc(m,i,5,k)+cc(m,ic,4,k)))) - &
                  wa1(i-1)* ((cc(m,i,1,k)+tr11* (cc(m,i,3, &
                  k)-cc(m,ic,2,k))+tr12* (cc(m,i,5,k)-cc(m, &
                  ic,4,k)))+ (ti11* (cc(m,i-1,3,k)-cc(m, &
                  ic-1,2,k))+ti12* (cc(m,i-1,5,k)-cc(m, &
                  ic-1,4,k))))
                  ch(m,i,k,2) = wa1(i-2)* ((cc(m,i,1,k)+tr11* (cc(m,i,3, &
                  k)-cc(m,ic,2,k))+tr12* (cc(m,i,5,k)-cc(m, &
                  ic,4,k)))+ (ti11* (cc(m,i-1,3,k)-cc(m,ic-1, &
                  &                        2,k))+ti12* (cc(m,i-1,5,k)-cc(m,ic-1,4, &
                  k)))) + wa1(i-1)* ((cc(m,i-1,1, &
                  k)+tr11* (cc(m,i-1,3,k)+cc(m,ic-1,2, &
                  k))+tr12* (cc(m,i-1,5,k)+cc(m,ic-1,4,k)))- &
                  (ti11* (cc(m,i,3,k)+cc(m,ic,2, &
                  k))+ti12* (cc(m,i,5,k)+cc(m,ic,4,k))))
                  ch(m,i-1,k,3) = wa2(i-2)* ((cc(m,i-1,1,k)+tr12* (cc(m, &
                  i-1,3,k)+cc(m,ic-1,2,k))+tr11* (cc(m,i-1, &
                  &                          5,k)+cc(m,ic-1,4,k)))- &
                  (ti12* (cc(m,i,3,k)+cc(m,ic,2, &
                  k))-ti11* (cc(m,i,5,k)+cc(m,ic,4,k)))) - &
                  wa2(i-1)* ((cc(m,i,1,k)+tr12* (cc(m,i,3, &
                  k)-cc(m,ic,2,k))+tr11* (cc(m,i,5,k)-cc(m, &
                  ic,4,k)))+ (ti12* (cc(m,i-1,3,k)-cc(m, &
                  ic-1,2,k))-ti11* (cc(m,i-1,5,k)-cc(m, &
                  ic-1,4,k))))
                  ch(m,i,k,3) = wa2(i-2)* ((cc(m,i,1,k)+tr12* (cc(m,i,3, &
                  k)-cc(m,ic,2,k))+tr11* (cc(m,i,5,k)-cc(m, &
                  ic,4,k)))+ (ti12* (cc(m,i-1,3,k)-cc(m,ic-1, &
                  &                        2,k))-ti11* (cc(m,i-1,5,k)-cc(m,ic-1,4, &
                  k)))) + wa2(i-1)* ((cc(m,i-1,1, &
                  k)+tr12* (cc(m,i-1,3,k)+cc(m,ic-1,2, &
                  k))+tr11* (cc(m,i-1,5,k)+cc(m,ic-1,4,k)))- &
                  (ti12* (cc(m,i,3,k)+cc(m,ic,2, &
                  k))-ti11* (cc(m,i,5,k)+cc(m,ic,4,k))))
                  ch(m,i-1,k,4) = wa3(i-2)* ((cc(m,i-1,1,k)+tr12* (cc(m, &
                  i-1,3,k)+cc(m,ic-1,2,k))+tr11* (cc(m,i-1, &
                  &                          5,k)+cc(m,ic-1,4,k)))+ &
                  (ti12* (cc(m,i,3,k)+cc(m,ic,2, &
                  k))-ti11* (cc(m,i,5,k)+cc(m,ic,4,k)))) - &
                  wa3(i-1)* ((cc(m,i,1,k)+tr12* (cc(m,i,3, &
                  k)-cc(m,ic,2,k))+tr11* (cc(m,i,5,k)-cc(m, &
                  ic,4,k)))- (ti12* (cc(m,i-1,3,k)-cc(m, &
                  ic-1,2,k))-ti11* (cc(m,i-1,5,k)-cc(m, &
                  ic-1,4,k))))
                  ch(m,i,k,4) = wa3(i-2)* ((cc(m,i,1,k)+tr12* (cc(m,i,3, &
                  k)-cc(m,ic,2,k))+tr11* (cc(m,i,5,k)-cc(m, &
                  ic,4,k)))- (ti12* (cc(m,i-1,3,k)-cc(m,ic-1, &
                  &                        2,k))-ti11* (cc(m,i-1,5,k)-cc(m,ic-1,4, &
                  k)))) + wa3(i-1)* ((cc(m,i-1,1, &
                  k)+tr12* (cc(m,i-1,3,k)+cc(m,ic-1,2, &
                  k))+tr11* (cc(m,i-1,5,k)+cc(m,ic-1,4,k)))+ &
                  (ti12* (cc(m,i,3,k)+cc(m,ic,2, &
                  k))-ti11* (cc(m,i,5,k)+cc(m,ic,4,k))))
                  ch(m,i-1,k,5) = wa4(i-2)* ((cc(m,i-1,1,k)+tr11* (cc(m, &
                  i-1,3,k)+cc(m,ic-1,2,k))+tr12* (cc(m,i-1, &
                  &                          5,k)+cc(m,ic-1,4,k)))+ &
                  (ti11* (cc(m,i,3,k)+cc(m,ic,2, &
                  k))+ti12* (cc(m,i,5,k)+cc(m,ic,4,k)))) - &
                  wa4(i-1)* ((cc(m,i,1,k)+tr11* (cc(m,i,3, &
                  k)-cc(m,ic,2,k))+tr12* (cc(m,i,5,k)-cc(m, &
                  ic,4,k)))- (ti11* (cc(m,i-1,3,k)-cc(m, &
                  ic-1,2,k))+ti12* (cc(m,i-1,5,k)-cc(m, &
                  ic-1,4,k))))
                  ch(m,i,k,5) = wa4(i-2)* ((cc(m,i,1,k)+tr11* (cc(m,i,3, &
                  k)-cc(m,ic,2,k))+tr12* (cc(m,i,5,k)-cc(m, &
                  ic,4,k)))- (ti11* (cc(m,i-1,3,k)-cc(m,ic-1, &
                  &                        2,k))+ti12* (cc(m,i-1,5,k)-cc(m,ic-1,4, &
                  k)))) + wa4(i-1)* ((cc(m,i-1,1, &
                  k)+tr11* (cc(m,i-1,3,k)+cc(m,ic-1,2, &
                  k))+tr12* (cc(m,i-1,5,k)+cc(m,ic-1,4,k)))+ &
                  (ti11* (cc(m,i,3,k)+cc(m,ic,2, &
                  k))+ti12* (cc(m,i,5,k)+cc(m,ic,4,k))))
30             END DO
40          END DO
50       END DO
         RETURN
      END SUBROUTINE vradb5

      SUBROUTINE vradbg(mp,ido,ip,l1,idl1,cc,c1,c2,ch,ch2, &
                        mdimc,wa)

!     VRFFTPK, VERSION 1, AUGUST 1985

         USE m_constants, ONLY : pimach
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: mp, ido, ip, l1, idl1, mdimc
         REAL, INTENT(IN)    :: cc(mdimc,ido,ip,l1), wa(ido)
         REAL, INTENT(INOUT) :: ch(mdimc,ido,l1,ip), ch2(mdimc,idl1,ip), &
                                c1(mdimc,ido,l1,ip), c2(mdimc,idl1,ip)

         INTEGER :: i, j, j2, k, l, m, ik, is, idp2, nbd, ipp2, ipph, &
                    ic, jc, lc, idij
         REAL :: ai1, ar1, ar2, ai2, ar1h, ar2h, tpi, arg, dcp, dsp, &
                 dc2, ds2

         tpi = 2.*pimach()
         arg = tpi/real(ip)
         dcp = cos(arg)
         dsp = sin(arg)
         idp2 = ido + 2
         nbd = (ido-1)/2
         ipp2 = ip + 2
         ipph = (ip+1)/2
         IF (ido < l1) GO TO 40
         DO 30 k = 1,l1
            DO 20 i = 1,ido
               DO 10 m = 1,mp
                  ch(m,i,k,1) = cc(m,i,1,k)
10             END DO
20          END DO
30       END DO
         GO TO 80
40       DO 70 i = 1,ido
            DO 60 k = 1,l1
               DO 50 m = 1,mp
                  ch(m,i,k,1) = cc(m,i,1,k)
50             END DO
60          END DO
70       END DO
80       DO 110 j = 2,ipph
            jc = ipp2 - j
            j2 = j + j
            DO 100 k = 1,l1
               DO 90 m = 1,mp
                  ch(m,1,k,j) = cc(m,ido,j2-2,k) + cc(m,ido,j2-2,k)
                  ch(m,1,k,jc) = cc(m,1,j2-1,k) + cc(m,1,j2-1,k)
90             END DO
100         END DO
110      END DO
         IF (ido == 1) GO TO 210
         IF (nbd < l1) GO TO 160
         DO 150 j = 2,ipph
            jc = ipp2 - j
            DO 140 k = 1,l1
               DO 130 i = 3,ido,2
                  ic = idp2 - i
                  DO 120 m = 1,mp
                     ch(m,i-1,k,j) = cc(m,i-1,2*j-1,k) + cc(m,ic-1,2*j-2,k)
                     ch(m,i-1,k,jc) = cc(m,i-1,2*j-1,k) - &
                                      cc(m,ic-1,2*j-2,k)
                     ch(m,i,k,j) = cc(m,i,2*j-1,k) - cc(m,ic,2*j-2,k)
                     ch(m,i,k,jc) = cc(m,i,2*j-1,k) + cc(m,ic,2*j-2,k)
120               END DO
130            END DO
140         END DO
150      END DO
         GO TO 210
160      DO 200 j = 2,ipph
            jc = ipp2 - j
            DO 190 i = 3,ido,2
               ic = idp2 - i
               DO 180 k = 1,l1
                  DO 170 m = 1,mp
                     ch(m,i-1,k,j) = cc(m,i-1,2*j-1,k) + cc(m,ic-1,2*j-2,k)
                     ch(m,i-1,k,jc) = cc(m,i-1,2*j-1,k) - &
                                      cc(m,ic-1,2*j-2,k)
                     ch(m,i,k,j) = cc(m,i,2*j-1,k) - cc(m,ic,2*j-2,k)
                     ch(m,i,k,jc) = cc(m,i,2*j-1,k) + cc(m,ic,2*j-2,k)
170               END DO
180            END DO
190         END DO
200      END DO
210      ar1 = 1.
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
220            END DO
230         END DO
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
240               END DO
250            END DO
260         END DO
270      END DO
         DO 300 j = 2,ipph
            DO 290 ik = 1,idl1
               DO 280 m = 1,mp
                  ch2(m,ik,1) = ch2(m,ik,1) + ch2(m,ik,j)
280            END DO
290         END DO
300      END DO
         DO 330 j = 2,ipph
            jc = ipp2 - j
            DO 320 k = 1,l1
               DO 310 m = 1,mp
                  ch(m,1,k,j) = c1(m,1,k,j) - c1(m,1,k,jc)
                  ch(m,1,k,jc) = c1(m,1,k,j) + c1(m,1,k,jc)
310            END DO
320         END DO
330      END DO
         IF (ido == 1) GO TO 430
         IF (nbd < l1) GO TO 380
         DO 370 j = 2,ipph
            jc = ipp2 - j
            DO 360 k = 1,l1
               DO 350 i = 3,ido,2
                  DO 340 m = 1,mp
                     ch(m,i-1,k,j) = c1(m,i-1,k,j) - c1(m,i,k,jc)
                     ch(m,i-1,k,jc) = c1(m,i-1,k,j) + c1(m,i,k,jc)
                     ch(m,i,k,j) = c1(m,i,k,j) + c1(m,i-1,k,jc)
                     ch(m,i,k,jc) = c1(m,i,k,j) - c1(m,i-1,k,jc)
340               END DO
350            END DO
360         END DO
370      END DO
         GO TO 430
380      DO 420 j = 2,ipph
            jc = ipp2 - j
            DO 410 i = 3,ido,2
               DO 400 k = 1,l1
                  DO 390 m = 1,mp
                     ch(m,i-1,k,j) = c1(m,i-1,k,j) - c1(m,i,k,jc)
                     ch(m,i-1,k,jc) = c1(m,i-1,k,j) + c1(m,i,k,jc)
                     ch(m,i,k,j) = c1(m,i,k,j) + c1(m,i-1,k,jc)
                     ch(m,i,k,jc) = c1(m,i,k,j) - c1(m,i-1,k,jc)
390               END DO
400            END DO
410         END DO
420      END DO
430      CONTINUE
         IF (ido == 1) RETURN
         DO 450 ik = 1,idl1
            DO 440 m = 1,mp
               c2(m,ik,1) = ch2(m,ik,1)
440         END DO
450      END DO
         DO 480 j = 2,ip
            DO 470 k = 1,l1
               DO 460 m = 1,mp
                  c1(m,1,k,j) = ch(m,1,k,j)
460            END DO
470         END DO
480      END DO
         IF (nbd > l1) GO TO 530
         is = -ido
         DO 520 j = 2,ip
            is = is + ido
            idij = is
            DO 510 i = 3,ido,2
               idij = idij + 2
               DO 500 k = 1,l1
                  DO 490 m = 1,mp
                     c1(m,i-1,k,j) = wa(idij-1)*ch(m,i-1,k,j) - &
                                     wa(idij)*ch(m,i,k,j)
                     c1(m,i,k,j) = wa(idij-1)*ch(m,i,k,j) + &
                                   wa(idij)*ch(m,i-1,k,j)
490               END DO
500            END DO
510         END DO
520      END DO
         GO TO 580
530      is = -ido
         DO 570 j = 2,ip
            is = is + ido
            DO 560 k = 1,l1
               idij = is
               DO 550 i = 3,ido,2
                  idij = idij + 2
                  DO 540 m = 1,mp
                     c1(m,i-1,k,j) = wa(idij-1)*ch(m,i-1,k,j) - &
                                     wa(idij)*ch(m,i,k,j)
                     c1(m,i,k,j) = wa(idij-1)*ch(m,i,k,j) + &
                                   wa(idij)*ch(m,i-1,k,j)
540               END DO
550            END DO
560         END DO
570      END DO
580      RETURN
      END SUBROUTINE vradbg

      SUBROUTINE vrftb1(m,n,c,ch,mdimc,wa,fac)

!     VRFFTPK, VERSION 1, AUGUST 1985

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
            IF (ip /= 4) GO TO 30
            ix2 = iw + ido
            ix3 = ix2 + ido
            IF (na /= 0) GO TO 10
            CALL vradb4(m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2),wa(ix3))
            GO TO 20
10          CALL vradb4(m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2),wa(ix3))
20          na = 1 - na
            GO TO 150
30          IF (ip /= 2) GO TO 60
            IF (na /= 0) GO TO 40
            CALL vradb2(m,ido,l1,c,ch,mdimc,wa(iw))
            GO TO 50
40          CALL vradb2(m,ido,l1,ch,c,mdimc,wa(iw))
50          na = 1 - na
            GO TO 150
60          IF (ip /= 3) GO TO 90
            ix2 = iw + ido
            IF (na /= 0) GO TO 70
            CALL vradb3(m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2))
            GO TO 80
70          CALL vradb3(m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2))
80          na = 1 - na
            GO TO 150
90          IF (ip /= 5) GO TO 120
            ix2 = iw + ido
            ix3 = ix2 + ido
            ix4 = ix3 + ido
            IF (na /= 0) GO TO 100
            CALL vradb5(m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2),wa(ix3),wa(ix4))
            GO TO 110
100         CALL vradb5(m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2),wa(ix3),wa(ix4))
110         na = 1 - na
            GO TO 150
120         IF (na /= 0) GO TO 130
            CALL vradbg(m,ido,ip,l1,idl1,c,c,c,ch,ch,mdimc,wa(iw))
            GO TO 140
130         CALL vradbg(m,ido,ip,l1,idl1,ch,ch,ch,c,c,mdimc,wa(iw))
140         IF (ido == 1) na = 1 - na
150         l1 = l2
            iw = iw + (ip-1)*ido
160      END DO
         scale = sqrt(1./n)
         IF (na == 0) GO TO 190
         DO 180 j = 1,n
            DO 170 i = 1,m
               c(i,j) = scale*ch(i,j)
170         END DO
180      END DO
         RETURN
190      DO 210 j = 1,n
            DO 200 i = 1,m
               c(i,j) = scale*c(i,j)
200         END DO
210      END DO
         RETURN
      END SUBROUTINE vrftb1

   END SUBROUTINE vrfftb

END MODULE m_rfft
