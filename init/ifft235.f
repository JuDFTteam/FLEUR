!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_ifft
      use m_juDFT
      CONTAINS
      INTEGER FUNCTION  ifft235 (iofile,ksfft,n,gmaxp)
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C   this function checks wether n can be expressed as :
C   n =  (2**P) * (3**Q) * (5**R) to match withe the MFFT - routines
C   used in this program. If close to n there is a number n' which
c   is solely expressable as n'= 2**p then this number sis choosen
c   and is outputted by the function. If n is not expressable as
C   n =  (2**P) * (3**Q) * (5**R) a number n' is choosen close to ns
c   is selected which fullfilles these requirements.
c
c
c                            Stefan Bl"ugel , kfa, Oct. 1993
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IMPLICIT NONE
c
c----> declaration part
c
      INTEGER iofile,ksfft,n
      REAL    gmaxp
c
c----> local variabels
c
      INTEGER np,nn,nl,iexp,is,itrial
      REAL rl
      REAL fac_1,two,fac_2,fac_3,one,fac_4
      PARAMETER ( two=2.0e0,fac_1=1.001e0,fac_2=1.7e0,fac_3=1.95e0 )
      PARAMETER ( one=1.0e0,fac_4=0.03e0 )
      INTRINSIC abs,log,mod,real
c
c
      ifft235=0
      IF ( n.eq.0 ) THEN
         WRITE (iofile,*) 'n should not be zero'
          CALL juDFT_error("n should not be zero",calledby="ifft235")
      ENDIF
c
c====> RADIX 2 CALCULATION ONLY
c
      IF ( ksfft .EQ. 0 ) THEN
          rl   = log(real(n))/log(two)
          nl   = rl
          IF ( abs( n - 2 ** nl ) .GT. abs(n - 2 ** (nl+1) ) ) THEN
               np = nl + 1
          ELSE
               np = nl
          ENDIF
          IF ( ( np .eq. nl ) .AND. ( gmaxp * real(nl)/rl .lt. fac_1 ) )
     >         np = nl + 1
          ifft235 = 2 ** np
          RETURN
       ENDIF
c
c====> RADIX 2 , 3, 5 CALCULATION
c
c
c----> since gmaxp is large enough, try also smaller n
c
      IF ( gmaxp .GT. fac_2 .AND. gmaxp .LT. fac_3 ) THEN
           is = -1
      ELSE
           is =  1
      ENDIF
      np   = n
      nn   = n
c
      DO 20 itrial = 1 , 200
c
c----> check whether there is a number close by, which is 2**NL
c      ( FFT very fast )
c
          rl   = log(real(nn))/log(two)
          nl   = rl
          IF (( abs( real(nl)/rl - one ) .LT. fac_4 ).AND.
     +        ( 2**nl.GE.n ) ) THEN
             IF (gmaxp * real(nl)/rl .GT. one) ifft235 = 2 ** nl
             RETURN
          ELSE IF (( abs( real(nl+1)/rl - one ) .lt. fac_4 ).AND.
     +             ( 2**(nl+1).GE.n ) ) THEN
             ifft235 = 2 ** (nl + 1)
             RETURN
          ELSE
c
c----> no , no binary number is arround, check wether number can be
c      divided by 2,3 or 5
c
             DO 10 iexp = 1 , nl+1
                IF (mod(nn,2) .EQ. 0) THEN
                   nn = nn / 2
                ELSE IF (mod(nn,3) .EQ. 0 ) THEN
                   nn = nn / 3
                ELSE IF (mod(nn,5) .EQ. 0 ) THEN
                   nn = nn / 5
                ENDIF
                IF ( nn.eq.2 .OR. nn.eq.3 .OR. nn.eq.5 ) THEN
c---> o.k.
                   ifft235 = np
                   RETURN
                ENDIF
  10         ENDDO
c
c----> no , n  cannot be expressed as (2**P) * (3**Q) * (5**R)
c      change number
             nn = n + (is ** (itrial)) * ((itrial + 1)/2)
             np = n + (is ** (itrial)) * ((itrial + 1)/2)
          ENDIF
  20    ENDDO
c
      WRITE (iofile,'('' to few trials '', i3)') itrial
      STOP
c
      END FUNCTION ifft235
C-----------------------------------------------------------------------
      INTEGER FUNCTION i2357(ii)

!
! simple setup to determine fft-length for ESSL calls
!
      IMPLICIT NONE
      INTEGER, INTENT (IN)  :: ii
      INTEGER i,j,k,h,m,n,nn

      nn = 12582912
      DO m = 0,1 
        DO k = 0,1 
          DO j = 0,1 
            DO i = 0,1 
              DO h = 1,25
                n = 2**h * 3**i * 5**j * 7**k * 11**m
                IF ( (n.GT.ii).AND.(n.LT.nn) ) nn = n
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      i2357 = nn

      END FUNCTION i2357
 
      END MODULE m_ifft
