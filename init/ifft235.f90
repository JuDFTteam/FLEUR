!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_ifft
   use m_juDFT
CONTAINS
   INTEGER FUNCTION ifft235(ksfft, n, gmaxp)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!   this function checks wether n can be expressed as :
!   n =  (2**P) * (3**Q) * (5**R) to match withe the MFFT - routines
!   used in this program. If close to n there is a number n' which
!   is solely expressable as n'= 2**p then this number sis choosen
!   and is outputted by the function. If n is not expressable as
!   n =  (2**P) * (3**Q) * (5**R) a number n' is choosen close to ns
!   is selected which fullfilles these requirements.

!                            Stefan Bl"ugel , kfa, Oct. 1993
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IMPLICIT NONE

!----> declaration part

      INTEGER :: ksfft, n
      REAL ::    gmaxp

!----> local variabels

      INTEGER :: np, nn, nl, iexp, is, itrial
      REAL :: rl
      REAL :: fac_1, two, fac_2, fac_3, one, fac_4
      PARAMETER(two=2.0e0, fac_1=1.001e0, fac_2=1.7e0, fac_3=1.95e0)
      PARAMETER(one=1.0e0, fac_4=0.03e0)
      INTRINSIC abs, log, mod, real

      ifft235 = 0
      IF (n == 0) THEN
         CALL juDFT_error("n should not be zero", calledby="ifft235")
      ENDIF

!====> RADIX 2 CALCULATION ONLY

      IF (ksfft == 0) THEN
         rl = log(real(n))/log(two)
         nl = rl
         IF (abs(n - 2**nl) > abs(n - 2**(nl + 1))) THEN
            np = nl + 1
         ELSE
            np = nl
         ENDIF
         IF ((np == nl) .AND. (gmaxp*real(nl)/rl < fac_1)) &
            np = nl + 1
         ifft235 = 2**np
         RETURN
      ENDIF

!====> RADIX 2 , 3, 5 CALCULATION

!----> since gmaxp is large enough, try also smaller n

      IF (gmaxp > fac_2 .AND. gmaxp < fac_3) THEN
         is = -1
      ELSE
         is = 1
      ENDIF
      np = n
      nn = n

      DO itrial = 1, 200

         !----> check whether there is a number close by, which is 2**NL
         !      ( FFT very fast )

         rl = log(real(nn))/log(two)
         nl = rl
         IF ((abs(real(nl)/rl - one) < fac_4) .AND. &
             (2**nl >= n)) THEN
            IF (gmaxp*real(nl)/rl > one) ifft235 = 2**nl
            RETURN
         ELSE IF ((abs(real(nl + 1)/rl - one) < fac_4) .AND. &
                  (2**(nl + 1) >= n)) THEN
            ifft235 = 2**(nl + 1)
            RETURN
         ELSE

            !----> no , no binary number is arround, check wether number can be
            !      divided by 2,3 or 5

            DO iexp = 1, nl + 1
               IF (mod(nn, 2) == 0) THEN
                  nn = nn/2
               ELSE IF (mod(nn, 3) == 0) THEN
                  nn = nn/3
               ELSE IF (mod(nn, 5) == 0) THEN
                  nn = nn/5
               ENDIF
               IF (nn == 2 .OR. nn == 3 .OR. nn == 5) THEN
                  !---> o.k.
                  ifft235 = np
                  RETURN
               ENDIF
            ENDDO

            !----> no , n  cannot be expressed as (2**P) * (3**Q) * (5**R)
            !      change number
            nn = n + (is**(itrial))*((itrial + 1)/2)
            np = n + (is**(itrial))*((itrial + 1)/2)
         ENDIF
      ENDDO

      call juDFT_error(' to few trials '//int2str(itrial))

   END FUNCTION ifft235
!-----------------------------------------------------------------------
   INTEGER FUNCTION i2357(ii)

! simple setup to determine fft-length for ESSL calls

      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: ii
      INTEGER :: i, j, k, h, m, n, nn

      nn = 12582912
      DO m = 0, 1
         DO k = 0, 1
            DO j = 0, 1
               DO i = 0, 1
                  DO h = 1, 25
                     n = 2**h*3**i*5**j*7**k*11**m
                     IF ((n >= ii) .AND. (n < nn)) nn = n
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      i2357 = nn

   END FUNCTION i2357


   function next_optimal_fft_size(n) result(fft_size)
      implicit none 
      integer, intent(in) :: n
      integer             :: fft_size 

      fft_size = n
      do while(.not. is_optimal_fft_size(fft_size))
         fft_size = fft_size + 1
      enddo
   end function next_optimal_fft_size

   logical function is_optimal_fft_size(n) 
      implicit none 
      integer, intent(in) :: n 
      integer :: i, cnt, divisor
      integer, parameter :: primes(4) = [2,3,5,7]

      i = n 
      do cnt =1,4 
         divisor = primes(cnt)
         do while(mod(i,divisor) == 0)
            i = i / divisor
         enddo
      enddo 

      is_optimal_fft_size = (i == 1)
   end function is_optimal_fft_size
END MODULE m_ifft
