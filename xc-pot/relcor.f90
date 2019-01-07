MODULE m_relcor
!************************************************************************

! calculate the relativistic  corrections for exchange as formulated by
!              A.H. MacDonald and S.H. Vosko, J. Phys. C12, 2977 (1979)
! either phi (for xc-energy) or psi (for xc-potential) are calculated
!             (if l_psi=.true. we call from vxc.. and psi is evaluated)

!************************************************************************
CONTAINS
   SUBROUTINE relcor( &
      mgrid,ngrid,jspins,krla,l_psi,rh, &
      phsi)
      IMPLICIT NONE

!     .. Scalar Arguments ..
      INTEGER, INTENT (IN)  :: mgrid,krla,ngrid,jspins
      LOGICAL, INTENT (IN)  :: l_psi

!     .. Array Arguments ..
      REAL,    INTENT (IN)  :: rh(mgrid,jspins)
      REAL,    INTENT (OUT) :: phsi(ngrid)

!     .. Local Parameters ..
      REAL, PARAMETER :: betac = 2.2576918e-2  ! alpha * (3 * pi)^(1/3)
      REAL, PARAMETER :: d_15 = 1.e-15 , d_3 = 1.e-3
      REAL, PARAMETER :: one = 1.0 , three = 3.0 , half = 0.5
      REAL, PARAMETER :: thrhalf = three * half , thrd = one/three
      REAL, PARAMETER :: bs1 = 0.75 , bs2 = 0.45 ,  bf2 = 0.4
      REAL, PARAMETER :: bf1 = 2*thrd

!     .. Locals ..
      INTEGER :: ir
      REAL :: beta              ! Fermi velocity devided by speed of light
      REAL :: rho,eta,xi,betasq

      INTRINSIC max,sqrt

      IF (krla == 1) THEN      !  evaluate relativistic corrections for exchange

         DO ir = 1,ngrid
            IF (jspins == 1) THEN
               rho = max(d_15,rh(ir,1))
            ELSE
               rho = max(d_15,rh(ir,1))+max(d_15,rh(ir,jspins))
            ENDIF
            beta = betac * rho**thrd
            betasq = beta*beta
            eta = sqrt(one+betasq)
            xi = beta + eta

            !----->    If beta.LT.10**(-3) use taylor series of psi,phi with respect to
            !          beta, because of accuracy considerations. Taylor series
            !          implemented is exact up to  beta**5 (see notes S.B.)

            IF (l_psi) THEN
               IF (beta < d_3) THEN
                  phsi(ir) = one - betasq + bs1*beta*betasq - bs2*betasq**2
               ELSE
                  phsi(ir) = half* (-one+three*alog(xi)/ (beta*eta))
               END IF
            ELSE
               IF (beta < d_3) THEN
                  phsi(ir) = one - bf1*betasq + bf2*betasq*betasq
               ELSE
                  phsi(ir) = one - thrhalf*((beta*eta-alog(xi))/betasq)**2
               END IF
            ENDIF
         ENDDO

      ELSE
         DO ir = 1,ngrid
            phsi(ir) = one
         ENDDO
      ENDIF

      RETURN
   END SUBROUTINE relcor
END MODULE m_relcor
