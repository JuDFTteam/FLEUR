MODULE m_rinvgj

   CONTAINS

   SUBROUTINE rinvgj(ainv,aMat,arraydim,n)
!   ********************************************************************
!   *                                                                  *
!   *                      AINV = A**(-1)                              *
!   *                                                                  *
!   *  invert aMat using the GAUSS-JORDAN - algorithm                     *
!   *  the 1- matrix is not set up and use is made of its structure    *
!   *                                                                  *
!   *                    REAL*8 VERSION                                *
!   *                                                                  *
!   ********************************************************************

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: arraydim,n

      REAL, INTENT (INOUT) :: aMat(arraydim,arraydim)
      REAL, INTENT (INOUT) :: ainv(arraydim,arraydim)

      REAL t, t1
      INTEGER icol, l, ll

      DO icol = 1,n
!                                               make A(ICOL,ICOL) = 1
         t1 = 1.0 / aMat(icol,icol)
         DO l = (icol+1),n
            aMat(icol,l) = aMat(icol,l)*t1
         END DO
         DO l = 1, (icol-1)
            ainv(icol,l) = ainv(icol,l)*t1
         END DO
         ainv(icol,icol) = t1
!                                    make A(LL,ICOL) = 0 for LL<>ICOL
         DO ll = 1,n
            IF (ll.NE.icol) THEN
               t = aMat(ll,icol)
               DO l = (icol+1),n
                  aMat(ll,l) = aMat(ll,l) - aMat(icol,l)*t
               END DO
               DO l = 1, (icol-1)
                  ainv(ll,l) = ainv(ll,l) - ainv(icol,l)*t
               END DO
               ainv(ll,icol) = -t1*t
            END IF
         END DO
      END DO

   END SUBROUTINE rinvgj

END MODULE m_rinvgj
