c............................................................rinvgj
      SUBROUTINE rinvgj(ainv,a,arraydim,n)
C   ********************************************************************
C   *                                                                  *
C   *                      AINV = A**(-1)                              *
C   *                                                                  *
C   *  invert A using the GAUSS-JORDAN - algorithm                     *
C   *  the 1- matrix is not set up and use is made of its structure    *
C   *                                                                  *
C   *                    REAL*8 VERSION                                *
C   *                                                                  *
C   ********************************************************************

      IMPLICIT NONE
c                                                        scan columns
C     .. Scalar Arguments ..
      INTEGER arraydim,n
C     ..
C     .. Array Arguments ..
      REAL a(arraydim,arraydim),ainv(arraydim,arraydim)
C     ..
C     .. Local Scalars ..
      REAL t,t1
      INTEGER icol,l,ll
C     ..
      DO 60 icol = 1,n
c                                               make A(ICOL,ICOL) = 1
         t1 = 1.0/a(icol,icol)
         DO 10 l = (icol+1),n
            a(icol,l) = a(icol,l)*t1
   10    CONTINUE
         DO 20 l = 1, (icol-1)
            ainv(icol,l) = ainv(icol,l)*t1
   20    CONTINUE
         ainv(icol,icol) = t1
c                                    make A(LL,ICOL) = 0 for LL<>ICOL
         DO 50 ll = 1,n
            IF (ll.NE.icol) THEN
               t = a(ll,icol)
               DO 30 l = (icol+1),n
                  a(ll,l) = a(ll,l) - a(icol,l)*t
   30          CONTINUE
               DO 40 l = 1, (icol-1)
                  ainv(ll,l) = ainv(ll,l) - ainv(icol,l)*t
   40          CONTINUE
               ainv(ll,icol) = -t1*t
            END IF
   50    CONTINUE
   60 CONTINUE
      RETURN
      END
