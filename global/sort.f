      MODULE m_sort
      CONTAINS
      SUBROUTINE sort(
     >                n,sk,
     <                js)
c********************************************************************
c     heapsort routine
c     input:   n     = number of objects
c              sk    = array of objects to be sorted
c     output:  js(i) = index of i'th smallest object       m.w.
c     modified to avoid numerical uncertainties for equal-lentgh
c     vectors and to make sorting machine independent  aug. 90, r.p.
c********************************************************************

      IMPLICIT NONE
c
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: n
C     ..
C     .. Array Arguments ..
      REAL,    INTENT (IN) :: sk(n)
      INTEGER, INTENT (OUT) :: js(n)
C     ..
C     .. Local Scalars ..
      REAL eps,q
      INTEGER i,ind,ir,j,l
C     ..
C     .. Data statements ..
      DATA eps/1.e-10/
C     ..
c
      IF (n == 0) RETURN ! Nothing to do
      IF (n == 1) THEN   ! Not much to do
          js(1) = 1
          RETURN
      END IF
      
      DO i = 1,n
         js(i) = i
      ENDDO
c
      l = n/2 + 1
      ir = n
   20 CONTINUE
      IF (l.GT.1) THEN
         l = l - 1
         ind = js(l)
         q = sk(ind)
      ELSE
         ind = js(ir)
         q = sk(ind)
         js(ir) = js(1)
         ir = ir - 1
         IF (ir.EQ.1) THEN
            js(1) = ind
            RETURN
         END IF
      END IF
      i = l
      j = l + l
   30 IF (j.LE.ir) THEN
         IF (j.LT.ir) THEN
c           if(sk(js(j)).lt.sk(js(j+1))) j=j+1
            IF ((sk(js(j+1))-sk(js(j))).GT.eps) j = j + 1
         END IF
c        if(q.lt.sk(js(j))) then
         IF ((sk(js(j))-q).GT.eps) THEN
            js(i) = js(j)
            i = j
            j = j + j
         ELSE
            j = ir + 1
         END IF
         GO TO 30
      END IF
      js(i) = ind
      GO TO 20
      END SUBROUTINE
      END
