!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_sort
   USE m_judft
CONTAINS

   SUBROUTINE sort(ind,lv,lv1)
      !********************************************************************
      !     heapsort routine
      !     input:   lv    = array of objects to be sorted
      !              lv1   = second array to use as secondary sort key
      !     output:  ind(i) = index of i'th smallest object
      !********************************************************************
      IMPLICIT NONE
      !
      REAL,    INTENT (IN) :: lv(:)
      INTEGER, INTENT (OUT) :: ind(:)
      REAL,INTENT(IN),OPTIONAL :: lv1(:)
      !     ..
      !     .. Local Scalars ..
      REAL eps,q,q1
      INTEGER i,idx,ir,j,l,n
      REAL,ALLOCATABLE :: llv(:)
      !     ..
      !     .. Data statements ..
      DATA eps/1.e-10/
      !     ..
      !
      n=SIZE(ind)
      IF (n>SIZE(lv)) CALL judft_error("BUG: incosistent dimensions")
      ALLOCATE(llv(n))
      IF (PRESENT(lv1)) THEN
         IF (n>SIZE(lv1)) CALL judft_error("BUG: incosistent dimensions")
         llv=lv1
      ELSE
         llv=(/(1.*i,i=1,n)/)
      END IF
      IF (n == 0) RETURN ! Nothing to do
      IF (n == 1) THEN   ! Not much to do
         ind(1) = 1
         RETURN
      END IF

      DO i = 1,n
         ind(i) = i
      ENDDO
      !
      l = n/2 + 1
      ir = n
      DO
         IF (l.GT.1) THEN
            l = l - 1
            idx = ind(l)
            q = lv(idx)
            q1= llv(idx)
         ELSE
            idx = ind(ir)
            q = lv(idx)
            q1= llv(idx)
            ind(ir) = ind(1)
            ir = ir - 1
            IF (ir.EQ.1) THEN
               ind(1) = idx
               RETURN
            END IF
         END IF
         i = l
         j = l + l
         DO WHILE(j.LE.ir)
            IF (j.LT.ir) THEN
               !           if(lv(ind(j)).lt.lv(ind(j+1))) j=j+1
               IF (((lv(ind(j+1))-lv(ind(j))).GE.eps).OR. &!Standard comparison
                   ((ABS((lv(ind(j+1))-lv(ind(j))))<eps).AND.&!Same length, check second key
                    ((llv(ind(j+1))-llv(ind(j))).GE.eps))) &
                  j=j+1
            END IF
            !        if(q.lt.lv(ind(j))) then
            IF ((lv(ind(j))-q).GE.eps.OR.&
                (ABS((lv(ind(j))-q))<eps.AND.(llv(ind(j))-q1).GE.eps))THEN
               ind(i) = ind(j)
               i = j
               j = j + j
            ELSE
               j = ir + 1
            END IF
         enddo
         ind(i) = idx
      ENDDO
   END SUBROUTINE sort

END MODULE m_sort
