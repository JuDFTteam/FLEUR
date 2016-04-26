      MODULE m_cnodes
      use m_juDFT
c...........................................................cnodes
c number of nodes
c
      CONTAINS
      SUBROUTINE cnodes(mrad,iflag,is,ec,l,xmj,nqn,vv,bb,rc,dx,
     +                  nmatch,nzero,gc,fc,pow,qow,piw,qiw,node)
c
      USE m_coredir
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: mrad
      REAL dx,ec,xmj
      INTEGER iflag,is,l,nmatch,node,nqn,nzero
C     ..
C     .. Array Arguments ..
      REAL bb(mrad),fc(2,2,mrad),gc(2,2,mrad),piw(2,2),pow(2,2),
     +     qiw(2,2),qow(2,2),rc(mrad),vv(mrad)
C     ..
C     .. Local Scalars ..
      INTEGER n
C     ..
c                    - outward solution -
      CALL coredir(mrad,ec,l,xmj,1,vv,bb,rc,dx,nmatch,nzero,
     +             gc,fc,pow,qow,piw,qiw)

      node = 0
      DO 10 n = 2,nmatch
         IF (gc(is,is,n)*gc(is,is,n-1).LT.0.0) node = node + 1
   10 CONTINUE
      IF (node.EQ. (nqn-l-1)) THEN
         IF ((gc(is,is,nmatch)/gc(is,is,nmatch-1).LE.0.0) .OR.
     +       (gc(is,is,nmatch)/gc(is,is,nmatch-1).GE.1.0)) THEN
            ec = 0.9*ec
            iflag = 1
!            write(*,*) '=',nmatch,is,node,ec
!            DO l = 1,nzero
!             write(*,*) l,gc(is,is,l)
!            ENDDO
!            stop
            IF (ec > -0.00000001) CALL juDFT_error("cnodes:1",calledby
     +           ="cnodes")
            GO TO 20
         END IF
      ELSE
         IF (node.GT. (nqn-l-1)) THEN
            ec = 1.2*ec
            write(*,*) '>',node,ec
         ELSE
            ec = 0.8*ec
            write(*,*) '<',node,ec
         END IF
         iflag = 1
         GO TO 20
      END IF
c                    - inward solution -
      CALL coredir(mrad,ec,l,xmj,2,vv,bb,rc,dx,nmatch,nzero,
     +             gc,fc,pow,qow,piw,qiw)
   20 CONTINUE

      END SUBROUTINE cnodes
      END MODULE m_cnodes
