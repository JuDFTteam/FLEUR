      SUBROUTINE nwrfst(mrad,nsol,is,it,nmatch,nzero,ferro,ec,rc,
     +                  pow,piw,gc,err,var,dv,varnew,errnew)

c..........................................................nwrfst
c starting values for E, A, A_out, A_inw
c
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: mrad
      REAL ec
      INTEGER is,it,nmatch,nsol,nzero
      LOGICAL ferro
C     ..
C     .. Array Arguments ..
      REAL dv(4),err(4),errnew(4),gc(2,2,mrad),piw(2,2),pow(2,2),
     +     rc(mrad),var(4),varnew(4)
C     ..
C     .. Local Scalars ..
      REAL ratt,rr,trymix
      INTEGER iv,j,n
C     ..
C     .. Local Arrays ..
      REAL niw(2),now(2)
C     ..
C     .. Data statements ..
C
      DATA trymix/0.010/
C     ..
C                                      --------------------
C                                       START VALUES FOR
C                                           PARAMETERS
C                                      --------------------
      var(1) = ec
      var(2) = pow(is,is)/piw(is,is)
c
      IF ((nsol.EQ.2) .AND. ferro) THEN
         DO 10 j = 1,nsol
            now(j) = 0.00
   10    CONTINUE
         DO 30 n = 1,nmatch - 1
            rr = rc(n)**3
            DO 20 j = 1,nsol
               now(j) = now(j) + gc(j,j,n)**2*rr
   20       CONTINUE
   30    CONTINUE
         DO 40 j = 1,nsol
            niw(j) = 0.00
   40    CONTINUE
         DO 60 n = nmatch,nzero - 1
            rr = rc(n)**3
            DO 50 j = 1,nsol
               niw(j) = niw(j) + gc(j,j,n)**2*rr
   50       CONTINUE
   60    CONTINUE
         ratt = pow(it,it)/piw(it,it)
         var(3) = trymix* (now(is)+niw(is)*var(2))/
     +            (now(it)+niw(it)*ratt)
         var(4) = ratt*var(3)/var(2)
      ELSE
         DO 70 iv = 3,4
            err(iv) = 0.00
            errnew(iv) = 0.00
            var(iv) = 0.00
            varnew(iv) = 0.00
            dv(iv) = 0.00
   70    CONTINUE
      END IF
c
      RETURN
      END
