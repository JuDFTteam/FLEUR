MODULE m_nwrfst

   CONTAINS
   
   SUBROUTINE nwrfst(mrad,nsol,is,it,nmatch,nzero,ferro,ec,rc,pow,piw,gc,err,var,dv,varnew,errnew)

!..........................................................nwrfst
! starting values for E, A, A_out, A_inw
!
      IMPLICIT NONE

      INTEGER, INTENT (IN) :: mrad
      INTEGER, INTENT (IN) :: is, it, nmatch, nsol, nzero
      LOGICAL, INTENT (IN) :: ferro
      REAL, INTENT    (IN) :: ec

      REAL, INTENT (INOUT) :: dv(4), err(4), errnew(4), var(4), varnew(4)
      REAL, INTENT (IN)    :: gc(2,2,mrad), piw(2,2), pow(2,2), rc(mrad)


      REAL ratt, rr, trymix
      INTEGER iv, j, n

      REAL niw(2), now(2)

      trymix = 0.010

!                                      --------------------
!                                       START VALUES FOR
!                                           PARAMETERS
!                                      --------------------
      var(1) = ec
      var(2) = pow(is,is)/piw(is,is)

      IF ((nsol.EQ.2) .AND. ferro) THEN
         DO j = 1,nsol
            now(j) = 0.00
         END DO
         DO n = 1,nmatch - 1
            rr = rc(n)**3
            DO j = 1,nsol
               now(j) = now(j) + gc(j,j,n)**2*rr
            END DO
         END DO
         DO j = 1,nsol
            niw(j) = 0.00
         END DO
         DO n = nmatch,nzero - 1
            rr = rc(n)**3
            DO j = 1,nsol
               niw(j) = niw(j) + gc(j,j,n)**2*rr
            END DO
         END DO
         ratt = pow(it,it)/piw(it,it)
         var(3) = trymix * (now(is)+niw(is)*var(2)) / (now(it)+niw(it)*ratt)
         var(4) = ratt * var(3) / var(2)
      ELSE
         DO iv = 3,4
            err(iv) = 0.00
            errnew(iv) = 0.00
            var(iv) = 0.00
            varnew(iv) = 0.00
            dv(iv) = 0.00
         END DO
      END IF

   END SUBROUTINE nwrfst
      
END MODULE m_nwrfst
