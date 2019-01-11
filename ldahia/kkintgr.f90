MODULE m_kkintgr

CONTAINS

   SUBROUTINE kkintgr(n_cut,sigma,del,im,g,bot)

      USE m_types
      USE m_constants

      IMPLICIT NONE

      COMPLEX,             INTENT(OUT)    :: g(:)
      REAL,                INTENT(IN)     :: im(:)

      INTEGER,             INTENT(IN)     :: n_cut
      REAL,                INTENT(IN)     :: sigma       
      REAL,                INTENT(IN)     :: del
      REAL,                INTENT(IN)     :: bot

      !Local declarations:
      INTEGER i, j
      REAL gfr, gfi
      COMPLEX clg


      !$OMP PARALLEL DEFAULT(none) &
      !$OMP SHARED(n_cut,sigma,del) &
      !$OMP SHARED(g,im) &
      !$OMP PRIVATE(gfr,gfi,j,clg)

      !$OMP DO
      DO i = 1, n_cut
         gfr = 0.0
         gfi = 0.0
         DO j = 1, i-2

            !value of the integral between j-1 and j if the imaginary part of G is assumed to be constant
            clg = LOG(REAL(i-j)/REAL(i-j-1))

            gfr = gfr + im(j) * REAL(clg) + (im(j+1)-im(j))*(REAL(clg)*(i-j)-1)
            
         ENDDO

         gfr = gfr + im(i-1)-im(i+1)

         DO j = i+2, n_cut

            !value of the integral between j-1 and j if the imaginary part of G is assumed to be constant
            clg = LOG(REAL(j-i-1)/REAL(j-i))

            gfr = gfr + im(j) * REAL(clg) + (im(j-1)-im(j))*(1+REAL(clg)*(j-i))
            
         ENDDO

         g(i) = 1/pi_const * gfr + ImagUnit * im(i)
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL


   ENDSUBROUTINE kkintgr

ENDMODULE m_kkintgr