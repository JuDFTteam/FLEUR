MODULE m_kkintgr

CONTAINS

   SUBROUTINE kkintgr(ne,sigma,del,im,g)

      USE m_types
      USE m_constants

      IMPLICIT NONE

      COMPLEX,             INTENT(OUT)    :: g(:)
      REAL,                INTENT(IN)     :: im(:)

      INTEGER,             INTENT(IN)     :: ne
      REAL,                INTENT(IN)     :: sigma       
      REAL,                INTENT(IN)     :: del

      !Local declarations:
      INTEGER i, j
      REAL gfr, gfi
      COMPLEX clg

      !$OMP PARALLEL DEFAULT(none) &
      !$OMP SHARED(ne,sigma,del) &
      !$OMP SHARED(g,im) &
      !$OMP PRIVATE(gfr,gfi,j,clg)

      !$OMP DO
      DO i = 1, ne
         gfr = 0.0
         gfi = 0.0
         DO j = 1, ne

            !value of the integral between j-1 and j if the imaginary part of G is assumed to be constant
            clg = LOG(((j-i)*del - ImagUnit * sigma)/((j-i-1)*del - ImagUnit * sigma))

            gfr = gfr + im(j) * REAL(clg)
            !gfi = gfi + im(j) * AIMAG(clg)
            
         ENDDO

         g(i) = (-1.0) * gfr/pi_const - ImagUnit * im(i)
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL


   ENDSUBROUTINE kkintgr

ENDMODULE m_kkintgr