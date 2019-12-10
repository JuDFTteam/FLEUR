MODULE m_tetsrt

   CONTAINS

   SUBROUTINE tetsrt(n,etetra,ind)

      IMPLICIT NONE

      INTEGER,          INTENT(IN)     :: n
      REAL,             INTENT(IN)     :: etetra(:)
      INTEGER,          INTENT(INOUT)  :: ind(:)

      INTEGER i,j
      INTEGER tmp


      DO i = 1, n
         ind(i) = i
      ENDDO

      !Sort the energies in the tetrahedron in ascending order
      DO i = 1, n-1
         DO j = i+1, n
            IF (etetra(ind(i)).GT.etetra(ind(j))) THEN
               tmp = ind(i)
               ind(i) = ind(j)
               ind(j) = tmp
            ENDIF
         ENDDO
      ENDDO

   END SUBROUTINE tetsrt

END MODULE m_tetsrt