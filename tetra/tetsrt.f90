MODULE m_tetsrt


   IMPLICIT NONE

   CONTAINS

   PURE FUNCTION tetsrt(etetra) Result(ind)

      REAL,             INTENT(IN)     :: etetra(:)

      INTEGER :: ind(SIZE(etetra))

      INTEGER i,j
      INTEGER tmp


      DO i = 1, SIZE(etetra)
         ind(i) = i
      ENDDO

      !Sort the energies in the tetrahedron in ascending order
      DO i = 1, SIZE(etetra)-1
         DO j = i+1, SIZE(etetra)
            IF (etetra(ind(i)).GT.etetra(ind(j))) THEN
               tmp = ind(i)
               ind(i) = ind(j)
               ind(j) = tmp
            ENDIF
         ENDDO
      ENDDO

   END FUNCTION tetsrt

END MODULE m_tetsrt