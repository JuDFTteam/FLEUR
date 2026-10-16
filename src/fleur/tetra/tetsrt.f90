!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_tetsrt


   IMPLICIT NONE
   PRIVATE
   PUBLIC :: tetsrt

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