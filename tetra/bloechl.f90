MODULE m_bloechl

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: bloechl

   CONTAINS

   PURE REAL FUNCTION bloechl(efermi,etetra,ind,vol,film)

      REAL,                INTENT(IN)     :: efermi
      REAL,                INTENT(IN)     :: etetra(:)
      INTEGER,             INTENT(IN)     :: ind
      REAL,                INTENT(IN)     :: vol
      LOGICAL,             INTENT(IN)     :: film

      REAL dos_ef

      bloechl = 0.0
      IF(.NOT. film) THEN
         dos_ef = dos_tetra(efermi,etetra,vol)
         bloechl = 1/40.0*dos_ef*sum(etetra-etetra(ind))
      ENDIF

   END FUNCTION bloechl

   PURE REAL FUNCTION dos_tetra(energy,etetra,vol)

      REAL,                INTENT(IN)     :: energy
      REAL,                INTENT(IN)     :: etetra(:)
      REAL,                INTENT(IN)     :: vol

      IF((energy.GT.etetra(4)).OR.(energy.LT.etetra(1))) THEN
         dos_tetra = 0.0
      ELSE IF(energy.GE.etetra(3)) THEN

         dos_tetra = 3.0 * vol * (etetra(4)-energy)**2/((etetra(4)-etetra(1))*(etetra(4)-etetra(2))*(etetra(4)-etetra(3)))

      ELSE IF(energy.GE.etetra(2)) THEN

         dos_tetra = vol * 3.0/((etetra(3)-etetra(1))*(etetra(4)-etetra(1))) * (etetra(2) - etetra(1) + 2*(energy - etetra(2)) &
               -(etetra(3)-etetra(1)+etetra(4)-etetra(2)) * (energy-etetra(2))**2/((etetra(3)-etetra(2))*(etetra(4)-etetra(2)))) 

      ELSE IF(energy.GE.etetra(1)) THEN

         dos_tetra = 3.0 * vol * (energy-etetra(1))**2/((etetra(2)-etetra(1))*(etetra(3)-etetra(1))*(etetra(4)-etetra(1)))

      END IF

   END FUNCTION dos_tetra

END MODULE m_bloechl