MODULE m_dosbin
CONTAINS
   SUBROUTINE dos_bin(jspins, wtkpt, e, eig, qal, g, energyShift)
      !! This subroutine generates the idos, the ldos, the partial
      !! ldos in the spheres and the z-dependent dos for the vacuum.

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: jspins
      REAL,    INTENT(IN)  :: wtkpt(:),e(:)
      REAL,    INTENT(IN)  :: eig(:,:,:),qal(:,:,:)
      REAL,    INTENT(OUT) :: g(:,:)

      REAL,    OPTIONAL, INTENT(IN) :: energyShift

      INTEGER :: nl, k, j, i, js
      REAL    :: de, wk, emin, shift

      de = abs(e(2)-e(1))
      g = 0.0
      shift = 0.0
      IF(PRESENT(energyShift)) shift = energyShift
      emin = minval(e)
      ! put weights in the right bins
      DO js=1, size(qal,3)
         DO k = 1 , size(qal,2)
            wk = wtkpt(k)/de
            DO j = 1 , size(eig,1)
               i = NINT((eig(j,k,js)-shift-emin)/de) + 1
               IF ( (i.LE.size(g,1)) .AND. (i.GE.1) ) THEN
                  g(i,js) = g(i,js) + wk*qal(j,k,js)* 2.0/jspins
               END IF
            END DO
         END DO
      END DO
   END SUBROUTINE dos_bin
END MODULE m_dosbin
