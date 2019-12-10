MODULE m_tetraWeight

   IMPLICIT NONE

   PUBLIC   :: tetraWeight
   PRIVATE  :: tetraWeightBulk,tetraWeightFilm

   CONTAINS

   SUBROUTINE tetraWeight(efermi,etetra,ind,vol,film,weight)

      IMPLICIT NONE

      REAL,       INTENT(IN)  :: efermi
      REAL,       INTENT(IN)  :: etetra(:)
      INTEGER,    INTENT(IN)  :: ind
      REAL,       INTENT(IN)  :: vol
      LOGICAL,    INTENT(IN)  :: film
      REAL,       INTENT(OUT) :: weight

      IF(film) THEN
         CALL tetraWeightFilm(efermi,etetra,ind,weight)
         weight = weight * vol/3.0
      ELSE
         CALL tetraWeightBulk(efermi,etetra,ind,weight)
         weight = weight * vol/4.0
      ENDIF

   END SUBROUTINE tetraWeight

   SUBROUTINE tetraWeightBulk(efermi,etetra,ind,weight)

      !-------------------------------------------------------------
      ! Returns the integration weights for a single tetrahedron
      ! at one corner for charge density integration with the
      ! given fermi efermi without the volume factor in Bulk
      !
      ! The energies at the corners etetra are ordered so that
      !     etetra(1) <= etetra(2) <= etetra(3) <= etetra(4)
      ! ind is the corner index of the kpoint we are interested in
      !-------------------------------------------------------------

      IMPLICIT NONE

      REAL,       INTENT(IN)  :: efermi
      REAL,       INTENT(IN)  :: etetra(:)
      INTEGER,    INTENT(IN)  :: ind
      REAL,       INTENT(OUT) :: weight

      REAL C1, C2, C3
      INTEGER i,j
      REAL f(4)
      REAL e(4,4)

      weight = 0.0

      !Prepare differences between eigenvalues at the corners e(i,j) = ei-ej
      !This is to improve readability in the formulas
      DO i = 1, 4
         DO j = 1, 4
            e(i,j) = etetra(i) - etetra(j)
         ENDDO
      ENDDO

      IF( efermi>=etetra(4) ) THEN

         weight = 1.0

      ELSE IF( efermi>=etetra(3) ) THEN

         f(1:3) = ( etetra(4)-efermi )/e(4,1:3)

         IF( ind==4 ) weight = 1.0 - f(1) * f(2) * f(3) * ( 4.0 - f(1) - f(2) - f(3))
         IF( ind/=4 ) weight = 1.0 - f(1) * f(2) * f(3) * f(ind)

      ELSE IF( efermi>=etetra(2) ) THEN

         C1     = ( efermi-etetra(1) )**2/( e(4,1)*e(3,1) )

         C2     = ( efermi-etetra(1) )*( efermi-etetra(2) )*( etetra(3)-efermi )/( e(3,1)*e(3,2)*e(4,1) )

         C3     = ( efermi-etetra(2) )**2*( etetra(4)-efermi )/( e(3,2)*e(4,1)*e(4,2) )

         IF( ind==1 ) weight = C1 + (C1 + C2) * (etetra(3)-efermi)/e(3,1) + &
                                 (C1 + C2 + C3) * (etetra(4)-efermi)/e(4,1)

         IF( ind==2 ) weight = C1 + C2 + C3 + (C2 + C3) * (etetra(3)-efermi)/e(3,2) + &
                                 C3 * (etetra(4)-efermi)/e(4,2)

         IF( ind==3 ) weight = (C1 + C2) * (efermi-etetra(1))/e(3,1)+&
                                 (C2 + C3) * (efermi-etetra(2))/e(3,2)

         IF( ind==4 ) weight = (C1 + C2 + C3) * (efermi-etetra(1))/e(4,1) +&
                                 C3 * (efermi-etetra(2))/e(4,2)

      ELSE IF( efermi>=etetra(1) ) THEN

         f(2:4) = ( efermi-etetra(1) )/e(2:4,1)

         IF( ind==1 ) weight = f(2) * f(3) * f(4) * (4.0 - f(2) - f(3) - f(4))
         IF( ind/=1 ) weight = f(2) * f(3) * f(4) * f(ind)

      END IF

   END SUBROUTINE tetraWeightBulk

   SUBROUTINE tetraWeightFilm(efermi,etetra,ind,weight)

      !-------------------------------------------------------------
      ! Returns the integration weights for a single tetrahedron
      ! at one corner for charge density integration with the
      ! given fermi efermi without the volume factor in Film
      !
      ! The energies at the corners etetra are ordered so that
      !     etetra(1) <= etetra(2) <= etetra(3)
      ! ind is the corner index of the kpoint we are interested in
      !-------------------------------------------------------------

      IMPLICIT NONE

      REAL,       INTENT(IN)  :: efermi
      REAL,       INTENT(IN)  :: etetra(:)
      INTEGER,    INTENT(IN)  :: ind
      REAL,       INTENT(OUT) :: weight

      REAL f(3)

      weight = 0.0

      IF( efermi>etetra(3) ) THEN

         weight = 1.0

      ELSE IF( efermi>=etetra(2) ) THEN

         f(1:2) = ( etetra(3)-efermi )/( etetra(3)-etetra(1:2) )

         IF( ind/=3 ) weight = 1.0 - f(1) * f(2) * f(ind)
         IF( ind==3 ) weight = 1.0 - f(1) * f(2) * ( 3.0 - f(1) - f(2) )

      ELSE IF( efermi>=etetra(1) ) THEN

         f(2:3) = ( efermi-etetra(1) )/( etetra(2:3)-etetra(1) )

         IF( ind==1 ) weight = f(2) * f(3) * ( 3.0 - f(2) - f(3) )
         IF( ind/=1 ) weight = f(2) * f(3) * f(ind)

      END IF

   END SUBROUTINE tetraWeightFilm

END MODULE m_tetraWeight