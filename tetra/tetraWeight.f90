MODULE m_tetraWeight

   IMPLICIT NONE

   PRIVATE
   PUBLIC   :: tetraWeight

   CONTAINS

   PURE REAL FUNCTION tetraWeight(efermi,etetra,ind,vol,film)

      REAL,       INTENT(IN)     :: efermi
      REAL,       INTENT(IN)     :: etetra(:)
      INTEGER,    INTENT(IN)     :: ind
      REAL,       INTENT(IN)     :: vol
      LOGICAL,    INTENT(IN)     :: film

      IF(film) THEN
         tetraWeight = tetraWeightFilm(efermi,etetra,ind) * vol/3.0
      ELSE
         tetraWeight = tetraWeightBulk(efermi,etetra,ind) * vol/4.0
      ENDIF

   END FUNCTION tetraWeight

   PURE REAL FUNCTION tetraWeightBulk(efermi,etetra,ind)

      !-------------------------------------------------------------
      ! Returns the integration weights for a single tetrahedron
      ! at one corner for charge density integration with the
      ! given fermi efermi without the volume factor in Bulk
      !
      ! The energies at the corners etetra are ordered so that
      !     etetra(1) <= etetra(2) <= etetra(3) <= etetra(4)
      ! ind is the corner index of the kpoint we are interested in
      !-------------------------------------------------------------

      REAL,       INTENT(IN)     :: efermi
      REAL,       INTENT(IN)     :: etetra(:)
      INTEGER,    INTENT(IN)     :: ind

      REAL C1, C2, C3
      INTEGER i,j
      REAL f(4)
      REAL e(4,4)

      tetraWeightBulk = 0.0

      !Prepare differences between eigenvalues at the corners e(i,j) = ei-ej
      !This is to improve readability in the formulas
      DO i = 1, 4
         DO j = 1, 4
            e(i,j) = etetra(i) - etetra(j)
         ENDDO
      ENDDO

      IF( efermi>=etetra(4) ) THEN

         tetraWeightBulk = 1.0

      ELSE IF( efermi>=etetra(3) ) THEN

         f(1:3) = ( etetra(4)-efermi )/e(4,1:3)

         SELECT CASE(ind)
         CASE(4)
            tetraWeightBulk = 1.0 - f(1) * f(2) * f(3) * ( 4.0 - f(1) - f(2) - f(3))
         CASE DEFAULT
            tetraWeightBulk = 1.0 - f(1) * f(2) * f(3) * f(ind)
         END SELECT

      ELSE IF( efermi>=etetra(2) ) THEN

         C1     = ( efermi-etetra(1) )**2/( e(4,1)*e(3,1) )

         C2     = ( efermi-etetra(1) )*( efermi-etetra(2) )*( etetra(3)-efermi )/( e(3,1)*e(3,2)*e(4,1) )

         C3     = ( efermi-etetra(2) )**2*( etetra(4)-efermi )/( e(3,2)*e(4,1)*e(4,2) )

         SELECT CASE(ind)
         CASE(1)
            tetraWeightBulk = C1 + (C1 + C2) * (etetra(3)-efermi)/e(3,1) + &
                             (C1 + C2 + C3) * (etetra(4)-efermi)/e(4,1)
         CASE(2)
            tetraWeightBulk = C1 + C2 + C3 + (C2 + C3) * (etetra(3)-efermi)/e(3,2) + &
                              C3 * (etetra(4)-efermi)/e(4,2)
         CASE(3)
            tetraWeightBulk = (C1 + C2) * (efermi-etetra(1))/e(3,1)+&
                              (C2 + C3) * (efermi-etetra(2))/e(3,2)
         CASE(4)
            tetraWeightBulk = (C1 + C2 + C3) * (efermi-etetra(1))/e(4,1) +&
                               C3 * (efermi-etetra(2))/e(4,2)
         CASE DEFAULT
         END SELECT

      ELSE IF( efermi>=etetra(1) ) THEN

         f(2:4) = ( efermi-etetra(1) )/e(2:4,1)

         SELECT CASE(ind)
         CASE(1)
            tetraWeightBulk = f(2) * f(3) * f(4) * (4.0 - f(2) - f(3) - f(4))
         CASE DEFAULT
            tetraWeightBulk = f(2) * f(3) * f(4) * f(ind)
         END SELECT

      END IF

   END FUNCTION tetraWeightBulk

   PURE REAL FUNCTION tetraWeightFilm(efermi,etetra,ind)

      !-------------------------------------------------------------
      ! Returns the integration weights for a single tetrahedron
      ! at one corner for charge density integration with the
      ! given fermi efermi without the volume factor in Film
      !
      ! The energies at the corners etetra are ordered so that
      !     etetra(1) <= etetra(2) <= etetra(3)
      ! ind is the corner index of the kpoint we are interested in
      !-------------------------------------------------------------

      REAL,       INTENT(IN)     :: efermi
      REAL,       INTENT(IN)     :: etetra(:)
      INTEGER,    INTENT(IN)     :: ind

      REAL f(3)

      tetraWeightFilm = 0.0

      IF( efermi>etetra(3) ) THEN

         tetraWeightFilm = 1.0

      ELSE IF( efermi>=etetra(2) ) THEN

         f(1:2) = ( etetra(3)-efermi )/( etetra(3)-etetra(1:2) )

         SELECT CASE(ind)
         CASE(3)
            tetraWeightFilm = 1.0 - f(1) * f(2) * ( 3.0 - f(1) - f(2) )
         CASE DEFAULT
            tetraWeightFilm = 1.0 - f(1) * f(2) * f(ind)
         END SELECT

      ELSE IF( efermi>=etetra(1) ) THEN

         f(2:3) = ( efermi-etetra(1) )/( etetra(2:3)-etetra(1) )

         SELECT CASE(ind)
         CASE(1)
            tetraWeightFilm = f(2) * f(3) * ( 3.0 - f(2) - f(3) )
         CASE DEFAULT
            tetraWeightFilm = f(2) * f(3) * f(ind)
         END SELECT

      END IF

   END FUNCTION tetraWeightFilm

END MODULE m_tetraWeight