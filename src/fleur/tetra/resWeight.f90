MODULE m_resWeight

   IMPLICIT NONE

   PRIVATE
   PUBLIC   :: resWeight

   CONTAINS

   PURE FUNCTION resWeight(eMesh,etetra,ind,vol,film) Result(weights)

      REAL,       INTENT(IN)     :: eMesh(:)
      REAL,       INTENT(IN)     :: etetra(:)
      INTEGER,    INTENT(IN)     :: ind
      REAL,       INTENT(IN)     :: vol
      LOGICAL,    INTENT(IN)     :: film

      REAL :: weights(SIZE(eMesh))
      INTEGER ie

      IF(.not.film) then
         weights = resWeightBulk(eMesh,etetra,ind) * vol
      ELSE
         DO ie = 1, SIZE(eMesh)
            IF(film) THEN
               weights(ie) = resWeightFilm(eMesh(ie),etetra,ind) * vol
            ENDIF
         ENDDO
      ENDIF

   END FUNCTION resWeight

   PURE FUNCTION resWeightBulk(eMesh,e,ind) Result(weights)

      !Calculates resolvent weights for point in tetrahedron (Compare PhysRevB.29.3430)

      REAL,          INTENT(IN)  :: eMesh(:)
      REAL,          INTENT(IN)  :: e(:)
      INTEGER,       INTENT(IN)  :: ind

      REAL :: weights(SIZE(eMesh))

      REAL, PARAMETER :: tol = 1e-7!tolerance for degeneracy
      REAL, PARAMETER :: fac = 5.0
      REAL, PARAMETER :: highEnergyCut = 0.1

      REAL a,b,cut,min,z,denom(4)
      REAL eShift
      REAL eDeg(4)
      INTEGER lowerCut,upperCut,ie
      INTEGER ndeg,i,j,k,l,m
      INTEGER ideg(6,2)

      REAL, ALLOCATABLE :: eMesh_shifted(:)
      REAL :: eig_shifted(4)

      eShift = 0.0
      if (MAXVAL(eMesh)*MINVAL(eMesh) < 0.0) then
         !0 is inside the interval
         !The asymptotic relations for energies far from the
         !eigenvalues at the corners has a pole for energies at 0
         !To prevent these we shift everything (eignevalues and energy grid so that zero)
         !is far away. This does not change the results for all other cases
         !Since there only energy differences contribute
         eShift = 2.0*MAX(ABS(MINVAL(eMesh)),ABS(MAXVAL(eMesh)))
      endif

      weights = 0.0
      cut = 0.0
      !
      DO i = 1, 3
         DO j=i+1,4
            IF(ABS(e(i)-e(j)).GT.cut) cut = MAX(ABS(e(i)-e(j)),tol)
         ENDDO
      ENDDO
      !Determine the cutoffs
      upperCut = SIZE(eMesh)
      DO
         IF(MINVAL(ABS(eMesh(upperCut)-e)).LT.fac*cut.OR.upperCut.EQ.1) THEN
            EXIT
         ELSE IF(MAXVAL(eMesh(upperCut)-e) < 0.0) THEN
            EXIT
         ELSE
            upperCut = upperCut - 1
         ENDIF
      ENDDO
      lowerCut = 1
      DO
         IF(MINVAL(ABS(eMesh(lowerCut)-e)).LT.fac*cut.OR.lowerCut.EQ.upperCut) THEN
            EXIT
         ELSE IF(MINVAL(eMesh(lowerCut)-e) > 0.0) THEN
            EXIT
         ELSE
            lowerCut = lowerCut + 1
         ENDIF
      ENDDO
      eig_shifted = e(:) - eShift

      ALLOCATE(eMesh_shifted,mold=eMesh)
      eMesh_shifted = eMesh(:) - eShift

      ndeg = 0
      ideg = 0
      !Search for degenerate eigenvalues
      eDeg = eig_shifted
      DO i = 1, 3
         DO j = i+1,4
            IF(ABS(eig_shifted(i)-eig_shifted(j)).LT.tol) THEN
               ndeg = ndeg + 1
               ideg(ndeg,1) = i
               ideg(ndeg,2) = j
               !Set the two values equal
               eDeg(i) = (eig_shifted(i)+eig_shifted(j))/2.0
               eDeg(j) = eDeg(i)
               DO k = 1, ndeg-1
                  IF(ideg(k,1).EQ.i.OR.ideg(k,1).EQ.j) THEN
                     eDeg(ideg(k,1)) = (eig_shifted(i)+eig_shifted(j)+eig_shifted(ideg(k,1)))/3.0
                     eDeg(j) = eDeg(ideg(k,1))
                     eDeg(i) = eDeg(ideg(k,1))
                  ELSE IF(ideg(k,2).EQ.i.OR.ideg(k,2).EQ.j) THEN
                     eDeg(ideg(k,2)) = (eig_shifted(i)+eig_shifted(j)+eig_shifted(ideg(k,2)))/3.0
                     eDeg(j) = eDeg(ideg(k,2))
                     eDeg(i) = eDeg(ideg(k,2))
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      ndeg = 0
      ideg = 0
      DO i = 1, 3
         DO j = i+1,4
            IF(ABS(eDeg(i)-eDeg(j)).LT.tol) THEN
               ndeg = ndeg + 1
               ideg(ndeg,1) = i
               ideg(ndeg,2) = j
            ENDIF
         ENDDO
      ENDDO

      IF(MAXVAL(MAXVAL(eMesh)-e)<-highEnergyCut) THEN
         lowerCut=SIZE(eMesh)
         upperCut=SIZE(eMesh)
      ENDIF

      IF(MINVAL(MINVAL(eMesh)-e)>highEnergyCut) THEN
         lowerCut=1
         upperCut=1
      ENDIF

      !asymptotic relation Eqs. 9-11
      !the formulas below are unstable for big z arguments
      a = (eDeg(ind) + SUM(eDeg(:)))/5.0
      b = 0.0
      DO j = 1, 4
         IF(j.EQ.ind) CYCLE
         b = b + 3*(eDeg(j)-eDeg(ind))**2
         DO k = 1, 4
            IF(k.EQ.j) CYCLE
            b = b + (eDeg(k)-eDeg(j))**2
         ENDDO
      ENDDO
      b = b/300.0

      weights(  :lowerCut) = 0.25/(eMesh_shifted(  :lowerCut)-a-b/eMesh_shifted(  :lowerCut))
      weights(upperCut:) = 0.25/(eMesh_shifted(upperCut:)-a-b/eMesh_shifted(upperCut:))

      IF(lowerCut==upperCut) RETURN

      IF(ndeg.EQ.0) THEN
         denom = 1.0
         DO i = 1, 4
            DO j = 1, 4
               IF(i.EQ.j) CYCLE
               denom(i) = denom(i)*(eDeg(j)-eDeg(i))
            ENDDO
            IF(i.EQ.ind) CYCLE
            denom(i) = denom(i)*(eDeg(ind)-eDeg(i))
         ENDDO
         DO ie = lowerCut, upperCut
            z = eMesh_shifted(ie)
            !all eigenvalues are non degenerate (EQ.7 from the paper)
            !First term
            weights(ie) = (z-eDeg(ind))**2/denom(ind)
            DO j = 1, 4
               IF(j.EQ.ind) CYCLE
               weights(ie) = weights(ie) + (z-eDeg(j))**3/denom(j)*LOG(ABS((z-eDeg(j)))/ABS((z-eDeg(ind))))
            ENDDO
         ENDDO

      ELSE IF(ndeg.EQ.1) THEN

         !Two eigenvalues are degenerate (EQs. A1 A2)

         l = ideg(1,1)
         m = ideg(1,2)

         IF(ind.EQ.l.OR.ind.EQ.m) THEN
            !EQ. A2

            !Get the two indices different from l,m and ind
            j = 0
            DO i = 1, 4
               IF(i.EQ.l.OR.i.EQ.m) CYCLE
               j = i
            ENDDO
            k = 0
            DO i = 1, 4
               IF(i.EQ.l.OR.i.EQ.m.OR.i.EQ.j) CYCLE
               k = i
            ENDDO
            DO ie = lowerCut, upperCut
               z = eMesh_shifted(ie)
               weights(ie) = (z-eDeg(k))**3/((eDeg(k)-eDeg(j))*(eDeg(k)-eDeg(m))**3)*LOG(ABS(z-eDeg(k))) &
                            +(z-eDeg(j))**3/((eDeg(j)-eDeg(k))*(eDeg(j)-eDeg(m))**3)*LOG(ABS(z-eDeg(j))) &
                            +(z-eDeg(m))/((eDeg(m)-eDeg(j))*(eDeg(m)-eDeg(k))) * (0.5 + (z-eDeg(j))/(eDeg(m)-eDeg(j))&
                           +(z-eDeg(k))/(eDeg(m)-eDeg(k)) + ((z-eDeg(j))**2/(eDeg(m)-eDeg(j))**2 &
                           +(z-eDeg(k))**2/(eDeg(m)-eDeg(k))**2 +(z-eDeg(j))/(eDeg(m)-eDeg(j))*(z-eDeg(k))/(eDeg(m)-eDeg(k))) &
                           *LOG(ABS(z-eDeg(m))))
            ENDDO
         ELSE
            !k is the one site not equal to ind, l or m
            k = 0
            DO i = 1, 4
               IF(i.EQ.ind.OR.i.EQ.l.OR.i.EQ.m) CYCLE
               k = i
            ENDDO
            !EQ. A1
            DO ie = lowerCut, upperCut
               z = eMesh_shifted(ie)
               weights(ie) = (z-eDeg(ind))**2/((eDeg(ind)-eDeg(m))**2*(eDeg(k)-eDeg(ind))) *&
                        (1 + (2*(z-eDeg(m))/(eDeg(ind)-eDeg(m))+(z-eDeg(k))/(eDeg(ind)-eDeg(k)))*LOG(ABS(z-eDeg(ind)))) &
                       +(z-eDeg(m))**2/((eDeg(m)-eDeg(ind))**2*(eDeg(k)-eDeg(m))) *&
                        (1 + (2*(z-eDeg(ind))/(eDeg(m)-eDeg(ind))+(z-eDeg(k))/(eDeg(m)-eDeg(k)))*LOG(ABS(z-eDeg(m)))) &
                       +(z-eDeg(k))**3/((eDeg(k)-eDeg(ind))**2*(eDeg(k)-eDeg(m))**2)*LOG(ABS(z-eDeg(k)))
            ENDDO
         ENDIF
      ELSE IF(ndeg.EQ.2) THEN
         !This is the case E1=E2<E3=E4 => A4
         IF(ind.LE.2) THEN
            DO ie = lowerCut, upperCut
               z = eMesh_shifted(ie)
               weights(ie) = 3.0*(z-eDeg(3))**2*(z-eDeg(2))/(eDeg(3)-eDeg(2))**4*LOG(ABS((z-eDeg(2)))/ABS((z-eDeg(3)))) &
                            - 3.0/2.0*(z-eDeg(2))*(2*(z-eDeg(3))-eDeg(3)+eDeg(2))/(eDeg(3)-eDeg(2))**3-1.0/(eDeg(3)-eDeg(2))
            ENDDO
         ELSE
            DO ie = lowerCut, upperCut
               z = eMesh_shifted(ie)
               weights(ie) = 3.0*(z-eDeg(2))**2*(z-eDeg(3))/(eDeg(3)-eDeg(2))**4*LOG(ABS((z-eDeg(3)))/ABS((z-eDeg(2)))) &
                            + 3.0/2.0*(z-eDeg(3))*(2*(z-eDeg(2))+eDeg(3)-eDeg(2))/(eDeg(3)-eDeg(2))**3+1.0/(eDeg(3)-eDeg(2))
            ENDDO
         ENDIF
      ELSE IF(ndeg.GE.3 .AND. ndeg.LT.6) THEN
         !EQs A3/A5 (here we explicitly write each weight)
         IF(ALL(ideg(:,:).NE.4)) THEN
            !A3
            IF(ind.NE.4) THEN
               DO ie = lowerCut, upperCut
                  z = eMesh_shifted(ie)
                  weights(ie) = (z-eDeg(4))**3/(eDeg(4)-eDeg(3))**4*LOG(ABS((z-eDeg(4))/(z-eDeg(3)))) &
                            + (6*(z-eDeg(4))**2-3*(eDeg(4)-eDeg(3))*(z-eDeg(4))+2*(eDeg(4)-eDeg(3))**2)/(6*(eDeg(4)-eDeg(3))**3)
               ENDDO
            ELSE
               DO ie = lowerCut, upperCut
                  z = eMesh_shifted(ie)
                  weights(ie) = 3.0*(z-eDeg(4))**2*(z-eDeg(3))/(eDeg(4)-eDeg(3))**4*LOG(ABS((z-eDeg(3))/(z-eDeg(4)))) &
                            - 3.0/2.0*(z-eDeg(3))*(2*(z-eDeg(4))-eDeg(4)+eDeg(3))/(eDeg(4)-eDeg(3))**3-1.0/(eDeg(4)-eDeg(3))
               ENDDO
            ENDIF
         ELSE IF(ALL(ideg(:,:).NE.1)) THEN
            !A5
            IF(ind.EQ.1) THEN
               DO ie = lowerCut, upperCut
                  z = eMesh_shifted(ie)
                  weights(ie) = 3.0*(z-eDeg(1))**2*(z-eDeg(2))/(eDeg(2)-eDeg(1))**4*LOG(ABS((z-eDeg(2))/(z-eDeg(1)))) &
                            + 3.0/2.0*(z-eDeg(2))*(2*(z-eDeg(1))-eDeg(1)+eDeg(2))/(eDeg(2)-eDeg(1))**3+1.0/(eDeg(2)-eDeg(1))
               ENDDO
            ELSE
               DO ie = lowerCut, upperCut
                  z = eMesh_shifted(ie)
                  weights(ie) = (z-eDeg(1))**3/(eDeg(2)-eDeg(1))**4*LOG(ABS((z-eDeg(1))/(z-eDeg(2)))) &
                            - (6*(z-eDeg(1))**2+3*(z-eDeg(1))*(eDeg(2)-eDeg(1))+2*(eDeg(2)-eDeg(1))**2)/(6*(eDeg(2)-eDeg(1))**3)
               ENDDO
            ENDIF
         ENDIF
      ELSE IF(ndeg.EQ.6) THEN
         !Eq. A6
         weights(lowerCut:upperCut) = 0.25/(eMesh_shifted(lowerCut:upperCut)-eDeg(1))
      ENDIF
   END FUNCTION resWeightBulk

   PURE REAL FUNCTION resWeightFilm(z,e,ind)

      !Calculates resolvent weights for point in triangle (Compare Solid State Phys. 20 (1987) 4823-4831.)

      REAL,          INTENT(IN)  :: z
      REAL,          INTENT(IN)  :: e(:)
      INTEGER,       INTENT(IN)  :: ind

      REAL, PARAMETER :: tol = 1e-8!tolerance for degeneracy
      REAL, PARAMETER :: fac = 4.0

      REAL denom,a,b,cut,min,prod
      INTEGER ndeg,i,j,k,l,m
      INTEGER ideg(2,2)

      min = 9e+20
      cut = 0.0
      DO i = 1, 2
         DO j=i+1,3
            IF(ABS(e(i)-e(j)).GT.cut) cut = MAX(ABS(e(i)-e(j)),tol)
         ENDDO
      ENDDO

      DO i =1, 3
         IF(ABS(z-e(i)).LT.min) min = ABS(z-e(i))
      ENDDO
      ndeg = 0
      ideg = 0
      IF(min.GT.fac*cut) THEN
         !asymptotic relation Eqs. 9-11
         !the formulas below are unstable for big z arguments
         a = (e(ind) + SUM(e(:)))/4.0
         b = 3*e(ind)**2
         prod = 1.0
         DO j = 1,3
            IF(j.EQ.ind) CYCLE
            b = b + 2*e(ind)*e(j)+e(j)**2
            prod = prod * e(j)
         ENDDO
         b = (b+prod)/300.0
         resWeightFilm = 1.0/3.0/(z-a-b/z)
         !The asymptotic relationship has a divergence at z=0
         !For now we drop the 1/z term in the denominator here
         IF(ABS(z).LT.2e-2) resWeightFilm = 1.0/3.0/(z-a)
      ELSE
         !Search for degenerate eigenvalues
         DO i = 1, 2
            DO j = i+1,3
               IF(ABS(e(i)-e(j)).LT.tol) THEN
                  ndeg = ndeg + 1
                  ideg(ndeg,1) = i
                  ideg(ndeg,2) = j
               ENDIF
            ENDDO
         ENDDO
         resWeightFilm = 0.0
         IF(ndeg.EQ.0) THEN

            !all eigenvalues are non degenerate
            denom = 1.0
            DO i = 1, 3
               IF(i.EQ.ind) CYCLE
               denom = denom*(e(i)-e(ind))
            ENDDO
            !First term
            resWeightFilm = (z-e(ind))/denom
            DO j = 1, 3
               IF(j.EQ.ind) CYCLE
               denom = 1.0
               DO i = 1, 3
                  IF(i.EQ.j) CYCLE
                  denom = denom*(e(i)-e(j))
               ENDDO
               resWeightFilm = resWeightFilm + (z-e(j))**2/denom*LOG(ABS((z-e(j))/(z-e(ind))))/(e(ind)-e(j))
            ENDDO
         ELSE IF(ndeg.EQ.1) THEN

            IF(ALL(ideg(:,:).NE.3)) THEN
               IF(ind.NE.3) THEN
                  resWeightFilm =  (z-e(3))**2/((e(1)-e(3))**3)*LOG(ABS((z-e(3))/(z-e(1)))) &
                                  -(z-e(1))/(e(3)-e(1))**2+1.5/(e(3)-e(1))
               ELSE
                  resWeightFilm =  (z-e(3))/((e(1)-e(3))*(e(2)-e(3))) &
                                  - 2*(z-e(1))*(z-e(3))/(e(3)-e(1))**3 * LOG(ABS((z-e(1))/(z-e(3)))) &
                                   +(z-e(1))/(e(3)-e(1))**2
               ENDIF
            ELSE
               IF(ind.EQ.1) THEN
                  resWeightFilm =  (z-e(1))/((e(2)-e(1))*(e(3)-e(1))) &
                                  - 2*(z-e(1))*(z-e(3))/(e(1)-e(3))**3 * LOG(ABS((z-e(3))/(z-e(1)))) &
                                   +(z-e(3))/(e(1)-e(3))**2
               ELSE
                  resWeightFilm =  (z-e(1))**2/((e(3)-e(1))**3)*LOG(ABS((z-e(1))/(z-e(3)))) &
                                  -(z-e(1))/(e(1)-e(3))**2+0.5/(e(1)-e(3))
               ENDIF
            ENDIF

         ELSE
            resWeightFilm = 1.0/3.0 * 1/(z-e(1))
         ENDIF
      ENDIF

   END FUNCTION resWeightFilm

END MODULE m_resWeight