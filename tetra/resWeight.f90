MODULE m_resWeight

   IMPLICIT NONE

   CONTAINS

   PURE REAL FUNCTION resWeight(z,e,vol,ind)

      !Calculates resolvent weights for point in tetrahedron (Compare PhysRevB.29.3430)

      REAL,          INTENT(IN)  :: z
      REAL,          INTENT(IN)  :: e(:)
      REAL,          INTENT(IN)  :: vol
      INTEGER,       INTENT(IN)  :: ind

      REAL, PARAMETER :: tol = 1e-8!tolerance for degeneracy
      REAL, PARAMETER :: fac = 9.0

      REAL denom,a,b,cut,min
      INTEGER ndeg,i,j,k,l,m
      INTEGER ideg(6,2)

      min = 9e+20
      cut = 0.0
      !
      DO i = 1, 3
         DO j=i+1,4
            IF(ABS(e(i)-e(j)).GT.cut) cut = MAX(ABS(e(i)-e(j)),tol)
         ENDDO
      ENDDO

      DO i =1, 4
         IF(ABS(z-e(i)).LT.min) min = ABS(z-e(i))
      ENDDO
      ndeg = 0
      ideg = 0

      IF(min.GT.fac*cut) THEN
         !asymptotic relation Eqs. 9-11
         !the formulas below are unstable for big z arguments
         a = (e(ind) + SUM(e(:)))/5.0
         b = 0.0
         DO j = 1, 4
            IF(j.EQ.ind) CYCLE
            b = b + 3*(e(j)-e(ind))**2
            DO k = 1, 4
               IF(k.EQ.j) CYCLE
               b = b + (e(k)-e(j))**2
            ENDDO
         ENDDO
         b = b/300.0
         resWeight = 0.25/(z-a-b/z)
      ELSE
         !Search for degenerate eigenvalues
         DO i = 1, 3
            DO j = i+1,4
               IF(ABS(e(i)-e(j)).LT.tol) THEN
                  ndeg = ndeg + 1
                  ideg(ndeg,1) = i
                  ideg(ndeg,2) = j
               ENDIF
            ENDDO
         ENDDO
         resWeight = 0.0
         IF(ndeg.EQ.0) THEN

            !all eigenvalues are non degenerate (EQ.7 from the paper)
            denom = 1.0
            DO i = 1, 4
               IF(i.EQ.ind) CYCLE
               denom = denom*(e(i)-e(ind))
            ENDDO
            !First term
            resWeight = (z-e(ind))**2/denom
            DO j = 1, 4
               IF(j.EQ.ind) CYCLE
               denom = 1.0
               DO i = 1, 4
                  IF(i.EQ.j) CYCLE
                  denom = denom*(e(i)-e(j))
               ENDDO
               resWeight = resWeight + (z-e(j))**3/denom*LOG(ABS((z-e(j))/(z-e(ind))))/(e(ind)-e(j))
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

               resWeight = (z-e(k))**3/((e(k)-e(j))*(e(k)-e(m))**3)*LOG(ABS(z-e(k))) &
                          +(z-e(j))**3/((e(j)-e(k))*(e(j)-e(m))**3)*LOG(ABS(z-e(j))) &
                          +(z-e(m))/((e(m)-e(j))*(e(m)-e(k))) * (0.5 + (z-e(j))/(e(m)-e(j))&
                           +(z-e(k))/(e(m)-e(k)) + ((z-e(j))**2/(e(m)-e(j))**2 &
                           +(z-e(k))**2/(e(m)-e(k))**2 +(z-e(j))/(e(m)-e(j))*(z-e(k))/(e(m)-e(k))) &
                           *LOG(ABS(z-e(m))))
            ELSE
               !k is the one site not equal to ind, l or m
               k = 0
               DO i = 1, 4
                  IF(i.EQ.ind.OR.i.EQ.l.OR.i.EQ.m) CYCLE
                  k = i
               ENDDO
               !EQ. A1
               resWeight = (z-e(ind))**2/((e(ind)-e(m))**2*(e(k)-e(ind))) *&
                           (1 + (2*(z-e(m))/(e(ind)-e(m))+(z-e(k))/(e(ind)-e(k)))*LOG(ABS(z-e(ind)))) &
                          +(z-e(m))**2/((e(m)-e(ind))**2*(e(k)-e(m))) *&
                           (1 + (2*(z-e(ind))/(e(m)-e(ind))+(z-e(k))/(e(m)-e(k)))*LOG(ABS(z-e(m)))) &
                          +(z-e(k))**3/((e(k)-e(ind))**2*(e(k)-e(m))**2)*LOG(ABS(z-e(k)))
            ENDIF
         ELSE IF(ndeg.EQ.2) THEN
            !This is the case E1=E2<E3=E4 => A4
            IF(ind.LE.2) THEN
               resWeight = 3.0*(z-e(3))**2*(z-e(2))/(e(3)-e(2))**4*LOG(ABS((z-e(2))/(z-e(3)))) &
                         - 3.0/2.0*(z-e(2))*(2*(z-e(3))-e(3)+e(2))/(e(3)-e(2))**3-1.0/(e(3)-e(2))
            ELSE
               resWeight = 3.0*(z-e(2))**2*(z-e(3))/(e(3)-e(2))**4*LOG(ABS((z-e(3))/(z-e(2)))) &
                         + 3.0/2.0*(z-e(3))*(2*(z-e(2))+e(3)-e(2))/(e(3)-e(2))**3+1.0/(e(3)-e(2))
            ENDIF
         ELSE IF(ndeg.GE.3 .AND. ndeg.LT.6) THEN
            !EQs A3/A5 (here we explicitly write each weight)
            IF(ALL(ideg(:,:).NE.4)) THEN
               !A3
               IF(ind.NE.4) THEN
                  resWeight = (z-e(4))**3/(e(4)-e(3))**4*LOG(ABS((z-e(4))/(z-e(3)))) &
                            + (6*(z-e(4))**2-3*(e(4)-e(3))*(z-e(4))+2*(e(4)-e(3))**2)/(6*(e(4)-e(3))**3)
               ELSE
                  resWeight = 3.0*(z-e(4))**2*(z-e(3))/(e(4)-e(3))**4*LOG(ABS((z-e(3))/(z-e(4)))) &
                            - 3.0/2.0*(z-e(3))*(2*(z-e(4))-e(4)+e(3))/(e(4)-e(3))**3-1.0/(e(4)-e(3))
               ENDIF
            ELSE IF(ALL(ideg(:,:).NE.1)) THEN
               !A5
               IF(ind.EQ.1) THEN
                  resWeight = 3.0*(z-e(1))**2*(z-e(2))/(e(2)-e(1))**4*LOG(ABS((z-e(2))/(z-e(1)))) &
                            + 3.0/2.0*(z-e(2))*(2*(z-e(1))-e(1)+e(2))/(e(2)-e(1))**3+1.0/(e(2)-e(1))
               ELSE
                  resWeight = (z-e(1))**3/(e(2)-e(1))**4*LOG(ABS((z-e(1))/(z-e(2)))) &
                            - (6*(z-e(1))**2+3*(z-e(1))*(e(2)-e(1))+2*(e(2)-e(1))**2)/(6*(e(2)-e(1))**3)
               ENDIF
            ENDIF
         ELSE IF(ndeg.EQ.6) THEN
            !Eq. A6
            resWeight = 0.25/(z-e(1))
         ENDIF
      ENDIF

      resWeight = resWeight * vol

   END FUNCTION resWeight


END MODULE m_resWeight