MODULE m_resWeights

   CONTAINS

   SUBROUTINE resWeightsCalc(ikpt,kpts,neig,eig,g,weights,boundInd)

      !Weights for analytical tetrahedron method for spectral functions 

      USE m_types
      USE m_juDFT

      IMPLICIT NONE

      INTEGER,          INTENT(IN)     :: ikpt
      TYPE(t_kpts),     INTENT(IN)     :: kpts
      INTEGER,          INTENT(IN)     :: neig
      REAL,             INTENT(IN)     :: eig(neig,kpts%nkpt)
      TYPE(t_greensfCoeffs), INTENT(IN) :: g
      REAL,          INTENT(INOUT)  :: weights(g%ne,neig)
      INTEGER,          INTENT(INOUT)  :: boundInd(neig,2)


      INTEGER itet,ib,i,j,iz,icorn,ind(4),k(4),tmp
      REAL e(4),vol,fac
      REAL weight

      !Here we do no truncation for now
      boundInd(:,1) = 1
      boundInd(:,2) = g%ne
      weights = 0.0

      fac = count(kpts%bkp(:).EQ.ikpt)

      DO itet = 1, kpts%ntet

         IF(ALL(kpts%ntetra(1:4,itet).NE.ikpt)) CYCLE
         IF(kpts%nkptf.NE.0) THEN
            DO i = 1, 4
               IF(kpts%ntetra(i,itet).GT.kpts%nkpt) THEN
                  k(i) = kpts%bkp(kpts%ntetra(i,itet))
               ELSE
                  k(i) = kpts%ntetra(i,itet)
               ENDIF
            ENDDO
         ELSE
            k(1:4) = kpts%ntetra(1:4,itet)
         ENDIF

         !$OMP PARALLEL DEFAULT(none) &
         !$OMP SHARED(ikpt,itet,neig,g,k,fac) &
         !$OMP SHARED(kpts,eig,weights) &
         !$OMP PRIVATE(ib,iz,i,j,icorn,tmp,vol) &
         !$OMP PRIVATE(ind,e,weight)

         !$OMP DO
         DO ib = 1, neig

            e(1:4) = eig(ib,k(1:4)) 
            ind=(/1,2,3,4/)
            !Sort the energies in the tetrahedron in ascending order
            DO i = 1, 3
               DO j = i+1, 4
                  IF (e(ind(i)).GT.e(ind(j))) THEN
                     tmp = ind(i)
                     ind(i) = ind(j)
                     ind(j) = tmp
                  ENDIF
               ENDDO
            ENDDO
            !search for the corner ikpt in the sorted array
            DO i = 1, 4
               IF(kpts%ntetra(ind(i),itet).EQ.ikpt) icorn = i
            ENDDO

            vol = kpts%voltet(itet)/kpts%ntet
            DO iz = 1, g%ne
               CALL resWeightTetra((iz-1)*g%del+g%e_bot,e(ind(1:4)),vol,weight,icorn)
               IF(ISNAN(weight)) CALL juDFT_error("Tetra Weight NaN",calledby="resWeights")
               weights(iz,ib) = weights(iz,ib) + weight * fac
            ENDDO
         ENDDO
         !$OMP END DO
         !$OMP END PARALLEL
      ENDDO

   END SUBROUTINE resWeightsCalc

   SUBROUTINE resWeightTetra(z,e,vol,weight,ind)

      !Calculates resolvent weights for point in tetrahedron (Compare PhysRevB.29.3430)

      USE m_juDFT

      IMPLICIT NONE 

      REAL,       INTENT(IN)  :: z
      REAL,          INTENT(IN)  :: e(4)
      REAL,          INTENT(IN)  :: vol
      REAL,       INTENT(OUT) :: weight 
      INTEGER,       INTENT(IN)  :: ind

      REAL tol,denom,a,b,cut,min,fac
      INTEGER ndeg,i,j,k,l,m
      INTEGER ideg(6,2)

      tol = 1e-9 !Tolerance for degeneracy
      fac = 50.0


      min = 9e+20
      cut = 9e+20
      DO i = 1, 3
         DO j=i+1,4
            IF(ABS(e(i)-e(j)).LT.cut.AND.ABS(e(i)-e(j)).GT.tol) cut = ABS(e(i)-e(j))
         ENDDO
      ENDDO

      DO i =1, 4
         IF(ABS(z-e(i)).LT.min) min = ABS(z-e(i)) 
      ENDDO
      ndeg = 0
      ideg = 0
      
      IF(min.GT.fac*cut) THEN
         !asymptotic relation Eqs. 9-11
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
         weight = vol*0.25/(z-a-b/z)
      ELSE
         DO i = 1, 3
            DO j = i+1,4
               IF(ABS(e(i)-e(j)).LT.tol) THEN
                  ndeg = ndeg + 1
                  IF(ndeg.GT.6) CALL juDFT_error("Dim Error ndeg>6",calledby="resWeightTetra")
                  ideg(ndeg,1) = i
                  ideg(ndeg,2) = j
               ENDIF
            ENDDO
         ENDDO
         weight = 0.0
         IF(ndeg.EQ.0) THEN

            !all eigenvalues are non degenerate (EQ.7 from the paper)
            denom = 1.0
            DO i = 1, 4
               IF(i.EQ.ind) CYCLE
               denom = denom*(e(i)-e(ind))
            ENDDO
            !First term
            weight = vol*(z-e(ind))**2/denom
            DO j = 1, 4
               IF(j.EQ.ind) CYCLE
               denom = 1.0
               DO i = 1, 4
                  IF(i.EQ.j) CYCLE
                  denom = denom*(e(i)-e(j))
               ENDDO
               weight = weight + vol*(z-e(j))**3/denom*LOG(ABS((z-e(j))/(z-e(ind))))/(e(ind)-e(j))
            ENDDO

         ELSE IF(ndeg.EQ.1) THEN

            !Two eigenvalues are degenerate (EQs. A1 A2)

            l = ideg(1,1)
            m = ideg(1,2)

            IF(ind.EQ.l.OR.ind.EQ.m) THEN
               !EQ. A2
               j = 0
               DO i = 1, 4
                  IF(i.EQ.l.OR.i.EQ.m) CYCLE
                  j = i
               ENDDO
               IF(j.EQ.0) CALL juDFT_error("j not found",calledby="resWeightTetra")
               k = 0
               DO i = 1, 4
                  IF(i.EQ.l.OR.i.EQ.m.OR.i.EQ.j) CYCLE
                  k = i
               ENDDO
               IF(k.EQ.0) CALL juDFT_error("k not found",calledby="resWeightTetra")
               weight = vol*(z-e(k))**3/((e(k)-e(j))*(e(k)-e(m))**3)*LOG(ABS(z-e(k))) & 
                        +vol*(z-e(j))**3/((e(j)-e(k))*(e(j)-e(m))**3)*LOG(ABS(z-e(j))) & 
                        +vol*(z-e(m))/((e(m)-e(j))*(e(m)-e(k))) * (0.5 + (z-e(j))/(e(m)-e(j))&
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
               IF(k.EQ.0) CALL juDFT_error("k not found",calledby="resWeightTetra")
               !EQ. A1
               weight = vol*(z-e(ind))**2/((e(ind)-e(m))**2*(e(k)-e(ind))) *&
                           (1 + (2*(z-e(m))/(e(ind)-e(m))+(z-e(k))/(e(ind)-e(k)))*LOG(ABS(z-e(ind)))) &
                        +vol*(z-e(m))**2/((e(m)-e(ind))**2*(e(k)-e(m))) *&
                           (1 + (2*(z-e(ind))/(e(m)-e(ind))+(z-e(k))/(e(m)-e(k)))*LOG(ABS(z-e(m)))) &
                        +vol*(z-e(k))**3/((e(k)-e(ind))**2*(e(k)-e(m))**2)*LOG(ABS(z-e(k)))
            ENDIF
         ELSE IF(ndeg.EQ.2) THEN
            !This is the case E1=E2<E3=E4 => A4 
            IF(ind.LE.2) THEN
               weight = vol*3*(z-e(3))**2*(z-e(2))/(e(3)-e(2))**4*LOG(ABS((z-e(2))/(z-e(3)))) &
                        - vol*3.0/2.0*(z-e(2))*(2*(z-e(3))-e(3)+e(2))/(e(3)-e(2))**3-vol/(e(3)-e(2))
            ELSE
               weight = vol*3*(z-e(2))**2*(z-e(3))/(e(3)-e(2))**4*LOG(ABS((z-e(3))/(z-e(2)))) &
                        + vol*3.0/2.0*(z-e(3))*(2*(z-e(2))+e(3)-e(2))/(e(3)-e(2))**3+vol/(e(3)-e(2))
            ENDIF
         ELSE IF(ndeg.GE.3.AND.ndeg.LT.6) THEN
            !EQs A3/A5 (here we explicitly write each weight)
            IF(ALL(ideg(:,:).NE.4)) THEN
               !A3
               IF(ind.NE.4) THEN
                  weight = vol*(z-e(4))**3/(e(4)-e(3))**4*LOG(ABS((z-e(4))/(z-e(3)))) &
                           + vol*(6*(z-e(4))**2-3*(e(4)-e(3))*(z-e(4))+2*(e(4)-e(3))**2)/(6*(e(4)-e(3))**3)
               ELSE
                  weight = vol*3*(z-e(4))**2*(z-e(3))/(e(4)-e(3))**4*LOG(ABS((z-e(3))/(z-e(4)))) &
                           - vol*3.0/2.0*(z-e(3))*(2*(z-e(4))-e(4)+e(3))/(e(4)-e(3))**3-vol/(e(4)-e(3))
               ENDIF
            ELSE IF(ALL(ideg(:,:).NE.1)) THEN
               !A5
               IF(ind.EQ.1) THEN
                  weight = vol*3*(z-e(1))**2*(z-e(2))/(e(2)-e(1))**4*LOG(ABS((z-e(2))/(z-e(1)))) &
                           + vol*3.0/2.0*(z-e(2))*(2*(z-e(1))-e(1)+e(2))/(e(2)-e(1))**3+vol/(e(2)-e(1))
               ELSE
                  weight = vol*(z-e(1))**3/(e(2)-e(1))**4*LOG(ABS((z-e(1))/(z-e(2)))) &
                           - vol*(6*(z-e(1))**2+3*(z-e(1))*(e(2)-e(1))+2*(e(2)-e(1))**2)/(6*(e(2)-e(1))**3)
               ENDIF
            ENDIF
         ELSE IF(ndeg.EQ.6) THEN
            !Eq. A6
            weight = vol/4.0 * 1.0/(z-e(1))
         ENDIF
      ENDIF

   END SUBROUTINE


END MODULE m_resWeights