MODULE m_tetra_weights

   CONTAINS

   SUBROUTINE tetra_weights(ikpt,kpts,neig,eig,g,weights,ef)

      USE m_types
      USE m_constants

      INTEGER,             INTENT(IN)     :: ikpt
      TYPE(t_kpts),        INTENT(IN)     :: kpts
      INTEGER,             INTENT(IN)     :: neig(:)
      REAL,                INTENT(IN)     :: eig(:,:)
      TYPE(t_greensf),     INTENT(IN)     :: g

      REAL,                INTENT(OUT)    :: weights(:,:)
      REAL,                INTENT(IN)     :: ef


      INTEGER icorn, itet, ib, j, k, l, nstart,corn_ind

      REAL e(4), weight

      INTEGER ind(4),tmp


      DO itet = 1, kpts%ntet

         IF(ALL(kpts%ntetra(1:4,itet).NE.ikpt)) CYCLE !search for the tetrahedra containing ikpt
         

         DO ib = 1, MINVAL(neig(kpts%ntetra(1:4,itet))) !Only go up to the band index which exists at all corners

            IF((MINVAL(eig(ib,kpts%ntetra(1:4,itet))).GT.g%e_top).OR.(MAXVAL(eig(ib,kpts%ntetra(1:4,itet))).LT.g%e_bot)) CYCLE !Maybe cancel the band loop completely if we go above top

            ind(1:4) = (/1,2,3,4/) !array for sorting the energies
            e(:) = eig(ib,kpts%ntetra(1:4,itet)) 


            DO l = 1, 3
               DO k = l+1, 4
                  IF (e(ind(l)).GT.e(ind(k))) THEN
                     tmp = ind(l)
                     ind(l) = ind(k)
                     ind(k) = tmp
                  ENDIF
               ENDDO
            ENDDO

            !check for energy degeneracies 
            DO l = 1, 3
               DO k = l+1, 4
                  IF(abs(e(ind(l))-e(ind(k))).LT.1.0E-9) THEN
                     e(ind(l)) = e(ind(l)) - l*1.0E-9
                     e(ind(k)) = e(ind(k)) + l*1.0E-9
                  ENDIF
               ENDDO
            ENDDO

            !weight = 0.0
            !CALL getTetraContrib(ef,kpts%voltet(itet)/kpts%ntet,e(:),weight(:),g)

            !search for the corner ikpt in the sorted array

            DO l = 1, 4
               IF(kpts%ntetra(ind(l),itet).EQ.ikpt) corn_ind = ind(l)
            ENDDO


            nstart = INT((e(1)-g%e_bot)/g%del)+1
            DO ie = MAX(1,nstart), g%ne

               CALL contrib_singletetra((ie-1)*g%del+g%e_bot,kpts%voltet(itet)/kpts%ntet,e(ind(:)),weight,corn_ind)
               weights(ie,ib) = weights(ie,ib) + weight

            ENDDO

         ENDDO
      ENDDO

   END SUBROUTINE tetra_weights

   SUBROUTINE contrib_singletetra(energy,vol,e,weight,ind)

      !Integration weights taken from 
      !PhysRevB.49.16223 

      IMPLICIT NONE

      REAL,                INTENT(IN)     :: energy
      REAL,                INTENT(IN)     :: vol
      REAL,                INTENT(IN)     :: e(4)
      REAL,                INTENT(OUT)    :: weight
      INTEGER,             INTENT(IN)     :: ind

      INTEGER j

      REAL C
      REAL C1, C2, C3

      IF(energy.GT.e(4)) THEN

         weight = 1/4.*vol

      ELSE IF(energy.GT.e(1)) THEN

         weight = 0.0

      ELSE IF((e(1).LE.energy).AND.(energy.LE.e(2))) THEN

         C     = 1./4. * vol * (energy-e(1))**3/((e(2)-e(1))*(e(3)-e(1))*(e(4)-e(1)))

         IF(ind.EQ.1) THEN
            weight = C * (4 - (energy-e(1)) * (1/(e(2)-e(1)) + 1/(e(3)-e(1)) + 1/(e(4)-e(1))))
         ELSE 
            weight = C * (energy-e(1))/(e(ind)-e(1)) 
         ENDIF
      ELSE IF((e(2).LE.energy).AND.(energy.LE.e(3))) THEN

         C1     = 1./4. * vol * (energy-e(1))**2/((e(3)-e(1))*(e(4)-e(1)))
      
         C2     = 1./4. * vol * (energy-e(1))*(energy-e(2))*(e(3)-energy)/((e(3)-e(1))*(e(3)-e(2))*(e(4)-e(1)))


         C3     = 1./4. * vol * (energy-e(2))**2*(e(4)-energy)/((e(3)-e(2))*(e(4)-e(1))*(e(4)-e(2)))

         IF(ind.EQ.1) THEN
            weight = C1 + (C1 + C2) * (e(3)-energy)/(e(3)-e(1))+&
                        (C1 + C2 + C3) * (e(4)-energy)/(e(4)-e(1)) 
         ELSE IF(ind.EQ.2) THEN
            weight = C1 + C2 + C3 + (C2 + C3) * (e(3)-energy)/(e(3)-e(2)) &
                        + C3 * (e(4)-energy)/(e(4)-e(2))
         ELSE IF(ind.EQ.3) THEN
            weight = (C1 + C2) * (energy-e(1))/(e(3)-e(1)) +&
                        (C2 + C3) * (energy-e(2))/(e(3)-e(2)) 
         ELSE IF(ind.EQ.4) THEN
            weight = (C1 + C2 + C3) * (energy-e(1))/(e(4)-e(1)) +&
                        C3 * (energy-e(2))/(e(4)-e(2))
         END IF
      ELSE IF((e(3).LE.energy).AND.(energy.LE.e(4))) THEN

         C     = 1./4. * vol * (e(4)-energy)**3/((e(4)-e(1))*(e(4)-e(2))*(e(4)-e(3)))

         IF(ind.LE.3) THEN
            weight = 1/4.*vol-C * (e(4)-energy)/(e(4)-e(ind))
         ELSE
            weight = 1/4.*vol-C * (4 - (e(4)-energy) * (1/(e(4)-e(1)) + 1/(e(4)-e(2)) + 1/(e(4)-e(3))))
         END IF

      END IF

      !IF(ISNAN(weight(1)).OR.ISNAN(weight(2)).OR.ISNAN(weight(3)).OR.ISNAN(weight(4))) THEN
      !   WRITE(*,*) "ENERGIES: ", e(1), e(2), e(3), e(4)
      !   WRITE(*,*) "WEIGHTS: ",weight(1), weight(2), weight(3), weight(4)
      !END IF

   END SUBROUTINE contrib_singletetra


   ! Not used at the moment
      SUBROUTINE getTetraContrib(energy,vol,ev,w,g)

      !Integration weights taken from 
      !PhysRevB.49.16223

      USE m_types

      IMPLICIT NONE

      REAL,                INTENT(IN)     :: energy
      REAL,                INTENT(IN)     :: vol
      REAL,                INTENT(IN)     :: ev(4)
      REAL,                INTENT(INOUT)    :: w(:)
      TYPE(t_greensf),     INTENT(IN)     :: g

      INTEGER i , j, n
      REAL emin, emax, e
      REAL weight 
      REAL a,b,c
      INTEGER ind(3)
      REAL deln

      REAL aw, bw
      INTEGER nstart, nend


      ind(1:3) = (/1,2,4/)

      DO i = 1, 3
         emin = ev(i)
         emax = ev(i+1)

         !determine the last index with a lower energy than emin/max
         nstart = INT((emin-g%e_bot)/g%del)+1
         nend = INT((emax-g%e_bot)/g%del)+1

         !distance between e(nstart+1) and emin 
         aw = g%del - (emin - ((nstart-1)*g%del+g%e_bot))
         !distance between e(nend) and emax
         bw = emax - ((nend-1)*g%del+g%e_bot)

         weight = 3.0 * vol
         DO j = 1, 4
            IF (ind(i).NE.j) weight = weight/ABS(ev(ind(i))-ev(j))
         ENDDO

         !calculate the right prefactors
         IF(i.EQ.2) THEN
            a = a - emin**2 * weight
            b = b + emin * weight
            c = c - 1./3. * weight
         ELSE
            IF(i.EQ.3) THEN
               a = emax**2 * weight
               b = (-1.0) * emax * weight
               c = 1./3. * weight
            ELSE 
               a = emin**2 * weight
               b = (-1.0) * emin * weight
               c = 1./3. * weight
            ENDIF
         ENDIF

         !If emin and emax are inside the same energy step we use a different formula (from where?)
         IF(nstart.EQ.nend) THEN
            deln = (emax-emin)/g%del * (a + b*(emax+emin) + c *((emax+emin)**2-emax*emin))
            IF((nstart.LE.g%ne).AND.(nstart.GE.1)) w(nstart) = w(nstart) + deln
         ELSE
            !Part of the integral between emin and the first grid point
            deln = aw/g%del * (a + b*(2*emin+aw) + c*((2*emin+aw)**2-emin*(emin+aw)))
            IF((nstart.LE.g%ne).AND.(nstart.GE.1)) w(nstart) = w(nstart) + deln
            e = emin + aw
            !loop over grid points between emin and emax
            DO n = nstart + 1, nend-1 
               deln = (a + b*(2*e+g%del) + c*((2*e+g%del)**2-e*(e+g%del)))
               IF((n.LE.g%ne).AND.(n.GE.1)) w(n) = w(n) + deln
               e = e + g%del
            ENDDO
            !Part of the integral between the last grid point before emax and emax
            deln = bw/g%del * (a + b*(2*e+bw) + c*((2*e+bw)**2-e*(e+bw)))
            IF((nend.LE.g%ne).AND.(nend.GE.1)) w(nend) = w(nend) + deln
         ENDIF
      ENDDO
   END SUBROUTINE getTetraContrib

END MODULE m_tetra_weights