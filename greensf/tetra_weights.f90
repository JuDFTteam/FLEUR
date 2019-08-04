MODULE m_tetra_weights

   !------------------------------------------------------------------------------
   !
   ! MODULE: m_tetra_weights
   !
   !> @author
   !> Henning Janßen
   !
   ! DESCRIPTION: 
   !>  calculates the weights for one k-point in the tetrahedron method 
   !>  the weights are calculated according to PhysRevB.49.16223 
   !>  When used to calculate a DOS we need to differentiate with respect to energy
   !
   !------------------------------------------------------------------------------

   CONTAINS

   SUBROUTINE tria_weights(ikpt,kpts,ntria,vol,itria,voltria,neig,eig,g,weights,e_ind,ef)

      USE m_types 
      USE m_constants
      USE m_juDFT


      TYPE(t_kpts),           INTENT(IN)     :: kpts
      INTEGER,                INTENT(IN)     :: ikpt
      INTEGER,                INTENT(IN)     :: ntria
      INTEGER,                INTENT(IN)     :: itria(3,2*kpts%nkpt)
      REAL,                   INTENT(IN)     :: vol 
      REAL,                   INTENT(IN)     :: voltria(2*kpts%nkpt)
      INTEGER,                INTENT(IN)     :: neig
      REAL,                   INTENT(IN)     :: eig(:,:)
      TYPE(t_greensfCoeffs),  INTENT(IN)     :: g

      REAL,                   INTENT(OUT)    :: weights(:,:)
      INTEGER,                INTENT(OUT)    :: e_ind(:,:)
      REAL,                   INTENT(IN)     :: ef

      !e_ind(:,1) = g%ne
      !e_ind(:,2) = 0
      !max_ib = 0
      !weights = 0.0
      !
      !DO itr = 1, ntria
!
      !   IF(ALL(itria(1:3,itr).NE.ikpt)) CYCLE 
!
      !   max_ib_tria = MINVAL(neig(itria(1:3,itr)))
      !   IF(max_ib_tet.GT.max_ib) max_ib = max_ib_tet
      !   !$OMP PARALLEL DEFAULT(none) &
      !   !$OMP SHARED(ikpt,itr,ef,max_ib_tria) &
      !   !$OMP SHARED(kpts,eig,g,weights,itria,ntria,atr,as) &
      !   !$OMP PRIVATE(ib,ie,i,j,nstart,nend,icorn,weight,tmp) &
      !   !$OMP PRIVATE(ind,e)
!
      !   !$OMP DO
      !   DO ib = 1, max_ib_tet!Only go up to the band index which exists at all 
      !      
      !      !check if the band is inside the energy window
      !      IF((MINVAL(eig(ib,itria(1:3,itr))).GT.g%e_top).OR.(MAXVAL(eig(ib,itria(1:3,itr))).LT.g%e_bot)) CYCLE !Maybe cancel the band loop completely if we go above top
!
      !      ind(1:3) = (/1,2,3/) !array for sorting the energies
      !      e(1:3) = eig(ib,itria(1:3,itr)) 
!
      !      !Sort the energies in the tetrahedron in ascending order
      !      DO i = 1, 2
      !         DO j = i+1, 3
      !            IF (e(ind(i)).GT.e(ind(j))) THEN
      !               tmp = ind(i)
      !               ind(i) = ind(j)
      !               ind(j) = tmp
      !            ENDIF
      !         ENDDO
      !      ENDDO
!
      !      !check for energy degeneracies 
      !      DO i = 1, 2
      !         DO j = i+1, 3
      !            IF(abs(e(ind(i))-e(ind(j))).LT.1.0E-9) THEN
      !               e(ind(i)) = e(ind(i)) - i*1.0E-9
      !               e(ind(j)) = e(ind(j)) + i*1.0E-9
      !            ENDIF
      !         ENDDO
      !      ENDDO
!
      !      !search for the corner ikpt in the sorted array
      !      DO i = 1, 3
      !         IF(itria(ind(itr),i).EQ.ikpt) icorn = i
      !      ENDDO
!
      !      !Below this index the weight is 0
      !      nstart = INT((e(ind(1))-g%e_bot)/g%del)+1
      !      nend = g%ne
!
      !      DO ie = MAX(1,nstart), g%ne
      !         CALL dos_tria
      !         weights(ie,ib) = weights(ie,ib) + weight
      !         IF(weight.EQ.1/4.0*vol) THEN
      !            !here we are above all eigenergies at the corners 
      !            !of the tetrahedron and we can simplify the rest of the loop
      !            nend = ie
      !            EXIT
      !         ENDIF
      !      ENDDO 
      !      IF(nend.NE.g%ne) weights(nend+1:g%ne,ib) = weights(nend+1:g%ne,ib) + 1/4.0*vol
!
      !   ENDDO
      !   !$OMP ENDDO
      !   !$OMP END PARALLEL
      !ENDDO    


   END SUBROUTINE tria_weights

   SUBROUTINE tetra_weights(ikpt,kpts,neig,eig,g,weights,e_ind,ef)

      USE m_types
      USE m_constants
      USE m_differentiate
      USE m_juDFT

      IMPLICIT NONE

      INTEGER,                INTENT(IN)     :: ikpt
      TYPE(t_kpts),           INTENT(IN)     :: kpts
      INTEGER,                INTENT(IN)     :: neig
      REAL,                   INTENT(IN)     :: eig(:,:)
      TYPE(t_greensfCoeffs),  INTENT(IN)     :: g

      REAL,                   INTENT(OUT)    :: weights(:,:)
      INTEGER,                INTENT(OUT)    :: e_ind(:,:)
      REAL,                   INTENT(IN)     :: ef

      !Local Scalars
      INTEGER itet, ie, ib,i, j,icorn, nstart,nend,tmp
      REAL    weight,dweight,tol,vol
      LOGICAL l_bloechl

      !Local Arrays
      REAL,ALLOCATABLE :: dos_weights(:)
      REAL             :: e(4)
      INTEGER          :: ind(4)

      l_bloechl = .TRUE.
      ALLOCATE( dos_weights(g%ne) )
      dos_weights = 0.0
      e_ind(:,1) = g%ne+1
      e_ind(:,2) = 0
      max_ib = 0
      weights = 0.0
      DO itet = 1, kpts%ntet

         IF(ALL(kpts%ntetra(1:4,itet).NE.ikpt)) CYCLE !search for the tetrahedra containing ikpt

         !$OMP PARALLEL DEFAULT(none) &
         !$OMP SHARED(ikpt,itet,ef,l_bloechl,neig) &
         !$OMP SHARED(kpts,eig,g,weights) &
         !$OMP PRIVATE(ib,ie,i,j,nstart,nend,icorn,weight,dweight,tmp,vol) &
         !$OMP PRIVATE(ind,e)

         !$OMP DO
         DO ib = 1, neig!Only go up to the band index which exists at all 
            
            !check if the band is inside the energy window
            IF((MINVAL(eig(ib,kpts%ntetra(1:4,itet))).GT.g%e_top).OR.(MAXVAL(eig(ib,kpts%ntetra(1:4,itet))).LT.g%e_bot)) CYCLE !Maybe cancel the band loop completely if we go above top

            ind(1:4) = (/1,2,3,4/) !array for sorting the energies
            e(1:4) = eig(ib,kpts%ntetra(1:4,itet)) 

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

            !check for energy degeneracies 
            DO i = 1, 3
               DO j = i+1, 4
                  IF(abs(e(ind(i))-e(ind(j))).LT.1.0E-9) THEN
                     e(ind(i)) = e(ind(i)) - i*1.0E-9
                     e(ind(j)) = e(ind(j)) + i*1.0E-9
                  ENDIF
               ENDDO
            ENDDO

            !search for the corner ikpt in the sorted array
            DO i = 1, 4
               IF(kpts%ntetra(ind(i),itet).EQ.ikpt) icorn = i
            ENDDO

            !Below this index the weight is 0
            nstart = INT((e(ind(1))-g%e_bot)/g%del)+1
            vol = kpts%voltet(itet)/kpts%ntet
            nend = g%ne
            IF(l_bloechl) CALL bloechl_corrections(ef,vol,e(ind(:)),dweight,icorn)
            DO ie = MAX(1,nstart), g%ne
               CALL contrib_singletetra((ie-1)*g%del+g%e_bot,vol,e(ind(:)),weight,icorn)
               weights(ie,ib) = weights(ie,ib) + weight
               IF(weight.EQ.1/4.0*vol) THEN
                  !here we are above all eigenergies at the corners 
                  !of the tetrahedron and we can simplify the rest of the loop
                  nend = ie
                  EXIT
               ENDIF
            ENDDO 
            IF(nend.NE.g%ne) weights(nend+1:g%ne,ib) = weights(nend+1:g%ne,ib) + 1/4.0*vol
         ENDDO
         !$OMP END DO
         !$OMP END PARALLEL
      ENDDO
      !
      !Differentiate the weights with respect to energy to get the correct weights for a DOS-calculation
      !
      tol = 1E-14
      !$OMP PARALLEL DEFAULT(none) &
      !$OMP SHARED(weights,g,e_ind,eig,ikpt,tol,neig) &
      !$OMP PRIVATE(ib,i,dos_weights) 

      !$OMP DO
      DO ib = 1, neig
         dos_weights = 0.0
         !Imaginary part of the weights
         CALL diff3(weights(:,ib),g%del,dos_weights(:))
 
         IF(.FALSE.) THEN
            !Weights for spin-offdiagonal part 
            !-----------------------------------------------------
            ! Here we have two things contributing to the imaginary part of the GF:
            !    1. Imaginary Part of the denominator * Real part of the coefficients
            !        This corresponds to the weights calculated above
            !    2. Real Part of the denominator * Imaginary part of the coefficients
            !------------------------------------------------------
            ! For the second part we need to multiply the integration weights with 
            !  1/(E-E_kv)
            !------------------------------------------------------
            DO ie = 1, g%ne
               weights(ie,ib) = weights(ie,ib) / ((ie-1)*g%del+g%e_bot-eig(ib,ikpt)) +ImagUnit * dos_weights(ie)
            ENDDO
         ENDIF

         !
         !Find the range where the weights are not equal to 0
         !
         !Start
         i = 1
         DO
            IF(dos_weights(i).GT.tol) THEN
               e_ind(ib,1) = i
               EXIT
            ELSE
               dos_weights(i) = 0.0
               i = i + 1
               IF(i.EQ.g%ne+1) THEN
                  e_ind(ib,1) = g%ne
                  EXIT
               ENDIF
            ENDIF
         ENDDO
         !End
         i = g%ne
         DO
            IF(dos_weights(i).GT.tol) THEN
               e_ind(ib,2) = i
               EXIT
            ELSE
               dos_weights(i) = 0.0
               i = i - 1
               IF(i.EQ.0) THEN
                  e_ind(ib,2) = 1
                  EXIT
               ENDIF
            ENDIF
         ENDDO
         IF(ANY(dos_weights(:).LT.0)) THEN
            CALL juDFT_error("Negative tetra weight",calledby="tetra_weights")
         ENDIF
         IF(e_ind(ib,1).GT.e_ind(ib,2)) THEN
            e_ind(ib,1) = 1
            e_ind(ib,2) = 1
         ENDIF
         weights(:,ib) = dos_weights(:)
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL
   END SUBROUTINE tetra_weights

   SUBROUTINE contrib_singletetra(energy,vol,e,weight,ind)

      USE m_juDFT
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

      ELSE IF(energy.LT.e(1)) THEN

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

   END SUBROUTINE contrib_singletetra

   SUBROUTINE bloechl_corrections(ef,vol,e,dweight,ind)

      REAL,                INTENT(IN)     :: ef
      REAL,                INTENT(IN)     :: vol
      REAL,                INTENT(IN)     :: e(4)
      REAL,                INTENT(OUT)    :: dweight
      INTEGER,             INTENT(IN)     :: ind

      INTEGER i
      REAL dos

      CALL dos_tetra(ef,vol,e,dos)

      dweight = 0.0
      DO i = 1, 4
         dweight = dweight + 1/40.0*dos*(e(i)-e(ind))
      ENDDO


   END SUBROUTINE bloechl_corrections

   SUBROUTINE dos_tetra(energy,vol,e,dos)

      IMPLICIT NONE

      REAL,                INTENT(IN)     :: energy
      REAL,                INTENT(IN)     :: vol
      REAL,                INTENT(IN)     :: e(4)
      REAL,                INTENT(OUT)    :: dos

      IF((energy.GT.e(4)).OR.(energy.LT.e(1))) THEN
         dos = 0.0
      ELSE IF((energy.GT.e(1)).AND.(energy.LT.e(2))) THEN

         dos = 3.0 * vol * (energy-e(1))**2/((e(2)-e(1))*(e(3)-e(1))*(e(4)-e(1)))

      ELSE IF((energy.GT.e(2)).AND.(energy.LT.e(3))) THEN

         dos = 3.0 * vol * 1./((e(3)-e(1))*(e(4)-e(1))) * (e(2) - e(1) + 2*(energy - e(2)) -&
               (e(3)-e(1)+e(4)-e(2)) * (energy-e(2))**2/((e(3)-e(2))*(e(4)-e(2))))

      ELSE IF((energy.GT.e(3)).AND.(energy.LT.e(4))) THEN

         dos = 3.0 * vol * (e(4)-energy)**2/((e(4)-e(1))*(e(4)-e(2))*(e(4)-e(3)))
      END IF
   END SUBROUTINE dos_tetra


   ! Not used at the moment
      SUBROUTINE getTetraContrib(energy,vol,ev,w,g)

      !Integration weights taken from 
      !PhysRevB.49.16223

      USE m_types

      IMPLICIT NONE

      REAL,                   INTENT(IN)     :: energy
      REAL,                   INTENT(IN)     :: vol
      REAL,                   INTENT(IN)     :: ev(4)
      REAL,                   INTENT(INOUT)  :: w(:)
      TYPE(t_greensfCoeffs),  INTENT(IN)     :: g

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