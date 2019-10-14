MODULE m_dosWeights

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

   SUBROUTINE sortEig(n,e,ind)

      IMPLICIT NONE

      INTEGER,          INTENT(IN)     :: n
      REAL,             INTENT(INOUT)  :: e(:)
      INTEGER,          INTENT(INOUT)  :: ind(:)

      INTEGER i,j
      INTEGER tmp


      DO i = 1, n
         ind(i) = i
      ENDDO

      !Sort the energies in the tetrahedron in ascending order
      DO i = 1, n-1
         DO j = i+1, n
            IF (e(ind(i)).GT.e(ind(j))) THEN
               tmp = ind(i)
               ind(i) = ind(j)
               ind(j) = tmp
            ENDIF
         ENDDO
      ENDDO
      !check for energy degeneracies
      DO i = 1, n-1
         DO j = i+1, n
            IF(abs(e(ind(i))-e(ind(j))).LT.1.0E-9) THEN
               e(ind(i)) = e(ind(i)) - i*1.0E-9
               e(ind(j)) = e(ind(j)) + i*1.0E-9
            ENDIF
         ENDDO
      ENDDO

   END SUBROUTINE sortEig



   SUBROUTINE dosWeightsCalc(ikpt,kpts,neig,eig,g,ef,weights,e_ind)

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
      INTEGER itet, ie, ib,i, icorn, nstart,nend
      REAL    weight,dweight,tol,vol,fac
      LOGICAL l_bloechl

      !Local Arrays
      REAL,ALLOCATABLE :: dos_weights(:)
      REAL             :: e(4)
      INTEGER          :: ind(4), k(4)

      l_bloechl = .FALSE.
      ALLOCATE( dos_weights(g%ne) )
      dos_weights = 0.0
      e_ind(:,1) = 0
      e_ind(:,2) = g%ne+1
      weights = 0.0

      fac = MERGE(1,count(kpts%bkp(:).EQ.ikpt),kpts%nkptf.EQ.0)

      DO itet = 1, kpts%ntet
         IF(ALL(kpts%ntetra(1:4,itet).NE.ikpt)) CYCLE !search for the tetrahedra containing ikpt

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
         !$OMP SHARED(ikpt,itet,ef,l_bloechl,neig,k,fac) &
         !$OMP SHARED(kpts,eig,g,weights) &
         !$OMP PRIVATE(ib,ie,i,nstart,nend,icorn,weight,dweight,vol) &
         !$OMP PRIVATE(ind,e)

         !$OMP DO
         DO ib = 1, neig

            !check if the band is inside the energy window
            IF((MINVAL(eig(ib,k(1:4))).GT.g%e_top).OR.(MAXVAL(eig(ib,k(1:4))).LT.g%e_bot)) CYCLE !Maybe cancel the band loop completely if we go above top


            e(1:4) = eig(ib,k(1:4))
            ind=(/1,2,3,4/)
            CALL sortEig(4,e,ind)
            !Sort the energies in the tetrahedron in ascending order
            !search for the corner ikpt in the sorted array
            DO i = 1, 4
               IF(kpts%ntetra(ind(i),itet).EQ.ikpt) icorn = i
            ENDDO

            !Below this index the weight is 0
            nstart = INT((e(ind(1))-g%e_bot)/g%del)+1
            vol = kpts%voltet(itet)/kpts%ntet
            nend = g%ne
            DO ie = MAX(1,nstart), g%ne
               CALL contribSingletetra((ie-1)*g%del+g%e_bot,vol,e(ind(1:4)),icorn,l_bloechl,weight)
               weights(ie,ib) = weights(ie,ib) + weight * fac
               IF(weight.EQ.1/4.0*vol) THEN
                  !here we are above all eigenergies at the corners
                  !of the tetrahedron and we can simplify the rest of the loop
                  nend = ie
                  EXIT
               ENDIF
            ENDDO
            IF(nend.NE.g%ne) weights(nend+1:g%ne,ib) = weights(nend+1:g%ne,ib) + fac/4.0*vol
         ENDDO
         !$OMP END DO
         !$OMP END PARALLEL
      ENDDO
      !
      !Differentiate the weights with respect to energy to get the correct weights for a DOS-calculation
      !
      tol = 1E-14
      !$OMP PARALLEL DEFAULT(none) &
      !$OMP SHARED(weights,g,e_ind,tol,neig) &
      !$OMP PRIVATE(ib,i,dos_weights)

      !$OMP DO
      DO ib = 1, neig
         dos_weights = 0.0
         !---------------------------------------------------
         ! Weights for DOS -> differentiate with respect to E
         !---------------------------------------------------
         CALL diff3(weights(:,ib),g%del,dos_weights)
         !----------------------------------------------------
         !Find the range where the weights are bigger than tol
         !----------------------------------------------------
         ! Lower bound
         !--------------------
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
         !--------------------
         ! Upper bound
         !--------------------
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
            CALL juDFT_error("Negative tetra weight",calledby="dosWeightsCalc")
         ENDIF
         !---------------------------------------------------------------------------------------
         ! All weights are smaller than the tolerance -> set the indices to equal to discard this
         !---------------------------------------------------------------------------------------
         IF(e_ind(ib,1).GT.e_ind(ib,2)) THEN
            e_ind(ib,1) = 1
            e_ind(ib,2) = 1
         ENDIF
         weights(:,ib) = dos_weights(:)
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL
   END SUBROUTINE dosWeightsCalc

   SUBROUTINE dosWeightsCalcTria(ikpt,kpts,neig,eig,g,ef,weights,e_ind)

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

      INTEGER itria,ib,nstart,nend,ie,icorn,i
      REAL vol,weight,tol

      !Local Arrays
      REAL,ALLOCATABLE :: dos_weights(:)
      REAL             :: e(3)
      INTEGER          :: ind(3)

      ALLOCATE( dos_weights(g%ne) )
      dos_weights = 0.0
      e_ind(:,1) = 0
      e_ind(:,2) = g%ne+1
      weights = 0.0

      DO itria = 1, kpts%ntria

         IF(ALL(kpts%itria(:,itria).NE.ikpt)) CYCLE

         !$OMP PARALLEL DEFAULT(none) &
         !$OMP SHARED(ikpt,itria,neig) &
         !$OMP SHARED(kpts,eig,g,weights) &
         !$OMP PRIVATE(ib,ie,i,nstart,nend,icorn,weight,vol) &
         !$OMP PRIVATE(ind,e)

         !$OMP DO
         DO ib = 1, neig
            IF((MINVAL(eig(ib,kpts%itria(1:3,itria))).GT.g%e_top).OR.(MAXVAL(eig(ib,kpts%itria(1:3,itria))).LT.g%e_bot)) CYCLE

            e(1:3) = eig(ib,kpts%itria(1:3,itria))
            ind(1:3) = (/1,2,3/)
            CALL sortEig(3,e,ind)
            !Sort the energies in the tetrahedron in ascending order
            !search for the corner ikpt in the sorted array
            DO i = 1, 3
               IF(kpts%itria(ind(i),itria).EQ.ikpt) icorn = i
            ENDDO

            !Below this index the weight is 0
            nstart = INT((e(ind(1))-g%e_bot)/g%del)+1
            vol = kpts%voltria(itria)
            nend = g%ne

            DO ie = MAX(1,nstart), g%ne
               CALL contribSingleTria((ie-1)*g%del+g%e_bot,vol,e(ind(1:3)),icorn,weight)
               weights(ie,ib) = weights(ie,ib) + weight
               IF(weight.EQ.1.0/3.0*vol) THEN
                  !here we are above all eigenergies at the corners
                  !of the tetrahedron and we can simplify the rest of the loop
                  nend = ie
                  EXIT
               ENDIF
            ENDDO
            IF(nend.NE.g%ne) weights(nend+1:g%ne,ib) = weights(nend+1:g%ne,ib) + 1.0/3.0*vol
         ENDDO
         !$OMP END DO
         !$OMP END PARALLEL
      ENDDO

      tol = 1E-14
      !$OMP PARALLEL DEFAULT(none) &
      !$OMP SHARED(weights,g,e_ind,tol,neig) &
      !$OMP PRIVATE(ib,i,dos_weights)

      !$OMP DO
      DO ib = 1, neig
         dos_weights = 0.0
         !---------------------------------------------------
         ! Weights for DOS -> differentiate with respect to E
         !---------------------------------------------------
         CALL diff3(weights(:,ib),g%del,dos_weights(:))
         !----------------------------------------------------
         !Find the range where the weights are bigger than tol
         !----------------------------------------------------
         ! Lower bound
         !--------------------
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
         !--------------------
         ! Upper bound
         !--------------------
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
            CALL juDFT_error("Negative tria weight",calledby="dosWeightsCalcTria")
         ENDIF
         !---------------------------------------------------------------------------------------
         ! All weights are smaller than the tolerance -> set the indices to equal to discard this
         !---------------------------------------------------------------------------------------
         IF(e_ind(ib,1).GT.e_ind(ib,2)) THEN
            e_ind(ib,1) = 1
            e_ind(ib,2) = 1
         ENDIF
         weights(:,ib) = dos_weights(:)
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL

   END SUBROUTINE dosWeightsCalcTria

   SUBROUTINE contribSingleTria(energy,vol,e,ind,weight)

      USE m_juDFT
      !Integration weights (equivalent to doswt)

      IMPLICIT NONE

      REAL,                INTENT(IN)     :: energy
      REAL,                INTENT(IN)     :: vol
      REAL,                INTENT(IN)     :: e(:)
      INTEGER,             INTENT(IN)     :: ind
      REAL,                INTENT(OUT)    :: weight


      weight = 0.0
      IF(energy.GT.e(3)) THEN

         weight = 1.0

      ELSE IF( energy.GT.e(2) ) THEN

         IF(ind.EQ.1) THEN
            weight = (1.0-(e(3)-energy)**3/((e(3)-e(1))**2*(e(3)-e(2))))
         ELSE IF(ind.EQ.2) THEN
            weight = (1.0-(e(3)-energy)**3/((e(3)-e(1))*(e(3)-e(2))**2))
         ELSE IF(ind.EQ.3) THEN
            weight = (1.0-(e(3)-energy)**2/((e(3)-e(1))*(e(3)-e(2)))*&
                     (3.0-(e(3)-energy)/(e(3)-e(1))-(e(3)-energy)/(e(3)-e(2))))
         ENDIF

      ELSE IF(energy.GT.e(1)) THEN

         IF(ind.EQ.1) THEN
            weight = (energy-e(1))**2/((e(3)-e(1))*(e(2)-e(1)))*&
                     (3.0-(energy-e(1))/(e(3)-e(1))-(energy-e(1))/(e(2)-e(1)))
         ELSE IF(ind.EQ.2) THEN
            weight = (energy-e(1))**3/((e(3)-e(1))*(e(2)-e(1))**2)
         ELSE IF(ind.EQ.3) THEN
            weight = (energy-e(1))**3/((e(3)-e(1))**2*(e(2)-e(1)))
         ENDIF

      ENDIF

      weight = weight * vol/3.0

   END SUBROUTINE contribSingleTria



   SUBROUTINE contribSingleTetra(energy,vol,e,ind,l_bloechl,weight)

      USE m_juDFT
      !Integration weights taken from
      !PhysRevB.49.16223

      IMPLICIT NONE

      REAL,                INTENT(IN)     :: energy
      REAL,                INTENT(IN)     :: vol
      REAL,                INTENT(IN)     :: e(:)
      INTEGER,             INTENT(IN)     :: ind
      LOGICAL,             INTENT(IN)     :: l_bloechl
      REAL,                INTENT(OUT)    :: weight

      INTEGER j

      REAL C, D_ef
      REAL C1, C2, C3

      weight = 0.0
      IF( energy.GT.e(4) ) THEN

         weight = 1.0

      ELSE IF( energy.GT.e(3) ) THEN

         C     =  (e(4)-energy)**3/((e(4)-e(1))*(e(4)-e(2))*(e(4)-e(3)))

         IF(ind.EQ.4) THEN
            weight = 1.0 - C * ( 4 - (e(4)-energy) * (1/(e(4)-e(1)) + 1/(e(4)-e(2)) + 1/(e(4)-e(3))) )
         ELSE
            weight = 1.0 - C * (e(4)-energy)/(e(4)-e(ind))
         END IF

      ELSE IF( energy.GT.e(2) ) THEN

         C1     = (energy-e(1))**2/((e(3)-e(1))*(e(4)-e(1)))

         C2     = (energy-e(1))*(energy-e(2))*(e(3)-energy)/((e(3)-e(1))*(e(3)-e(2))*(e(4)-e(1)))

         C3     = (energy-e(2))**2*(e(4)-energy)/((e(3)-e(2))*(e(4)-e(1))*(e(4)-e(2)))

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

      ELSE IF( energy.GT.e(1) ) THEN

         C     = (energy-e(1))**3/((e(2)-e(1))*(e(3)-e(1))*(e(4)-e(1)))

         IF(ind.EQ.1) THEN
            weight = C * (4 - (energy-e(1)) * (1/(e(2)-e(1)) + 1/(e(3)-e(1)) + 1/(e(4)-e(1))))
         ELSE
            weight = C * (energy-e(1))/(e(ind)-e(1))
         ENDIF

      END IF

      IF(l_bloechl) THEN
         IF( energy.GT.e(4) ) THEN
            D_ef = 0.0
         ELSE IF( energy.LT.e(1) ) THEN
            D_ef = 0.0
         ELSE IF( energy.GT.e(3) ) THEN
            D_ef = 3.0 * (e(4)-energy)**2/( (e(4)-e(1))*(e(4)-e(2))*(e(4)-e(3)) )
         ELSE IF( energy.GT.e(2) ) THEN
            D_ef = 3.0 * 1./( (e(3)-e(1))*(e(4)-e(2)) ) *&
                       ( (energy-e(2))*(e(3)-energy)/(e(3)-e(2)) &
                        +(energy-e(1))*(e(4)-energy)/(e(4)-e(1)) )
         ELSE IF( energy.GT.e(1) ) THEN
            D_ef = 3.0 * (energy-e(1))**2/( (e(2)-e(1))*(e(3)-e(1))*(e(4)-e(1)) )
         ENDIF

         weight = weight + 1/10.0 * D_ef * SUM(e(:)-e(ind))
      ENDIF

      weight = weight * vol/4.0

   END SUBROUTINE contribSingleTetra


   ! Not used at the moment
      SUBROUTINE getTetraContrib(energy,vol,ev,w,g)

      !Integration weights taken from
      !PhysRevB.49.16223

      USE m_types

      IMPLICIT NONE

      REAL,                   INTENT(IN)     :: energy
      REAL,                   INTENT(IN)     :: vol
      REAL,                   INTENT(IN)     :: ev(:)
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

END MODULE m_dosWeights
