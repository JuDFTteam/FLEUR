MODULE m_tetrahedronInit

   !------------------------------------------------------------------------------------
   !This module provides the weights for the linear tetrahedron/triangular method
   !for Brillouin-zone integration with the step function as the weight
   !
   !This module should be used by calling the interface tetrahedronInit
   !
   !Structure:
   !     tetrahedronInit:
   !        getWeightKpoints:    gives weight at one (fermi) energy for all kpoints
   !        getWeightEnergyMesh: gives weight on an energy grid for one kpoint.
   !                             By providing the dos switch as .TRUE. the weights are
   !                             numerically differentiated for use for DOS calculations
   !
   !These subroutines all call getWeightSingleBand which handles the actual calculation
   !of weights. These differ in the order of the loops and the arguments provided.
   !------------------------------------------------------------------------------------

   USE m_types
   USE m_juDFT

   IMPLICIT NONE

   PRIVATE
   PUBLIC   :: tetrahedronInit

   INTERFACE tetrahedronInit
      PROCEDURE getWeightKpoints, getWeightEnergyMesh
   END INTERFACE

   REAL, PARAMETER :: weightCutoff=1e-14

   CONTAINS

   SUBROUTINE getWeightKpoints(kpts,eig,neig,efermi,film,weightSum,weights)

      TYPE(t_kpts),  INTENT(IN)    :: kpts
      REAL,          INTENT(IN)    :: eig(:,:)
      INTEGER,       INTENT(IN)    :: neig
      REAL,          INTENT(IN)    :: efermi
      LOGICAL,       INTENT(IN)    :: film
      REAL, OPTIONAL,INTENT(INOUT) :: weightSum
      REAL, OPTIONAL,INTENT(INOUT) :: weights(:,:)

      INTEGER :: ikpt,ncorn,itet,icorn,iband,k(4)
      REAL    :: w(1),etetra(4),fac,vol
      REAL    :: weightSum_Band

      IF(.NOT.PRESENT(weightSum).AND..NOT.PRESENT(weights)) THEN
         CALL juDFT_error("No output variable provided (either weightSum or weights)",&
                           calledby="getWeightKpoints")
      ENDIF


      !Tetrahedra or Triangles?
      ncorn = MERGE(3,4,film)

      IF(PRESENT(weights)) weights = 0.0
      IF(PRESENT(weightSum)) weightSum = 0.0

      !More efficient to just loop through all tetrahedra
      DO itet = 1, kpts%ntet
         IF(kpts%nkptf.NE.0) THEN
            DO icorn = 1, ncorn
               IF(kpts%ntetra(icorn,itet).GT.kpts%nkpt) THEN
                  k(icorn) = kpts%bkp(kpts%ntetra(icorn,itet))
               ELSE
                  k(icorn) = kpts%ntetra(icorn,itet)
               ENDIF
            ENDDO
         ELSE
            k(:ncorn) = kpts%ntetra(:ncorn,itet)
         ENDIF
         DO icorn = 1, ncorn
            ikpt = kpts%ntetra(icorn,itet)
            IF(ikpt.GT.kpts%nkpt) CYCLE
            fac = REAL(MERGE(1,COUNT(kpts%bkp(:).EQ.ikpt),kpts%nkptf.EQ.0))
            vol = kpts%voltet(itet)/kpts%ntet*fac
            !$OMP PARALLEL DEFAULT(none) &
            !$OMP SHARED(itet,neig,ikpt,film,ncorn,k,vol) &
            !$OMP SHARED(kpts,eig,weights,efermi,weightSum) &
            !$OMP PRIVATE(iband,etetra,w,weightSum_Band)
            weightSum_Band = 0.0
            !$OMP DO SCHEDULE(DYNAMIC,1)
            DO iband = 1, neig

               etetra(:ncorn) = eig(iband,k(:ncorn))

               IF( ALL(etetra(:ncorn)>efermi) ) CYCLE

               w = 0.0
               CALL getWeightSingleBand([efermi],etetra(:ncorn),ikpt,kpts%ntetra(:,itet),&
                                        vol,film,.FALSE.,w)

               IF(PRESENT(weights)) weights(iband,ikpt) = weights(iband,ikpt) + w(1)
               weightSum_Band = weightSum_Band + w(1)
            ENDDO
            !$OMP END DO
            IF(PRESENT(weightSum)) THEN
               !$OMP CRITICAL
               weightSum = weightSum + weightSum_Band
               !$OMP END CRITICAL
            ENDIF
            !$OMP END PARALLEL
         ENDDO
      ENDDO

   END SUBROUTINE getWeightKpoints

   SUBROUTINE getWeightEnergyMesh(kpts,ikpt,eig,neig,eMesh,ne,film,weights,bounds,dos)

      USE m_differentiate

      TYPE(t_kpts),     INTENT(IN)    :: kpts
      REAL,             INTENT(IN)    :: eig(:,:)
      REAL,             INTENT(INOUT) :: weights(:,:)
      INTEGER,OPTIONAL, INTENT(INOUT) :: bounds(:,:)

      INTEGER,          INTENT(IN)  :: neig
      INTEGER,          INTENT(IN)  :: ikpt,ne
      REAL,             INTENT(IN)  :: eMesh(:)
      LOGICAL,          INTENT(IN)  :: film
      LOGICAL,OPTIONAL, INTENT(IN)  :: dos

      INTEGER :: itet,iband,ncorn,ie,icorn,k(4)
      LOGICAL :: l_dos
      REAL    :: etetra(4),del,fac,vol
      REAL, ALLOCATABLE :: dos_weights(:)
      !Temporary Arrays to include end points
      !to avoid numerical trouble with differentiation
      REAL, ALLOCATABLE :: calc_weights(:,:)
      REAL, ALLOCATABLE :: calc_eMesh(:)

      !Tetrahedra or Triangles?
      ncorn = MERGE(3,4,film)
      fac = REAL(MERGE(1,COUNT(kpts%bkp(:).EQ.ikpt),kpts%nkptf.EQ.0))

      l_dos = PRESENT(dos)
      IF(PRESENT(dos))THEN
         l_dos = dos.AND.ne>1
      ENDIF

      IF(l_dos) THEN
         ALLOCATE(calc_weights(ne+2,neig),source=0.0)
         ALLOCATE(calc_eMesh(ne+2),source=0.0)
         del = eMesh(2)-eMesh(1)
         calc_eMesh = [eMesh(1)-del,eMesh(:),eMesh(ne)+del]
      ELSE
         ALLOCATE(calc_weights(ne,neig),source=0.0)
         ALLOCATE(calc_eMesh(ne),source=0.0)
         calc_eMesh = eMesh
      ENDIF


      weights = 0.0
      DO itet = 1, kpts%ntet
         IF(ALL(kpts%ntetra(:ncorn,itet).NE.ikpt)) CYCLE
         vol = kpts%voltet(itet)/kpts%ntet*fac
         !The array k is only for getting the right indices in the eigenvalues
         IF(kpts%nkptf.NE.0) THEN
            DO icorn = 1, ncorn
               IF(kpts%ntetra(icorn,itet).GT.kpts%nkpt) THEN
                  k(icorn) = kpts%bkp(kpts%ntetra(icorn,itet))
               ELSE
                  k(icorn) = kpts%ntetra(icorn,itet)
               ENDIF
            ENDDO
         ELSE
            k(:ncorn) = kpts%ntetra(:ncorn,itet)
         ENDIF
         !$OMP PARALLEL DEFAULT(none) &
         !$OMP SHARED(itet,neig,ikpt,film,ncorn,vol,k) &
         !$OMP SHARED(kpts,eig,calc_weights,calc_eMesh) &
         !$OMP PRIVATE(iband,etetra)
         !$OMP DO SCHEDULE(DYNAMIC,1)
         DO iband = 1, neig

            etetra(:ncorn) = eig(iband,k(:ncorn))
            IF( ALL(etetra(:ncorn)>MAXVAL(calc_eMesh)) ) CYCLE

            CALL getWeightSingleBand(calc_eMesh,etetra(:ncorn),ikpt,kpts%ntetra(:,itet),&
                                     vol,film,.FALSE.,calc_weights(:,iband))

         ENDDO
         !$OMP END DO
         !$OMP END PARALLEL
      ENDDO


      !-------------------------------------
      ! PostProcess weights
      !-------------------------------------
      !$OMP PARALLEL DEFAULT(none) &
      !$OMP SHARED(neig,l_dos,ne,del) &
      !$OMP SHARED(calc_weights,weights,bounds) &
      !$OMP PRIVATE(iband,dos_weights,ie)
      IF(l_dos) ALLOCATE(dos_weights(ne+2),source=0.0)
      !$OMP DO SCHEDULE(DYNAMIC,1)
      DO iband = 1, neig
         !---------------------------------------------------
         ! Weights for DOS -> differentiate with respect to E
         !---------------------------------------------------
         IF(l_dos) THEN
            CALL diff3(calc_weights(:,iband),del,dos_weights)
            weights(1:ne,iband) = dos_weights(2:ne+1)
         ELSE
            weights(:,:neig) = calc_weights
         ENDIF
         !--------------------------------------------------------------
         ! Find the range where the weights are bigger than weightCutoff
         !--------------------------------------------------------------
         IF(PRESENT(bounds)) THEN
            !--------------------
            ! Lower bound
            !--------------------
            ie = 1
            DO
               IF(ABS(weights(ie,iband)).GT.weightCutoff) THEN
                  bounds(iband,1) = ie
                  EXIT
               ELSE
                  weights(ie,iband) = 0.0
                  ie = ie + 1
                  IF(ie.EQ.ne) THEN
                     bounds(iband,1) = ne
                     EXIT
                  ENDIF
               ENDIF
            ENDDO
            !--------------------
            ! Upper bound
            !--------------------
            ie = ne
            DO
               IF(ABS(weights(ie,iband)).GT.weightCutoff) THEN
                  bounds(iband,2) = ie
                  EXIT
               ELSE
                  weights(ie,iband) = 0.0
                  ie = ie - 1
                  IF(ie.EQ.0) THEN
                     bounds(iband,2) = 1
                     EXIT
                  ENDIF
               ENDIF
            ENDDO
            !For safety
            IF(bounds(iband,1).GT.bounds(iband,2)) THEN
               bounds(iband,1) = 1
               bounds(iband,2) = 1
            ENDIF
            IF(ANY(weights(bounds(iband,1):bounds(iband,2),iband)<0.0)) THEN
               CALL juDFT_error("TetraWeight error: Unexpected negative weight", calledby="getWeightEnergyMesh")
            ENDIF
         ENDIF
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL


   END SUBROUTINE getWeightEnergyMesh

   SUBROUTINE getWeightSingleBand(eMesh,etetra,ikpt,ntetra,vol,film,bloechl,weight)

      !--------------------------------------------------------------
      ! This is the core routine calculating the integration
      ! weight for one kpoint in one tetrahedron for one band index
      ! for the provided energy points
      ! All other routines wrap this procedure so that it is most
      ! efficient in all cases
      !--------------------------------------------------------------

      USE m_tetsrt
      USE m_tetraWeight

      REAL,             INTENT(IN)     :: eMesh(:)    !Energy points, where the weights are calculated
      REAL,             INTENT(IN)     :: etetra(:)   !Eigenvalues at the corners of the tetrahedron
      INTEGER,          INTENT(IN)     :: ntetra(:)   !k-point indices at the corners of the tetrahedron
      REAL,             INTENT(INOUT)  :: weight(:)   !Output array for the weights

      INTEGER,          INTENT(IN)     :: ikpt        !Current k-point index (We calculate the weights for this corner)
      REAL,             INTENT(IN)     :: vol         !Volume of the tetrahedron
      LOGICAL,          INTENT(IN)     :: film        !Switch controls wether tetrahedron/triangular method is used
      LOGICAL,          INTENT(IN)     :: bloechl     !Controls bloechl corrections (not atm)

      INTEGER :: icorn,i,ie,nstart,nend,ind(SIZE(etetra))
      REAL    :: w,eb,et,del,lastWeight

      !Sort the eigenvalues at the corners (ascending order in ind)
      CALL tetsrt(SIZE(etetra),etetra,ind)

      !search for the corner ikpt in the sorted array
      DO icorn = 1, SIZE(etetra)
         IF(ntetra(ind(icorn)).EQ.ikpt) EXIT
      ENDDO

      IF(ALL(ntetra.NE.ikpt)) THEN
         CALL juDFT_error("kpoint not found in ntetra"&
                          ,hint="This is a bug in FLEUR, please report"&
                          ,calledby="getWeightSingleBand")
      ENDIF

      !Find the range in the (equidistant) energy mesh where the weights are changing
      IF( SIZE(eMesh)>1 ) THEN
         !Extract basic parameters of the equidistant eMesh
         eb = MINVAL(eMesh)
         et = MAXVAL(eMesh)
         del = eMesh(2)-eMesh(1)

         !Get last point to the left of the lowest eigenvalue in the tetrahedron
         nstart = INT((etetra(ind(1))-eb)/del)+1
         nstart = MAX(1,nstart)

         !Get first point to the right of the highest eigenvalue in the tetrahedron
         nend   = INT((etetra(ind(SIZE(etetra)))-eb)/del)+2
         nend   = MIN(SIZE(eMesh),nend)
         nend   = MAX(1,nend)
      ELSE
         nstart = 1
         nend   = SIZE(eMesh)
      ENDIF

      !Calculate the weights
      DO ie = nstart, nend
         CALL tetraWeight(eMesh(ie),etetra(ind),icorn,vol,film,w)
         lastWeight = w
         weight(ie) = weight(ie) + w
      ENDDO

      !The loop terminates if the energy is larger than
      !all eigenvalues at the tetrahedron/triangle corners (nend)
      !For all consecutive values the weight is constant and stored in lastWeight
      IF(nend.NE.SIZE(eMesh)) weight(nend+1:) = weight(nend+1:) + lastWeight

   END SUBROUTINE getWeightSingleBand

END MODULE m_tetrahedronInit