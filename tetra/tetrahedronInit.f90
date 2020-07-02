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

      INTEGER :: ikpt,ncorn,itet,icorn,iband,k(SIZE(kpts%ntetra,1))
      REAL    :: w(1),etetra(SIZE(kpts%ntetra,1)),fac,vol
      logical :: l_weights_pres, l_weightsum_pres

      IF(.NOT.PRESENT(weightSum).AND..NOT.PRESENT(weights)) THEN
         CALL juDFT_error("No output variable provided (either weightSum or weights)",&
                           calledby="getWeightKpoints")
      ENDIF

      l_weights_pres = PRESENT(weights)
      l_weightsum_pres = PRESENT(weightSum)
      !Tetrahedra or Triangles?
      ncorn = MERGE(3,4,film)

      IF(PRESENT(weights)) weights = 0.0
      IF(PRESENT(weightSum)) weightSum = 0.0

      !More efficient to just loop through all tetrahedra
      DO itet = 1, kpts%ntet
         k = kpts%ntetra(:,itet)
         where(k.GT.kpts%nkpt) k = kpts%bkp(k)

         DO icorn = 1, ncorn
            ikpt = kpts%ntetra(icorn,itet)
            IF(ikpt.GT.kpts%nkpt) CYCLE
            fac = REAL(COUNT(kpts%bkp(:).EQ.ikpt))
            vol = kpts%voltet(itet)/kpts%ntet*fac
            !$OMP parallel do default(none) &
            !$OMP shared(itet,neig,ikpt,icorn,film,ncorn,k,vol,l_weights_pres,l_weightsum_pres) &
            !$OMP shared(eig,weights,efermi) &
            !$OMP private(iband,etetra,w) &
            !$OMP reduction(+:weightSum) schedule(dynamic,1)
            DO iband = 1, neig

               etetra = eig(iband,k)

               IF( ALL(etetra>efermi) ) CYCLE
               IF( ALL(etetra<efermi) ) THEN
                  w = vol/REAL(ncorn)
               ELSE
                  w  = getWeightSingleBand([efermi],etetra,icorn,vol,film,.FALSE.)
               ENDIF

               IF(l_weights_pres)  weights(iband,ikpt) = weights(iband,ikpt) + w(1)
               IF(l_weightsum_pres) weightSum = weightSum + w(1)
            ENDDO
            !$OMP end parallel do
         ENDDO
      ENDDO

   END SUBROUTINE getWeightKpoints

   SUBROUTINE getWeightEnergyMesh(kpts,ikpt,eig,neig,eMesh,ne,film,weights,resWeights,bounds,dos)

      USE m_differentiate

      TYPE(t_kpts),     INTENT(IN)    :: kpts
      REAL,             INTENT(IN)    :: eig(:,:)
      REAL,             INTENT(INOUT) :: weights(:,:)
      INTEGER,OPTIONAL, INTENT(INOUT) :: bounds(:,:)
      REAL, OPTIONAL,   INTENT(INOUT) :: resWeights(:,:)

      INTEGER,          INTENT(IN)  :: neig
      INTEGER,          INTENT(IN)  :: ikpt,ne
      REAL,             INTENT(IN)  :: eMesh(:)
      LOGICAL,          INTENT(IN)  :: film
      LOGICAL,OPTIONAL, INTENT(IN)  :: dos

      INTEGER :: itet,iband,ncorn,ie,icorn,i,k(SIZE(kpts%ntetra,1))
      LOGICAL :: l_dos, l_bounds_pres, l_resWeights_pres
      REAL    :: etetra(SIZE(kpts%ntetra,1)),del,fac,vol
      REAL, ALLOCATABLE :: dos_weights(:)
      !Temporary Arrays to include end points
      !to avoid numerical trouble with differentiation
      REAL, ALLOCATABLE :: calc_weights(:,:)
      REAL, ALLOCATABLE :: calc_eMesh(:)

      fac = REAL(COUNT(kpts%bkp(:).EQ.ikpt))

      l_dos = PRESENT(dos)
      IF(PRESENT(dos))THEN
         l_dos = dos.AND.ne>1
      ENDIF

      IF(PRESENT(resWeights)) resWeights = 0.0

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

      l_resWeights_pres = PRESENT(resWeights)
      DO itet = 1, kpts%ntet
         IF(ALL(kpts%ntetra(:,itet).NE.ikpt)) CYCLE
         vol = kpts%voltet(itet)/kpts%ntet*fac
         !The array k is only for getting the right indices in the eigenvalues
         k = kpts%ntetra(:,itet)
         DO i = 1, SIZE(k)
            IF(i==ikpt) icorn = i
         ENDDO
         where(k.GT.kpts%nkpt) k = kpts%bkp(k)
         !$OMP parallel do default(none) schedule(dynamic,1) &
         !$OMP shared(icorn,neig,film,vol,k,l_resWeights_pres) &
         !$OMP shared(eig,calc_weights,calc_eMesh,eMesh,resWeights) &
         !$OMP private(iband,etetra)
         DO iband = 1, neig

            etetra = eig(iband,k)

            IF(l_resWeights_pres) THEN
               resWeights(:,iband) = resWeights(:,iband) + getWeightSingleBand(eMesh,etetra,icorn,&
                                                                               vol,film,.FALSE.,l_res=.TRUE.)
            ENDIF

            IF( ALL(etetra>MAXVAL(calc_eMesh)) ) CYCLE

            calc_weights(:,iband) = calc_weights(:,iband) + getWeightSingleBand(calc_eMesh,etetra,icorn,&
                                                                                vol,film,.FALSE.)

         ENDDO
         !$OMP end parallel do
      ENDDO

      weights = 0.0
      l_bounds_pres = PRESENT(bounds)
      !-------------------------------------
      ! PostProcess weights
      !-------------------------------------
      !$OMP parallel default(none) &
      !$OMP shared(neig,l_dos,ne,del,eMesh) &
      !$OMP shared(calc_weights,weights,bounds, l_bounds_pres,l_resWeights_pres) &
      !$OMP private(iband,dos_weights,ie)
      IF(l_dos) ALLOCATE(dos_weights(ne+2),source=0.0)
      !$OMP do schedule(dynamic,1)
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
         IF(l_bounds_pres.AND..NOT.l_resWeights_pres) THEN
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
         ELSE IF(l_bounds_pres) THEN
            bounds(iband,1) = 1
            bounds(iband,2) = SIZE(eMesh)
         ENDIF
      ENDDO
      !$OMP end do
      !$OMP end parallel

      IF(ANY(weights(:,:neig)<0.0)) THEN
         CALL juDFT_error("TetraWeight error: Unexpected negative weight", calledby="getWeightEnergyMesh")
      ENDIF

   END SUBROUTINE getWeightEnergyMesh


   PURE FUNCTION getWeightSingleBand(eMesh,etetra,icorn,vol,film,bloechl,l_res) Result(weight)

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
      INTEGER,          INTENT(IN)     :: icorn       !Current k-point index (We calculate the weights for this corner)
      REAL,             INTENT(IN)     :: vol         !Volume of the tetrahedron
      LOGICAL,          INTENT(IN)     :: film        !Switch controls wether tetrahedron/triangular method is used
      LOGICAL,          INTENT(IN)     :: bloechl     !Controls bloechl corrections (not atm)
      LOGICAL, OPTIONAL,INTENT(IN)     :: l_res

      REAL    :: weight(SIZE(eMesh))
      INTEGER :: i,ie,nstart,nend,ind(SIZE(etetra))
      REAL    :: eb,et,del
      LOGICAL :: l_calcres


      !Sort the eigenvalues at the corners (ascending order in ind)
      ind = tetsrt(etetra)

      l_calcres = .FALSE.
      IF(PRESENT(l_res)) l_calcres = l_res

      !Find the range in the (equidistant) energy mesh where the weights are changing
      IF( SIZE(eMesh)>1 .AND..NOT.l_calcres) THEN
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

      IF(nstart /= 1) weight(:nstart-1) = 0.0
      !Calculate the weights
      DO ie = nstart, nend
         weight(ie) = tetraWeight(eMesh(ie),etetra(ind),ind(icorn),vol,film)
      ENDDO

      !The loop terminates if the energy is larger than
      !all eigenvalues at the tetrahedron/triangle corners (nend)
      !For all consecutive values the weight is constant
      IF(nend.NE.SIZE(eMesh)) weight(nend+1:) = weight(nend)

   END FUNCTION getWeightSingleBand

END MODULE m_tetrahedronInit