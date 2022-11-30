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

   USE m_types_kpts
   USE m_types_input
   USE m_juDFT

   IMPLICIT NONE

   PRIVATE
   PUBLIC   :: tetrahedronInit

   INTERFACE tetrahedronInit
      PROCEDURE getWeightKpoints, getWeightEnergyMesh
   END INTERFACE

   REAL, PARAMETER :: weightCutoff=1e-14
   real, parameter :: tol_degeneracy=1e-9

   CONTAINS

   SUBROUTINE getWeightKpoints(kpts,input,eig,neig,efermi,weightSum,weights)

      TYPE(t_kpts),  INTENT(IN)    :: kpts
      TYPE(t_input), INTENT(IN)    :: input
      REAL,          INTENT(IN)    :: eig(:,:)
      INTEGER,       INTENT(IN)    :: neig
      REAL,          INTENT(IN)    :: efermi
      REAL, OPTIONAL,INTENT(INOUT) :: weightSum
      REAL, OPTIONAL,INTENT(INOUT) :: weights(:,:)

      INTEGER :: ikpt,ncorn,itet,icorn,iband,num_degenerate_states,last_band
      REAL    :: w(1),etetra(SIZE(kpts%ntetra,1)),vol(kpts%ntet)
      logical :: l_weights_pres, l_weightsum_pres

      IF(.NOT.PRESENT(weightSum).AND..NOT.PRESENT(weights)) THEN
         CALL juDFT_error("No output variable provided (either weightSum or weights)",&
                           calledby="getWeightKpoints")
      ENDIF

      l_weights_pres = PRESENT(weights)
      l_weightsum_pres = PRESENT(weightSum)
      !Tetrahedra or Triangles?
      ncorn = MERGE(3,4,input%film)

      IF(PRESENT(weights)) weights = 0.0
      IF(PRESENT(weightSum)) weightSum = 0.0

      vol = kpts%voltet(:)/kpts%ntet

      !More efficient to just loop through all tetrahedra
      !$OMP parallel do default(none) &
      !$OMP shared(neig,ncorn,vol,l_weights_pres,l_weightsum_pres) &
      !$OMP shared(kpts,input,eig,weights,efermi) &
      !$OMP private(itet,ikpt,icorn,iband,etetra,w) &
      !$OMP reduction(+:weightSum) collapse(3)
      DO itet = 1, kpts%ntet
         DO icorn = 1, ncorn
            DO iband = 1, neig

               ikpt = kpts%ntetra(icorn,itet)
               etetra = eig(iband,kpts%ntetra(:,itet))

               IF( ALL(etetra>efermi + 1E-8) .AND. .NOT.input%l_bloechl ) CYCLE

               w  = getWeightSingleBand([efermi],etetra,icorn,vol(itet),input%film,input%l_bloechl)

               IF(l_weightsum_pres) weightSum = weightSum + w(1)

               IF(l_weights_pres) THEN
                  !$OMP critical
                  weights(iband,ikpt) = weights(iband,ikpt) + w(1)
                  !$OMP end critical
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      !$OMP end parallel do

      if (l_weights_pres) then
         !find degenerate states
         do ikpt=1,kpts%nkpt
            iband = 1
            do while(iband<neig)
               num_degenerate_states=1
               do while (abs(eig(iband,ikpt)-eig(iband+num_degenerate_states,ikpt))<tol_degeneracy)
                  num_degenerate_states=num_degenerate_states+1
                  if (iband+num_degenerate_states>neig) exit
               enddo

               if (num_degenerate_states>1) THEN
                  last_band=iband+num_degenerate_states-1
                  !Make sure all weights are equal
                  weights(iband:last_band,ikpt)=sum(weights(iband:last_band,ikpt))/num_degenerate_states
                  iband=iband+num_degenerate_states
               endif
               iband=iband+1
            enddo
         enddo
      endif

   END SUBROUTINE getWeightKpoints

   SUBROUTINE getWeightEnergyMesh(kpts,input,ikpt,eig,neig,eMesh,weights,resWeights,bounds,dos)

      USE m_differentiate

      TYPE(t_kpts),     INTENT(IN)    :: kpts
      TYPE(t_input),    INTENT(IN)    :: input
      REAL,             INTENT(IN)    :: eig(:,:)
      REAL,             INTENT(INOUT) :: weights(:,:)
      INTEGER,OPTIONAL, INTENT(INOUT) :: bounds(:,:)
      REAL, OPTIONAL,   INTENT(INOUT) :: resWeights(:,:)

      INTEGER,          INTENT(IN)  :: neig
      INTEGER,          INTENT(IN)  :: ikpt
      REAL,             INTENT(IN)  :: eMesh(:)
      LOGICAL,OPTIONAL, INTENT(IN)  :: dos

      INTEGER :: itet,iband,ncorn,ie,icorn,ntet,num_degenerate_states, last_band
      LOGICAL :: l_dos, l_bounds_pres, l_resWeights_pres
      REAL    :: etetra(SIZE(kpts%ntetra,1),neig,SIZE(kpts%tetraList,1)),del,vol
      REAL, ALLOCATABLE :: dos_weights(:)
      !Temporary Arrays to include end points
      !to avoid numerical trouble with differentiation
      REAL, ALLOCATABLE :: calc_weights(:,:)
      REAL, ALLOCATABLE :: calc_weights_thread(:,:),resWeights_thread(:,:)
      REAL, ALLOCATABLE :: calc_eMesh(:)

      !Extract the right eigenvalues
      IF(kpts%tetraList(1,ikpt)==0) CALL juDFT_error("No tetrahedrons in tetraList", calledby="tetrahedronInit")
      ntet = 1
      do
         itet = kpts%tetraList(ntet,ikpt)

         do iband = 1, neig
            etetra(:,iband,ntet) = eig(iband,kpts%ntetra(:,itet))
         enddo
         if(ntet < SIZE(kpts%tetraList,1)) THEN
            if(kpts%tetraList(ntet+1,ikpt)>0) THEN
               ntet = ntet + 1
            else
               exit
            endif
         else
            exit
         endif
      enddo

      l_dos = PRESENT(dos)
      IF(PRESENT(dos))THEN
         l_dos = dos.AND.SIZE(eMesh)>1
      ENDIF

      IF(PRESENT(resWeights)) resWeights = 0.0

      IF(l_dos) THEN
         ALLOCATE(calc_weights(SIZE(eMesh)+2,neig),source=0.0)
         ALLOCATE(calc_eMesh(SIZE(eMesh)+2),source=0.0)
         del = eMesh(2)-eMesh(1)
         calc_eMesh = [eMesh(1)-del,eMesh(:),eMesh(SIZE(eMesh))+del]
      ELSE
         ALLOCATE(calc_weights(SIZE(eMesh),neig),source=0.0)
         ALLOCATE(calc_eMesh(SIZE(eMesh)),source=0.0)
         calc_eMesh = eMesh
      ENDIF

      l_resWeights_pres = PRESENT(resWeights)
      !$OMP parallel default(none)&
      !$OMP shared(neig,ntet,l_resWeights_pres,kpts,input,ikpt,resWeights,calc_weights) &
      !$OMP shared(etetra,calc_eMesh,eMesh) &
      !$OMP private(itet,iband,vol,icorn,resWeights_thread,calc_weights_thread)
      IF(l_resWeights_pres) ALLOCATE(resWeights_thread(SIZE(resWeights,1),SIZE(resWeights,2)),source = 0.0)
      ALLOCATE(calc_weights_thread(SIZE(calc_weights,1),SIZE(calc_weights,2)),source = 0.0)
      !$OMP do collapse(3)
      DO itet = 1, ntet
         DO iband = 1, neig
            DO icorn = 1, SIZE(kpts%ntetra,1)

               vol = kpts%voltet(kpts%tetraList(itet,ikpt))/kpts%ntet
               IF(kpts%ntetra(icorn,kpts%tetraList(itet,ikpt)).NE.ikpt) CYCLE

               IF(l_resWeights_pres) THEN
                  resWeights_thread(:,iband) = resWeights_thread(:,iband) &
                                              + getWeightSingleBand(eMesh,etetra(:,iband,itet),icorn,vol,&
                                                                    input%film,input%l_bloechl,l_res=.TRUE.)
               ENDIF

               IF( ALL(etetra(:,iband,itet)>MAXVAL(calc_eMesh) + 1E-8) .AND. .NOT.input%l_bloechl ) CYCLE

               calc_weights_thread(:,iband) = calc_weights_thread(:,iband) &
                                             + getWeightSingleBand(calc_eMesh,etetra(:,iband,itet),icorn,vol,&
                                                                   input%film,input%l_bloechl)
            ENDDO
         ENDDO
      ENDDO
      !$OMP end do
      !$OMP critical
      IF(l_resWeights_pres) resWeights = resWeights + resWeights_thread
      calc_weights = calc_weights + calc_weights_thread
      !$OMP end critical
      DEALLOCATE(calc_weights_thread)
      IF(l_resWeights_pres) DEALLOCATE(resWeights_thread)
      !$OMP end parallel

      !find degenerate states
      iband = 1
      do while(iband<neig)
         num_degenerate_states=1
         do while (abs(eig(iband,ikpt)-eig(iband+num_degenerate_states,ikpt))<tol_degeneracy)
            num_degenerate_states=num_degenerate_states+1
            if (iband+num_degenerate_states>neig) exit
         enddo

         if (num_degenerate_states>1) THEN
            last_band=iband+num_degenerate_states-1
            !Make sure all weights are equal
            calc_weights(:,iband:last_band)=sum(calc_weights(:,iband:last_band))/num_degenerate_states
            iband=iband+num_degenerate_states
         endif
         iband=iband+1
      enddo

      weights = 0.0
      l_bounds_pres = PRESENT(bounds)
      !-------------------------------------
      ! PostProcess weights
      !-------------------------------------
      !$OMP parallel default(none) &
      !$OMP shared(neig,l_dos,del,eMesh) &
      !$OMP shared(calc_weights,weights,bounds,l_bounds_pres,l_resWeights_pres) &
      !$OMP private(iband,dos_weights,ie)
      IF(l_dos) ALLOCATE(dos_weights(SIZE(eMesh)+2),source=0.0)
      !$OMP do schedule(dynamic,1)
      DO iband = 1, neig
         !---------------------------------------------------
         ! Weights for DOS -> differentiate with respect to E
         !---------------------------------------------------
         IF(l_dos) THEN
            CALL diff3(calc_weights(:,iband),del,dos_weights)
            weights(1:SIZE(eMesh),iband) = dos_weights(2:SIZE(eMesh)+1)
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
                  IF(ie.EQ.SIZE(eMesh)) THEN
                     bounds(iband,1) = SIZE(eMesh)
                     EXIT
                  ENDIF
               ENDIF
            ENDDO
            !--------------------
            ! Upper bound
            !--------------------
            ie = SIZE(eMesh)
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
      IF(l_dos) DEALLOCATE(dos_weights)
      !$OMP end parallel

      IF(ANY(weights(:,:neig)<0.0).AND..NOT.input%l_bloechl) THEN
         CALL juDFT_error("TetraWeight error: Unexpected negative weight", calledby="getWeightEnergyMesh")
      ENDIF

   END SUBROUTINE getWeightEnergyMesh


   PURE FUNCTION getWeightSingleBand(eMesh,etetra,icorn,vol,film,l_bloechl,l_res) Result(weight)

      !--------------------------------------------------------------
      ! This is the core routine calculating the integration
      ! weight for one kpoint in one tetrahedron for one band index
      ! for the provided energy points
      ! All other routines wrap this procedure so that it is most
      ! efficient in all cases
      !--------------------------------------------------------------

      USE m_tetsrt
      USE m_tetraWeight
      USE m_resWeight
      USE m_bloechl

      REAL,             INTENT(IN)     :: eMesh(:)    !Energy points, where the weights are calculated
      REAL,             INTENT(IN)     :: etetra(:)   !Eigenvalues at the corners of the tetrahedron
      INTEGER,          INTENT(IN)     :: icorn       !Current k-point index (We calculate the weights for this corner)
      REAL,             INTENT(IN)     :: vol         !Volume of the tetrahedron
      LOGICAL,          INTENT(IN)     :: film        !Switch controls wether tetrahedron/triangular method is used
      LOGICAL,          INTENT(IN)     :: l_bloechl   !Controls bloechl corrections
      LOGICAL, OPTIONAL,INTENT(IN)     :: l_res

      REAL    :: weight(SIZE(eMesh))
      INTEGER :: i,ie,nstart,nend,ind(SIZE(etetra)),icornSorted
      REAL    :: eb,et,del
      LOGICAL :: l_calcres


      !Sort the eigenvalues at the corners (ascending order in ind)
      ind = tetsrt(etetra)

      !Find out where the corner went
      DO i = 1, SIZE(ind)
         IF(ind(i)==icorn) icornSorted = i
      ENDDO

      l_calcres = .FALSE.
      IF(PRESENT(l_res)) l_calcres = l_res

      IF(l_calcres) THEN
         weight = resWeight(eMesh,etetra(ind),icornSorted,vol,film)
      ELSE
         !Find the range in the (equidistant) energy mesh where the weights are changing
         IF( SIZE(eMesh)>1 .AND. .NOT.l_bloechl) THEN !With bloechl corrections we cannot cut as clean
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
            weight(ie) = tetraWeight(eMesh(ie),etetra(ind),icornSorted,vol,film)
            IF(l_bloechl) weight(ie) = weight(ie) + bloechl(eMesh(ie),etetra(ind),icornSorted,vol,film)
         ENDDO

         !The loop terminates if the energy is larger than
         !all eigenvalues at the tetrahedron/triangle corners (nend)
         !For all consecutive values the weight is constant
         IF(nend /= SIZE(eMesh)) weight(nend+1:) = weight(nend)
      ENDIF

   END FUNCTION getWeightSingleBand

END MODULE m_tetrahedronInit
