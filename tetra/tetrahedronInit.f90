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

   IMPLICIT NONE

   PUBLIC   :: tetrahedronInit
   PRIVATE  :: getWeightKpoints, getWeightEnergyMesh, getWeightSingleBand

   INTERFACE tetrahedronInit
      PROCEDURE getWeightKpoints, getWeightEnergyMesh
   END INTERFACE

   REAL, PARAMETER :: weightCutoff=1e-14

   CONTAINS

   SUBROUTINE getWeightKpoints(kpts,eig,neig,efermi,film,weights)

      USE m_types
      USE m_juDFT

      IMPLICIT NONE

      TYPE(t_kpts),  INTENT(IN)  :: kpts
      REAL,          INTENT(IN)  :: eig(:,:)
      REAL,          INTENT(OUT) :: weights(:,:)

      INTEGER,       INTENT(IN)  :: neig
      REAL,          INTENT(IN)  :: efermi
      LOGICAL,       INTENT(IN)  :: film

      INTEGER :: ikpt,ncorn,itet,icorn,iband,k(4),fac
      REAL    :: eMesh(1),weight_tmp(1),etetra(4)

      eMesh(1) = efermi  !Only a single energy point but getWeightSingleBand takes an array

      !Tetrahedra or Triangles?
      ncorn = MERGE(3,4,film)

      weights = 0.0
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
            fac = MERGE(1,count(kpts%bkp(:).EQ.ikpt),kpts%nkptf.EQ.0)
            !$OMP PARALLEL DEFAULT(none) &
            !$OMP SHARED(itet,neig,ikpt,film,ncorn,k,fac) &
            !$OMP SHARED(kpts,eig,weights,eMesh) &
            !$OMP PRIVATE(iband,etetra,weight_tmp)
            !$OMP DO
            DO iband = 1, neig

               etetra(:ncorn) = eig(iband,k(:ncorn))
               IF( ALL(MAXVAL(eMesh)<=etetra) ) CYCLE

               weight_tmp = 0.0
               CALL getWeightSingleBand(eMesh,etetra,1,ncorn,ikpt,kpts%ntetra(:,itet),&
                                        kpts%voltet(itet)*fac,film,.FALSE.,weight_tmp(:))

               weights(iband,ikpt) = weights(iband,ikpt) + weight_tmp(1)
            ENDDO
            !$OMP END DO
            !$OMP END PARALLEL
         ENDDO
      ENDDO

   END SUBROUTINE getWeightKpoints

   SUBROUTINE getWeightEnergyMesh(kpts,ikpt,eig,neig,eMesh,ne,film,weights,bounds,dos)

      USE m_types
      USE m_juDFT
      USE m_differentiate

      IMPLICIT NONE

      TYPE(t_kpts),     INTENT(IN)  :: kpts
      REAL,             INTENT(IN)  :: eig(:,:)
      REAL,             INTENT(OUT) :: weights(:,:)
      INTEGER,OPTIONAL, INTENT(OUT) :: bounds(:,:)

      INTEGER,          INTENT(IN)  :: ikpt,neig,ne
      REAL,             INTENT(IN)  :: eMesh(:)
      LOGICAL,          INTENT(IN)  :: film
      LOGICAL,OPTIONAL, INTENT(IN)  :: dos

      INTEGER :: itet,iband,ncorn,ie,icorn,k(4),fac
      LOGICAL :: l_dos
      REAL    :: etetra(4),del
      REAL, ALLOCATABLE :: dos_weights(:), end_weights(:,:), occ_weights(:)

      !Tetrahedra or Triangles?
      ncorn = MERGE(3,4,film)
      IF(kpts%nkptf.NE.0) THEN
         fac = count(kpts%bkp(:).EQ.ikpt)
      ELSE
         fac = 1
      ENDIF

      l_dos = PRESENT(dos)
      IF(PRESENT(dos))THEN
         l_dos = dos.AND.ne>1
         IF(l_dos) THEN
            ALLOCATE(end_weights(2,neig))
            end_weights = 0.0
         ENDIF
      ENDIF
      IF(ne>1) del = eMesh(2)-eMesh(1)

      weights = 0.0
      DO itet = 1, kpts%ntet
         IF(ALL(kpts%ntetra(:ncorn,itet).NE.ikpt)) CYCLE
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
         !$OMP SHARED(itet,neig,ikpt,film,l_dos,ncorn,ne,k,fac) &
         !$OMP SHARED(kpts,eig,weights,eMesh,end_weights,del) &
         !$OMP PRIVATE(iband,etetra)
         !$OMP DO
         DO iband = 1, neig

            etetra(:ncorn) = eig(iband,k(:ncorn))
            IF( ALL(MAXVAL(eMesh)<etetra) ) CYCLE

            CALL getWeightSingleBand(eMesh,etetra,ne,ncorn,ikpt,kpts%ntetra(:,itet),&
                                     kpts%voltet(itet)*fac,film,.FALSE.,weights(:,iband))
            !---------------------------------------------------------------------------
            ! For the dos we want to avoid the formulas for the numerical
            ! derivatives at the endpoints so we calculate additional energy points here
            !---------------------------------------------------------------------------

            IF(l_dos) THEN
               !Energy point one step to the left of the grid
               CALL getWeightSingleBand((/MINVAL(eMesh)-del/),etetra,1,ncorn,ikpt,kpts%ntetra(:,itet),&
                                     kpts%voltet(itet)*fac,film,.FALSE.,end_weights(1:1,iband))
               !Energy point one step to the right of the grid
               CALL getWeightSingleBand((/MAXVAL(eMesh)+del/),etetra,1,ncorn,ikpt,kpts%ntetra(:,itet),&
                                     kpts%voltet(itet)*fac,film,.FALSE.,end_weights(2:2,iband))
            ENDIF

         ENDDO
         !$OMP END DO
         !$OMP END PARALLEL
      ENDDO

      !---------------------------------------------------
      ! Weights for DOS -> differentiate with respect to E
      !---------------------------------------------------
      IF(l_dos) THEN
         ALLOCATE(dos_weights(ne+2),occ_weights(ne+2))
         DO iband = 1, neig
            dos_weights = 0.0
            occ_weights(2:ne+1) = weights(1:ne,iband)
            occ_weights(1) = end_weights(1,iband)
            occ_weights(ne+2) = end_weights(2,iband)
            CALL diff3(occ_weights,del,dos_weights)
            weights(1:ne,iband) = dos_weights(2:ne+1)
         ENDDO
      ENDIF

      IF(PRESENT(bounds)) THEN
         DO iband = 1, neig
            !--------------------------------------------------------------
            !Find the range where the weights are bigger than weightCutoff
            !--------------------------------------------------------------
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
            IF(bounds(iband,1).GT.bounds(iband,2)) THEN
               bounds(iband,1) = 1
               bounds(iband,2) = 1
            ENDIF
            IF(ANY(weights(bounds(iband,1):bounds(iband,2),iband)<0.0)) THEN
               CALL juDFT_error("TetraWeight error: Unexpected negative weight", calledby="getWeightEnergyMesh")
            ENDIF
         ENDDO
      ENDIF

   END SUBROUTINE getWeightEnergyMesh

   SUBROUTINE getWeightSingleBand(eMesh,etetra,ne,ncorn,ikpt,ntetra,vol,film,bloechl,weight)

      !--------------------------------------------------------------
      ! This is the core routine calculating the integration
      ! weight for one kpoint in one tetrahedron for one band index
      ! for the provided energy points
      ! All other routines wrap this procedure so that it is most
      ! efficient in all cases
      !--------------------------------------------------------------

      USE m_juDFT
      USE m_tetsrt
      USE m_tetraWeight

      IMPLICIT NONE

      REAL,             INTENT(IN)     :: eMesh(:)
      REAL,             INTENT(IN)     :: etetra(:)
      INTEGER,          INTENT(IN)     :: ntetra(:)
      REAL,             INTENT(INOUT)  :: weight(:)

      INTEGER,          INTENT(IN)     :: ne,ncorn,ikpt
      REAL,             INTENT(IN)     :: vol
      LOGICAL,          INTENT(IN)     :: film,bloechl

      INTEGER :: icorn,i,ie,nstart,nend,ind(ncorn)
      REAL    :: w,eb,et,del,lastWeight

      !Sort the eigenvalues at the corners (ascending order in ind)
      CALL tetsrt(ncorn,etetra,ind)

      !search for the corner ikpt in the sorted array
      DO i = 1, ncorn
         IF(ntetra(ind(i)).EQ.ikpt) icorn = i
      ENDDO

      IF(ALL(ntetra.NE.ikpt)) CALL juDFT_error("kpoint not found in ntetra",&
                                               hint="This is a bug in FLEUR, please report"&
                                               ,calledby="getWeightSingleBand")

      !Find the range in the (equidistant) energy mesh where the weights are changing
      IF( ne>1 ) THEN
         eb = MINVAL(eMesh)
         et = MAXVAL(eMesh)
         del = eMesh(2)-eMesh(1)
         nstart = INT((etetra(ind(1))-eb)/del)+1
         nstart = MAX(1,nstart)
         nend   = INT((etetra(ind(ncorn))-eb)/del)+2
         nend   = MIN(ne,nend)
         !We need this to catch the case where all tetrahedron energies lie below the energy mesh
         nend   = MAX(1,nend)
      ELSE
         nstart = 1
         nend   = ne
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
      IF(nend.NE.ne) weight(nend+1:) = weight(nend+1:) + lastWeight

   END SUBROUTINE getWeightSingleBand

END MODULE m_tetrahedronInit