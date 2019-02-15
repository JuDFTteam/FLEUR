!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_chkmt
      USE m_juDFT
      private
      public chkmt
!---------------------------------------------------------------------
!  Check muffin tin radii and determine a reasonable choice for MTRs.
!  Derive also other parameters for the input file, to provide some
!  help in the out-file.
!                         GM'16
!---------------------------------------------------------------------
      CONTAINS
      SUBROUTINE chkmt(&
     &                 atoms,input,vacuum,cell,oneD,&
     &                 l_gga,noel,l_test,&
     &                 kmax,dtild,dvac1,lmax1,jri1,rmt1,dx1)

      USE m_types
      USE m_sort
      USE m_inv3
      USE m_juDFT
      IMPLICIT NONE
!     ..
!     .. Scalar Arguments ..
      TYPE(t_atoms),INTENT(IN) :: atoms
      TYPE(t_input),INTENT(IN) :: input
      TYPE(t_vacuum),INTENT(IN):: vacuum
      TYPE(t_cell),INTENT(IN)  :: cell
      TYPE(t_oneD),INTENT(IN)  :: oneD
      CHARACTER*3, INTENT (IN) :: noel(atoms%ntype)
      LOGICAL, INTENT (IN)     :: l_gga,l_test
      REAL,    INTENT (OUT)    :: kmax,dtild,dvac1
!     ..
!     .. Array Arguments ..
      INTEGER, INTENT (OUT)    :: lmax1(atoms%ntype),jri1(atoms%ntype)
      REAL,    INTENT (OUT)    :: rmt1(atoms%ntype),dx1(atoms%ntype)
!     ..
!     .. Local Scalars ..
      INTEGER na,n
      INTEGER i,j,k,l,jri11,lmax11
      INTEGER maxCubeAtoms, iAtom, numAtoms, iNeighborAtom, identicalAtoms
      INTEGER typeA, typeB
      REAL    dx11,rkm,sum_r,facA,facB
      REAL    rmtMax, rmtMin, rmtMaxDefault, rmtDelta
      REAL    rmtFac, cubeLength, amatAuxDet
      REAL    maxSqrDist, dist, currentDist
      LOGICAL error, outOfBounds
!     ..
!     .. Local Arrays ..
      REAL    t_rmt(0:103), minRmts(0:103)
      REAL    amatAux(3,3), invAmatAux(3,3)
      REAL    taualAux(3,atoms%nat), posAux(3,atoms%nat)
      REAL    minPos(3), maxPos(3), pos(3), point(3), realCellPos(3)
      REAL    offsetPos(3)
      REAL    nearestAtomDists(atoms%ntype)
      INTEGER nearestAtoms(atoms%ntype)
      INTEGER sortedDistList(atoms%ntype)
      INTEGER minCubeIndex(3), maxCubeIndex(3), cubeIndex(3)
      INTEGER minCellIndices(3), maxCellIndices(3)

      INTEGER, ALLOCATABLE :: numAtomsInCubes(:,:,:)
      INTEGER, ALLOCATABLE :: atomRefsInCubes(:,:,:,:)
      INTEGER, ALLOCATABLE :: refCubes(:,:)
      INTEGER, ALLOCATABLE :: nearestNeighbors(:,:)
      INTEGER, ALLOCATABLE :: numNearestNeighbors(:)
      INTEGER, ALLOCATABLE :: neighborAtoms(:)
      INTEGER, ALLOCATABLE :: distIndexList(:)
      REAL,    ALLOCATABLE :: posInCubes(:,:,:,:,:)
      REAL,    ALLOCATABLE :: refPos(:,:)
      REAL,    ALLOCATABLE :: nearestNeighborDists(:,:)
      REAL,    ALLOCATABLE :: sqrDistances(:)

!     Plan for this routine:
!     0. Do initializations and set some constants
!     1. Locally replace unit cell by an auxiliary unit cell with:
!        a) all atoms within the unit cell
!        b) basis vectors obtained by lattice reduction of the original cell. 
!           [not in 1st (this) version of routine. Can be implemented with the LLL algorithm when needed.]
!     2. Get minimal and maximal coordinates within auxiliary unit cell
!     3. Construct mesh of cubes covering the auxiliary unit cell and a boundary of width 2*rmtMax + rmtDelta
!     4. Fill mesh of cubes with atoms
!        a) Store atoms in cubes and representative cube for each atom type
!     5. For each atom in auxiliary unit cell select cube and collect shortest distances to other atoms in neighborhood
!        a) Sort distances and set MT radii for the atoms
!     6. Correct bad choices and set missing MT radii, vacuum distances, and other parameters
!     7. Test old MT radii


!     0. Do initializations and set some constants

      rmtMaxDefault = 2.8
      rmtMax = rmtMaxDefault
      rmtMin = 1.0
      IF (l_test) THEN
         rmtMax = MAX(rmtMax,MAXVAL(atoms%rmt(:)))
         rmtMin = MIN(rmtMin,MINVAL(atoms%rmt(:)))
      END IF
      rmtDelta = 0.3
      IF (input%film) THEN
        rmtFac = 0.95
      ELSE
        rmtFac = 0.975
      ENDIF
      t_rmt(0:103) = 2.3 ! default value
      t_rmt(1) = 1.0 ; t_rmt(5:9) = 1.3 ; t_rmt(16:17) = 1.8
      cubeLength = 2*rmtMax+rmtDelta
      maxCubeAtoms = (FLOOR(cubeLength / (0.7*2.0*rmtMin)) + 1)**3
      error  = .FALSE.

!     1. For the 1st version the auxiliary unit cell is just a copy of the original unit cell with
!        all atoms within the cell.

      DO i = 1, 3
         DO j = 1, 3
            amatAux(i,j) = cell%amat(i,j)
         END DO
      END DO

      DO i = 1, atoms%nat
         taualAux(1,i) = atoms%taual(1,i) - FLOOR(atoms%taual(1,i))
         taualAux(2,i) = atoms%taual(2,i) - FLOOR(atoms%taual(2,i))
         taualAux(3,i) = atoms%taual(3,i) - FLOOR(atoms%taual(3,i))
         posAux(:,i) = MATMUL(amatAux,taualAux(:,i))
      END DO

!     2. Get minimal and maximal coordinates for auxiliary unit cell

      minPos = 0.0
      maxPos = 0.0

      DO i = 0, 1
         DO j = 0, 1
            DO k = 0, 1
               DO l = 1, 3
                  pos(l) = i*amatAux(l,1) + j*amatAux(l,2) + k*amatAux(l,3)
                  IF (pos(l).GT.maxPos(l)) maxPos(l) = pos(l)
                  IF (pos(l).LT.minPos(l)) minPos(l) = pos(l)
               END DO
            END DO
         END DO
      END DO
      
!     3. Construct cube mesh:
!        In each dimension cube i covers the interval from i*cubeLength to (i+1)*cubeLength
!        Each cube may cover up to maxCubeAtoms atoms. This should be set to a save size.

      DO i = 1, 3
         minPos(i) = minPos(i) - cubeLength
         maxPos(i) = maxPos(i) + cubeLength
         minCubeIndex(i) = FLOOR(minPos(i)/cubeLength)
         maxCubeIndex(i) = CEILING(maxPos(i)/cubeLength)
      END DO

      ALLOCATE (numAtomsInCubes(minCubeIndex(1):maxCubeIndex(1),&
                                minCubeIndex(2):maxCubeIndex(2),&
                                minCubeIndex(3):maxCubeIndex(3)))
      ALLOCATE (atomRefsInCubes(maxCubeAtoms,minCubeIndex(1):maxCubeIndex(1),&
                                             minCubeIndex(2):maxCubeIndex(2),&
                                             minCubeIndex(3):maxCubeIndex(3)))
      ALLOCATE (posInCubes(3,maxCubeAtoms,minCubeIndex(1):maxCubeIndex(1),&
                                          minCubeIndex(2):maxCubeIndex(2),&
                                          minCubeIndex(3):maxCubeIndex(3)))
      ALLOCATE (refCubes(3,atoms%ntype),refPos(3,atoms%ntype))
      ALLOCATE (nearestNeighbors(maxCubeAtoms,atoms%ntype),numNearestNeighbors(atoms%ntype))
      ALLOCATE (nearestNeighborDists(maxCubeAtoms,atoms%ntype))

      numAtomsInCubes = 0

!     4. Fill mesh of cubes with atoms
!        First obtain minimal and maximal indices for relevant unit cells

      minCellIndices = 0
      maxCellIndices = 0

      CALL inv3(amatAux,invAmatAux,amatAuxDet)

      DO i = 0, 1
         DO j = 0, 1
            DO k = 0, 1
               point(:) = minPos(:)
               IF(i.EQ.1) point(1) = maxPos(1)
               IF(j.EQ.1) point(2) = maxPos(2)
               IF(k.EQ.1) point(3) = maxPos(3)
               realCellPos(:) = matmul(invAmatAux,point(:))
               DO l = 1, 3
                  IF(minCellIndices(l).GT.realCellPos(l)) THEN
                     minCellIndices(l) = FLOOR(realCellPos(l))
                  END IF
                  IF(maxCellIndices(l).LT.realCellPos(l)) THEN
                     maxCellIndices(l) = FLOOR(realCellPos(l)) ! Is 'FLOOR' enough?
                  END IF
               END DO
            END DO
         END DO
      END DO

!        Store atoms in cubes and representative cube for each atom type

      DO i = minCellIndices(1), maxCellIndices(1)
         DO j = minCellIndices(2), maxCellIndices(2)
            DO k = minCellIndices(3), maxCellIndices(3)
               DO l = 1, 3
                  offsetPos(l) = i*amatAux(l,1) + j*amatAux(l,2) + k*amatAux(l,3)
               END DO
               iAtom = 0
               DO n = 1, atoms%ntype
                  DO na = 1, atoms%neq(n)
                     iAtom = iAtom + 1
                     pos(:) = posAux(:,iAtom) + offsetPos(:)
                     outOfBounds = .FALSE.
                     DO l = 1, 3
                        cubeIndex(l) = FLOOR(pos(l)/cubeLength)
                        IF(cubeIndex(l).LT.minCubeIndex(l)) outOfBounds = .TRUE.
                        IF(cubeIndex(l).GT.maxCubeIndex(l)) outOfBounds = .TRUE.
                     END DO
                     IF(.NOT.outOfBounds) THEN
                        numAtomsInCubes(cubeIndex(1),cubeIndex(2),cubeIndex(3)) = &
                           numAtomsInCubes(cubeIndex(1),cubeIndex(2),cubeIndex(3)) + 1
                        numAtoms = numAtomsInCubes(cubeIndex(1),cubeIndex(2),cubeIndex(3))
                        IF(numAtoms.GT.maxCubeAtoms) THEN
                           STOP 'ERROR: maxCubeAtoms is not large enough in chkmt.'
                        END IF
                        atomRefsInCubes(numAtoms,cubeIndex(1),cubeIndex(2),cubeIndex(3)) = n
                        posInCubes(:,numAtoms,cubeIndex(1),cubeIndex(2),cubeIndex(3)) = pos(:)
                        IF((i.EQ.0).AND.(j.EQ.0).AND.(k.EQ.0).AND.(na.EQ.1)) THEN
                           refCubes(:,n) = cubeIndex(:)
                           refPos(:,n) = pos(:)
                        END IF
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO

!     5. For each atom type in auxiliary unit cell select cube and collect shortest distances 
!        to other atoms in neighborhood

      maxSqrDist = cubeLength**2
      ALLOCATE(sqrDistances(8*maxCubeAtoms)) ! Formally 27, but 8 should be enough due to maxSqrDist
      ALLOCATE(neighborAtoms(8*maxCubeAtoms))
      ALLOCATE(distIndexList(8*maxCubeAtoms))

      DO n = 1, atoms%ntype
         cubeIndex(:) = refCubes(:,n)
         neighborAtoms = 0
         iNeighborAtom = 0
         identicalAtoms = 0
         DO i = cubeIndex(1) - 1, cubeIndex(1) + 1
            DO j = cubeIndex(2) - 1, cubeIndex(2) + 1
               DO k = cubeIndex(3) - 1, cubeIndex(3) + 1
                  DO iAtom = 1, numAtomsInCubes(i,j,k)
                     currentDist = (refPos(1,n) - posInCubes(1,iAtom,i,j,k))**2 + &
                                   (refPos(2,n) - posInCubes(2,iAtom,i,j,k))**2 + &
                                   (refPos(3,n) - posInCubes(3,iAtom,i,j,k))**2
                     IF (currentDist.LT.0.000001) THEN
                        identicalAtoms = identicalAtoms + 1
                     ELSE IF (currentDist.LT.maxSqrDist) THEN
                        iNeighborAtom = iNeighborAtom + 1
                        neighborAtoms(iNeighborAtom) = atomRefsInCubes(iAtom,i,j,k)
                        sqrDistances(iNeighborAtom) = currentDist
                     END IF
                  END DO
               END DO
            END DO
         END DO
         IF (identicalAtoms.GT.1) THEN
            WRITE(*,*) 'Position: ', refPos(:,n)
            CALL juDFT_error("Too many atoms at same position.",calledby ="chkmt")
         END IF
         numNearestNeighbors(n) = MIN(maxCubeAtoms,iNeighborAtom)
         CALL sort(distIndexList(:iNeighborAtom),sqrDistances(:iNeighborAtom))
         DO i = 1, numNearestNeighbors(n)
            nearestNeighbors(i,n) = neighborAtoms(distIndexList(i))
            nearestNeighborDists(i,n) = SQRT(sqrDistances(distIndexList(i)))
         END DO
      END DO

      DO i = 1, atoms%ntype
         IF(numNearestNeighbors(i).GE.1) THEN
            nearestAtoms(i) = nearestNeighbors(1,i)
            nearestAtomDists(i) = nearestNeighborDists(1,i)
         ELSE
            nearestAtoms(i) = -1
            nearestAtomDists(i) = 5000.0 * cubeLength
         END IF
      END DO

!        Sort distances and set MT radii for the atoms

      CALL sort(sortedDistList,nearestAtomDists)
      rmt1 = -1.0
      minRmts = -1.0
      DO i = 1, atoms%ntype
         typeA = sortedDistList(i)
         typeB = nearestAtoms(typeA)
         IF(typeB.LT.0) CYCLE
         dist = nearestAtomDists(typeA)
         IF (dist.LT.0.5) THEN
            WRITE (*,*) "Distance between atoms too small!"
            WRITE (*,*) "atom type A: ", typeA
            WRITE (*,*) "atom type B: ", typeB
            WRITE (*,*) "distance: ", dist
            CALL juDFT_error("Distance between atoms too small!",calledby ="chkmt")
         END IF
         sum_r = 1.0 / ( t_rmt(atoms%nz(typeA)) + t_rmt(atoms%nz(typeB)) )
         facA = t_rmt(atoms%nz(typeA)) * sum_r
         facB = t_rmt(atoms%nz(typeB)) * sum_r
         ! Note: The result of this section may be slightly different from the old version
         !       iff the nearest atom is another atom of the same type
         IF (minRmts(atoms%nz(typeA)).LT.0.0) THEN
            IF (minRmts(atoms%nz(typeB)).LT.0.0) THEN
               minRmts(atoms%nz(typeA)) = rmtFac * dist * facA
               minRmts(atoms%nz(typeB)) = rmtFac * dist * facB
            ELSE
               minRmts(atoms%nz(typeA)) = rmtFac * (dist - minRmts(atoms%nz(typeB)))
            END IF
         ELSE IF (minRmts(atoms%nz(typeB)).LT.0.0) THEN
           minRmts(atoms%nz(typeB)) = rmtFac * (dist - minRmts(atoms%nz(typeA)))
         END IF
      END DO

!     6. Correct bad choices and set missing MT radii, vacuum distances, and other parameters

      DO i = 1, atoms%ntype
         IF((minRmts(atoms%nz(i)).LT.0.0).OR.(minRmts(atoms%nz(i)).GE.rmtMaxDefault)) THEN
            minRmts(atoms%nz(i)) = rmtMaxDefault
         END IF
         rmt1(i) = minRmts(atoms%nz(i))
      END DO

      ! NOTE: The result of this section may be slightly different from the old version
      !       iff the old version would enlarge a MT sphere at this point.
      !       Also the old version does not propagate the changes of the MT radii to all
      !       atoms with the same atomic number
      DO i = 1, atoms%ntype
         DO j = 1, numNearestNeighbors(i)
            k = nearestNeighbors(j,i)
            IF (rmt1(i)+rmt1(k).GE.nearestNeighborDists(j,i)) THEN
               minRmts(atoms%nz(i)) = MIN(rmtFac*nearestNeighborDists(j,i)/2.0,MIN(rmt1(i),minRmts(atoms%nz(i))))
               minRmts(atoms%nz(k)) = MIN(rmtFac*nearestNeighborDists(j,i)/2.0,MIN(rmt1(k),minRmts(atoms%nz(k))))
            END IF
         END DO
      END DO

      DO i = 1, atoms%ntype
         rmt1(i) = minRmts(atoms%nz(i))
      END DO

      WRITE (6,*) '----------------------------------------------------'
      WRITE (6,*) 'Suggested values for input: '
      WRITE (6,*) 

      dvac1 = 0.0
      IF (input%film) THEN
         iAtom = 0
         DO i = 1, atoms%ntype
            DO na = 1, atoms%neq(i)
               iAtom = iAtom + 1
               IF (oneD%odd%d1) THEN
                  dvac1 = MAX(dvac1, SQRT(atoms%pos(1,iAtom)**2+atoms%pos(2,iAtom)**2)+rmt1(i))
               ELSE
                  dvac1 = MAX(dvac1, ABS(atoms%pos(3,iAtom))+rmt1(i))
               END IF
            END DO
         END DO
         dvac1 = 2.0 * (dvac1+0.3)
         dtild = dvac1 + 1.5 * MAXVAL(rmt1(:))
         WRITE (6,'("vacuum distance dvac =",f10.5)') dvac1
         WRITE (6,'("extra vac.dist. dtild=",f10.5)') dtild
      END IF

      rkm = 0.0
      WRITE (6,*) 'Atom Z  lmax jri    rmt         dx'
      DO n = 1, atoms%ntype
         IF (rmt1(n).LT.1.8) THEN
            lmax11 = 6
         ELSE IF (rmt1(n).LT.2.4) THEN
            lmax11 = 8
         ELSE 
            lmax11 = 10
         END IF
         IF (l_gga) THEN
            jri11 = NINT(330*rmt1(n)) 
         ELSE
            jri11 = NINT(220*rmt1(n)) 
         END IF
         jri11 = NINT(jri11*0.5) * 2 + 1
         IF (atoms%nz(n) > 0) THEN
           dx11 = LOG(3200*atoms%nz(n)*rmt1(n))/(jri11-1)
         ELSE
           dx11 = LOG(3200*rmt1(n))/(jri11-1)
         ENDIF
         rkm = MAX(rkm, lmax11/rmt1(n))
         WRITE (6,'(a3,i3,2i5,2f10.6)') noel(n),atoms%nz(n),lmax11,jri11,rmt1(n),dx11
         dx1(n) = dx11
         lmax1(n) = lmax11
         jri1(n) = jri11
      END DO
      WRITE (6,'("k_max =",f8.5)') rkm
      WRITE (6,'("G_max =",f8.5)') 3*rkm
      kmax = rkm

!     7. Test old MT radii

      IF (l_test) THEN
         iAtom = 0
         DO i = 1, atoms%ntype
            DO j = 1, numNearestNeighbors(i)
               k = nearestNeighbors(j,i)
               IF (atoms%rmt(i)+atoms%rmt(k).GE.nearestNeighborDists(j,i)) THEN
                  error = .TRUE.
                  WRITE(6,240) i,k,nearestNeighborDists(j,i),atoms%rmt(i),atoms%rmt(k)
               END IF
            END DO
            IF (input%film) THEN
               DO na = 1, atoms%neq(i)
                  iAtom = iAtom + 1
                  IF (oneD%odd%d1) THEN
                     IF ((sqrt(atoms%pos(1,iAtom)**2+atoms%pos(2,iAtom)**2)+&
                         atoms%rmt(i)).GT.vacuum%dvac/2.) THEN
                        error=.TRUE.
                        WRITE(6,241) i ,na
                        WRITE(6,*) sqrt(atoms%pos(1,iAtom)**2+atoms%pos(2,iAtom)**2),&
                                   atoms%rmt(i),vacuum%dvac/2.
                     END IF
                  ELSE
                     IF (((atoms%pos(3,iAtom)+atoms%rmt(i)).GT. vacuum%dvac/2.).OR.&
                         ((atoms%pos(3,iAtom)-atoms%rmt(i)).LT.-vacuum%dvac/2.)) THEN
                        error=.TRUE.
                        WRITE(6,241) i ,na
                        WRITE(6,*) atoms%pos(3,iAtom),atoms%rmt(i),vacuum%dvac/2.
                     ENDIF
                  ENDIF
               END DO
            END IF
         END DO
         IF (error) CALL juDFT_error("Error checking M.T. radii",calledby ="chkmt")
      END IF

      DEALLOCATE(nearestNeighbors,numNearestNeighbors,nearestNeighborDists)
      DEALLOCATE(distIndexList,neighborAtoms,sqrDistances)
      DEALLOCATE(numAtomsInCubes,atomRefsInCubes,posInCubes,refCubes,refPos)

  240 FORMAT('Error in muffin tin radii pair (',i5,',',i5,'):',3f10.5)
  241 FORMAT('   error: atom ',i3,' # ',i3,'reaches out into vaccuum')

      END SUBROUTINE chkmt
      END MODULE m_chkmt
