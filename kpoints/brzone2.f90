!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_brzone2
USE m_juDFT

CONTAINS
SUBROUTINE brzone2(rcmt,nsym,idrot,mface,nbsz,nv48,&
                   cpoint,xvec,ncorn,nedge,nface,fnorm,fdist)

   USE m_constants, ONLY : pimach

   IMPLICIT NONE

   ! This subroutine constructs the irreducible wedge of the Brillouin
   ! zone (IBZ). The subroutine signature is based on an earlier
   ! implementation of this functionality in brzone but the algorithm 
   ! is a new design as the old routine featured too many bugs. Some 
   ! code parts are taken from the old routine.
   !
   !                                  GM, 2016

   INTEGER, INTENT (IN) :: mface,nbsz,nv48
   INTEGER, INTENT (IN) :: nsym               ! number of symmetry elements
   REAL,    INTENT (IN) :: rcmt(3,3)          ! reciprocal lattice basis (2\pi/a.u.)
   REAL,    INTENT (IN) :: idrot(3,3,48)      ! rotation matrices in cartesian repr.

   INTEGER, INTENT (OUT) :: ncorn,nedge,nface ! number of corners, faces and edges of the IBZ
   REAL,    INTENT (OUT) :: fnorm(3,mface)    ! normal vector of the planes bordering the IBZ
   REAL,    INTENT (OUT) :: fdist(mface)      ! distance vector of the planes bordering the IBZ
   REAL,    INTENT (OUT) :: cpoint(3,mface)   ! cartesian coordinates of corner points of IBZ
   REAL,    INTENT (OUT) :: xvec(3)           ! arbitrary vector lying in the IBZ

   INTEGER               :: info, n1, n2, n3, n4, i, j, k, iPlane, iSym
   INTEGER               :: nPlanes, nOuterPlanes, nCorners, numIBZPlanes
   INTEGER               :: maxNumPlanes, nUniqueCorners, duplicateNum
   INTEGER               :: nbszLocal, planesDim
   INTEGER               :: ipiv(3), krecip(3)
   INTEGER               :: cornerPlanes(3,5000)
   INTEGER               :: numCornerPlanes(5000)
   LOGICAL               :: foundDuplicate, filterOut
   LOGICAL               :: isIBZCorner(5000)
   REAL, PARAMETER       :: eps09 = 1.0e-9
   REAL                  :: pi, maxRecDist, vecDist, recScale, scalarProduct
   REAL                  :: norm, normalScalarProduct
   REAL                  :: denominator, vec1Fac, vec2Fac, edgeDist
   REAL                  :: testVec(3), sk(3), vecA(3), vecB(3), vecC(3)
   REAL                  :: corners(3,5000)
   REAL                  :: planesMatrix(3,3), solutions(3,1)
   REAL                  :: equationSystem(3,3)
   REAL                  :: edgeDirec(3), edgeDistVec(3), distVec(3)

   INTEGER, ALLOCATABLE  :: cornerPlaneList(:,:)
   INTEGER, ALLOCATABLE  :: planeCorners(:,:)
   INTEGER, ALLOCATABLE  :: nPlaneCorners(:)
   LOGICAL, ALLOCATABLE  :: isIBZPlane(:)

   REAL, ALLOCATABLE     :: dvec(:,:), ddist(:)

   ! Subroutine plan:
   ! 0. Initializations
   ! 1. Construct all possible boundary planes. Remove duplicate planes.
   ! 1.1. Construct boundary planes according to bisections of the lines
   !      to reciprocal lattice vectors <> 0.
   ! 1.2. Construct all boundary planes due to symmetry relations.
   ! 2. Determine all corner points of the IBZ.
   ! 2.1. Determine all possibly relevant intersection points of each 
   !      subset of 3 boundary planes, respectively.
   ! 2.2. Select those points that are on the correct side of all of the
   !      possible boundary planes.
   ! 3. Remove all boundary planes that don't feature at least three 
   !    corner points. These planes are not part of the real boundary.
   !    Also remove the corner points associated with each removed plane.
   ! 4. Construct edges from the set of boundary points. Actually
   !    only the number of edges is needed. I might use the Euler 
   !    characteristic to calculate it.

   ! 0. Initializations

   pi = pimach()

   maxRecDist = 0.0
   DO n1 = -1,1,2
      DO n2 = -1,1,2
         DO n3 = -1,1,2
            testVec(:) = n1*rcmt(:,1)+n2*rcmt(:,2)+n3*rcmt(:,3)
            vecDist = SQRT(testVec(1)**2.0+testVec(2)**2.0+testVec(3)**2.0)
            IF(vecDist.GT.maxRecDist) maxRecDist = vecDist
         END DO
      END DO
   END DO
   maxRecDist = 1.01 * maxRecDist

   ! We ignore the nbsz (the number of neighboring cells to consider) passed to the routine
   ! and calculate an own nbszLocal on the basis of maxRecDist.

   nbszLocal = 1
   DO i = -1, 1, 2
      DO j = -1, 1, 2
         DO k = -1, 1, 2
            solutions(1,1) = i * maxRecDist
            solutions(2,1) = j * maxRecDist
            solutions(3,1) = k * maxRecDist
            equationSystem(:,:) = rcmt(:,:)
            ipiv = 0
            info = 0
            CALL DGETRF(3,3, equationSystem,3,ipiv,info)
            CALL DGETRS('N',3,1,equationSystem,3,ipiv,solutions,3,info)
            ! I assume that info == 0: The reciprocal lattice vectors should be linearly independent.
            if(nbszLocal.LT.ABS(solutions(1,1))) nbszLocal = CEILING(ABS(solutions(1,1)))
            if(nbszLocal.LT.ABS(solutions(2,1))) nbszLocal = CEILING(ABS(solutions(2,1)))
            if(nbszLocal.LT.ABS(solutions(3,1))) nbszLocal = CEILING(ABS(solutions(3,1)))
         END DO
      END DO
   END DO

   maxRecDist = 0.5 * maxRecDist

   ! 1. Construct all possible boundary planes. Remove duplicate planes.
   ! 1.1. Construct boundary planes according to bisections of the lines
   !      to reciprocal lattice vectors <> 0

   ! 1.1.1. Construct an arbitrary point xvec in the IBZ
   ! 1.1.1.1. Determine a "scale" for the point on the basis of the
   !          reciprocal lattice basis.
   DO i = 1,3
      sk(i) = 0.0
      DO j = 1,3
         sk(i)=sk(i)+rcmt(j,i)**2
      END DO
   END DO
   recScale = SQRT(MIN(sk(1),sk(2),sk(3)))*0.1

   ! 1.1.1.2. Use the scale to construct the arbitrary point xvec.
   xvec(1) = recScale
   xvec(2) = recScale/SQRT(pi)
   xvec(3) = recScale/pi

   ! loops over all neighboring lattice vectors. The size of the 
   ! neigborhood is defined by nbszLocal.

   planesDim = 2*nbszLocal+1
   planesDim = planesDim*planesDim*planesDim + nsym
   ALLOCATE(dvec(3,planesDim), ddist(planesDim))
   ALLOCATE(nPlaneCorners(planesDim),isIBZPlane(planesDim))

   nPlanes = 0
   DO n1 = -nbszLocal, nbszLocal
      krecip(1) = n1
      DO n2 = -nbszLocal, nbszLocal
         krecip(2) = n2
         DO n3 = -nbszLocal, nbszLocal
            IF ( .NOT.(n1.EQ.0.AND.n2.EQ.0.AND.n3.EQ.0) ) THEN
               krecip(3) = n3
               nPlanes = nPlanes + 1

               ! determine distance vector dvec to the selected reciprocal
               ! lattice vector
               DO i = 1,3
                  dvec(i,nPlanes) = 0.0
                  DO j = 1,3
                     dvec(i,nPlanes) = dvec(i,nPlanes) + rcmt(i,j)*krecip(j)
                  END DO
               END DO

               ! Determine the norm of the distance vector
               norm = 0.0
               DO i = 1,3
                  norm = norm + dvec(i,nPlanes)**2
               END DO
               norm = SQRT(norm)

               ! Deterime the distance for the intersection
               ddist(nPlanes) = 0.5*norm

               IF(ddist(nPlanes).GT.maxRecDist) THEN
                  nPlanes = nPlanes - 1
                  CYCLE
               END IF

               ! Normalize distance vector. Note that this is now the
               ! normal to the intersecting plane. The intersecting plane 
               ! is defined by the normal and the distance.
               ! Note: With the current sign the vector points outwards.
               dvec(:,nPlanes) = dvec(:,nPlanes) / norm

               ! Remove plane again if it is a duplicate of another plane
               ! or only use the plane nearer to the origin if it is parallel 
               ! to another plane and has the same direction of the normal.
               foundDuplicate = .FALSE.
               DO iPlane = 1, nPlanes-1
                  IF (ABS(dvec(1,iPlane)-dvec(1,nPlanes)).GT.eps09) CYCLE
                  IF (ABS(dvec(2,iPlane)-dvec(2,nPlanes)).GT.eps09) CYCLE
                  IF (ABS(dvec(3,iPlane)-dvec(3,nPlanes)).GT.eps09) CYCLE
                  IF (ddist(nPlanes).LT.ddist(iPlane)) ddist(iPlane) = ddist(nPlanes)
                  foundDuplicate = .TRUE.
                  EXIT
               END DO
               IF(foundDuplicate) nPlanes = nPlanes - 1
            END IF
         END DO
      END DO
   END DO

   nOuterPlanes = nPlanes

   ! 1.2. Construct all boundary planes due to symmetry relations.
   !      These are the planes bisecting the line connecting xvec
   !      with an element of the star of xvec ( <> xvec )
   !      Note that the loop over the symmetry operations actually
   !      is a loop over all group element of all the relevant 
   !      symmetry groups. These are not only the generating elements
   !      and therefore it is not neccessary to invert the matrices to
   !      also obtain boundary planes on the other side.

   DO iSym = 2, nsym ! leave out identity symmetry operation
      nPlanes = nPlanes + 1

      ! The origin is part of all the planes determined in this way.
      ! -> distance = 0.
      ddist(nPlanes) = 0.0

      ! Determine the vector connecting xvec and its image after performing
      ! the symmetry operation.
      DO i = 1,3
         dvec(i,nPlanes)=-xvec(i)
         DO j=1,3
            dvec(i,nPlanes) = dvec(i,nPlanes) + idrot(i,j,iSym)*xvec(j)
         END DO
      END DO

      ! Normalize the vector (it is normal to the possible boundary plane).
      ! Note that the vector points away from xvec.
      norm = 0.0
      DO i = 1,3
         norm = norm + dvec(i,nPlanes)**2
      END DO
      norm = SQRT(norm)
      dvec(:,nPlanes)=dvec(:,nPlanes) / norm

      ! Remove plane again if it is a duplicate of another plane
      foundDuplicate = .FALSE.
      DO iPlane = nOuterPlanes+1, nPlanes-1
         IF (ABS(ddist(iPlane)-ddist(nPlanes)).GT.eps09) CYCLE
         IF (ABS(dvec(1,iPlane)-dvec(1,nPlanes)).GT.eps09) CYCLE
         IF (ABS(dvec(2,iPlane)-dvec(2,nPlanes)).GT.eps09) CYCLE
         IF (ABS(dvec(3,iPlane)-dvec(3,nPlanes)).GT.eps09) CYCLE
         foundDuplicate = .TRUE.
         EXIT
      END DO
      IF(foundDuplicate) nPlanes = nPlanes - 1
   END DO

   ! 2. Determine all corner points of the IBZ.
   ! 2.1. Determine all possibly relevant intersection points of each 
   !      subset of 3 boundary planes, respectively.


   nCorners = 0
   nPlaneCorners = 0
   DO n1 = 1, nOuterPlanes
      DO n2 = n1+1, nPlanes
         ! 2.1.1. Calculate intersection edge between planes n1 and n2
         !        and only consider those n1 and n2 where the cuttig 
         !        edge is not too far from the origin. (This would mean
         !        that the crossing point of these two planes and a third 
         !        plane is never relevant for the construction of the
         !        IBZ.)

         ! The direction of the intersection edge is the cross product
         ! of the normals on planes n1 and n2:
         edgeDirec(1) = dvec(2,n1)*dvec(3,n2) - dvec(3,n1)*dvec(2,n2)
         edgeDirec(2) = dvec(3,n1)*dvec(1,n2) - dvec(1,n1)*dvec(3,n2)
         edgeDirec(3) = dvec(1,n1)*dvec(2,n2) - dvec(2,n1)*dvec(1,n2)

         ! Ignore parallel planes
         IF ((ABS(edgeDirec(1)).LT.eps09).AND.&
             (ABS(edgeDirec(2)).LT.eps09).AND.&
             (ABS(edgeDirec(3)).LT.eps09)) CYCLE

         ! The distance vector of the intersection edge to the origin is given
         ! by (since dvec is normalized):
         !
         !       ddist(n1) - ddist(n2) * <dvec(:,n1)|dvec(:,n2)>              ddist(n2) - ddist(n1) * <dvec(:,n1)|dvec(:,n2)>
         ! d_e = ----------------------------------------------- dvec(:,n1) + ----------------------------------------------- dvec(:,n2)
         !              1.0 - <dvec(:,n1)|dvec(:,n2)>^2                               1.0 - <dvec(:,n1)|dvec(:,n2)>^2

         normalScalarProduct = dvec(1,n1)*dvec(1,n2)+&
                               dvec(2,n1)*dvec(2,n2)+&
                               dvec(3,n1)*dvec(3,n2)
         denominator = 1.0 - (normalScalarProduct**2.0)
         vec1Fac = ddist(n1) - ddist(n2) * normalScalarProduct
         vec2Fac = ddist(n2) - ddist(n1) * normalScalarProduct
         edgeDistVec(:) = vec1Fac*dvec(:,n1) + vec2Fac*dvec(:,n2)
         edgeDistVec(:) = edgeDistVec(:) / denominator
         edgeDist = SQRT(edgeDistVec(1)**2.0+edgeDistVec(2)**2.0+edgeDistVec(3)**2.0)

         ! Ignore planes if intersection edge is too far from origin
         ! (is this criterion ok?)
         IF (edgeDist.GT.maxRecDist) CYCLE

         ! Now calculate the possibly relevant crossing points
         innerPlaneLoop: DO n3 = n2+1, nPlanes
            ! Set up system of linear equations to determine the crossing point
            DO i = 1, 3
               planesMatrix(1,i) = dvec(i,n1)
               planesMatrix(2,i) = dvec(i,n2)
               planesMatrix(3,i) = dvec(i,n3)
            END DO
            solutions(1,1) = ddist(n1)
            solutions(2,1) = ddist(n2)
            solutions(3,1) = ddist(n3)
            ! Solve system of linear equations and cycle if no solution is found
            ! or solution is not in relevant region.
            ipiv = 0
            info = 0
            CALL DGETRF(3,3, planesMatrix,3,ipiv,info)
            IF(info.NE.0) CYCLE
            CALL DGETRS('N',3,1,planesMatrix,3,ipiv,solutions,3,info)
            IF(info.NE.0) CYCLE
            vecDist = SQRT(solutions(1,1)**2.0+solutions(2,1)**2.0+solutions(3,1)**2.0)
            IF(vecDist.GT.maxRecDist) CYCLE

   ! 2.2. Select those points that are on the correct side of all of the
   !      possible boundary planes.

            DO n4 = 1, nPlanes
               vecDist = dvec(1,n4)*solutions(1,1)+dvec(2,n4)*solutions(2,1)+dvec(3,n4)*solutions(3,1) - ddist(n4)
               IF(vecDist*10..GT.eps09) THEN
                  CYCLE innerPlaneLoop
               END IF
            END DO

            !Add new crossing point to list and associate it with the 3 planes

            nCorners = nCorners + 1
            corners(:,nCorners) = solutions(:,1)
            cornerPlanes(1,nCorners) = n1
            cornerPlanes(2,nCorners) = n2
            cornerPlanes(3,nCorners) = n3
            nPlaneCorners(n1) = nPlaneCorners(n1) + 1
            nPlaneCorners(n2) = nPlaneCorners(n2) + 1
            nPlaneCorners(n3) = nPlaneCorners(n3) + 1
         END DO innerPlaneLoop
      END DO
   END DO

   ! 3. Remove all boundary planes that don't feature at least three 
   !    corner points. These planes are not part of the real boundary.
   !    Also remove the corner points associated with each removed plane.
   !    ...and also those planes and corners that are irrelevant due to 
   !    other reasons.

   ! Corners might be present multiple times (for different triples of planes). 
   ! Remove duplicates. Collect planes intersecting at the respective point.
   maxNumPlanes = MAX(nPlanes-nOuterPlanes,3*nCorners) ! Just a save guess for the array dimension in the next line
   ALLOCATE(cornerPlaneList(maxNumPlanes,nCorners+1)) ! The +1 is because of the origin that is added later.
   cornerPlaneList = 0
   numCornerPlanes = 0

   nUniqueCorners = 0
   DO i = 1, nCorners
      duplicateNum = -1
      compareLoop: DO j = 1, i-1
         distVec(:) = corners(:,i) - corners(:,j)
         vecDist = SQRT(distVec(1)**2+distVec(2)**2+distVec(3)**2)
         IF (vecDist.LT.eps09) THEN
            duplicateNum = j
            EXIT compareLoop
         END IF 
      END DO compareLoop
      IF (duplicateNum.EQ.-1) THEN
         nUniqueCorners = nUniqueCorners + 1
         numCornerPlanes(nUniqueCorners) = 3
         corners(:,nUniqueCorners) = corners(:,i)
         cornerPlaneList(1,nUniqueCorners) = cornerPlanes(1,i)
         cornerPlaneList(2,nUniqueCorners) = cornerPlanes(2,i)
         cornerPlaneList(3,nUniqueCorners) = cornerPlanes(3,i)
      ELSE
         DO j = 1, 3
            foundDuplicate = .FALSE.
            DO k = 1, numCornerPlanes(duplicateNum)
               IF (cornerPlaneList(k,duplicateNum).EQ.cornerPlanes(j,i)) THEN
                  foundDuplicate = .TRUE.
               END IF
            END DO
            IF (.NOT.foundDuplicate) THEN
               numCornerPlanes(duplicateNum) = numCornerPlanes(duplicateNum) + 1
               cornerPlaneList(numCornerPlanes(duplicateNum),duplicateNum) = cornerPlanes(j,i)
            END IF
         END DO
      END IF
   END DO
   nCorners = nUniqueCorners

   ! Add origin to corner points
   IF ((nPlanes-nOuterPlanes).GE.3) THEN ! The origin is only a corner if at least 3 planes feature this point
      nCorners = nCorners + 1
      corners(:,nCorners) = 0.0
      DO i = nOuterPlanes + 1, nPlanes
         cornerPlaneList(i-nOuterPlanes,nCorners) = i
      END DO
      numCornerPlanes(nCorners) = nPlanes - nOuterPlanes
   END IF

   ! Filter out "corners" found for sets of planes that do not meet in a single
   ! point but have a common intersection edge.
   ! (This absurd case actually occured so it has to be treated. LAPACK might
   ! find a "corner point" in this situation. Also the origin can be such a point.)
   nUniqueCorners = 0
   DO i = 1, nCorners
      filterOut = .FALSE.
      IF(numCornerPlanes(i).GE.3) THEN
         filterOut = .TRUE.
         vecA(:) = dvec(:,cornerPlaneList(1,i))
         vecB(:) = dvec(:,cornerPlaneList(2,i))
         cornerPlaneLoop: DO n3 = 3, numCornerPlanes(i)
            vecC(:) = dvec(:,cornerPlaneList(n3,i))
            testVec(1) = vecA(2)*vecB(3)-vecA(3)*vecB(2)
            testVec(2) = vecA(3)*vecB(1)-vecA(1)*vecB(3)
            testVec(3) = vecA(1)*vecB(2)-vecA(2)*vecB(1)
            scalarProduct = testVec(1)*vecC(1) + testVec(2)*vecC(2) + testVec(3)*vecC(3)
            IF (ABS(scalarProduct).GT.eps09) THEN
               filterOut = .FALSE.
               EXIT cornerPlaneLoop
            END IF
         END DO cornerPlaneLoop
      END IF
      IF(filterOut) THEN
         CYCLE
      END IF
      nUniqueCorners = nUniqueCorners + 1
      IF(nUniqueCorners.NE.i) THEN
         numCornerPlanes(nUniqueCorners) = numCornerPlanes(i)
         DO j = 1, numCornerPlanes(i)
            cornerPlaneList(j,nUniqueCorners) = cornerPlaneList(j,i)
         END DO
         corners(:,nUniqueCorners) = corners(:,i)
      END IF
   END DO
   nCorners = nUniqueCorners

   ! Count the number of corners for each plane
   nPlaneCorners = 0
   ALLOCATE(planeCorners(nPlanes,nCorners))
   planeCorners = 0
   DO i = 1, nCorners
      DO j = 1, numCornerPlanes(i)
         nPlaneCorners(cornerPlaneList(j,i)) = nPlaneCorners(cornerPlaneList(j,i)) + 1
         planeCorners(cornerPlaneList(j,i),nPlaneCorners(cornerPlaneList(j,i))) = i
      END DO
   END DO

   ! Remove irrelevant planes:
   nface = 0
   isIBZPlane(:) = .TRUE.
   DO n1 = 1, nPlanes
      IF (nPlaneCorners(n1).LE.2) THEN
         isIBZPlane(n1) = .FALSE.
         CYCLE
      END IF
!      WRITE(*,*) 'plane ', n1
!      WRITE(*,'(4f20.13)') dvec(:,n1), ddist(n1)
!      WRITE(*,*) 'corners:'
!      DO i = 1, nPlaneCorners(n1)
!         WRITE(*,'(i5,3f20.13)') planeCorners(n1,i), corners(:,planeCorners(n1,i))
!      END DO
      nface = nface + 1
   END DO
   
   ! Remove irrelevant corners:
   ncorn = 0
   isIBZCorner(:) = .TRUE.
!   WRITE(*,*) 'IBZ corners:'
   DO i = 1, nCorners
      numIBZPlanes = 0
      DO j = 1, numCornerPlanes(i)
         IF(isIBZPlane(cornerPlaneList(j,i))) THEN
            numIBZPlanes = numIBZPlanes + 1
         END IF
      END DO
      IF(numIBZPlanes.LE.2) isIBZCorner(i) = .FALSE.
      IF(.NOT.isIBZCorner(i)) CYCLE
      ncorn = ncorn + 1
!      WRITE(*,'(i5,3f20.13)') i, corners(:,i)
   END DO

   DEALLOCATE(cornerPlaneList)

   ! 4. Construct edges from the set of boundary points. Actually
   !    only the number of edges is needed.

   ! Count number of edges based on number of corners for each face
   ! Each face has the same number of edges as it has corners.
   ! Counting the number of edges in this way counts each edge twice.
   nedge = 0
   DO n1 = 1, nPlanes
      IF(.NOT.isIBZPlane(n1)) CYCLE
      nedge = nedge + nPlaneCorners(n1)
   END DO
   nedge = nedge / 2

   ! 5. Fill the output arrays.

   IF((mface.LT.ncorn).OR.(mface.LT.nface)) THEN
      WRITE(*,*) "mface: ", mface
      WRITE(*,*) "ncorn: ", ncorn
      WRITE(*,*) "nface: ", nface
      WRITE(*,*) "mface has to be larger or equal to nface and ncorn."
      WRITE(*,*) "Its value is hardcoded. Adapt it."
      CALL juDFT_error("mface was chosen too small.",calledby ="brzone2")
   END IF

   cpoint = 0.0
   j = 0
   DO i = 1, nCorners
      IF (.NOT.isIBZCorner(i)) CYCLE
      j = j + 1
      cpoint(:,j) = corners(:,i)
   END DO

   fdist = 0.0
   fnorm = 0.0
   j = 0
   DO i = 1, nPlanes
      IF(.NOT.isIBZPlane(i)) CYCLE
      j = j + 1
      fdist(j) = ddist(i)
      fnorm(:,j) = dvec(:,i)
   END DO

!   WRITE(*,*) 'ncorn', ncorn
!   WRITE(*,*) 'nedge', nedge
!   WRITE(*,*) 'nface', nface
!   WRITE(*,*) 'faces:'
!   DO i  =1,nface
!      WRITE(*,'(4f20.13)') fnorm(:,i), fdist(i)
!   END DO
!   WRITE(*,*) 'coners:'
!   DO i = 1,ncorn
!      WRITE(*,'(3f20.13)') cpoint(:,i)
!   END DO

   ! 5.1. Check Euler characteristic

   IF ((ncorn + nface - nedge).NE.2) THEN
      WRITE(*,*) "ncorn: ", ncorn
      WRITE(*,*) "nface: ", nface
      WRITE(*,*) "nedge: ", nedge
      WRITE(*,*) "corners: "
      DO i = 1, ncorn
         WRITE(*,'(3f20.13)') cpoint(:,i)
      END DO
      WRITE(*,*) "faces: "
      DO i = 1, nface
         WRITE(*,'(4f20.13)') fnorm(:,i), fdist(i)
      END DO
      CALL juDFT_error("Brillouin zone does not fulfill Euler characteristic.",calledby ="brzone2")
   END IF

   DEALLOCATE (planeCorners,dvec,ddist,nPlaneCorners,isIBZPlane)

END SUBROUTINE brzone2

END MODULE m_brzone2
