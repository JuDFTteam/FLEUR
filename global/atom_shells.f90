MODULE m_atom_shells

   !Module for calculating the shells for Jij exchange constants
   !Author: HJ 2021

   USE m_types_atoms
   USE m_types_sym
   USE m_types_cell
   USE m_juDFT
   USE m_constants
   USE m_sort
   USE m_inv3

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: construct_atom_shells

   CONTAINS

   SUBROUTINE construct_atom_shells(referenceAtom, nshells, atoms, cell, sym, film, shellDistances,&
                                    shellDiffs, shellAtoms, shellOps, numAtomsShell, generatedShells)

      !----------------------------------------------------------------------------------------------
      !Construct neighbour shells around a given reference atom
      !Reduces the found shells according to symmetry
      !
      !Takes the reference atom type and the number of shells as input (plus fleur input types)
      !
      ! Returns:
      !    - generatedShells: How many shells were generated (not necessarily at different distances)
      !    - shellDistances: Array of the distances of each shell
      !    - shellDiffs: Array of connecting vectors (lattice coordinates) in the shells
      !    - shellAtoms: Which Atoms (not atom types) correspond to the vectors in the shell
      !    - shellOps: The symmetry operations, which can be used to construct the other elements
      !                in the shell from the first element
      !    - numAtomsShell: How many elements are in each shell
      !----------------------------------------------------------------------------------------------


      INTEGER,             INTENT(IN)   :: referenceAtom !which is the reference atom
      INTEGER,             INTENT(IN)   :: nshells !How many nearest neighbour shells are requested
      TYPE(t_atoms),       INTENT(IN)   :: atoms
      TYPE(t_cell),        INTENT(IN)   :: cell
      TYPE(t_sym),         INTENT(IN)   :: sym
      LOGICAL,             INTENT(IN)   :: film
      REAL, ALLOCATABLE,   INTENT(OUT)  :: shellDistances(:)
      REAL, ALLOCATABLE,   INTENT(OUT)  :: shellDiffs(:,:,:)
      INTEGER, ALLOCATABLE,INTENT(OUT)  :: shellAtoms(:,:,:)
      INTEGER, ALLOCATABLE,INTENT(OUT)  :: shellOps(:,:)
      INTEGER, ALLOCATABLE,INTENT(OUT)  :: numAtomsShell(:)
      INTEGER,             INTENT(OUT)  :: generatedShells

      REAL, PARAMETER :: eps = 1e-5

      INTEGER :: newNeighbours,atomShells,atomShells1,actualShells,num_cells
      INTEGER :: ishell,i,iAtom,atomTypep,refAt
      LOGICAL :: l_unfinished_shell,l_found_shell,l_add
      REAL :: min_distance_to_border(3), distance,distance_to_border, lastDist
      REAL :: leftBorder(3),rightBorder(3),refPos(3)

      INTEGER, ALLOCATABLE :: newAtoms(:,:), distIndexList(:)
      REAL,    ALLOCATABLE :: newDiffs(:,:), newDistances(:)
      INTEGER, ALLOCATABLE :: sortedAtoms(:,:)
      REAL,    ALLOCATABLE :: sortedDiffs(:,:), sortedDistances(:)

      CALL timestart('Atom shells')

      CALL alloc_shells(shellAtoms, shellDiffs, shellDistances, numAtomsShell,&
                        array_size = atoms%nat * atoms%neq(referenceAtom) * 27)

      WRITE(oUnit,'(/,/,A,I0,A,I0)') "Generating ", nshells, " shells for atomType: ", referenceAtom

      atomShells = 0
      actualShells = 0
      num_cells = 0
      l_unfinished_shell = .TRUE.
      DO WHILE(l_unfinished_shell)

         WRITE(oUnit,'(A,I0)') "Number of unit cells in each direction: ", num_cells + 1
         !Calculate the vectors and distances to neighbours in the next
         !extension of unit cells
         CALL calculate_next_neighbours(referenceAtom, atoms, cell%amat, film, newNeighbours, newAtoms,&
                                        newDiffs, newDistances, num_cells)

         WRITE(oUnit,'(A,I0)') "New neighbours found: ", newNeighbours

         CALL timestart('Sorting')
         IF(ALLOCATED(distIndexList)) DEALLOCATE(distIndexList)
         IF(ALLOCATED(sortedDiffs)) DEALLOCATE(sortedDiffs)
         IF(ALLOCATED(sortedDistances)) DEALLOCATE(sortedDistances)
         IF(ALLOCATED(sortedAtoms)) DEALLOCATE(sortedAtoms)
         ALLOCATE(distIndexList(newNeighbours))
         ALLOCATE(sortedDiffs(3,newNeighbours),source=0.0)
         ALLOCATE(sortedDistances(newNeighbours),source=0.0)
         ALLOCATE(sortedAtoms(2,newNeighbours),source=0)
         !Sort the atoms according to distance
         CALL sort(distIndexList(:newNeighbours),newDistances(:newNeighbours))
         DO i = 1, newNeighbours
            sortedDiffs(:,i) = newDiffs(:,distIndexList(i))
            sortedDistances(i) = newDistances(distIndexList(i))
            sortedAtoms(:,i) = newAtoms(:,distIndexList(i))
         END DO
         CALL timestop('Sorting')

         CALL timestart('Grouping Elements into Shells')
         !Sort the nearestNeighbours into shells
         DO iAtom = 1, newNeighbours
            IF(ABS(sortedDistances(iAtom))<1e-12) CYCLE

            l_found_shell = .FALSE.
            !Search for the same shell
            DO ishell = 1, actualShells
               IF(ABS(shellDistances(ishell)-sortedDistances(iAtom)).GT.eps) CYCLE
               atomTypep = atoms%itype(sortedAtoms(2,iAtom))
               IF(atoms%itype(shellAtoms(2,1,ishell))/=atomTypep) CYCLE

               IF(numAtomsShell(ishell)+1>SIZE(shellDistances)) THEN
                  CALL alloc_shells(shellAtoms, shellDiffs, shellDistances, numAtomsShell)
               ENDIF

               numAtomsShell(ishell) = numAtomsShell(ishell) + 1
               shellAtoms(:,numAtomsShell(ishell),ishell) = sortedAtoms(:,iAtom)
               shellDiffs(:,numAtomsShell(ishell),ishell) = sortedDiffs(:,iAtom)
               l_found_shell = .TRUE.
            ENDDO

            !Search for shells with the same distance
            IF(.NOT.l_found_shell) THEN
               l_add = .FALSE.
               DO ishell = 1, actualShells
                  IF(shellDistances(ishell)-sortedDistances(iAtom).GT.eps) THEN
                     !Here we are at the first shell with a larger distance
                     l_add = .TRUE.
                     EXIT
                  ENDIF
               ENDDO

               IF(l_add) THEN
                  CALL insert_new_shell(ishell, actualShells, shellAtoms, shellDiffs, shellDistances, numAtomsShell)

                  IF(.NOT.ANY(ABS(shellDistances-sortedDistances(iAtom))<eps)) atomShells = atomShells + 1
                  numAtomsShell(ishell) = numAtomsShell(ishell) + 1
                  shellAtoms(:,numAtomsShell(ishell),ishell) = sortedAtoms(:,iAtom)
                  shellDistances(ishell) = sortedDistances(iAtom)
                  shellDiffs(:,numAtomsShell(ishell),ishell) = sortedDiffs(:,iAtom)
                  l_found_shell = .TRUE.
               ENDIF
            ENDIF

            !Add shell with new distance
            IF(.NOT.l_found_shell) THEN
               ishell = actualShells + 1
               CALL insert_new_shell(ishell, actualShells, shellAtoms, shellDiffs, shellDistances, numAtomsShell)
               
               atomShells = atomShells + 1
               numAtomsShell(ishell) = numAtomsShell(ishell) + 1
               shellAtoms(:,numAtomsShell(ishell),ishell) = sortedAtoms(:,iAtom)
               shellDistances(ishell) = sortedDistances(iAtom)
               shellDiffs(:,numAtomsShell(ishell),ishell) = sortedDiffs(:,iAtom)
            ENDIF
         ENDDO
         CALL timestop('Grouping Elements into Shells')

         WRITE(oUnit,'(A,I0)') "Shells found: ", atomShells

         IF(atomShells>=nshells) THEN
            CALL timestart('Checking completeness of shell')
            !Calculate how many shells correspond to the requested nshells
            !We look if there can possibly be more elements outside the currently
            !chosen quadrant of cells
            lastDist = 0.0
            atomShells1 = 0
            DO ishell = 1, actualShells
               IF(shellDistances(ishell)-lastDist > eps) atomShells1 = atomShells1 + 1
               lastDist = shellDistances(ishell)
               IF(atomShells1==nshells) THEN
                  distance = SQRT(shellDistances(ishell))

                  min_distance_to_border = 9e99
                  DO refAt = SUM(atoms%neq(:referenceAtom-1)) + 1, SUM(atoms%neq(:referenceAtom))
                     refPos(:) = atoms%taual(:,refAt)

                     DO i = 1, 3
                        !Distance to border in direction of lattice vector
                        leftBorder  = cell%amat(:,i) * (num_cells + refPos(i))
                        rightBorder = cell%amat(:,i) * (num_cells + 1 - refPos(i))
                        distance_to_border = MIN(norm2(leftBorder),norm2(rightBorder))
                        IF(distance_to_border<min_distance_to_border(i)) THEN
                           min_distance_to_border(i) = distance_to_border
                        ENDIF
                     ENDDO
                  ENDDO

                  !-1e-12 to avoid uneccesary calculations where both are equla to numerical precision
                  IF(ALL(min_distance_to_border(:)-distance > -eps)) THEN
                     WRITE(oUnit,'(A)') "Shells finished."
                     l_unfinished_shell = .FALSE.
                  ENDIF
                  EXIT
               ENDIF
            ENDDO
            CALL timestop('Checking completeness of shell')
         ENDIF
      ENDDO
      shellDistances(:actualShells) = SQRT(shellDistances(:actualShells))

      !Calculate how many shells correspond to the requested nshells
      lastDist = 0.0
      atomShells1 = 0
      DO ishell = 1, actualShells
         IF(shellDistances(ishell)-lastDist > eps) atomShells1 = atomShells1 + 1
         lastDist = shellDistances(ishell)
         IF(atomShells1>nshells) THEN
            generatedShells = ishell - 1
            EXIT
         ELSE IF (ishell==actualShells) THEN
            generatedShells = actualShells
            EXIT
         ENDIF
      ENDDO

      IF(ANY(numAtomsShell(:generatedShells)==0)) THEN
         CALL juDFT_error("Empty shells were generated",&
                          hint="This is bug in FLEUR,please report",&
                          calledby="atom_shells")
      ENDIF

      ALLOCATE(shellOps(SIZE(shellDistances),SIZE(shellDistances)), source=0)
      !Symmetry reduction (modernized and modified version of nshell.f from v26)
      CALL apply_sym_to_shell(generatedShells, atoms, sym, shellAtoms, shellDiffs, shellDistances, numAtomsShell, shellOps)

      CALL timestop('Atom shells')

   END SUBROUTINE construct_atom_shells


   SUBROUTINE calculate_next_neighbours(referenceAtom, atoms, amat, film, neighboursFound, neighbourAtoms,&
                                        neighbourDiffs, neighbourDistances, lastBorder)

      !Calculate the distances and vectors to neighbour atoms to a reference atom in a
      !supercell
      !This is done in steps to avoid too much unnecessary calculations
      !The first call will calculate all neighbours in a 3x3x3 supercell and after that
      !only the neighbours in the next layer of supercells is calculated
      !all connecting vectors are returned in lattice coordinates and the distances
      !are not squared, but the root is already taken

      INTEGER,              INTENT(IN)    :: referenceAtom
      TYPE(t_atoms),        INTENT(IN)    :: atoms
      REAL,                 INTENT(IN)    :: amat(:,:)
      LOGICAL,              INTENT(IN)    :: film
      INTEGER,              INTENT(OUT)   :: neighboursFound
      INTEGER, ALLOCATABLE, INTENT(OUT)   :: neighbourAtoms(:,:)
      REAL,    ALLOCATABLE, INTENT(OUT)   :: neighbourDiffs(:,:)
      REAL,    ALLOCATABLE, INTENT(OUT)   :: neighbourDistances(:)
      INTEGER,              INTENT(INOUT) :: lastBorder

      INTEGER :: maxNeighbours,iAtom,refAt,identicalAtoms,i,j,k,n,na,zmax
      REAL :: amatDet, currentDist
      REAL :: tau(3),refPos(3),offsetPos(3),currentDiff(3),pos(3)
      REAL :: invAmat(3,3),posCart(3,atoms%nat)

      CALL timestart('Atom shells: Calculate neighbours')

      CALL inv3(amat,invAmat,amatDet)

      IF(ALLOCATED(neighbourAtoms)) DEALLOCATE(neighbourAtoms)
      IF(ALLOCATED(neighbourDiffs)) DEALLOCATE(neighbourDiffs)
      IF(ALLOCATED(neighbourDistances)) DEALLOCATE(neighbourDistances)

      lastBorder = lastBorder + 1
      IF(lastBorder == 1) THEN
         maxNeighbours = atoms%nat * atoms%neq(referenceAtom) * 27 !For the first one we calculate two in one go
      ELSE
         maxNeighbours = atoms%nat * atoms%neq(referenceAtom) * ((2*lastBorder+1)**3 - (2*lastBorder-1)**3)
      ENDIF

      zmax = lastBorder
      IF(film) zmax = 0

      ALLOCATE(neighbourAtoms(2,maxNeighbours), source=0)
      ALLOCATE(neighbourDiffs(3,maxNeighbours), source=0.0)
      ALLOCATE(neighbourDistances(maxNeighbours), source=0.0)

      !Calculate Atom positions in cartesian coordinates
      DO iAtom = 1, atoms%nat
         tau = atoms%taual(:,iAtom) - FLOOR(atoms%taual(:,iAtom))
         posCart(:,iAtom) = MATMUL(amat,tau)
      END DO

      neighboursFound = 0
      DO refAt = SUM(atoms%neq(:referenceAtom-1)) + 1, SUM(atoms%neq(:referenceAtom))
         refPos(:) = posCart(:,refAt)
         identicalAtoms = 0
         DO i = -lastBorder, lastBorder
            DO j = -lastBorder, lastBorder
               DO k = -zmax, zmax

                  IF(ALL(ABS([i,j,k]) < lastBorder).AND.lastBorder/=1) CYCLE

                  offsetPos = matmul(amat, [i,j,k])

                  iAtom = 0
                  DO n = 1, atoms%ntype
                     DO na = 1, atoms%neq(n)
                        iAtom = iAtom + 1
                        pos(:) = posCart(:,iAtom) + offsetPos(:)
                        currentDist = (refPos(1) - pos(1))**2 + &
                                      (refPos(2) - pos(2))**2 + &
                                      (refPos(3) - pos(3))**2
                        currentDiff = refPos(:) - pos(:)
                        IF (currentDist.LT.0.000001) THEN
                           identicalAtoms = identicalAtoms + 1
                        ELSE
                           neighboursFound = neighboursFound + 1
                           neighbourAtoms(1,neighboursFound) = refAt
                           neighbourAtoms(2,neighboursFound) = iAtom
                           neighbourDiffs(:,neighboursFound) = MATMUL(invAmat,currentDiff(:))
                           neighbourDistances(neighboursFound) = currentDist
                        END IF
                     ENDDO
                  END DO
               END DO
            END DO
         END DO
         IF (identicalAtoms.GT.1) THEN
            WRITE(*,*) 'Position: ', refPos(:)
            CALL juDFT_error("Too many atoms at same position.",calledby ="calculate_next_neighbours")
         END IF
      ENDDO

      CALL timestop('Atom shells: Calculate neighbours')

   END SUBROUTINE calculate_next_neighbours

   SUBROUTINE apply_sym_to_shell(actualShells, atoms, sym, shellAtoms, shellDiffs, shellDistances, numAtomsShell, shellOps)

      !Reduce the calculated shells according to the given symmetry
      !Each shell is attempted to be fully recovered by applying symmerty operations
      !to the first element in the shell
      !If it cannot be fully constructed it is split up

      INTEGER,              INTENT(INOUT) :: actualShells
      TYPE(t_atoms),        INTENT(IN)    :: atoms
      TYPE(t_sym),          INTENT(IN)    :: sym
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: shellAtoms(:,:,:)
      REAL,    ALLOCATABLE, INTENT(INOUT) :: shellDiffs(:,:,:)
      REAL,    ALLOCATABLE, INTENT(INOUT) :: shellDistances(:)
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: numAtomsShell(:)
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: shellOps(:,:)

      REAL,    PARAMETER :: sym_tol = 1e-5

      INTEGER :: current_shell,refAtom,refAtomp,ishellAtom,nshellAtom,nshellAtom1
      INTEGER :: iop, iAtom,found_index,atom_alpha,atom_beta
      LOGICAL :: l_diff_in_shell

      REAL :: diffRot(3)
      REAL :: shellDiffAux(3,SIZE(shellDistances))
      REAL :: shellDiffAux1(3,SIZE(shellDistances))
      INTEGER :: shellOpsAux(SIZE(shellDistances))
      INTEGER :: shellAtomsAux(2,SIZE(shellDistances))
      INTEGER :: shellAtomsAux1(2,SIZE(shellDistances))

      CALL timestart('Atom shells: Apply Symmetries')

      current_shell = 1
      DO WHILE(current_shell<=actualShells)

         !Take the representative element of the shell
         shellDiffAux = 0.0
         refAtom = SUM(atoms%neq(:atoms%itype(shellAtoms(1,1,current_shell))-1)) + 1
         refAtomp = SUM(atoms%neq(:atoms%itype(shellAtoms(2,1,current_shell))-1)) + 1
         DO ishellAtom = 1, numAtomsShell(current_shell)
            IF(shellAtoms(1,ishellAtom,current_shell) == refAtom.AND.&
               shellAtoms(2,ishellAtom,current_shell) == refAtomp) THEN
               shellDiffAux(:,1) = shellDiffs(:,ishellAtom,current_shell)
               shellOpsAux(1) = 1 !Identity operation
               shellAtomsAux(:,1) = shellAtoms(:,ishellAtom,current_shell)
               EXIT
            ELSE IF(ishellAtom == numAtomsShell(current_shell)) THEN
               shellDiffAux(:,1) = shellDiffs(:,1,current_shell)
               shellOpsAux(1) = 1 !Identity operation
               shellAtomsAux(:,1) = shellAtoms(:,1,current_shell)
               !CALL juDFT_error("No representative element found", calledby="addNearestNeighbours_gfelem")
            ENDIF
         ENDDO

         nshellAtom = 1
         symLoop: DO iop = 1, sym%nop
            diffRot = matmul(sym%mrot(:,:,iop),shellDiffAux(:,1))
            atom_alpha = sym%mapped_atom(iop,shellAtomsAux(1,1))
            atom_beta = sym%mapped_atom(iop,shellAtomsAux(2,1))

            IF(atom_alpha == 0 .OR. atom_beta == 0) CALL juDFT_error("No mapped atom available",&
                                                                     hint="This is a bug in FLEUR, please report",&
                                                                     calledby='apply_sym_to_shell')

            DO ishellAtom = 1, nshellAtom
               !Is the atom equivalent to another atom already in the shell
               IF(ALL(ABS(diffRot-shellDiffAux(:,ishellAtom)).LT.sym_tol).AND.&
                  atom_alpha == shellAtomsAux(1,ishellAtom).AND.&
                  atom_beta == shellAtomsAux(2,ishellAtom)) CYCLE symLoop
            ENDDO

            l_diff_in_shell = .FALSE.
            DO ishellAtom = 1, numAtomsShell(current_shell)
               !Is the atom equivalent to any other atom in the previously found shell
               IF(ALL(ABS(diffRot-shellDiffs(:,ishellAtom,current_shell)).LT.sym_tol).AND.&
                  atom_alpha == shellAtoms(1,ishellAtom,current_shell).AND.&
                  atom_beta == shellAtoms(2,ishellAtom,current_shell)) THEN
                  l_diff_in_shell = .TRUE.
                  found_index = ishellAtom
                  EXIT
               ENDIF
            ENDDO
            IF(.NOT.l_diff_in_shell) CYCLE symLoop

            nshellAtom = nshellAtom + 1
            shellDiffAux(:,nshellAtom) = diffRot
            shellOpsAux(nshellAtom) = iop
            shellAtomsAux(:,nshellAtom) = shellAtoms(:,found_index,current_shell)
         ENDDO symLoop

         IF(nshellAtom < numAtomsShell(current_shell)) THEN
            !Not all elements can be constructed from the representative element

            shellDiffAux1 = 0.0
            shellAtomsAux1 = 0
            !Find the atoms which are not represented
            nshellAtom1 = 0
            atomLoop: DO iAtom = 1, numAtomsShell(current_shell)
               DO ishellAtom = 1, nshellAtom
                  IF(ALL(ABS(shellDiffs(:,iAtom,current_shell)-shellDiffAux(:,ishellAtom)).LT.sym_tol).AND.&
                     shellAtomsAux(1,ishellAtom) == shellAtoms(1,iAtom,current_shell).AND.&
                     shellAtomsAux(2,ishellAtom) == shellAtoms(2,iAtom,current_shell)) CYCLE atomLoop!.OR.&
                     !ALL(ABS(shellDiffs(:,iAtom,current_shell)+shellAux(:,ishellAtom)).LT.sym_tol)) CYCLE atomLoop
               ENDDO

               nshellAtom1 = nshellAtom1 + 1
               shellDiffAux1(:,nshellAtom1) = shellDiffs(:,iAtom,current_shell)
               shellAtomsAux1(:,nshellAtom1) = shellAtoms(:,iAtom,current_shell)
            ENDDO atomLoop

            CALL insert_new_shell(current_shell+1, actualShells, shellAtoms, shellDiffs, shellDistances,&
                                  numAtomsShell, shellOps=shellOps)

            !Modify current_shell (fewer atoms)
            numAtomsShell(current_shell) = nshellAtom
            DO ishellAtom = 1, nshellAtom
               shellDiffs(:,ishellAtom,current_shell) = shellDiffAux(:,ishellAtom)
               shellOps(ishellAtom,current_shell) = shellOpsAux(ishellAtom)
               shellAtoms(:,ishellAtom,current_shell) = shellAtomsAux(:,ishellAtom)
            ENDDO

            !Insert Element at current_shell+1 (This way it will be the next element in the
            !loop if it needs to be deconstructed further)
            shellDistances(current_shell+1) = shellDistances(current_shell)
            numAtomsShell(current_shell+1) = nshellAtom1
            DO ishellAtom = 1, nshellAtom1
               shellDiffs(:,ishellAtom,current_shell+1) = shellDiffAux1(:,ishellAtom)
               shellAtoms(:,ishellAtom,current_shell+1) = shellAtomsAux1(:,ishellAtom)
            ENDDO

         ELSE
            !The shell could be fully constructed from the first element
            shellOps(:,current_shell) = shellOpsAux(:)
            shellDiffs(:,:,current_shell) = shellDiffAux(:,:)
            shellAtoms(:,:,current_shell) = shellAtomsAux(:,:)
         ENDIF
         current_shell = current_shell + 1

      ENDDO

      CALL timestop('Atom shells: Apply Symmetries')

   END SUBROUTINE apply_sym_to_shell

   SUBROUTINE insert_new_shell(index, actualShells, shellAtoms, shellDiffs, shellDistances, numAtomsShell, shellOps)

      !Insert a shell at the given index (either appending or in the middle)
      !Takes care of reallocation if the arrays are too short

      INTEGER,              INTENT(IN)    :: index
      INTEGER,              INTENT(INOUT) :: actualShells
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: shellAtoms(:,:,:)
      REAL,    ALLOCATABLE, INTENT(INOUT) :: shellDiffs(:,:,:)
      REAL,    ALLOCATABLE, INTENT(INOUT) :: shellDistances(:)
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: numAtomsShell(:)

      INTEGER, ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: shellOps(:,:)

      INTEGER :: ishell

      IF(actualShells+1>SIZE(shellDistances)) THEN
         CALL alloc_shells(shellAtoms, shellDiffs, shellDistances, numAtomsShell, shellOps=shellOps)
      ENDIF
      IF(index>actualShells+1) THEN
         CALL juDFT_error('Too large index for inserting new shell',&
                          hint='This is a Bug in FLEUR , please report',&
                          calledby='insert_new_shell')
      ENDIF

      actualShells = actualShells + 1
      !Move everything above index one index up
      DO ishell = actualShells, index + 1, -1
         shellAtoms(:,:,ishell)  = shellAtoms(:,:,ishell-1)
         shellDistances(ishell) = shellDistances(ishell-1)
         numAtomsShell(ishell) = numAtomsShell(ishell-1)
         shellDiffs(:,:,ishell) = shellDiffs(:,:,ishell-1)
         IF(PRESENT(shellOps)) shellOps(:,ishell) = shellOps(:,ishell-1)
      ENDDO

      !Reset the shell index
      shellAtoms(:,:,index) = 0
      shellDistances(index) = 0.0
      numAtomsShell(index) = 0
      shellDiffs(:,:,index) = 0.0
      IF(PRESENT(shellOps)) shellOps(:,index) = 0

   END SUBROUTINE insert_new_shell

   SUBROUTINE alloc_shells(shellAtoms, shellDiffs, shellDistances, numAtomsShell, shellOps, array_size)

      !Allocate more space in the shell arrays (Or allocate initial sizes)
      !The size is always doubled

      INTEGER, OPTIONAL, INTENT(IN) :: array_size

      INTEGER, ALLOCATABLE, INTENT(INOUT) :: shellAtoms(:,:,:)
      REAL,    ALLOCATABLE, INTENT(INOUT) :: shellDiffs(:,:,:)
      REAL,    ALLOCATABLE, INTENT(INOUT) :: shellDistances(:)
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: numAtomsShell(:)

      INTEGER, ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: shellOps(:,:)

      INTEGER :: previousShells

      INTEGER, ALLOCATABLE :: int_tmp1(:), int_tmp2(:,:),int_tmp3(:,:,:)
      REAL, ALLOCATABLE :: real_tmp1(:), real_tmp3(:,:,:)

      IF(ALLOCATED(shellDistances)) THEN
         previousShells = SIZE(shellDistances)
      ELSE IF(PRESENT(array_size)) THEN
         ALLOCATE(shellAtoms(2,array_size,array_size), source = 0)
         ALLOCATE(shellDiffs(3,array_size,array_size), source = 0.0)
         ALLOCATE(shellDistances(array_size), source = 0.0)
         ALLOCATE(numAtomsShell(array_size), source=0)
         RETURN
      ELSE
         CALL juDFT_error("If the shell arrays are not allocated the array_size argument is needed",&
                          hint='This is a bug in FLEUR, please report', calledby='alloc_shells')
      ENDIF

      ALLOCATE(int_tmp1(2*previousShells), source = 0)
      int_tmp1(:previousShells) = numAtomsShell(:)
      CALL move_alloc(int_tmp1,numAtomsShell)

      ALLOCATE(real_tmp1(2*previousShells), source = 0.0)
      real_tmp1(:previousShells) = shellDistances(:)
      CALL move_alloc(real_tmp1,shellDistances)

      ALLOCATE(real_tmp3(3,2*previousShells,2*previousShells), source = 0.0)
      real_tmp3(:,:previousShells,:previousShells) = shellDiffs(:,:,:)
      CALL move_alloc(real_tmp3,shellDiffs)

      ALLOCATE(int_tmp3(2,2*previousShells,2*previousShells), source = 0)
      int_tmp3(:,:previousShells,:previousShells) = shellAtoms(:,:,:)
      CALL move_alloc(int_tmp3,shellAtoms)

      IF(PRESENT(shellOps)) THEN
         ALLOCATE(int_tmp2(2*previousShells,2*previousShells), source = 0)
         int_tmp2(:previousShells,:previousShells) = shellOps(:,:)
         CALL move_alloc(int_tmp2,shellOps)
      ENDIF

   END SUBROUTINE alloc_shells

END MODULE m_atom_shells