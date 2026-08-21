!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_local_coordination
   USE m_local_coordination_core, ONLY: resolve_periodic_neighbours_cartesian, &
                                        t_coordination_diagnostics, t_periodic_neighbor
   USE m_types_atoms, ONLY: t_atoms
   USE m_types_cell, ONLY: t_cell
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: resolve_fleur_periodic_neighbours
   PUBLIC :: select_fleur_physical_atoms_by_type
   PUBLIC :: t_coordination_diagnostics, t_periodic_neighbor

CONTAINS

   SUBROUTINE resolve_fleur_periodic_neighbours(cell, atoms, film, central_atom, candidate_atoms, nneighbors, &
                                                neighbors, diagnostics, success, message, maximum_image_range)
      ! central_atom and candidate_atoms are physical atom indices in
      ! atoms%pos, not atom-type indices. Equivalent atoms therefore remain
      ! distinct sites. A periodic image is identified by
      ! (neighbor%atom_id,neighbor%translation).
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_atoms), INTENT(IN) :: atoms
      LOGICAL, INTENT(IN) :: film
      INTEGER, INTENT(IN) :: central_atom, candidate_atoms(:), nneighbors
      TYPE(t_periodic_neighbor), ALLOCATABLE, INTENT(OUT) :: neighbors(:)
      TYPE(t_coordination_diagnostics), INTENT(OUT) :: diagnostics
      LOGICAL, INTENT(OUT) :: success
      CHARACTER(LEN=*), INTENT(OUT) :: message
      INTEGER, OPTIONAL, INTENT(IN) :: maximum_image_range

      REAL, ALLOCATABLE :: candidate_positions(:,:)
      LOGICAL :: periodic(3)

      success = .FALSE.
      message = ""
      diagnostics = t_coordination_diagnostics()
      ALLOCATE(neighbors(0))
      IF (.NOT.ALLOCATED(atoms%pos)) THEN
         message = "atoms%pos is not allocated"
         RETURN
      END IF
      IF (SIZE(atoms%pos,1) /= 3 .OR. SIZE(atoms%pos,2) /= atoms%nat) THEN
         message = "atoms%pos is inconsistent with atoms%nat"
         RETURN
      END IF
      IF (central_atom < 1 .OR. central_atom > atoms%nat) THEN
         message = "central physical-atom index is out of range"
         RETURN
      END IF
      IF (ANY(candidate_atoms < 1) .OR. ANY(candidate_atoms > atoms%nat)) THEN
         message = "candidate physical-atom index is out of range"
         RETURN
      END IF

      ALLOCATE(candidate_positions(3,SIZE(candidate_atoms)))
      candidate_positions = atoms%pos(:,candidate_atoms)
      periodic = [.TRUE.,.TRUE.,.NOT.film]
      DEALLOCATE(neighbors)
      CALL resolve_periodic_neighbours_cartesian(cell%amat,periodic,atoms%pos(:,central_atom), &
                                                 candidate_positions,candidate_atoms,nneighbors,neighbors, &
                                                 diagnostics,success,message,maximum_image_range)
   END SUBROUTINE resolve_fleur_periodic_neighbours

   SUBROUTINE select_fleur_physical_atoms_by_type(atoms, allowed_types, physical_atoms, success, message)
      ! Expand explicitly selected FLEUR atom-group/type indices to the physical
      ! atoms already stored in atoms%pos/taual. Output is canonically ordered
      ! by physical atom index and contains no duplicates.
      TYPE(t_atoms), INTENT(IN) :: atoms
      INTEGER, INTENT(IN) :: allowed_types(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: physical_atoms(:)
      LOGICAL, INTENT(OUT) :: success
      CHARACTER(LEN=*), INTENT(OUT) :: message

      INTEGER :: atom, count_atoms

      success = .FALSE.
      message = ""
      ALLOCATE(physical_atoms(0))
      IF (.NOT.ALLOCATED(atoms%itype) .OR. SIZE(atoms%itype) /= atoms%nat) THEN
         message = "atoms%itype is not initialized consistently with atoms%nat"
         RETURN
      END IF
      IF (SIZE(allowed_types) == 0) THEN
         message = "allowed candidate atom-type set is empty"
         RETURN
      END IF
      IF (ANY(allowed_types < 1) .OR. ANY(allowed_types > atoms%ntype)) THEN
         message = "allowed candidate atom-type index is out of range"
         RETURN
      END IF

      count_atoms = COUNT([(ANY(atoms%itype(atom) == allowed_types),atom=1,atoms%nat)])
      DEALLOCATE(physical_atoms)
      ALLOCATE(physical_atoms(count_atoms))
      count_atoms = 0
      DO atom = 1, atoms%nat
         IF (.NOT.ANY(atoms%itype(atom) == allowed_types)) CYCLE
         count_atoms = count_atoms+1
         physical_atoms(count_atoms) = atom
      END DO
      success = .TRUE.
      message = "OK"
   END SUBROUTINE select_fleur_physical_atoms_by_type

END MODULE m_local_coordination
