!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_local_structural_frame
   IMPLICIT NONE
   PRIVATE

   INTEGER, PARAMETER :: n_ligands = 6
   INTEGER, PARAMETER :: n_axes = 3
   REAL, PARAMETER :: pi = 4.0*ATAN(1.0)

   TYPE, PUBLIC :: t_structural_frame_diagnostics
      REAL :: bond_lengths(n_ligands) = 0.0
      INTEGER :: opposite_pairs(2, n_axes) = 0
      REAL :: opposite_pair_angles_deg(n_axes) = 0.0
      REAL :: opposite_pair_costs(n_axes) = 0.0
      REAL :: total_pair_cost = 0.0
      REAL :: raw_axes(3, n_axes) = 0.0
      REAL :: raw_mutual_angles_deg(n_axes) = 0.0
      REAL :: raw_determinant = 0.0
      REAL :: singular_values(n_axes) = 0.0
      REAL :: condition_number = 0.0
      REAL :: raw_to_orthonormal_angles_deg(n_axes) = 0.0
      REAL :: reference_alignment_score = 0.0
      REAL :: reference_alignment_gap = 0.0
      LOGICAL :: svd_reflection_corrected = .FALSE.
   END TYPE t_structural_frame_diagnostics

   PUBLIC :: construct_structural_frame_from_displacements
   PUBLIC :: construct_structural_frame_from_positions
   PUBLIC :: closest_proper_rotation

CONTAINS

   SUBROUTINE construct_structural_frame_from_positions(center, ligand_positions, reference_frame, &
                                                        local_to_global, diagnostics, success, message)
      ! The caller supplies the physically correct ligand images. Displacements
      ! are defined strictly as ligand_position-center; no PBC search is done.
      REAL, INTENT(IN) :: center(3), ligand_positions(3, n_ligands), reference_frame(3, 3)
      REAL, INTENT(OUT) :: local_to_global(3, 3)
      TYPE(t_structural_frame_diagnostics), INTENT(OUT) :: diagnostics
      LOGICAL, INTENT(OUT) :: success
      CHARACTER(LEN=*), INTENT(OUT) :: message

      REAL :: displacements(3, n_ligands)
      INTEGER :: ligand

      DO ligand = 1, n_ligands
         displacements(:, ligand) = ligand_positions(:, ligand) - center
      END DO
      CALL construct_structural_frame_from_displacements(displacements, reference_frame, local_to_global, &
                                                         diagnostics, success, message)
   END SUBROUTINE construct_structural_frame_from_positions

   SUBROUTINE construct_structural_frame_from_displacements(displacements, reference_frame, &
                                                            local_to_global, diagnostics, success, message)
      ! local_to_global stores x_local,y_local,z_local as columns in global
      ! Cartesian coordinates: v_global=local_to_global*v_local.
      REAL, INTENT(IN) :: displacements(3, n_ligands), reference_frame(3, 3)
      REAL, INTENT(OUT) :: local_to_global(3, 3)
      TYPE(t_structural_frame_diagnostics), INTENT(OUT) :: diagnostics
      LOGICAL, INTENT(OUT) :: success
      CHARACTER(LEN=*), INTENT(OUT) :: message

      INTEGER :: pairing_table(2, n_axes, 15)
      INTEGER :: raw_pairs(2, n_axes), assigned_pairs(2, n_axes)
      INTEGER :: pairing, axis, i, j, best_pairing
      INTEGER :: permutation(n_axes), signs(n_axes), best_permutation(n_axes), best_signs(n_axes)
      REAL :: directions(3, n_ligands), raw_unassigned(3, n_axes), assigned_axes(3, n_axes)
      REAL :: candidate_axes(3, n_axes), pair_costs(n_axes), candidate_cost
      REAL :: best_cost, candidate_score, candidate_gap, best_score, best_gap
      REAL :: norm_value, scale, determinant_reference
      LOGICAL :: assignment_success, rotation_success
      CHARACTER(LEN=256) :: rotation_message

      diagnostics = t_structural_frame_diagnostics()
      local_to_global = 0.0
      success = .FALSE.
      message = ""

      determinant_reference = determinant_3x3(reference_frame)
      IF (MAXVAL(ABS(MATMUL(TRANSPOSE(reference_frame), reference_frame) - identity_3())) > 1.0e-8 .OR. &
          determinant_reference <= 0.0) THEN
         message = "Reference frame must be orthonormal and right handed"
         RETURN
      END IF

      DO i = 1, n_ligands
         diagnostics%bond_lengths(i) = vector_norm(displacements(:, i))
      END DO
      scale = MAXVAL(diagnostics%bond_lengths)
      IF (scale <= 0.0 .OR. MINVAL(diagnostics%bond_lengths) <= 100.0*EPSILON(scale)*scale) THEN
         message = "All six ligand displacement vectors must be nonzero"
         RETURN
      END IF
      DO i = 1, n_ligands
         directions(:, i) = displacements(:, i)/diagnostics%bond_lengths(i)
      END DO

      CALL initialize_pairing_table(pairing_table)
      best_cost = HUGE(1.0)
      best_score = -HUGE(1.0)
      best_gap = 0.0
      best_pairing = 0
      assigned_axes = 0.0
      DO pairing = 1, SIZE(pairing_table, 3)
         raw_pairs = pairing_table(:, :, pairing)
         candidate_cost = 0.0
         DO axis = 1, n_axes
            i = raw_pairs(1, axis)
            j = raw_pairs(2, axis)
            pair_costs(axis) = 1.0 + DOT_PRODUCT(directions(:, i), directions(:, j))
            candidate_cost = candidate_cost + pair_costs(axis)
            raw_unassigned(:, axis) = directions(:, i) - directions(:, j)
            norm_value = vector_norm(raw_unassigned(:, axis))
            IF (norm_value <= 1.0e-10) EXIT
            raw_unassigned(:, axis) = raw_unassigned(:, axis)/norm_value
         END DO
         IF (axis <= n_axes) CYCLE

         CALL choose_reference_assignment(raw_unassigned, reference_frame, candidate_axes, permutation, signs, &
                                          candidate_score, candidate_gap, assignment_success)
         IF (.NOT. assignment_success) CYCLE
         IF (candidate_cost < best_cost - 1.0e-12 .OR. &
             (ABS(candidate_cost - best_cost) <= 1.0e-12 .AND. candidate_score > best_score + 1.0e-12) .OR. &
             (ABS(candidate_cost - best_cost) <= 1.0e-12 .AND. ABS(candidate_score - best_score) <= 1.0e-12 &
              .AND. lexicographically_precedes(candidate_axes, assigned_axes))) THEN
            best_cost = candidate_cost
            best_score = candidate_score
            best_gap = candidate_gap
            best_pairing = pairing
            assigned_axes = candidate_axes
            best_permutation = permutation
            best_signs = signs
         END IF
      END DO

      IF (best_pairing == 0) THEN
         message = "Ligand directions do not define three independent opposite-pair axes"
         RETURN
      END IF

      raw_pairs = pairing_table(:, :, best_pairing)
      DO axis = 1, n_axes
         assigned_pairs(:, axis) = raw_pairs(:, best_permutation(axis))
         IF (best_signs(axis) < 0) assigned_pairs(:, axis) = assigned_pairs(2:1:-1, axis)
      END DO
      diagnostics%opposite_pairs = assigned_pairs
      diagnostics%raw_axes = assigned_axes
      diagnostics%total_pair_cost = best_cost
      diagnostics%reference_alignment_score = best_score
      diagnostics%reference_alignment_gap = best_gap
      DO axis = 1, n_axes
         i = assigned_pairs(1, axis)
         j = assigned_pairs(2, axis)
         diagnostics%opposite_pair_costs(axis) = 1.0 + DOT_PRODUCT(directions(:, i), directions(:, j))
         diagnostics%opposite_pair_angles_deg(axis) = angle_degrees(directions(:, i), directions(:, j))
      END DO
      diagnostics%raw_mutual_angles_deg(1) = angle_degrees(assigned_axes(:, 1), assigned_axes(:, 2))
      diagnostics%raw_mutual_angles_deg(2) = angle_degrees(assigned_axes(:, 1), assigned_axes(:, 3))
      diagnostics%raw_mutual_angles_deg(3) = angle_degrees(assigned_axes(:, 2), assigned_axes(:, 3))
      diagnostics%raw_determinant = determinant_3x3(assigned_axes)

      CALL closest_proper_rotation(assigned_axes, local_to_global, diagnostics%singular_values, &
                                   diagnostics%condition_number, diagnostics%svd_reflection_corrected, &
                                   rotation_success, rotation_message)
      IF (.NOT. rotation_success) THEN
         message = "Pathological octahedral geometry: "//TRIM(rotation_message)
         RETURN
      END IF
      DO axis = 1, n_axes
         diagnostics%raw_to_orthonormal_angles_deg(axis) = &
            angle_degrees(assigned_axes(:, axis), local_to_global(:, axis))
      END DO
      success = .TRUE.
      message = "OK"
   END SUBROUTINE construct_structural_frame_from_displacements

   SUBROUTINE closest_proper_rotation(raw_axes, rotation, singular_values, condition_number, &
                                      reflection_corrected, success, message)
      ! Orthogonal Procrustes/polar factor from A=U Sigma V^T. If U V^T is
      ! improper, the singular vector associated with the smallest singular
      ! value is flipped, giving the closest proper rotation.
      REAL, INTENT(IN) :: raw_axes(3, 3)
      REAL, INTENT(OUT) :: rotation(3, 3), singular_values(3), condition_number
      LOGICAL, INTENT(OUT) :: reflection_corrected, success
      CHARACTER(LEN=*), INTENT(OUT) :: message

      REAL :: a_work(3, 3), u(3, 3), vt(3, 3), work(64)
      REAL :: improper_rotation(3, 3), minimum_projector(3, 3), reflection(3, 3)
      REAL :: reflection_normal(3), norm_value, degeneracy_tolerance
      INTEGER :: axis, info, singular_index

      a_work = raw_axes
      rotation = 0.0
      singular_values = 0.0
      condition_number = HUGE(1.0)
      reflection_corrected = .FALSE.
      success = .FALSE.
      message = ""
      CALL DGESVD('A', 'A', 3, 3, a_work, 3, singular_values, u, 3, vt, 3, work, SIZE(work), info)
      IF (info /= 0) THEN
         WRITE (message, '(a,i0)') "SVD failed with INFO=", info
         RETURN
      END IF
      IF (singular_values(1) <= 0.0 .OR. singular_values(3) <= 1.0e-8*singular_values(1)) THEN
         WRITE (message, '(a,3es13.5)') "raw axes are nearly linearly dependent; singular values=", singular_values
         RETURN
      END IF
      condition_number = singular_values(1)/singular_values(3)
      rotation = MATMUL(u, vt)
      IF (determinant_3x3(rotation) < 0.0) THEN
         ! R=Q(I-2vv^T), where v belongs to the smallest right-singular
         ! subspace, is the closest proper factor. If that singular value is
         ! degenerate, select v deterministically by projecting the first
         ! available global Cartesian basis vector into the full degenerate
         ! subspace. This avoids depending on LAPACK's basis inside it.
         improper_rotation = rotation
         minimum_projector = 0.0
         degeneracy_tolerance = 1000.0*EPSILON(singular_values(1))*singular_values(1)
         DO singular_index = 1, 3
            IF (ABS(singular_values(singular_index) - singular_values(3)) <= degeneracy_tolerance) THEN
               minimum_projector = minimum_projector &
                  + outer_product(vt(singular_index, :), vt(singular_index, :))
            END IF
         END DO
         reflection_normal = 0.0
         DO axis = 1, 3
            reflection_normal = minimum_projector(:, axis)
            norm_value = vector_norm(reflection_normal)
            IF (norm_value > 100.0*EPSILON(norm_value)) EXIT
         END DO
         reflection_normal = reflection_normal/norm_value
         reflection = identity_3() - 2.0*outer_product(reflection_normal, reflection_normal)
         rotation = MATMUL(improper_rotation, reflection)
         reflection_corrected = .TRUE.
      END IF
      IF (determinant_3x3(rotation) <= 0.0) THEN
         message = "SVD reflection correction did not produce a proper rotation"
         RETURN
      END IF
      success = .TRUE.
      message = "OK"
   END SUBROUTINE closest_proper_rotation

   SUBROUTINE choose_reference_assignment(raw_axes, reference_frame, selected_axes, selected_permutation, &
                                          selected_signs, best_score, score_gap, success)
      REAL, INTENT(IN) :: raw_axes(3, 3), reference_frame(3, 3)
      REAL, INTENT(OUT) :: selected_axes(3, 3), best_score, score_gap
      INTEGER, INTENT(OUT) :: selected_permutation(3), selected_signs(3)
      LOGICAL, INTENT(OUT) :: success

      INTEGER :: permutations(3, 6), permutation, sx, sy, sz, signs(3), axis
      REAL :: candidate(3, 3), score, second_score

      CALL initialize_permutations(permutations)
      selected_axes = 0.0
      selected_permutation = 0
      selected_signs = 0
      best_score = -HUGE(1.0)
      second_score = -HUGE(1.0)
      success = .FALSE.
      DO permutation = 1, 6
         DO sx = -1, 1, 2
            DO sy = -1, 1, 2
               DO sz = -1, 1, 2
                  signs = [sx, sy, sz]
                  DO axis = 1, 3
                     candidate(:, axis) = REAL(signs(axis))*raw_axes(:, permutations(axis, permutation))
                  END DO
                  IF (determinant_3x3(candidate) <= 1.0e-10) CYCLE
                  score = SUM(candidate*reference_frame)
                  IF (score > best_score + 1.0e-12 .OR. &
                      (ABS(score - best_score) <= 1.0e-12 .AND. &
                       lexicographically_precedes(candidate, selected_axes))) THEN
                     second_score = best_score
                     best_score = score
                     selected_axes = candidate
                     selected_permutation = permutations(:, permutation)
                     selected_signs = signs
                     success = .TRUE.
                  ELSE IF (score > second_score) THEN
                     second_score = score
                  END IF
               END DO
            END DO
         END DO
      END DO
      IF (success) THEN
         score_gap = best_score - second_score
      ELSE
         score_gap = 0.0
      END IF
   END SUBROUTINE choose_reference_assignment

   SUBROUTINE initialize_pairing_table(pairings)
      INTEGER, INTENT(OUT) :: pairings(2, 3, 15)

      pairings(:, :,  1) = RESHAPE([1,2, 3,4, 5,6], [2,3])
      pairings(:, :,  2) = RESHAPE([1,2, 3,5, 4,6], [2,3])
      pairings(:, :,  3) = RESHAPE([1,2, 3,6, 4,5], [2,3])
      pairings(:, :,  4) = RESHAPE([1,3, 2,4, 5,6], [2,3])
      pairings(:, :,  5) = RESHAPE([1,3, 2,5, 4,6], [2,3])
      pairings(:, :,  6) = RESHAPE([1,3, 2,6, 4,5], [2,3])
      pairings(:, :,  7) = RESHAPE([1,4, 2,3, 5,6], [2,3])
      pairings(:, :,  8) = RESHAPE([1,4, 2,5, 3,6], [2,3])
      pairings(:, :,  9) = RESHAPE([1,4, 2,6, 3,5], [2,3])
      pairings(:, :, 10) = RESHAPE([1,5, 2,3, 4,6], [2,3])
      pairings(:, :, 11) = RESHAPE([1,5, 2,4, 3,6], [2,3])
      pairings(:, :, 12) = RESHAPE([1,5, 2,6, 3,4], [2,3])
      pairings(:, :, 13) = RESHAPE([1,6, 2,3, 4,5], [2,3])
      pairings(:, :, 14) = RESHAPE([1,6, 2,4, 3,5], [2,3])
      pairings(:, :, 15) = RESHAPE([1,6, 2,5, 3,4], [2,3])
   END SUBROUTINE initialize_pairing_table

   SUBROUTINE initialize_permutations(permutations)
      INTEGER, INTENT(OUT) :: permutations(3, 6)

      permutations(:, 1) = [1, 2, 3]
      permutations(:, 2) = [1, 3, 2]
      permutations(:, 3) = [2, 1, 3]
      permutations(:, 4) = [2, 3, 1]
      permutations(:, 5) = [3, 1, 2]
      permutations(:, 6) = [3, 2, 1]
   END SUBROUTINE initialize_permutations

   LOGICAL FUNCTION lexicographically_precedes(left, right)
      REAL, INTENT(IN) :: left(3, 3), right(3, 3)
      INTEGER :: column, row

      lexicographically_precedes = .FALSE.
      DO column = 1, 3
         DO row = 1, 3
            IF (left(row, column) < right(row, column) - 1.0e-12) THEN
               lexicographically_precedes = .TRUE.
               RETURN
            ELSE IF (left(row, column) > right(row, column) + 1.0e-12) THEN
               RETURN
            END IF
         END DO
      END DO
   END FUNCTION lexicographically_precedes

   REAL FUNCTION vector_norm(vector)
      REAL, INTENT(IN) :: vector(3)
      vector_norm = SQRT(DOT_PRODUCT(vector, vector))
   END FUNCTION vector_norm

   FUNCTION outer_product(left, right) RESULT(product)
      REAL, INTENT(IN) :: left(3), right(3)
      REAL :: product(3, 3)
      INTEGER :: i, j

      DO i = 1, 3
         DO j = 1, 3
            product(i, j) = left(i)*right(j)
         END DO
      END DO
   END FUNCTION outer_product

   REAL FUNCTION determinant_3x3(matrix)
      REAL, INTENT(IN) :: matrix(3, 3)
      determinant_3x3 = matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2)) &
                      - matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1)) &
                      + matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
   END FUNCTION determinant_3x3

   REAL FUNCTION angle_degrees(left, right)
      REAL, INTENT(IN) :: left(3), right(3)
      REAL :: cosine

      cosine = DOT_PRODUCT(left, right)/(vector_norm(left)*vector_norm(right))
      cosine = MAX(-1.0, MIN(1.0, cosine))
      angle_degrees = ACOS(cosine)*180.0/pi
   END FUNCTION angle_degrees

   FUNCTION identity_3() RESULT(identity)
      REAL :: identity(3, 3)
      INTEGER :: i

      identity = 0.0
      DO i = 1, 3
         identity(i, i) = 1.0
      END DO
   END FUNCTION identity_3

END MODULE m_local_structural_frame
