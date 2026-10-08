!--------------------------------------------------------------------------------
! Standalone validation of distorted-octahedron frames and common orbital-spin
! density-matrix rotations. No production XAS/RIXS or jDOS code is used.
!--------------------------------------------------------------------------------

PROGRAM local_d_spin_frame_test
   USE m_local_d_spin_analysis, ONLY: local_d_complex_to_real_unitary, local_d_extract_t2g, &
                                      local_d_jeff_half_minus, local_d_jeff_three_plus_three, &
                                      local_d_jeff_unitary, local_d_jeff_weights, &
                                      local_d_n_jeff, local_d_n_t2g, local_d_transform_to_real
   USE m_local_d_spin_frame_transform, ONLY: local_d_orbital_global_to_local, &
                                              local_d_spin_native_to_structural, &
                                              spinor_local_to_global, transform_local_d_spin_density
   USE m_local_d_spin_projector_core, ONLY: local_d_l, local_d_n_orbitals, local_d_n_spins
   USE m_local_structural_frame, ONLY: closest_proper_rotation, &
                                       construct_structural_frame_from_displacements, &
                                       construct_structural_frame_from_positions, &
                                       t_structural_frame_diagnostics
   IMPLICIT NONE

   REAL, PARAMETER :: tolerance = 2.0e-11
   INTEGER :: failures
   REAL :: maximum_error

   failures = 0
   maximum_error = 0.0
   CALL test_geometry(failures, maximum_error)
   CALL test_orbital_rotations(failures, maximum_error)
   CALL test_spin_rotations(failures, maximum_error)
   CALL test_combined_covariance(failures, maximum_error)

   WRITE (*, '(a,es24.16)') "Maximum selected absolute error = ", maximum_error
   IF (failures /= 0) THEN
      WRITE (*, '(a,i0)') "LOCAL D-SPIN FRAME TEST: FAIL; failures = ", failures
      ERROR STOP 1
   END IF
   WRITE (*, '(a)') "LOCAL D-SPIN FRAME TEST: PASS"

CONTAINS

   SUBROUTINE test_geometry(n_failures, largest_error)
      INTEGER, INTENT(INOUT) :: n_failures
      REAL, INTENT(INOUT) :: largest_error
      TYPE(t_structural_frame_diagnostics) :: diagnostics, distorted_diagnostics
      REAL :: ideal(3,6), rotated(3,6), distorted(3,6), permuted(3,6), positions(3,6)
      REAL :: reference(3,3), frame(3,3), expected(3,3), distorted_frame(3,3), permuted_frame(3,3)
      REAL :: center(3), singular_values(3), condition_number, polar_product(3,3)
      REAL :: rx(3,3), ry(3,3), rz(3,3), generic(3,3), reflected(3,3)
      REAL :: pathological(3,6)
      INTEGER :: permutation(6), i
      LOGICAL :: success, reflection_corrected
      CHARACTER(LEN=256) :: message

      reference = identity_3()
      CALL initialize_ideal_octahedron(ideal)
      CALL construct_structural_frame_from_displacements(ideal, reference, frame, diagnostics, success, message)
      CALL check_logical("ideal octahedron accepted", success, message, n_failures)
      IF (success) THEN
         CALL check_bound("ideal unequal-bond frame", MAXVAL(ABS(frame - reference)), &
                          tolerance, n_failures, largest_error)
         CALL check_bound("ideal opposite-pair cost", diagnostics%total_pair_cost, &
                          tolerance, n_failures, largest_error)
         CALL check_bound("ideal opposite-pair angles", &
                          MAXVAL(ABS(diagnostics%opposite_pair_angles_deg - 180.0)), &
                          1.0e-10, n_failures, largest_error)
         CALL check_bound("ideal raw mutual angles", &
                          MAXVAL(ABS(diagnostics%raw_mutual_angles_deg - 90.0)), &
                          1.0e-10, n_failures, largest_error)
         CALL check_vector("ideal bond lengths", diagnostics%bond_lengths, &
                           [1.1,1.1,1.3,1.3,1.7,1.7], tolerance, n_failures, largest_error)
      END IF

      rx = axis_rotation([1.0,0.0,0.0], 0.5*ACOS(-1.0))
      ry = axis_rotation([0.0,1.0,0.0], 0.5*ACOS(-1.0))
      rz = axis_rotation([0.0,0.0,1.0], 0.5*ACOS(-1.0))
      generic = axis_rotation(normalized([1.0,2.0,-1.5]), 0.713)
      CALL check_rigid_geometry("90-degree x rotation", ideal, rx, n_failures, largest_error)
      CALL check_rigid_geometry("90-degree y rotation", ideal, ry, n_failures, largest_error)
      CALL check_rigid_geometry("90-degree z rotation", ideal, rz, n_failures, largest_error)
      CALL check_rigid_geometry("generic rigid rotation", ideal, generic, n_failures, largest_error)
      CALL check_vector("local-global vector convention", &
                        MATMUL(TRANSPOSE(generic), MATMUL(generic, [0.2,-0.4,0.7])), &
                        [0.2,-0.4,0.7], tolerance, n_failures, largest_error)

      CALL initialize_distorted_octahedron(distorted)
      CALL construct_structural_frame_from_displacements(distorted, reference, distorted_frame, &
                                                         distorted_diagnostics, success, message)
      CALL check_logical("angularly distorted octahedron accepted", success, message, n_failures)
      IF (success) THEN
         CALL check_bound("distorted-frame orthonormality", orthogonality_error(distorted_frame), &
                          tolerance, n_failures, largest_error)
         CALL check_real("distorted-frame determinant", determinant_3x3(distorted_frame), 1.0, &
                         tolerance, n_failures, largest_error)
         CALL check_real("distorted pair-cost diagnostic sum", &
                         SUM(distorted_diagnostics%opposite_pair_costs), &
                         distorted_diagnostics%total_pair_cost, tolerance, n_failures, largest_error)
         polar_product = MATMUL(TRANSPOSE(distorted_frame), distorted_diagnostics%raw_axes)
         CALL check_bound("closest-frame polar symmetry", MAXVAL(ABS(polar_product - TRANSPOSE(polar_product))), &
                          5.0e-11, n_failures, largest_error)
         IF (MINVAL(distorted_diagnostics%opposite_pair_angles_deg) >= 179.999 .OR. &
             MAXVAL(ABS(distorted_diagnostics%raw_mutual_angles_deg - 90.0)) <= 1.0e-3) THEN
            WRITE (*, '(a)') "FAIL: distorted fixture did not generate intended angular diagnostics"
            n_failures = n_failures + 1
         ELSE
            WRITE (*, '(a,3f11.5)') "PASS: distorted opposite-pair angles = ", &
                                    distorted_diagnostics%opposite_pair_angles_deg
            WRITE (*, '(a,3f11.5)') "PASS: distorted raw mutual angles = ", &
                                    distorted_diagnostics%raw_mutual_angles_deg
         END IF
      END IF

      center = [0.025, -0.018, 0.022]
      positions = distorted
      CALL construct_structural_frame_from_positions(center, positions, reference, frame, diagnostics, success, message)
      CALL check_logical("off-centered cage accepted", success, message, n_failures)
      IF (success) THEN
         CALL check_physical_bound("off-centered normalized-pair frame change", &
                                   MAXVAL(ABS(frame - distorted_frame)), 5.0e-2, n_failures)
      END IF

      permutation = [4,1,6,3,2,5]
      DO i = 1, 6
         permuted(:, i) = distorted(:, permutation(i))
      END DO
      CALL construct_structural_frame_from_displacements(permuted, reference, permuted_frame, diagnostics, &
                                                         success, message)
      CALL check_logical("permuted ligands accepted", success, message, n_failures)
      IF (success) CALL check_bound("ligand-permutation invariance", &
                                    MAXVAL(ABS(permuted_frame - distorted_frame)), &
                                    tolerance, n_failures, largest_error)

      permutation = [2,1,4,3,6,5]
      DO i = 1, 6
         permuted(:, i) = distorted(:, permutation(i))
      END DO
      CALL construct_structural_frame_from_displacements(permuted, reference, permuted_frame, diagnostics, &
                                                         success, message)
      CALL check_logical("pair-member reversal accepted", success, message, n_failures)
      IF (success) CALL check_bound("pair-member sign invariance", &
                                    MAXVAL(ABS(permuted_frame - distorted_frame)), &
                                    tolerance, n_failures, largest_error)

      reflected = identity_3()
      reflected(:,3) = -reflected(:,3)
      CALL closest_proper_rotation(reflected, frame, singular_values, condition_number, reflection_corrected, &
                                   success, message)
      CALL check_logical("reflection polar factor accepted", success, message, n_failures)
      IF (success) THEN
         CALL check_true("SVD reflection correction reported", reflection_corrected, n_failures)
         CALL check_real("reflection corrected determinant", determinant_3x3(frame), 1.0, &
                         tolerance, n_failures, largest_error)
         expected = 0.0
         expected(1,1) = -1.0
         expected(2,2) = 1.0
         expected(3,3) = -1.0
         CALL check_bound("degenerate-reflection deterministic choice", MAXVAL(ABS(frame-expected)), &
                          tolerance, n_failures, largest_error)
      END IF

      pathological(:,1) = [ 1.0, 0.0, 0.0]
      pathological(:,2) = [-1.0, 0.0, 0.0]
      pathological(:,3) = [ 1.0, 1.0e-10, 0.0]
      pathological(:,4) = [-1.0,-1.0e-10, 0.0]
      pathological(:,5) = [ 1.0, 0.0, 1.0e-10]
      pathological(:,6) = [-1.0, 0.0,-1.0e-10]
      CALL construct_structural_frame_from_displacements(pathological, reference, frame, diagnostics, success, message)
      IF (success) THEN
         WRITE (*, '(a)') "FAIL: nearly dependent ligand axes were accepted"
         n_failures = n_failures + 1
      ELSE IF (LEN_TRIM(message) == 0) THEN
         WRITE (*, '(a)') "FAIL: pathological geometry returned no diagnostic"
         n_failures = n_failures + 1
      ELSE
         WRITE (*, '(a,a)') "PASS: pathological geometry rejected: ", TRIM(message)
      END IF
   END SUBROUTINE test_geometry

   SUBROUTINE check_rigid_geometry(name, ideal, rotation, n_failures, largest_error)
      CHARACTER(LEN=*), INTENT(IN) :: name
      REAL, INTENT(IN) :: ideal(3,6), rotation(3,3)
      INTEGER, INTENT(INOUT) :: n_failures
      REAL, INTENT(INOUT) :: largest_error
      TYPE(t_structural_frame_diagnostics) :: diagnostics
      REAL :: rotated(3,6), recovered(3,3)
      INTEGER :: ligand
      LOGICAL :: success
      CHARACTER(LEN=256) :: message

      DO ligand = 1, 6
         rotated(:,ligand) = MATMUL(rotation, ideal(:,ligand))
      END DO
      ! A perfect octahedron is invariant under cubic 90-degree rotations.
      ! The correspondingly oriented reference supplies the otherwise absent
      ! axis labels and signs, as required by the frame definition.
      CALL construct_structural_frame_from_displacements(rotated, rotation, recovered, diagnostics, success, message)
      CALL check_logical(TRIM(name)//" accepted", success, message, n_failures)
      IF (success) CALL check_bound(name, MAXVAL(ABS(recovered - rotation)), &
                                    tolerance, n_failures, largest_error)
   END SUBROUTINE check_rigid_geometry

   SUBROUTINE test_orbital_rotations(n_failures, largest_error)
      INTEGER, INTENT(INOUT) :: n_failures
      REAL, INTENT(INOUT) :: largest_error
      COMPLEX :: d_matrix(-local_d_l:local_d_l, -local_d_l:local_d_l)
      COMPLEX :: d_inverse(-local_d_l:local_d_l, -local_d_l:local_d_l)
      COMPLEX :: complex_to_real(-local_d_l:local_d_l, local_d_n_orbitals)
      COMPLEX :: global_state(-local_d_l:local_d_l), local_state(-local_d_l:local_d_l)
      COMPLEX :: local_real(local_d_n_orbitals), expected_real(local_d_n_orbitals)
      REAL :: rz(3,3), generic(3,3), expected_sign(5)
      INTEGER :: expected_orbital(5), orbital

      rz = axis_rotation([0.0,0.0,1.0], 0.5*ACOS(-1.0))
      generic = axis_rotation(normalized([1.0,-2.0,0.7]), 0.619)
      CALL local_d_complex_to_real_unitary(complex_to_real)
      CALL local_d_orbital_global_to_local(rz, d_matrix)
      CALL check_bound("l=2 orbital transformation unitarity", unitary_error(d_matrix), &
                       tolerance, n_failures, largest_error)
      expected_orbital = [1,3,2,4,5]
      expected_sign = [-1.0,1.0,-1.0,-1.0,1.0]
      DO orbital = 1, local_d_n_orbitals
         global_state = complex_to_real(:,orbital)
         local_state = MATMUL(d_matrix, global_state)
         local_real = MATMUL(CONJG(TRANSPOSE(complex_to_real)), local_state)
         expected_real = CMPLX(0.0,0.0)
         expected_real(expected_orbital(orbital)) = CMPLX(expected_sign(orbital),0.0)
         CALL check_bound("analytical 90-degree real-d mapping", MAXVAL(ABS(local_real - expected_real)), &
                          tolerance, n_failures, largest_error)
      END DO
      CALL local_d_orbital_global_to_local(generic, d_matrix)
      CALL local_d_orbital_global_to_local(TRANSPOSE(generic), d_inverse)
      CALL check_bound("l=2 inverse rotation", MAXVAL(ABS(d_inverse - CONJG(TRANSPOSE(d_matrix)))), &
                       tolerance, n_failures, largest_error)
      CALL test_all_orbital_density_recovery(generic, n_failures, largest_error)
   END SUBROUTINE test_orbital_rotations

   SUBROUTINE test_all_orbital_density_recovery(rotation, n_failures, largest_error)
      REAL, INTENT(IN) :: rotation(3,3)
      INTEGER, INTENT(INOUT) :: n_failures
      REAL, INTENT(INOUT) :: largest_error
      COMPLEX :: complex_to_real(-local_d_l:local_d_l, local_d_n_orbitals)
      COMPLEX :: d_matrix(-local_d_l:local_d_l, -local_d_l:local_d_l)
      COMPLEX :: spin_identity(2,2), local_state(10), global_state(10)
      COMPLEX :: rho_local(-local_d_l:local_d_l,2,-local_d_l:local_d_l,2)
      COMPLEX :: rho_global(-local_d_l:local_d_l,2,-local_d_l:local_d_l,2)
      COMPLEX :: rho_recovered(-local_d_l:local_d_l,2,-local_d_l:local_d_l,2)
      COMPLEX :: transform(10,10), spin_transform(2,2), orbital_transform(-2:2,-2:2)
      REAL :: trace_before, trace_after
      INTEGER :: orbital, m

      CALL local_d_complex_to_real_unitary(complex_to_real)
      spin_identity = identity_2_complex()
      CALL local_d_spin_native_to_structural(rotation, spin_identity, orbital_transform, spin_transform, transform)
      DO orbital = 1, local_d_n_orbitals
         local_state = CMPLX(0.0,0.0)
         DO m = -local_d_l, local_d_l
            local_state((m+local_d_l)*2+1) = complex_to_real(m,orbital)
         END DO
         global_state = MATMUL(CONJG(TRANSPOSE(transform)), local_state)
         CALL flat_state_density(global_state, rho_global)
         CALL flat_state_density(local_state, rho_local)
         CALL transform_local_d_spin_density(rho_global, rotation, spin_identity, rho_recovered)
         CALL check_bound("pure real-d density recovery", MAXVAL(ABS(rho_recovered-rho_local)), &
                          tolerance, n_failures, largest_error)
         trace_before = density_trace(rho_global)
         trace_after = density_trace(rho_recovered)
         CALL check_real("orbital density trace conservation", trace_after, trace_before, &
                         tolerance, n_failures, largest_error)
         CALL check_bound("orbital density Hermiticity", hermiticity_error(rho_recovered), &
                          tolerance, n_failures, largest_error)
         CALL check_minimum_eigenvalue("orbital density positivity", rho_recovered, &
                                       -tolerance, n_failures, largest_error)
      END DO
   END SUBROUTINE test_all_orbital_density_recovery

   SUBROUTINE test_spin_rotations(n_failures, largest_error)
      INTEGER, INTENT(INOUT) :: n_failures
      REAL, INTENT(INOUT) :: largest_error
      REAL :: structural_rotation(3,3), native_rotation(3,3), global_direction(3), expected_local(3)
      REAL :: directions(3,6), local_direction(3)
      COMPLEX :: structural_spin(2,2), native_to_global(2,2), orbital_transform(-2:2,-2:2)
      COMPLEX :: native_to_local(2,2), global_density(2,2), local_density(2,2), native_density(2,2)
      COMPLEX :: pauli(2,2,3), covariance(2,2), expected_matrix(2,2)
      INTEGER :: axis, state

      structural_rotation = axis_rotation(normalized([1.0,2.0,-0.5]), 0.827)
      CALL spinor_local_to_global(structural_rotation, structural_spin)
      CALL initialize_pauli(pauli)
      CALL check_bound("spin-1/2 transformation unitarity", unitary_error(structural_spin), &
                       tolerance, n_failures, largest_error)
      DO axis = 1, 3
         covariance = MATMUL(structural_spin, MATMUL(pauli(:,:,axis), CONJG(TRANSPOSE(structural_spin))))
         expected_matrix = CMPLX(0.0,0.0)
         expected_matrix = structural_rotation(1,axis)*pauli(:,:,1) &
                         + structural_rotation(2,axis)*pauli(:,:,2) &
                         + structural_rotation(3,axis)*pauli(:,:,3)
         CALL check_bound("SU(2)-SO(3) covariance", MAXVAL(ABS(covariance-expected_matrix)), &
                          tolerance, n_failures, largest_error)
      END DO
      directions(:,1)=[1.0,0.0,0.0]; directions(:,2)=[-1.0,0.0,0.0]
      directions(:,3)=[0.0,1.0,0.0]; directions(:,4)=[0.0,-1.0,0.0]
      directions(:,5)=[0.0,0.0,1.0]; directions(:,6)=[0.0,0.0,-1.0]
      DO state = 1, 6
         global_direction = directions(:,state)
         CALL bloch_density(global_direction, global_density)
         local_density = MATMUL(CONJG(TRANSPOSE(structural_spin)), &
                                MATMUL(global_density, structural_spin))
         CALL density_bloch(local_density, local_direction)
         expected_local = MATMUL(TRANSPOSE(structural_rotation),global_direction)
         CALL check_vector("Cartesian spin global-to-structural", local_direction, expected_local, &
                           tolerance, n_failures, largest_error)
      END DO

      native_rotation = axis_rotation(normalized([-0.3,0.8,1.2]), 1.137)
      CALL spinor_local_to_global(native_rotation, native_to_global)
      CALL local_d_spin_native_to_structural(structural_rotation, native_to_global, &
                                             orbital_transform, native_to_local)
      CALL bloch_density(normalized([0.4,-0.2,0.9]), native_density)
      local_density = MATMUL(native_to_local, MATMUL(native_density, CONJG(TRANSPOSE(native_to_local))))
      CALL density_bloch(local_density, local_direction)
      expected_local = MATMUL(TRANSPOSE(structural_rotation), &
                              MATMUL(native_rotation,normalized([0.4,-0.2,0.9])))
      CALL check_vector("native-MT to global to structural spin route", local_direction, expected_local, &
                        tolerance, n_failures, largest_error)
      CALL check_bound("SU(2) double-cover density invariance", &
         MAXVAL(ABS(MATMUL(native_to_local,MATMUL(native_density,CONJG(TRANSPOSE(native_to_local)))) &
                   -MATMUL(-native_to_local,MATMUL(native_density,CONJG(TRANSPOSE(-native_to_local)))))), &
         tolerance, n_failures, largest_error)
   END SUBROUTINE test_spin_rotations

   SUBROUTINE test_combined_covariance(n_failures, largest_error)
      INTEGER, INTENT(INOUT) :: n_failures
      REAL, INTENT(INOUT) :: largest_error
      TYPE(t_structural_frame_diagnostics) :: diagnostics
      REAL :: ideal(3,6), rotated_ligands(3,6), physical_rotation(3,3), recovered_frame(3,3)
      REAL :: mj_weights(local_d_n_jeff), j_weights(2), normalized_mj(local_d_n_jeff), normalized_j(2)
      COMPLEX :: jeff(6,6), complex_to_real(-2:2,5), spin_identity(2,2), transform(10,10)
      COMPLEX :: orbital_transform(-2:2,-2:2), spin_transform(2,2)
      COMPLEX :: local_state(10), global_state(10), mixed_state(10)
      COMPLEX :: rho_global(-2:2,2,-2:2,2), rho_recovered(-2:2,2,-2:2,2)
      COMPLEX :: rho_real(5,2,5,2), rho_t2g(3,2,3,2)
      INTEGER :: ligand
      LOGICAL :: success, normalization_defined
      CHARACTER(LEN=256) :: message

      CALL initialize_ideal_octahedron(ideal)
      physical_rotation = axis_rotation(normalized([0.7,-1.1,0.4]), 0.893)
      DO ligand=1,6
         rotated_ligands(:,ligand)=MATMUL(physical_rotation,ideal(:,ligand))
      END DO
      CALL construct_structural_frame_from_displacements(rotated_ligands, physical_rotation, recovered_frame, &
                                                         diagnostics, success, message)
      CALL check_logical("combined-covariance frame accepted", success, message, n_failures)
      IF (.NOT.success) RETURN
      CALL local_d_jeff_unitary(jeff)
      CALL local_d_complex_to_real_unitary(complex_to_real)
      spin_identity=identity_2_complex()
      CALL local_d_spin_native_to_structural(recovered_frame,spin_identity,orbital_transform,spin_transform,transform)

      CALL embed_t2g_real_spinor(jeff(:,local_d_jeff_half_minus),complex_to_real,local_state)
      global_state=MATMUL(CONJG(TRANSPOSE(transform)),local_state)
      CALL flat_state_density(global_state,rho_global)
      CALL transform_local_d_spin_density(rho_global,recovered_frame,spin_identity,rho_recovered)
      CALL check_bound("entangled j_eff state covariance", &
                       MAXVAL(ABS(rho_recovered-flat_density(local_state))), &
                       tolerance,n_failures,largest_error)
      CALL analyze_jeff(rho_recovered,mj_weights,j_weights,normalized_mj,normalized_j,normalization_defined)
      CALL check_vector2("rotated ideal j_eff weights",j_weights,[1.0,0.0], &
                         tolerance,n_failures,largest_error)

      CALL embed_t2g_real_spinor(jeff(:,local_d_jeff_half_minus),complex_to_real,local_state)
      CALL embed_t2g_real_spinor(jeff(:,local_d_jeff_three_plus_three),complex_to_real,mixed_state)
      mixed_state=SQRT(0.60)*local_state+SQRT(0.30)*mixed_state
      mixed_state=mixed_state+SQRT(0.10)*eg_spinor(complex_to_real)
      global_state=MATMUL(CONJG(TRANSPOSE(transform)),mixed_state)
      CALL flat_state_density(global_state,rho_global)
      CALL transform_local_d_spin_density(rho_global,recovered_frame,spin_identity,rho_recovered)
      CALL analyze_jeff(rho_recovered,mj_weights,j_weights,normalized_mj,normalized_j,normalization_defined)
      CALL check_vector2("physical j_eff-sector mixing retained",j_weights,[0.60,0.30], &
                         tolerance,n_failures,largest_error)
      CALL local_d_transform_to_real(rho_recovered,rho_real)
      CALL check_real("physical t2g-eg mixing retained", &
                      SUM([(REAL(rho_real(5,ligand,5,ligand)),ligand=1,2)]),0.10, &
                      tolerance,n_failures,largest_error)
      CALL check_real("mixed-state trace",density_trace(rho_recovered),1.0, &
                      tolerance,n_failures,largest_error)
   END SUBROUTINE test_combined_covariance

   SUBROUTINE analyze_jeff(rho_complex,mj_weights,j_weights,normalized_mj,normalized_j,defined)
      COMPLEX,INTENT(IN)::rho_complex(-2:2,2,-2:2,2)
      REAL,INTENT(OUT)::mj_weights(6),j_weights(2),normalized_mj(6),normalized_j(2)
      LOGICAL,INTENT(OUT)::defined
      COMPLEX::rho_real(5,2,5,2),rho_t2g(3,2,3,2)
      CALL local_d_transform_to_real(rho_complex,rho_real)
      CALL local_d_extract_t2g(rho_real,rho_t2g)
      CALL local_d_jeff_weights(rho_t2g,mj_weights,j_weights,normalized_mj,normalized_j,defined)
   END SUBROUTINE analyze_jeff

   SUBROUTINE embed_t2g_real_spinor(t2g_state,complex_to_real,complex_state)
      COMPLEX,INTENT(IN)::t2g_state(6),complex_to_real(-2:2,5)
      COMPLEX,INTENT(OUT)::complex_state(10)
      COMPLEX::real_state(5,2)
      INTEGER::orbital,spin,m,row
      real_state=CMPLX(0.0,0.0)
      DO orbital=1,3
         DO spin=1,2
            real_state(orbital,spin)=t2g_state((orbital-1)*2+spin)
         END DO
      END DO
      complex_state=CMPLX(0.0,0.0)
      DO m=-2,2
         DO spin=1,2
            row=(m+2)*2+spin
            DO orbital=1,5
               complex_state(row)=complex_state(row)+complex_to_real(m,orbital)*real_state(orbital,spin)
            END DO
         END DO
      END DO
   END SUBROUTINE embed_t2g_real_spinor

   FUNCTION eg_spinor(complex_to_real) RESULT(state)
      COMPLEX,INTENT(IN)::complex_to_real(-2:2,5)
      COMPLEX::state(10)
      INTEGER::m
      state=CMPLX(0.0,0.0)
      DO m=-2,2
         state((m+2)*2+1)=complex_to_real(m,5)
      END DO
   END FUNCTION eg_spinor

   SUBROUTINE initialize_ideal_octahedron(displacements)
      REAL, INTENT(OUT) :: displacements(3,6)
      displacements(:,1)=[ 1.1,0.0,0.0];displacements(:,2)=[-1.1,0.0,0.0]
      displacements(:,3)=[0.0, 1.3,0.0];displacements(:,4)=[0.0,-1.3,0.0]
      displacements(:,5)=[0.0,0.0, 1.7];displacements(:,6)=[0.0,0.0,-1.7]
   END SUBROUTINE initialize_ideal_octahedron

   SUBROUTINE initialize_distorted_octahedron(displacements)
      REAL, INTENT(OUT) :: displacements(3,6)
      displacements(:,1)=1.10*normalized([ 1.00, 0.08, 0.02])
      displacements(:,2)=1.04*normalized([-1.00, 0.03,-0.04])
      displacements(:,3)=1.28*normalized([ 0.05, 1.00, 0.10])
      displacements(:,4)=1.34*normalized([-0.02,-1.00, 0.06])
      displacements(:,5)=1.63*normalized([ 0.04,-0.07, 1.00])
      displacements(:,6)=1.72*normalized([-0.06, 0.02,-1.00])
   END SUBROUTINE initialize_distorted_octahedron

   FUNCTION axis_rotation(axis,angle) RESULT(rotation)
      REAL,INTENT(IN)::axis(3),angle
      REAL::rotation(3,3),n(3),c,s
      n=normalized(axis);c=COS(angle);s=SIN(angle)
      rotation=0.0
      rotation(1,1)=c+n(1)*n(1)*(1-c);rotation(1,2)=n(1)*n(2)*(1-c)-n(3)*s
      rotation(1,3)=n(1)*n(3)*(1-c)+n(2)*s
      rotation(2,1)=n(2)*n(1)*(1-c)+n(3)*s;rotation(2,2)=c+n(2)*n(2)*(1-c)
      rotation(2,3)=n(2)*n(3)*(1-c)-n(1)*s
      rotation(3,1)=n(3)*n(1)*(1-c)-n(2)*s;rotation(3,2)=n(3)*n(2)*(1-c)+n(1)*s
      rotation(3,3)=c+n(3)*n(3)*(1-c)
   END FUNCTION axis_rotation

   FUNCTION normalized(vector) RESULT(unit)
      REAL,INTENT(IN)::vector(3)
      REAL::unit(3)
      unit=vector/SQRT(DOT_PRODUCT(vector,vector))
   END FUNCTION normalized

   FUNCTION identity_3() RESULT(identity)
      REAL::identity(3,3);INTEGER::i
      identity=0.0;DO i=1,3;identity(i,i)=1.0;END DO
   END FUNCTION identity_3

   FUNCTION identity_2_complex() RESULT(identity)
      COMPLEX::identity(2,2)
      identity=CMPLX(0.0,0.0);identity(1,1)=1.0;identity(2,2)=1.0
   END FUNCTION identity_2_complex

   REAL FUNCTION determinant_3x3(matrix)
      REAL,INTENT(IN)::matrix(3,3)
      determinant_3x3=matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2)) &
         -matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1)) &
         +matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
   END FUNCTION determinant_3x3

   REAL FUNCTION orthogonality_error(matrix)
      REAL,INTENT(IN)::matrix(3,3)
      orthogonality_error=MAXVAL(ABS(MATMUL(TRANSPOSE(matrix),matrix)-identity_3()))
   END FUNCTION orthogonality_error

   REAL FUNCTION unitary_error(matrix)
      COMPLEX,INTENT(IN)::matrix(:,:)
      COMPLEX::identity(SIZE(matrix,2),SIZE(matrix,2));INTEGER::i
      identity=CMPLX(0.0,0.0);DO i=1,SIZE(identity,1);identity(i,i)=1.0;END DO
      unitary_error=MAXVAL(ABS(MATMUL(CONJG(TRANSPOSE(matrix)),matrix)-identity))
   END FUNCTION unitary_error

   SUBROUTINE initialize_pauli(pauli)
      COMPLEX,INTENT(OUT)::pauli(2,2,3)
      pauli=CMPLX(0.0,0.0)
      pauli(1,2,1)=1.0;pauli(2,1,1)=1.0
      pauli(1,2,2)=CMPLX(0.0,-1.0);pauli(2,1,2)=CMPLX(0.0,1.0)
      pauli(1,1,3)=1.0;pauli(2,2,3)=-1.0
   END SUBROUTINE initialize_pauli

   SUBROUTINE bloch_density(direction,density)
      REAL,INTENT(IN)::direction(3);COMPLEX,INTENT(OUT)::density(2,2)
      COMPLEX::pauli(2,2,3)
      CALL initialize_pauli(pauli)
      density=0.5*(identity_2_complex()+direction(1)*pauli(:,:,1)+direction(2)*pauli(:,:,2)+direction(3)*pauli(:,:,3))
   END SUBROUTINE bloch_density

   SUBROUTINE density_bloch(density,direction)
      COMPLEX,INTENT(IN)::density(2,2);REAL,INTENT(OUT)::direction(3)
      direction(1)=2.0*REAL(density(2,1));direction(2)=2.0*AIMAG(density(2,1))
      direction(3)=REAL(density(1,1)-density(2,2))
   END SUBROUTINE density_bloch

   SUBROUTINE flat_state_density(state,density)
      COMPLEX,INTENT(IN)::state(10);COMPLEX,INTENT(OUT)::density(-2:2,2,-2:2,2)
      INTEGER::m,mp,s,sp,row,column
      DO m=-2,2;DO s=1,2;row=(m+2)*2+s
         DO mp=-2,2;DO sp=1,2;column=(mp+2)*2+sp
            density(m,s,mp,sp)=state(row)*CONJG(state(column))
         END DO;END DO
      END DO;END DO
   END SUBROUTINE flat_state_density

   FUNCTION flat_density(state) RESULT(density)
      COMPLEX,INTENT(IN)::state(10);COMPLEX::density(-2:2,2,-2:2,2)
      CALL flat_state_density(state,density)
   END FUNCTION flat_density

   REAL FUNCTION density_trace(density)
      COMPLEX,INTENT(IN)::density(-2:2,2,-2:2,2);INTEGER::m,s
      density_trace=0.0;DO m=-2,2;DO s=1,2;density_trace=density_trace+REAL(density(m,s,m,s));END DO;END DO
   END FUNCTION density_trace

   REAL FUNCTION hermiticity_error(density)
      COMPLEX,INTENT(IN)::density(-2:2,2,-2:2,2);INTEGER::m,mp,s,sp
      hermiticity_error=0.0;DO m=-2,2;DO s=1,2;DO mp=-2,2;DO sp=1,2
         hermiticity_error=MAX(hermiticity_error,ABS(density(m,s,mp,sp)-CONJG(density(mp,sp,m,s))))
      END DO;END DO;END DO;END DO
   END FUNCTION hermiticity_error

   SUBROUTINE check_minimum_eigenvalue(name,density,minimum_allowed,n_failures,largest_error)
      CHARACTER(LEN=*),INTENT(IN)::name;COMPLEX,INTENT(IN)::density(-2:2,2,-2:2,2)
      REAL,INTENT(IN)::minimum_allowed;INTEGER,INTENT(INOUT)::n_failures;REAL,INTENT(INOUT)::largest_error
      COMPLEX::matrix(10,10),work(40);REAL::eigenvalues(10),rwork(40);INTEGER::m,mp,s,sp,row,column,info
      DO m=-2,2;DO s=1,2;row=(m+2)*2+s;DO mp=-2,2;DO sp=1,2;column=(mp+2)*2+sp
         matrix(row,column)=density(m,s,mp,sp)
      END DO;END DO;END DO;END DO
      CALL ZHEEV('N','U',10,matrix,10,eigenvalues,work,SIZE(work),rwork,info)
      IF(info/=0)THEN;WRITE(*,'(a,a,i0)')'FAIL: ',TRIM(name)//' ZHEEV info=',info;n_failures=n_failures+1;RETURN;END IF
      largest_error=MAX(largest_error,MAX(0.0,-MINVAL(eigenvalues)))
      IF(MINVAL(eigenvalues)<minimum_allowed)THEN
         WRITE(*,'(a,a,es24.16)')'FAIL: ',TRIM(name)//' minimum=',MINVAL(eigenvalues);n_failures=n_failures+1
      ELSE;WRITE(*,'(a,a,es24.16)')'PASS: ',TRIM(name)//' minimum=',MINVAL(eigenvalues);END IF
   END SUBROUTINE check_minimum_eigenvalue

   SUBROUTINE check_bound(name,value,bound,n_failures,largest_error)
      CHARACTER(LEN=*),INTENT(IN)::name;REAL,INTENT(IN)::value,bound
      INTEGER,INTENT(INOUT)::n_failures;REAL,INTENT(INOUT)::largest_error
      largest_error=MAX(largest_error,value)
      IF(value>bound)THEN;WRITE(*,'(a,a,a,es24.16,a,es24.16)')'FAIL: ',TRIM(name),' error ',value,' > ',bound
         n_failures=n_failures+1
      ELSE;WRITE(*,'(a,a,a,es24.16)')'PASS: ',TRIM(name),' error = ',value;END IF
   END SUBROUTINE check_bound

   SUBROUTINE check_real(name,value,expected,bound,n_failures,largest_error)
      CHARACTER(LEN=*),INTENT(IN)::name;REAL,INTENT(IN)::value,expected,bound
      INTEGER,INTENT(INOUT)::n_failures;REAL,INTENT(INOUT)::largest_error
      CALL check_bound(name,ABS(value-expected),bound,n_failures,largest_error)
   END SUBROUTINE check_real

   SUBROUTINE check_physical_bound(name,value,bound,n_failures)
      CHARACTER(LEN=*),INTENT(IN)::name;REAL,INTENT(IN)::value,bound;INTEGER,INTENT(INOUT)::n_failures
      IF(value>bound)THEN;WRITE(*,'(a,a,a,es24.16,a,es24.16)')'FAIL: ',TRIM(name),' ',value,' > ',bound
         n_failures=n_failures+1
      ELSE;WRITE(*,'(a,a,a,es24.16)')'PASS: ',TRIM(name),' = ',value;END IF
   END SUBROUTINE check_physical_bound

   SUBROUTINE check_vector(name,value,expected,bound,n_failures,largest_error)
      CHARACTER(LEN=*),INTENT(IN)::name;REAL,INTENT(IN)::value(:),expected(:),bound
      INTEGER,INTENT(INOUT)::n_failures;REAL,INTENT(INOUT)::largest_error
      CALL check_bound(name,MAXVAL(ABS(value-expected)),bound,n_failures,largest_error)
   END SUBROUTINE check_vector

   SUBROUTINE check_vector2(name,value,expected,bound,n_failures,largest_error)
      CHARACTER(LEN=*),INTENT(IN)::name;REAL,INTENT(IN)::value(2),expected(2),bound
      INTEGER,INTENT(INOUT)::n_failures;REAL,INTENT(INOUT)::largest_error
      CALL check_bound(name,MAXVAL(ABS(value-expected)),bound,n_failures,largest_error)
   END SUBROUTINE check_vector2

   SUBROUTINE check_logical(name,value,message,n_failures)
      CHARACTER(LEN=*),INTENT(IN)::name,message;LOGICAL,INTENT(IN)::value;INTEGER,INTENT(INOUT)::n_failures
      IF(.NOT.value)THEN;WRITE(*,'(a,a,a,a)')'FAIL: ',TRIM(name),': ',TRIM(message);n_failures=n_failures+1
      ELSE;WRITE(*,'(a,a)')'PASS: ',TRIM(name);END IF
   END SUBROUTINE check_logical

   SUBROUTINE check_true(name,value,n_failures)
      CHARACTER(LEN=*),INTENT(IN)::name;LOGICAL,INTENT(IN)::value;INTEGER,INTENT(INOUT)::n_failures
      IF(.NOT.value)THEN;WRITE(*,'(a,a)')'FAIL: ',TRIM(name);n_failures=n_failures+1
      ELSE;WRITE(*,'(a,a)')'PASS: ',TRIM(name);END IF
   END SUBROUTINE check_true

END PROGRAM local_d_spin_frame_test
