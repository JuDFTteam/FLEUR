!--------------------------------------------------------------------------------
! Standalone validation of periodic physical-neighbour resolution and its bridge
! to the distorted-octahedron structural-frame constructor.
!--------------------------------------------------------------------------------

PROGRAM local_coordination_test
   USE m_local_coordination, ONLY: resolve_fleur_periodic_neighbours, &
                                   select_fleur_physical_atoms_by_type
   USE m_local_coordination_core, ONLY: resolve_periodic_neighbours_cartesian, &
                                        t_coordination_diagnostics, t_periodic_neighbor
   USE m_local_structural_frame, ONLY: construct_structural_frame_from_displacements, &
                                       t_structural_frame_diagnostics
   USE m_types_atoms, ONLY: t_atoms
   USE m_types_cell, ONLY: t_cell
   IMPLICIT NONE

   REAL, PARAMETER :: tolerance = 5.0e-11
   INTEGER :: failures
   REAL :: maximum_error

   failures = 0
   maximum_error = 0.0
   CALL test_cubic_boundaries_filtering_and_permutation(failures,maximum_error)
   CALL test_skew_cell(failures,maximum_error)
   CALL test_distinct_images_of_one_atom(failures,maximum_error)
   CALL test_distorted_shell_and_ambiguity(failures,maximum_error)
   CALL test_film_periodicity(failures,maximum_error)
   CALL test_structural_frame_bridge(failures,maximum_error)

   WRITE(*,'(a,es24.16)') "Maximum selected absolute error = ", maximum_error
   IF (failures /= 0) THEN
      WRITE(*,'(a,i0)') "LOCAL COORDINATION TEST: FAIL; failures = ", failures
      ERROR STOP 1
   END IF
   WRITE(*,'(a)') "LOCAL COORDINATION TEST: PASS"

CONTAINS

   SUBROUTINE test_cubic_boundaries_filtering_and_permutation(n_failures,largest_error)
      INTEGER, INTENT(INOUT) :: n_failures
      REAL, INTENT(INOUT) :: largest_error
      TYPE(t_periodic_neighbor), ALLOCATABLE :: neighbors(:), permuted_neighbors(:), adapter_neighbors(:)
      TYPE(t_coordination_diagnostics) :: diagnostics, permuted_diagnostics, adapter_diagnostics
      TYPE(t_cell) :: cell
      TYPE(t_atoms) :: atoms
      REAL :: lattice(3,3), center(3), positions(3,7), permuted_positions(3,7)
      INTEGER :: ids(7), permutation(7), permuted_ids(7), expected_translation(3,6)
      INTEGER, ALLOCATABLE :: ligand_atoms(:)
      INTEGER :: i, index
      LOGICAL :: success
      CHARACTER(LEN=256) :: message

      lattice = identity_3()
      center = [0.95,0.95,0.95]
      positions(:,1) = [0.05,0.95,0.95]
      positions(:,2) = [0.85,0.95,0.95]
      positions(:,3) = [0.95,0.05,0.95]
      positions(:,4) = [0.95,0.85,0.95]
      positions(:,5) = [0.95,0.95,0.05]
      positions(:,6) = [0.95,0.95,0.85]
      positions(:,7) = [0.25,0.95,0.95]
      ids = [2,3,4,5,6,7,8]
      expected_translation(:,1) = [1,0,0]
      expected_translation(:,2) = [0,0,0]
      expected_translation(:,3) = [0,1,0]
      expected_translation(:,4) = [0,0,0]
      expected_translation(:,5) = [0,0,1]
      expected_translation(:,6) = [0,0,0]

      CALL resolve_periodic_neighbours_cartesian(lattice,[.TRUE.,.TRUE.,.TRUE.],center,positions,ids,6, &
                                                 neighbors,diagnostics,success,message)
      CALL check_success("cubic boundary image resolution",success,message,n_failures)
      IF (success) THEN
         CALL check_real("cubic sixth distance",neighbors(6)%distance,0.1,tolerance,n_failures,largest_error)
         CALL check_real("cubic seventh distance",diagnostics%next_neighbor_distance,0.3,tolerance, &
                         n_failures,largest_error)
         CALL check_real("cubic sixth-seventh gap",diagnostics%shell_gap,0.2,tolerance,n_failures,largest_error)
         CALL check_real("cubic seventh-to-sixth ratio",diagnostics%next_to_last_ratio,3.0,tolerance, &
                         n_failures,largest_error)
         DO i = 1, 6
            index = find_atom_id(neighbors,ids(i))
            CALL check_true("selected cubic ligand id",index > 0,n_failures)
            IF (index > 0) THEN
               CALL check_integer_vector("cubic image translation",neighbors(index)%translation, &
                                         expected_translation(:,i),n_failures)
            END IF
         END DO
      END IF

      permutation = [5,2,7,1,6,3,4]
      DO i = 1, 7
         permuted_positions(:,i) = positions(:,permutation(i))
         permuted_ids(i) = ids(permutation(i))
      END DO
      CALL resolve_periodic_neighbours_cartesian(lattice,[.TRUE.,.TRUE.,.TRUE.],center,permuted_positions, &
                                                 permuted_ids,6,permuted_neighbors,permuted_diagnostics, &
                                                 success,message)
      CALL check_success("permuted candidate resolution",success,message,n_failures)
      IF (success) CALL compare_neighbor_lists("candidate input-order invariance",neighbors,permuted_neighbors, &
                                               n_failures,largest_error)

      ! Exercise the thin FLEUR adapter and type-based candidate filtering.
      cell%amat = lattice
      atoms%nat = 8
      atoms%ntype = 3
      ALLOCATE(atoms%pos(3,8),atoms%itype(8))
      atoms%pos(:,1) = center
      atoms%pos(:,2:8) = positions
      atoms%itype = [1,2,2,2,2,2,2,3]
      CALL select_fleur_physical_atoms_by_type(atoms,[2],ligand_atoms,success,message)
      CALL check_success("FLEUR atom-type candidate filtering",success,message,n_failures)
      IF (success) THEN
         CALL check_integer("filtered physical-atom count",SIZE(ligand_atoms),6,n_failures)
         CALL resolve_fleur_periodic_neighbours(cell,atoms,.FALSE.,1,ligand_atoms,6,adapter_neighbors, &
                                                adapter_diagnostics,success,message)
         CALL check_success("FLEUR physical-atom adapter",success,message,n_failures)
         IF (success) THEN
            CALL check_true("non-ligand candidate excluded",find_atom_id(adapter_neighbors,8) == 0,n_failures)
            CALL check_real("adapter sixth distance",adapter_neighbors(6)%distance,0.1,tolerance, &
                            n_failures,largest_error)
         END IF
      END IF
   END SUBROUTINE test_cubic_boundaries_filtering_and_permutation

   SUBROUTINE test_skew_cell(n_failures,largest_error)
      INTEGER, INTENT(INOUT) :: n_failures
      REAL, INTENT(INOUT) :: largest_error
      TYPE(t_periodic_neighbor), ALLOCATABLE :: neighbors(:)
      TYPE(t_coordination_diagnostics) :: diagnostics
      REAL :: lattice(3,3), center(3), positions(3,1), naive_displacement(3), expected(3)
      INTEGER :: ids(1)
      LOGICAL :: success
      CHARACTER(LEN=256) :: message

      lattice(:,1) = [1.0,0.0,0.0]
      lattice(:,2) = [0.9,0.2,0.0]
      lattice(:,3) = [0.0,0.0,2.0]
      center = 0.0
      naive_displacement = MATMUL(lattice,[0.49,0.49,0.0])
      positions(:,1) = naive_displacement
      ids = [11]
      expected = MATMUL(lattice,[0.49,-0.51,0.0])

      CALL resolve_periodic_neighbours_cartesian(lattice,[.TRUE.,.TRUE.,.TRUE.],center,positions,ids,1, &
                                                 neighbors,diagnostics,success,message)
      CALL check_success("skew-cell nearest image",success,message,n_failures)
      IF (success) THEN
         CALL check_integer_vector("skew-cell translation",neighbors(1)%translation,[0,-1,0],n_failures)
         CALL check_vector("skew-cell Cartesian displacement",neighbors(1)%displacement,expected,tolerance, &
                           n_failures,largest_error)
         CALL check_true("skew-cell result beats component wrapping", &
                         neighbors(1)%distance < vector_norm(naive_displacement),n_failures)
         CALL check_true("skew-cell adaptive search proved completeness",diagnostics%search_range > 1,n_failures)
      END IF
   END SUBROUTINE test_skew_cell

   SUBROUTINE test_distinct_images_of_one_atom(n_failures,largest_error)
      INTEGER, INTENT(INOUT) :: n_failures
      REAL, INTENT(INOUT) :: largest_error
      TYPE(t_periodic_neighbor), ALLOCATABLE :: neighbors(:)
      TYPE(t_coordination_diagnostics) :: diagnostics
      REAL :: position(3,1)
      LOGICAL :: success
      CHARACTER(LEN=256) :: message

      position(:,1) = [0.5,0.0,0.0]
      CALL resolve_periodic_neighbours_cartesian(identity_3(),[.TRUE.,.TRUE.,.TRUE.],[0.0,0.0,0.0], &
                                                 position,[21],2,neighbors,diagnostics,success,message)
      CALL check_success("distinct images of one representative atom",success,message,n_failures)
      IF (success) THEN
         CALL check_integer("same representative retained twice",COUNT(neighbors%atom_id == 21),2,n_failures)
         CALL check_integer_vector("negative periodic image",neighbors(1)%translation,[-1,0,0],n_failures)
         CALL check_integer_vector("positive periodic image",neighbors(2)%translation,[0,0,0],n_failures)
         CALL check_real("first representative-image distance",neighbors(1)%distance,0.5,tolerance, &
                         n_failures,largest_error)
         CALL check_real("second representative-image distance",neighbors(2)%distance,0.5,tolerance, &
                         n_failures,largest_error)
      END IF
   END SUBROUTINE test_distinct_images_of_one_atom

   SUBROUTINE test_distorted_shell_and_ambiguity(n_failures,largest_error)
      INTEGER, INTENT(INOUT) :: n_failures
      REAL, INTENT(INOUT) :: largest_error
      TYPE(t_periodic_neighbor), ALLOCATABLE :: neighbors(:)
      TYPE(t_coordination_diagnostics) :: diagnostics
      REAL :: positions(3,9), ambiguous_positions(3,7), distances(7)
      INTEGER :: ids(9), ambiguous_ids(7), i
      LOGICAL :: success
      CHARACTER(LEN=256) :: message

      positions(:,1) = [ 0.93, 0.02, 0.01]
      positions(:,2) = [-1.04, 0.01,-0.02]
      positions(:,3) = [ 0.03, 1.08, 0.02]
      positions(:,4) = [-0.02,-0.97, 0.01]
      positions(:,5) = [ 0.01, 0.04, 1.13]
      positions(:,6) = [-0.02,-0.03,-1.01]
      positions(:,7) = [ 1.55, 0.20, 0.00]
      positions(:,8) = [ 0.10,-1.73, 0.10]
      positions(:,9) = [ 0.00, 0.00, 1.91]
      ids = [(30+i,i=1,9)]
      CALL resolve_periodic_neighbours_cartesian(4.0*identity_3(),[.FALSE.,.FALSE.,.FALSE.],[0.0,0.0,0.0], &
                                                 positions,ids,6,neighbors,diagnostics,success,message)
      CALL check_success("distorted first-shell selection",success,message,n_failures)
      IF (success) THEN
         CALL check_true("distorted shell has positive sixth-seventh gap",diagnostics%shell_gap > 0.3,n_failures)
         CALL check_integer("distorted shell selected six intended ids",COUNT(neighbors%atom_id <= 36),6,n_failures)
      END IF

      distances = [0.91,0.93,0.95,0.97,0.99,1.0,1.0+2.0e-7]
      ambiguous_positions = 0.0
      DO i = 1, 7
         ambiguous_positions(1,i) = distances(i)
      END DO
      ambiguous_ids = [(50+i,i=1,7)]
      CALL resolve_periodic_neighbours_cartesian(4.0*identity_3(),[.FALSE.,.FALSE.,.FALSE.],[0.0,0.0,0.0], &
                                                 ambiguous_positions,ambiguous_ids,6,neighbors,diagnostics, &
                                                 success,message)
      CALL check_success("sixth-seventh ambiguity retained",success,message,n_failures)
      IF (success) THEN
         CALL check_true("seventh-neighbor diagnostic available",diagnostics%has_next_neighbor,n_failures)
         CALL check_real("small sixth-seventh gap",diagnostics%shell_gap,2.0e-7,5.0e-12,n_failures,largest_error)
      END IF
   END SUBROUTINE test_distorted_shell_and_ambiguity

   SUBROUTINE test_film_periodicity(n_failures,largest_error)
      INTEGER, INTENT(INOUT) :: n_failures
      REAL, INTENT(INOUT) :: largest_error
      TYPE(t_periodic_neighbor), ALLOCATABLE :: neighbors(:)
      TYPE(t_coordination_diagnostics) :: diagnostics
      REAL :: lattice(3,3), position(3,1)
      LOGICAL :: success
      CHARACTER(LEN=256) :: message

      lattice = 0.0
      lattice(1,1) = 1.0
      lattice(2,2) = 1.0
      lattice(3,3) = 10.0
      position(:,1) = [0.2,0.2,-4.9]
      CALL resolve_periodic_neighbours_cartesian(lattice,[.TRUE.,.TRUE.,.FALSE.],[0.2,0.2,4.9],position,[71],1, &
                                                 neighbors,diagnostics,success,message)
      CALL check_success("film excludes vacuum-direction images",success,message,n_failures)
      IF (success) THEN
         CALL check_integer("film z translation",neighbors(1)%translation(3),0,n_failures)
         CALL check_real("film nonperiodic z distance",neighbors(1)%distance,9.8,tolerance, &
                         n_failures,largest_error)
      END IF
   END SUBROUTINE test_film_periodicity

   SUBROUTINE test_structural_frame_bridge(n_failures,largest_error)
      INTEGER, INTENT(INOUT) :: n_failures
      REAL, INTENT(INOUT) :: largest_error
      TYPE(t_periodic_neighbor), ALLOCATABLE :: neighbors(:)
      TYPE(t_coordination_diagnostics) :: coordination_diagnostics
      TYPE(t_structural_frame_diagnostics) :: frame_diagnostics
      REAL :: lattice(3,3), center(3), rotation(3,3), symmetric_distortion(3,3), raw_axes(3,3)
      REAL :: displacements(3,6), positions(3,8), resolved_displacements(3,6), frame(3,3)
      REAL :: lengths(6), target(3), fractional(3)
      INTEGER :: ids(8), ligand, axis, translation(3)
      LOGICAL :: success
      CHARACTER(LEN=256) :: message

      lattice = 5.0*identity_3()
      center = [4.82,0.18,4.75]
      rotation = axis_rotation(normalized([1.0,-2.0,0.7]),0.619)
      symmetric_distortion = 0.0
      symmetric_distortion(1,1) = SQRT(1.0-0.04**2-0.02**2)
      symmetric_distortion(2,2) = SQRT(1.0-0.04**2-0.03**2)
      symmetric_distortion(3,3) = SQRT(1.0-0.02**2-0.03**2)
      symmetric_distortion(1,2) = 0.04
      symmetric_distortion(2,1) = 0.04
      symmetric_distortion(1,3) = -0.02
      symmetric_distortion(3,1) = -0.02
      symmetric_distortion(2,3) = 0.03
      symmetric_distortion(3,2) = 0.03
      raw_axes = MATMUL(rotation,symmetric_distortion)
      DO axis = 1, 3
         raw_axes(:,axis) = raw_axes(:,axis)/vector_norm(raw_axes(:,axis))
      END DO
      lengths = [0.92,1.04,1.08,0.97,1.15,1.01]
      DO axis = 1, 3
         displacements(:,2*axis-1) = lengths(2*axis-1)*raw_axes(:,axis)
         displacements(:,2*axis) = -lengths(2*axis)*raw_axes(:,axis)
      END DO
      DO ligand = 1, 6
         target = center+displacements(:,ligand)
         fractional = target/5.0
         translation = FLOOR(fractional)
         positions(:,ligand) = target-MATMUL(lattice,REAL(translation))
      END DO
      positions(:,7) = center+[1.7,0.2,0.1]
      positions(:,8) = center+[-0.1,-1.9,0.2]
      DO ligand = 7, 8
         fractional = positions(:,ligand)/5.0
         translation = FLOOR(fractional)
         positions(:,ligand) = positions(:,ligand)-MATMUL(lattice,REAL(translation))
      END DO
      ids = [(80+ligand,ligand=1,8)]

      CALL resolve_periodic_neighbours_cartesian(lattice,[.TRUE.,.TRUE.,.TRUE.],center,positions,ids,6,neighbors, &
                                                 coordination_diagnostics,success,message)
      CALL check_success("periodic distorted-octahedron resolution",success,message,n_failures)
      IF (.NOT.success) RETURN
      DO ligand = 1, 6
         resolved_displacements(:,ligand) = neighbors(ligand)%displacement
      END DO
      CALL construct_structural_frame_from_displacements(resolved_displacements,rotation,frame,frame_diagnostics, &
                                                         success,message)
      CALL check_success("coordination-to-structural-frame bridge",success,message,n_failures)
      IF (success) THEN
         CALL check_bound("known distorted structural frame",MAXVAL(ABS(frame-rotation)),2.0e-10, &
                          n_failures,largest_error)
         CALL check_true("distorted bridge retained nonorthogonal raw axes", &
                         MAXVAL(ABS(frame_diagnostics%raw_mutual_angles_deg-90.0)) > 0.1,n_failures)
      END IF
   END SUBROUTINE test_structural_frame_bridge

   SUBROUTINE compare_neighbor_lists(name,first,second,n_failures,largest_error)
      CHARACTER(LEN=*), INTENT(IN) :: name
      TYPE(t_periodic_neighbor), INTENT(IN) :: first(:), second(:)
      INTEGER, INTENT(INOUT) :: n_failures
      REAL, INTENT(INOUT) :: largest_error
      INTEGER :: i

      IF (SIZE(first) /= SIZE(second)) THEN
         WRITE(*,'(a,a)') "FAIL: ",TRIM(name)//" (different sizes)"
         n_failures = n_failures+1
         RETURN
      END IF
      DO i = 1, SIZE(first)
         IF (first(i)%atom_id /= second(i)%atom_id .OR. &
             ANY(first(i)%translation /= second(i)%translation)) THEN
            WRITE(*,'(a,a)') "FAIL: ",TRIM(name)//" (different identities)"
            n_failures = n_failures+1
            RETURN
         END IF
         CALL check_vector(TRIM(name)//" displacement",first(i)%displacement,second(i)%displacement, &
                           tolerance,n_failures,largest_error)
      END DO
   END SUBROUTINE compare_neighbor_lists

   INTEGER FUNCTION find_atom_id(neighbors,atom_id)
      TYPE(t_periodic_neighbor), INTENT(IN) :: neighbors(:)
      INTEGER, INTENT(IN) :: atom_id
      INTEGER :: i
      find_atom_id = 0
      DO i = 1, SIZE(neighbors)
         IF (neighbors(i)%atom_id == atom_id) THEN
            find_atom_id = i
            RETURN
         END IF
      END DO
   END FUNCTION find_atom_id

   SUBROUTINE check_success(name,success,message,n_failures)
      CHARACTER(LEN=*), INTENT(IN) :: name, message
      LOGICAL, INTENT(IN) :: success
      INTEGER, INTENT(INOUT) :: n_failures
      IF (success) THEN
         WRITE(*,'(a,a)') "PASS: ",TRIM(name)
      ELSE
         WRITE(*,'(a,a,a)') "FAIL: ",TRIM(name),": "//TRIM(message)
         n_failures = n_failures+1
      END IF
   END SUBROUTINE check_success

   SUBROUTINE check_true(name,condition,n_failures)
      CHARACTER(LEN=*), INTENT(IN) :: name
      LOGICAL, INTENT(IN) :: condition
      INTEGER, INTENT(INOUT) :: n_failures
      IF (condition) THEN
         WRITE(*,'(a,a)') "PASS: ",TRIM(name)
      ELSE
         WRITE(*,'(a,a)') "FAIL: ",TRIM(name)
         n_failures = n_failures+1
      END IF
   END SUBROUTINE check_true

   SUBROUTINE check_integer(name,value,expected,n_failures)
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER, INTENT(IN) :: value, expected
      INTEGER, INTENT(INOUT) :: n_failures
      IF (value == expected) THEN
         WRITE(*,'(a,a)') "PASS: ",TRIM(name)
      ELSE
         WRITE(*,'(a,a,2(a,i0))') "FAIL: ",TRIM(name)," value=",value," expected=",expected
         n_failures = n_failures+1
      END IF
   END SUBROUTINE check_integer

   SUBROUTINE check_integer_vector(name,value,expected,n_failures)
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER, INTENT(IN) :: value(3), expected(3)
      INTEGER, INTENT(INOUT) :: n_failures
      IF (ALL(value == expected)) THEN
         WRITE(*,'(a,a)') "PASS: ",TRIM(name)
      ELSE
         WRITE(*,'(a,a,2(a,3i5))') "FAIL: ",TRIM(name)," value=",value," expected=",expected
         n_failures = n_failures+1
      END IF
   END SUBROUTINE check_integer_vector

   SUBROUTINE check_real(name,value,expected,bound,n_failures,largest_error)
      CHARACTER(LEN=*), INTENT(IN) :: name
      REAL, INTENT(IN) :: value, expected, bound
      INTEGER, INTENT(INOUT) :: n_failures
      REAL, INTENT(INOUT) :: largest_error
      REAL :: error
      error = ABS(value-expected)
      largest_error = MAX(largest_error,error)
      IF (error <= bound) THEN
         WRITE(*,'(a,a,a,es12.4)') "PASS: ",TRIM(name)," error=",error
      ELSE
         WRITE(*,'(a,a,3(a,es12.4))') "FAIL: ",TRIM(name)," value=",value," expected=",expected," error=",error
         n_failures = n_failures+1
      END IF
   END SUBROUTINE check_real

   SUBROUTINE check_vector(name,value,expected,bound,n_failures,largest_error)
      CHARACTER(LEN=*), INTENT(IN) :: name
      REAL, INTENT(IN) :: value(3), expected(3), bound
      INTEGER, INTENT(INOUT) :: n_failures
      REAL, INTENT(INOUT) :: largest_error
      CALL check_bound(name,MAXVAL(ABS(value-expected)),bound,n_failures,largest_error)
   END SUBROUTINE check_vector

   SUBROUTINE check_bound(name,error,bound,n_failures,largest_error)
      CHARACTER(LEN=*), INTENT(IN) :: name
      REAL, INTENT(IN) :: error, bound
      INTEGER, INTENT(INOUT) :: n_failures
      REAL, INTENT(INOUT) :: largest_error
      largest_error = MAX(largest_error,error)
      IF (error <= bound) THEN
         WRITE(*,'(a,a,a,es12.4)') "PASS: ",TRIM(name)," error=",error
      ELSE
         WRITE(*,'(a,a,2(a,es12.4))') "FAIL: ",TRIM(name)," error=",error," bound=",bound
         n_failures = n_failures+1
      END IF
   END SUBROUTINE check_bound

   PURE FUNCTION identity_3() RESULT(identity)
      REAL :: identity(3,3)
      identity = 0.0
      identity(1,1) = 1.0
      identity(2,2) = 1.0
      identity(3,3) = 1.0
   END FUNCTION identity_3

   PURE REAL FUNCTION vector_norm(vector)
      REAL, INTENT(IN) :: vector(3)
      vector_norm = SQRT(MAX(0.0,DOT_PRODUCT(vector,vector)))
   END FUNCTION vector_norm

   PURE FUNCTION normalized(vector) RESULT(unit_vector)
      REAL, INTENT(IN) :: vector(3)
      REAL :: unit_vector(3)
      unit_vector = vector/vector_norm(vector)
   END FUNCTION normalized

   PURE FUNCTION axis_rotation(axis,angle) RESULT(rotation)
      REAL, INTENT(IN) :: axis(3), angle
      REAL :: rotation(3,3), unit_axis(3), cosine, sine
      unit_axis = normalized(axis)
      cosine = COS(angle)
      sine = SIN(angle)
      rotation = 0.0
      rotation(1,1) = cosine+unit_axis(1)**2*(1.0-cosine)
      rotation(2,2) = cosine+unit_axis(2)**2*(1.0-cosine)
      rotation(3,3) = cosine+unit_axis(3)**2*(1.0-cosine)
      rotation(1,2) = unit_axis(1)*unit_axis(2)*(1.0-cosine)-unit_axis(3)*sine
      rotation(1,3) = unit_axis(1)*unit_axis(3)*(1.0-cosine)+unit_axis(2)*sine
      rotation(2,1) = unit_axis(2)*unit_axis(1)*(1.0-cosine)+unit_axis(3)*sine
      rotation(2,3) = unit_axis(2)*unit_axis(3)*(1.0-cosine)-unit_axis(1)*sine
      rotation(3,1) = unit_axis(3)*unit_axis(1)*(1.0-cosine)-unit_axis(2)*sine
      rotation(3,2) = unit_axis(3)*unit_axis(2)*(1.0-cosine)+unit_axis(1)*sine
   END FUNCTION axis_rotation

END PROGRAM local_coordination_test
