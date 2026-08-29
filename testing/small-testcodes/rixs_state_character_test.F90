!--------------------------------------------------------------------------------
! Standalone validation of RIXS state-character selection, identity caches, and
! the HDF5 sidecar schema. It does not invoke RIXS or any spectral calculation.
!--------------------------------------------------------------------------------

PROGRAM rixs_state_character_test
   USE m_rixs_state_character_core, ONLY: t_rixs_state_character_context, t_rixs_state_record, &
                                           t_rixs_state_site, rixs_state_band_roles, &
                                           rixs_state_is_owner, rixs_state_physical_atom, &
                                           rixs_state_select_atoms_by_z, rixs_state_shard_filename
   IMPLICIT NONE

   TYPE(t_rixs_state_character_context) :: context
   TYPE(t_rixs_state_site) :: site
   TYPE(t_rixs_state_record) :: record
   INTEGER, ALLOCATABLE :: selected_atoms(:)
   INTEGER :: roles(8), failures, site_index, row, column
   LOGICAL :: added
   CHARACTER(LEN=512) :: filename

   failures = 0
   CALL rixs_state_band_roles(8,2,5,5,7,roles)
   CALL check_integer_array("band role-mask union",roles,[0,1,1,1,3,2,2,0],failures)
   CALL check_integer("equivalent-to-physical mapping", &
      rixs_state_physical_atom([1,3],[2,3],[1,1,2,2,2],2,2),4,failures)
   CALL rixs_state_select_atoms_by_z([77,8,8],[1,2,3,2,1,3],8,selected_atoms)
   CALL check_integer_array("multiple ligand-type expansion",selected_atoms,[2,3,4,6],failures)
   CALL check_true("subgroup root owns output",rixs_state_is_owner(0),failures)
   CALL check_true("subgroup collaborator does not own output",.NOT.rixs_state_is_owner(1),failures)
   CALL check_true("rank shard naming",TRIM(rixs_state_shard_filename("test","L3",7)) == &
                   "test_L3_state_character_rank0007.hdf",failures)

   CALL get_command_argument(1,filename)
   IF (LEN_TRIM(filename) == 0) filename = "rixs_state_character_test.hdf"
   CALL context%initialize(.TRUE.,8,7,TRIM(filename))
   site%physical_atom = 4
   site%atom_type = 2
   site%iatom_l = 2
   site%atomic_number = 77
   site%local_to_global = identity_3()
   site%reference_frame = identity_3()
   site%native_mt_to_global(1,1) = CMPLX(1.0,0.0)
   site%native_mt_to_global(2,2) = CMPLX(1.0,0.0)
   site%orbital_global_to_local = identity_complex(5)
   site%spin_native_to_local = identity_complex(2)
   site%combined_native_to_local = identity_complex(10)
   site%ligand_atom_ids = [11,12,13,14,15,16]
   site%ligand_translations = RESHAPE([(row,row=1,18)],[3,6])
   site%ligand_displacements = 0.1*REAL(site%ligand_translations)
   site%bond_distances = [1.0,1.1,1.2,1.3,1.4,1.5]
   site%next_neighbor_distance = 1.8
   site%shell_gap = 0.3
   site%next_to_last_ratio = 1.2
   site%opposite_pairs = RESHAPE([1,2,3,4,5,6],[2,3])
   site%opposite_pair_angles = [178.0,176.0,174.0]
   site%opposite_pair_costs = [0.01,0.02,0.03]
   site%total_pair_cost = 0.06
   site%raw_axis_angles = [89.0,91.0,88.0]
   site%raw_axes = identity_3()
   site%raw_to_orthonormal_angles = [0.1,0.2,0.3]
   site%raw_condition_number = 1.1
   site%reference_alignment_score = 2.9
   site%reference_alignment_gap = 0.4
   CALL context%add_site(site,site_index,added)
   CALL check_true("first site cache insertion",added .AND. site_index == 1,failures)
   CALL context%add_site(site,site_index,added)
   CALL check_true("duplicate site cache suppression",.NOT.added .AND. context%n_sites == 1,failures)

   record%ikpt = 3
   record%band = 5
   record%physical_atom = 4
   record%atom_type = 2
   record%iatom_l = 2
   record%role_mask = 3
   record%k_vector = [0.1,0.2,0.3]
   record%k_weight = 0.04
   record%band_energy = 1.25
   record%occupation = 0.75
   record%d_weight = 0.9
   record%orbital_weights = [0.1,0.2,0.15,0.25,0.2]
   record%t2g_weight = 0.45
   record%eg_weight = 0.45
   record%jeff_weights = [0.3,0.15]
   record%jeff_mj_weights = [0.1,0.2,0.01,0.02,0.03,0.09]
   record%spin_expectation = [0.01,0.02,0.03]
   record%orbital_expectation = [0.04,0.05,0.06]
   DO column = 1, 10
      DO row = 1, 10
         record%rho_native(row,column) = CMPLX(0.01*REAL(row+10*column),0.001*REAL(row-column))
      END DO
   END DO
   CALL context%add_record(record,added)
   CALL check_true("first state cache insertion",added .AND. context%n_records == 1,failures)
   CALL context%add_record(record,added)
   CALL check_true("duplicate state-key suppression",.NOT.added .AND. context%n_records == 1,failures)
   CALL context%write_hdf()

   IF (failures /= 0) THEN
      WRITE(*,'(a,i0)') "RIXS STATE CHARACTER TEST: FAIL; failures = ",failures
      ERROR STOP 1
   END IF
   WRITE(*,'(a)') "RIXS STATE CHARACTER TEST: PASS"

CONTAINS

   SUBROUTINE check_integer(name,value,expected,n_failures)
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER, INTENT(IN) :: value, expected
      INTEGER, INTENT(INOUT) :: n_failures
      IF (value /= expected) THEN
         WRITE(*,'(a,a,2(1x,i0))') "FAIL: ",TRIM(name),value,expected
         n_failures = n_failures+1
      END IF
   END SUBROUTINE check_integer

   SUBROUTINE check_integer_array(name,value,expected,n_failures)
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER, INTENT(IN) :: value(:), expected(:)
      INTEGER, INTENT(INOUT) :: n_failures
      IF (SIZE(value) /= SIZE(expected) .OR. ANY(value /= expected)) THEN
         WRITE(*,'(a,a)') "FAIL: ",TRIM(name)
         n_failures = n_failures+1
      END IF
   END SUBROUTINE check_integer_array

   SUBROUTINE check_true(name,value,n_failures)
      CHARACTER(LEN=*), INTENT(IN) :: name
      LOGICAL, INTENT(IN) :: value
      INTEGER, INTENT(INOUT) :: n_failures
      IF (.NOT.value) THEN
         WRITE(*,'(a,a)') "FAIL: ",TRIM(name)
         n_failures = n_failures+1
      END IF
   END SUBROUTINE check_true

   PURE FUNCTION identity_3() RESULT(matrix)
      REAL :: matrix(3,3)
      matrix = 0.0
      matrix(1,1) = 1.0
      matrix(2,2) = 1.0
      matrix(3,3) = 1.0
   END FUNCTION identity_3

   PURE FUNCTION identity_complex(size_matrix) RESULT(matrix)
      INTEGER, INTENT(IN) :: size_matrix
      COMPLEX :: matrix(size_matrix,size_matrix)
      INTEGER :: index
      matrix = CMPLX(0.0,0.0)
      DO index = 1, size_matrix
         matrix(index,index) = CMPLX(1.0,0.0)
      END DO
   END FUNCTION identity_complex

END PROGRAM rixs_state_character_test
