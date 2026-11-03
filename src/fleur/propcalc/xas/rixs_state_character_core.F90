!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rixs_state_character_core
   USE m_juDFT, ONLY: juDFT_error
#ifdef CPP_HDF
   USE hdf5
   USE m_hdf_tools, ONLY: io_write_att
#endif
   IMPLICIT NONE
   PRIVATE

   INTEGER, PARAMETER, PUBLIC :: rixs_state_schema_version = 1
   INTEGER, PARAMETER, PUBLIC :: rixs_state_n_orbitals = 5
   INTEGER, PARAMETER, PUBLIC :: rixs_state_n_spins = 2
   INTEGER, PARAMETER, PUBLIC :: rixs_state_matrix_size = 10
   INTEGER, PARAMETER, PUBLIC :: rixs_state_n_ligands = 6
   INTEGER, PARAMETER, PUBLIC :: rixs_state_n_mj = 6
   INTEGER, PARAMETER, PUBLIC :: rixs_state_role_valence = 1
   INTEGER, PARAMETER, PUBLIC :: rixs_state_role_intermediate = 2
#ifdef CPP_HDF
   LOGICAL, PARAMETER, PUBLIC :: rixs_state_hdf_available = .TRUE.
#else
   LOGICAL, PARAMETER, PUBLIC :: rixs_state_hdf_available = .FALSE.
#endif

   TYPE, PUBLIC :: t_rixs_state_site
      INTEGER :: physical_atom = 0
      INTEGER :: atom_type = 0
      INTEGER :: iatom_l = 0
      INTEGER :: atomic_number = 0
      REAL :: local_to_global(3,3) = 0.0
      REAL :: reference_frame(3,3) = 0.0
      COMPLEX :: native_mt_to_global(2,2) = CMPLX(0.0,0.0)
      COMPLEX :: orbital_global_to_local(5,5) = CMPLX(0.0,0.0)
      COMPLEX :: spin_native_to_local(2,2) = CMPLX(0.0,0.0)
      COMPLEX :: combined_native_to_local(10,10) = CMPLX(0.0,0.0)
      INTEGER :: ligand_atom_ids(rixs_state_n_ligands) = 0
      INTEGER :: ligand_translations(3,rixs_state_n_ligands) = 0
      REAL :: ligand_displacements(3,rixs_state_n_ligands) = 0.0
      REAL :: bond_distances(rixs_state_n_ligands) = 0.0
      REAL :: next_neighbor_distance = HUGE(1.0)
      REAL :: shell_gap = HUGE(1.0)
      REAL :: next_to_last_ratio = HUGE(1.0)
      INTEGER :: opposite_pairs(2,3) = 0
      REAL :: opposite_pair_angles(3) = 0.0
      REAL :: opposite_pair_costs(3) = 0.0
      REAL :: total_pair_cost = 0.0
      REAL :: raw_axis_angles(3) = 0.0
      REAL :: raw_axes(3,3) = 0.0
      REAL :: raw_to_orthonormal_angles(3) = 0.0
      REAL :: raw_condition_number = 0.0
      REAL :: reference_alignment_score = 0.0
      REAL :: reference_alignment_gap = 0.0
   END TYPE t_rixs_state_site

   TYPE, PUBLIC :: t_rixs_state_record
      INTEGER :: ikpt = 0
      INTEGER :: band = 0
      INTEGER :: physical_atom = 0
      INTEGER :: atom_type = 0
      INTEGER :: iatom_l = 0
      INTEGER :: role_mask = 0
      REAL :: k_vector(3) = 0.0
      REAL :: k_weight = 0.0
      REAL :: band_energy = 0.0
      REAL :: occupation = 0.0
      REAL :: d_weight = 0.0
      REAL :: orbital_weights(rixs_state_n_orbitals) = 0.0
      REAL :: t2g_weight = 0.0
      REAL :: eg_weight = 0.0
      REAL :: jeff_weights(2) = 0.0
      REAL :: jeff_mj_weights(rixs_state_n_mj) = 0.0
      REAL :: spin_expectation(3) = 0.0
      REAL :: orbital_expectation(3) = 0.0
      COMPLEX :: rho_native(rixs_state_matrix_size,rixs_state_matrix_size) = CMPLX(0.0,0.0)
   END TYPE t_rixs_state_record

   TYPE, PUBLIC :: t_rixs_state_character_context
      LOGICAL :: enabled = .FALSE.
      INTEGER :: ligand_z = -1
      INTEGER :: global_rank = -1
      CHARACTER(LEN=512) :: filename = ""
      TYPE(t_rixs_state_site), ALLOCATABLE :: sites(:)
      TYPE(t_rixs_state_record), ALLOCATABLE :: records(:)
      INTEGER :: n_sites = 0
      INTEGER :: n_records = 0
   CONTAINS
      PROCEDURE :: initialize => rixs_state_context_initialize
      PROCEDURE :: add_site => rixs_state_context_add_site
      PROCEDURE :: add_record => rixs_state_context_add_record
      PROCEDURE :: has_record => rixs_state_context_has_record
      PROCEDURE :: find_site => rixs_state_context_find_site
      PROCEDURE :: write_hdf => rixs_state_context_write_hdf
      PROCEDURE :: clear => rixs_state_context_clear
   END TYPE t_rixs_state_character_context

   PUBLIC :: rixs_state_band_roles
   PUBLIC :: rixs_state_physical_atom
   PUBLIC :: rixs_state_select_atoms_by_z
   PUBLIC :: rixs_state_is_owner
   PUBLIC :: rixs_state_shard_filename

CONTAINS

   SUBROUTINE rixs_state_band_roles(nbands, valence_min, valence_max, intermediate_min, intermediate_max, roles)
      INTEGER, INTENT(IN) :: nbands, valence_min, valence_max, intermediate_min, intermediate_max
      INTEGER, INTENT(OUT) :: roles(nbands)
      INTEGER :: band

      roles = 0
      DO band = MAX(1,valence_min), MIN(nbands,valence_max)
         roles(band) = IOR(roles(band),rixs_state_role_valence)
      END DO
      DO band = MAX(1,intermediate_min), MIN(nbands,intermediate_max)
         roles(band) = IOR(roles(band),rixs_state_role_intermediate)
      END DO
   END SUBROUTINE rixs_state_band_roles

   INTEGER FUNCTION rixs_state_physical_atom(first_atom, neq, atom_types, itype, iatom_l) RESULT(iatom)
      INTEGER, INTENT(IN) :: first_atom(:), neq(:), atom_types(:), itype, iatom_l

      IF (itype < 1 .OR. itype > SIZE(first_atom) .OR. itype > SIZE(neq)) THEN
         CALL juDFT_error("RIXS state-character atom type is out of range", calledby="m_rixs_state_character_core")
      END IF
      IF (iatom_l < 1 .OR. iatom_l > neq(itype)) THEN
         CALL juDFT_error("RIXS state-character equivalent-atom index is out of range", &
                          calledby="m_rixs_state_character_core")
      END IF
      iatom = first_atom(itype)+iatom_l-1
      IF (iatom < 1 .OR. iatom > SIZE(atom_types) .OR. atom_types(iatom) /= itype) THEN
         CALL juDFT_error("RIXS state-character physical/equivalent atom mapping is inconsistent", &
                          calledby="m_rixs_state_character_core")
      END IF
   END FUNCTION rixs_state_physical_atom

   SUBROUTINE rixs_state_select_atoms_by_z(type_atomic_numbers, atom_types, requested_z, physical_atoms)
      INTEGER, INTENT(IN) :: type_atomic_numbers(:), atom_types(:), requested_z
      INTEGER, ALLOCATABLE, INTENT(OUT) :: physical_atoms(:)
      INTEGER :: atom, count_atoms

      IF (requested_z <= 0) THEN
         CALL juDFT_error("RIXS state-character ligand atomic number must be positive", &
                          calledby="m_rixs_state_character_core")
      END IF
      IF (ANY(atom_types < 1) .OR. ANY(atom_types > SIZE(type_atomic_numbers))) THEN
         CALL juDFT_error("RIXS state-character atom-type map is inconsistent", &
                          calledby="m_rixs_state_character_core")
      END IF
      count_atoms = 0
      DO atom = 1, SIZE(atom_types)
         IF (type_atomic_numbers(atom_types(atom)) == requested_z) count_atoms = count_atoms+1
      END DO
      ALLOCATE(physical_atoms(count_atoms))
      count_atoms = 0
      DO atom = 1, SIZE(atom_types)
         IF (type_atomic_numbers(atom_types(atom)) /= requested_z) CYCLE
         count_atoms = count_atoms+1
         physical_atoms(count_atoms) = atom
      END DO
   END SUBROUTINE rixs_state_select_atoms_by_z

   LOGICAL PURE FUNCTION rixs_state_is_owner(n_rank)
      INTEGER, INTENT(IN) :: n_rank
      rixs_state_is_owner = n_rank == 0
   END FUNCTION rixs_state_is_owner

   FUNCTION rixs_state_shard_filename(prefix, edge, global_rank) RESULT(filename)
      CHARACTER(LEN=*), INTENT(IN) :: prefix, edge
      INTEGER, INTENT(IN) :: global_rank
      CHARACTER(LEN=512) :: filename
      CHARACTER(LEN=16) :: rank_label

      WRITE(rank_label,'(i4.4)') global_rank
      filename = TRIM(prefix)//"_"//TRIM(edge)//"_state_character_rank"//TRIM(ADJUSTL(rank_label))//".hdf"
   END FUNCTION rixs_state_shard_filename

   SUBROUTINE rixs_state_context_initialize(this, enabled, ligand_z, global_rank, filename)
      CLASS(t_rixs_state_character_context), INTENT(INOUT) :: this
      LOGICAL, INTENT(IN) :: enabled
      INTEGER, INTENT(IN) :: ligand_z, global_rank
      CHARACTER(LEN=*), INTENT(IN) :: filename

      CALL this%clear()
      this%enabled = enabled
      this%ligand_z = ligand_z
      this%global_rank = global_rank
      this%filename = TRIM(filename)
      ALLOCATE(this%sites(0),this%records(0))
   END SUBROUTINE rixs_state_context_initialize

   SUBROUTINE rixs_state_context_add_site(this, site, site_index, added)
      CLASS(t_rixs_state_character_context), INTENT(INOUT) :: this
      TYPE(t_rixs_state_site), INTENT(IN) :: site
      INTEGER, INTENT(OUT) :: site_index
      LOGICAL, INTENT(OUT) :: added
      TYPE(t_rixs_state_site), ALLOCATABLE :: temporary(:)

      site_index = this%find_site(site%physical_atom)
      IF (site_index > 0) THEN
         added = .FALSE.
         RETURN
      END IF
      ALLOCATE(temporary(this%n_sites+1))
      IF (this%n_sites > 0) temporary(1:this%n_sites) = this%sites(1:this%n_sites)
      temporary(this%n_sites+1) = site
      CALL MOVE_ALLOC(temporary,this%sites)
      this%n_sites = this%n_sites+1
      site_index = this%n_sites
      added = .TRUE.
   END SUBROUTINE rixs_state_context_add_site

   INTEGER FUNCTION rixs_state_context_find_site(this, physical_atom) RESULT(site_index)
      CLASS(t_rixs_state_character_context), INTENT(IN) :: this
      INTEGER, INTENT(IN) :: physical_atom
      INTEGER :: site

      site_index = 0
      DO site = 1, this%n_sites
         IF (this%sites(site)%physical_atom == physical_atom) THEN
            site_index = site
            RETURN
         END IF
      END DO
   END FUNCTION rixs_state_context_find_site

   LOGICAL FUNCTION rixs_state_context_has_record(this, ikpt, band, physical_atom)
      CLASS(t_rixs_state_character_context), INTENT(IN) :: this
      INTEGER, INTENT(IN) :: ikpt, band, physical_atom
      INTEGER :: record

      rixs_state_context_has_record = .FALSE.
      DO record = 1, this%n_records
         IF (this%records(record)%ikpt == ikpt .AND. this%records(record)%band == band .AND. &
             this%records(record)%physical_atom == physical_atom) THEN
            rixs_state_context_has_record = .TRUE.
            RETURN
         END IF
      END DO
   END FUNCTION rixs_state_context_has_record

   SUBROUTINE rixs_state_context_add_record(this, record, added)
      CLASS(t_rixs_state_character_context), INTENT(INOUT) :: this
      TYPE(t_rixs_state_record), INTENT(IN) :: record
      LOGICAL, INTENT(OUT) :: added
      TYPE(t_rixs_state_record), ALLOCATABLE :: temporary(:)
      INTEGER :: old_capacity, new_capacity

      IF (this%has_record(record%ikpt,record%band,record%physical_atom)) THEN
         added = .FALSE.
         RETURN
      END IF
      old_capacity = SIZE(this%records)
      IF (this%n_records == old_capacity) THEN
         new_capacity = MAX(64,2*old_capacity)
         ALLOCATE(temporary(new_capacity))
         IF (this%n_records > 0) temporary(1:this%n_records) = this%records(1:this%n_records)
         CALL MOVE_ALLOC(temporary,this%records)
      END IF
      this%n_records = this%n_records+1
      this%records(this%n_records) = record
      added = .TRUE.
   END SUBROUTINE rixs_state_context_add_record

   SUBROUTINE rixs_state_context_write_hdf(this)
      CLASS(t_rixs_state_character_context), INTENT(IN) :: this
#ifdef CPP_HDF
      INTEGER(HID_T) :: file_id, state_group, site_group
      INTEGER :: hdf_error, record, site
      INTEGER, ALLOCATABLE :: state_identity(:,:), site_identity(:,:), ligand_atoms(:,:), ligand_translations(:,:,:)
      INTEGER, ALLOCATABLE :: opposite_pairs(:,:,:)
      REAL, ALLOCATABLE :: k_vectors(:,:), scalar_state(:,:), orbital_weights(:,:), jeff_weights(:,:), mj_weights(:,:)
      REAL, ALLOCATABLE :: spin_expectation(:,:), orbital_expectation(:,:), rho_real(:,:,:), rho_imag(:,:,:)
      REAL, ALLOCATABLE :: site_rotations(:,:,:), reference_frames(:,:,:), spin_real(:,:,:), spin_imag(:,:,:)
      REAL, ALLOCATABLE :: orbital_transform_real(:,:,:), orbital_transform_imag(:,:,:)
      REAL, ALLOCATABLE :: spin_to_local_real(:,:,:), spin_to_local_imag(:,:,:)
      REAL, ALLOCATABLE :: ligand_displacements(:,:,:), bond_distances(:,:), shell_data(:,:)
      REAL, ALLOCATABLE :: frame_angles(:,:), raw_axes(:,:,:), frame_quality(:,:), pair_costs(:,:)

      IF (.NOT.this%enabled) RETURN
      IF (this%n_sites <= 0) CALL juDFT_error("RIXS state-character output has no prepared absorber sites", &
                                              calledby="m_rixs_state_character_core")
      CALL h5open_f(hdf_error)
      CALL check_hdf("initialize HDF5",hdf_error)
      CALL h5fcreate_f(TRIM(this%filename),H5F_ACC_TRUNC_F,file_id,hdf_error)
      CALL check_hdf("create state-character HDF5 file",hdf_error)
      CALL io_write_att(file_id,"schema_version",rixs_state_schema_version)
      CALL io_write_att(file_id,"description", &
         "Band-state annotations only; not an incoherent decomposition of RIXS intensity")
      CALL io_write_att(file_id,"native_basis", &
         "complex Y_2m m=-2..2, each with native-MT spin components 1,2")
      CALL io_write_att(file_id,"rho_convention","unnormalized native local MT d-spin reduced density matrix")
      CALL io_write_att(file_id,"state_key","ikpt band physical_atom")
      CALL io_write_att(file_id,"geometry_units","bohr; Cartesian vectors use the FLEUR cell frame")
      CALL io_write_att(file_id,"global_rank",this%global_rank)
      CALL io_write_att(file_id,"ligand_atomic_number",this%ligand_z)
      CALL io_write_att(file_id,"number_of_states",this%n_records)
      CALL io_write_att(file_id,"number_of_sites",this%n_sites)
      CALL h5gcreate_f(file_id,"states",state_group,hdf_error)
      CALL check_hdf("create states group",hdf_error)
      CALL h5gcreate_f(file_id,"sites",site_group,hdf_error)
      CALL check_hdf("create sites group",hdf_error)

      ALLOCATE(state_identity(6,this%n_records),k_vectors(3,this%n_records),scalar_state(7,this%n_records))
      ALLOCATE(orbital_weights(5,this%n_records),jeff_weights(2,this%n_records),mj_weights(6,this%n_records))
      ALLOCATE(spin_expectation(3,this%n_records),orbital_expectation(3,this%n_records))
      ALLOCATE(rho_real(10,10,this%n_records),rho_imag(10,10,this%n_records))
      DO record = 1, this%n_records
         state_identity(:,record) = [this%records(record)%ikpt,this%records(record)%band, &
            this%records(record)%physical_atom,this%records(record)%atom_type, &
            this%records(record)%iatom_l,this%records(record)%role_mask]
         k_vectors(:,record) = this%records(record)%k_vector
         scalar_state(:,record) = [this%records(record)%k_weight,this%records(record)%band_energy, &
            this%records(record)%occupation,this%records(record)%d_weight,this%records(record)%t2g_weight, &
            this%records(record)%eg_weight,SUM(this%records(record)%jeff_weights)]
         orbital_weights(:,record) = this%records(record)%orbital_weights
         jeff_weights(:,record) = this%records(record)%jeff_weights
         mj_weights(:,record) = this%records(record)%jeff_mj_weights
         spin_expectation(:,record) = this%records(record)%spin_expectation
         orbital_expectation(:,record) = this%records(record)%orbital_expectation
         ! HDF5 reverses Fortran dimension order for C/Python readers. Store
         ! matrix transposes here so h5py exposes [state,row,column].
         rho_real(:,:,record) = TRANSPOSE(REAL(this%records(record)%rho_native))
         rho_imag(:,:,record) = TRANSPOSE(AIMAG(this%records(record)%rho_native))
      END DO
      CALL write_integer_2d(state_group,"identity",state_identity)
      CALL write_real_2d(state_group,"k_vector",k_vectors)
      CALL write_real_2d(state_group,"scalars",scalar_state)
      CALL write_real_2d(state_group,"real_orbital_weights",orbital_weights)
      CALL write_real_2d(state_group,"jeff_weights",jeff_weights)
      CALL write_real_2d(state_group,"jeff_mj_weights",mj_weights)
      CALL write_real_2d(state_group,"spin_expectation",spin_expectation)
      CALL write_real_2d(state_group,"orbital_expectation",orbital_expectation)
      CALL write_real_3d(state_group,"rho_native_real",rho_real)
      CALL write_real_3d(state_group,"rho_native_imag",rho_imag)
      CALL io_write_att(state_group,"identity_columns","ikpt band physical_atom atom_type iatom_l role_mask")
      CALL io_write_att(state_group,"k_vector_convention","fractional reciprocal-lattice coordinates")
      CALL io_write_att(state_group,"scalar_rows","k_weight energy_Ha occupation d_weight t2g_weight eg_weight t2g_check")
      CALL io_write_att(state_group,"orbital_rows","d_xy d_yz d_xz d_x2-y2 d_z2")
      CALL io_write_att(state_group,"jeff_rows","j_eff_1/2 j_eff_3/2")
      CALL io_write_att(state_group,"mj_rows","1/2:-1/2,+1/2; 3/2:-3/2,-1/2,+1/2,+3/2")
      CALL io_write_att(state_group,"expectation_units","hbar=1; components use the local structural frame")

      ALLOCATE(site_identity(4,this%n_sites),site_rotations(3,3,this%n_sites),reference_frames(3,3,this%n_sites))
      ALLOCATE(spin_real(2,2,this%n_sites),spin_imag(2,2,this%n_sites),ligand_atoms(6,this%n_sites))
      ALLOCATE(orbital_transform_real(5,5,this%n_sites),orbital_transform_imag(5,5,this%n_sites))
      ALLOCATE(spin_to_local_real(2,2,this%n_sites),spin_to_local_imag(2,2,this%n_sites))
      ALLOCATE(ligand_translations(3,6,this%n_sites),ligand_displacements(3,6,this%n_sites))
      ALLOCATE(bond_distances(6,this%n_sites),shell_data(3,this%n_sites),frame_angles(9,this%n_sites))
      ALLOCATE(raw_axes(3,3,this%n_sites),frame_quality(3,this%n_sites),opposite_pairs(2,3,this%n_sites))
      ALLOCATE(pair_costs(4,this%n_sites))
      DO site = 1, this%n_sites
         site_identity(:,site) = [this%sites(site)%physical_atom,this%sites(site)%atom_type, &
                                  this%sites(site)%iatom_l,this%sites(site)%atomic_number]
         site_rotations(:,:,site) = TRANSPOSE(this%sites(site)%local_to_global)
         reference_frames(:,:,site) = TRANSPOSE(this%sites(site)%reference_frame)
         spin_real(:,:,site) = TRANSPOSE(REAL(this%sites(site)%native_mt_to_global))
         spin_imag(:,:,site) = TRANSPOSE(AIMAG(this%sites(site)%native_mt_to_global))
         orbital_transform_real(:,:,site) = TRANSPOSE(REAL(this%sites(site)%orbital_global_to_local))
         orbital_transform_imag(:,:,site) = TRANSPOSE(AIMAG(this%sites(site)%orbital_global_to_local))
         spin_to_local_real(:,:,site) = TRANSPOSE(REAL(this%sites(site)%spin_native_to_local))
         spin_to_local_imag(:,:,site) = TRANSPOSE(AIMAG(this%sites(site)%spin_native_to_local))
         ligand_atoms(:,site) = this%sites(site)%ligand_atom_ids
         ligand_translations(:,:,site) = this%sites(site)%ligand_translations
         ligand_displacements(:,:,site) = this%sites(site)%ligand_displacements
         bond_distances(:,site) = this%sites(site)%bond_distances
         shell_data(:,site) = [this%sites(site)%next_neighbor_distance,this%sites(site)%shell_gap, &
                               this%sites(site)%next_to_last_ratio]
         opposite_pairs(:,:,site) = this%sites(site)%opposite_pairs
         pair_costs(:,site) = [this%sites(site)%opposite_pair_costs,this%sites(site)%total_pair_cost]
         frame_angles(:,site) = [this%sites(site)%opposite_pair_angles,this%sites(site)%raw_axis_angles, &
                                 this%sites(site)%raw_to_orthonormal_angles]
         raw_axes(:,:,site) = TRANSPOSE(this%sites(site)%raw_axes)
         frame_quality(:,site) = [this%sites(site)%raw_condition_number, &
            this%sites(site)%reference_alignment_score,this%sites(site)%reference_alignment_gap]
      END DO
      CALL write_integer_2d(site_group,"identity",site_identity)
      CALL write_real_3d(site_group,"local_to_global",site_rotations)
      CALL write_real_3d(site_group,"reference_frame",reference_frames)
      CALL write_real_3d(site_group,"native_mt_to_global_real",spin_real)
      CALL write_real_3d(site_group,"native_mt_to_global_imag",spin_imag)
      CALL write_real_3d(site_group,"orbital_global_to_local_real",orbital_transform_real)
      CALL write_real_3d(site_group,"orbital_global_to_local_imag",orbital_transform_imag)
      CALL write_real_3d(site_group,"spin_native_to_local_real",spin_to_local_real)
      CALL write_real_3d(site_group,"spin_native_to_local_imag",spin_to_local_imag)
      CALL write_integer_2d(site_group,"ligand_atom_ids",ligand_atoms)
      CALL write_integer_3d(site_group,"ligand_translations",ligand_translations)
      CALL write_real_3d(site_group,"ligand_displacements",ligand_displacements)
      CALL write_real_2d(site_group,"bond_distances",bond_distances)
      CALL write_real_2d(site_group,"shell_diagnostics",shell_data)
      CALL write_integer_3d(site_group,"opposite_pairs",opposite_pairs)
      CALL write_real_2d(site_group,"opposite_pair_costs",pair_costs)
      CALL write_real_2d(site_group,"frame_angles_degrees",frame_angles)
      CALL write_real_3d(site_group,"raw_axes",raw_axes)
      CALL write_real_2d(site_group,"frame_quality",frame_quality)
      CALL io_write_att(site_group,"identity_columns","physical_atom atom_type iatom_l atomic_number")
      CALL io_write_att(site_group,"translation_convention", &
         "integer coefficients of the cell lattice-vector columns applied to each ligand atom")
      CALL io_write_att(site_group,"shell_rows","d7 d7_minus_d6 d7_over_d6")
      CALL io_write_att(site_group,"pair_cost_rows","pair_x pair_y pair_z total")
      CALL io_write_att(site_group,"frame_angle_rows", &
         "opposite_pair_x,y,z raw_xy,xz,yz raw_to_final_x,y,z")
      CALL io_write_att(site_group,"frame_quality_rows","condition_number alignment_score alignment_gap")
      CALL io_write_att(site_group,"frame_convention", &
         "R_loc_to_glob columns are local structural axes in global Cartesian coordinates")

      CALL h5gclose_f(state_group,hdf_error)
      CALL check_hdf("close states group",hdf_error)
      CALL h5gclose_f(site_group,hdf_error)
      CALL check_hdf("close sites group",hdf_error)
      CALL h5fclose_f(file_id,hdf_error)
      CALL check_hdf("close state-character HDF5 file",hdf_error)
      CALL h5close_f(hdf_error)
      CALL check_hdf("finalize HDF5",hdf_error)
#else
      IF (this%enabled) CALL juDFT_error("RIXS state-character output requires an HDF5-enabled FLEUR build", &
                                         calledby="m_rixs_state_character_core")
#endif
   END SUBROUTINE rixs_state_context_write_hdf

   SUBROUTINE rixs_state_context_clear(this)
      CLASS(t_rixs_state_character_context), INTENT(INOUT) :: this

      IF (ALLOCATED(this%sites)) DEALLOCATE(this%sites)
      IF (ALLOCATED(this%records)) DEALLOCATE(this%records)
      this%enabled = .FALSE.
      this%ligand_z = -1
      this%global_rank = -1
      this%filename = ""
      this%n_sites = 0
      this%n_records = 0
   END SUBROUTINE rixs_state_context_clear

#ifdef CPP_HDF
   SUBROUTINE check_hdf(operation,error)
      CHARACTER(LEN=*), INTENT(IN) :: operation
      INTEGER, INTENT(IN) :: error
      IF (error /= 0) CALL juDFT_error("Cannot "//TRIM(operation)//" for RIXS state-character sidecar", &
                                       calledby="m_rixs_state_character_core")
   END SUBROUTINE check_hdf

   SUBROUTINE write_integer_2d(group_id,name,data)
      INTEGER(HID_T), INTENT(IN) :: group_id
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER, INTENT(IN) :: data(:,:)
      INTEGER(HID_T) :: space_id, dataset_id
      INTEGER(HSIZE_T) :: dimensions(2)
      INTEGER :: error
      dimensions = SHAPE(data)
      CALL h5screate_simple_f(2,dimensions,space_id,error); CALL check_hdf("create "//name//" dataspace",error)
      CALL h5dcreate_f(group_id,name,H5T_NATIVE_INTEGER,space_id,dataset_id,error); CALL check_hdf("create "//name,error)
      CALL h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER,data,dimensions,error); CALL check_hdf("write "//name,error)
      CALL h5dclose_f(dataset_id,error); CALL check_hdf("close "//name,error)
      CALL h5sclose_f(space_id,error); CALL check_hdf("close "//name//" dataspace",error)
   END SUBROUTINE write_integer_2d

   SUBROUTINE write_integer_3d(group_id,name,data)
      INTEGER(HID_T), INTENT(IN) :: group_id
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER, INTENT(IN) :: data(:,:,:)
      INTEGER(HID_T) :: space_id, dataset_id
      INTEGER(HSIZE_T) :: dimensions(3)
      INTEGER :: error
      dimensions = SHAPE(data)
      CALL h5screate_simple_f(3,dimensions,space_id,error); CALL check_hdf("create "//name//" dataspace",error)
      CALL h5dcreate_f(group_id,name,H5T_NATIVE_INTEGER,space_id,dataset_id,error); CALL check_hdf("create "//name,error)
      CALL h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER,data,dimensions,error); CALL check_hdf("write "//name,error)
      CALL h5dclose_f(dataset_id,error); CALL check_hdf("close "//name,error)
      CALL h5sclose_f(space_id,error); CALL check_hdf("close "//name//" dataspace",error)
   END SUBROUTINE write_integer_3d

   SUBROUTINE write_real_2d(group_id,name,data)
      INTEGER(HID_T), INTENT(IN) :: group_id
      CHARACTER(LEN=*), INTENT(IN) :: name
      REAL, INTENT(IN) :: data(:,:)
      INTEGER(HID_T) :: space_id, dataset_id
      INTEGER(HSIZE_T) :: dimensions(2)
      INTEGER :: error
      dimensions = SHAPE(data)
      CALL h5screate_simple_f(2,dimensions,space_id,error); CALL check_hdf("create "//name//" dataspace",error)
      CALL h5dcreate_f(group_id,name,H5T_NATIVE_DOUBLE,space_id,dataset_id,error); CALL check_hdf("create "//name,error)
      CALL h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,data,dimensions,error); CALL check_hdf("write "//name,error)
      CALL h5dclose_f(dataset_id,error); CALL check_hdf("close "//name,error)
      CALL h5sclose_f(space_id,error); CALL check_hdf("close "//name//" dataspace",error)
   END SUBROUTINE write_real_2d

   SUBROUTINE write_real_3d(group_id,name,data)
      INTEGER(HID_T), INTENT(IN) :: group_id
      CHARACTER(LEN=*), INTENT(IN) :: name
      REAL, INTENT(IN) :: data(:,:,:)
      INTEGER(HID_T) :: space_id, dataset_id
      INTEGER(HSIZE_T) :: dimensions(3)
      INTEGER :: error
      dimensions = SHAPE(data)
      CALL h5screate_simple_f(3,dimensions,space_id,error); CALL check_hdf("create "//name//" dataspace",error)
      CALL h5dcreate_f(group_id,name,H5T_NATIVE_DOUBLE,space_id,dataset_id,error); CALL check_hdf("create "//name,error)
      CALL h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,data,dimensions,error); CALL check_hdf("write "//name,error)
      CALL h5dclose_f(dataset_id,error); CALL check_hdf("close "//name,error)
      CALL h5sclose_f(space_id,error); CALL check_hdf("close "//name//" dataspace",error)
   END SUBROUTINE write_real_3d
#endif

END MODULE m_rixs_state_character_core
