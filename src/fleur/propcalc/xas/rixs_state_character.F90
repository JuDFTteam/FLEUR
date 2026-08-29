!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rixs_state_character
   USE m_juDFT, ONLY: juDFT_error
   USE m_local_coordination, ONLY: resolve_fleur_periodic_neighbours, t_coordination_diagnostics, &
                                   t_periodic_neighbor
   USE m_local_d_spin_analysis, ONLY: local_d_character_weights, local_d_extract_t2g, &
                                      local_d_jeff_weights, local_d_orbital_expectation, &
                                      local_d_real_angular_momentum, local_d_spin_expectation, &
                                      local_d_transform_to_real
   USE m_local_d_spin_frame_transform, ONLY: local_d_spin_native_to_structural, &
                                              transform_local_d_spin_density_with_matrix
   USE m_local_d_spin_projector, ONLY: local_d_spin_density_matrix
   USE m_local_d_spin_projector_core, ONLY: local_d_l
   USE m_local_structural_frame, ONLY: closest_proper_rotation, construct_structural_frame_from_displacements, &
                                       t_structural_frame_diagnostics
   USE m_rixs_state_character_core, ONLY: t_rixs_state_character_context, t_rixs_state_record, &
                                           t_rixs_state_site, rixs_state_band_roles, &
                                           rixs_state_matrix_size, rixs_state_n_ligands, &
                                           rixs_state_physical_atom, rixs_state_select_atoms_by_z, &
                                           rixs_state_shard_filename
   USE m_types_abc, ONLY: t_abc
   USE m_types_atoms, ONLY: t_atoms
   USE m_types_cell, ONLY: t_cell
   USE m_types_nococonv, ONLY: t_nococonv
   USE m_types_radfun, ONLY: t_radfun
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: t_rixs_state_character_context
   PUBLIC :: rixs_prepare_state_character_context
   PUBLIC :: rixs_characterize_site_bands
   PUBLIC :: rixs_finalize_state_character_context

CONTAINS

   SUBROUTINE rixs_prepare_state_character_context(context, ligand_z, output_prefix, edge, global_rank, &
                                                   absorber_z, atoms, cell, film, nococonv)
      TYPE(t_rixs_state_character_context), INTENT(INOUT) :: context
      INTEGER, INTENT(IN) :: ligand_z, global_rank, absorber_z
      CHARACTER(LEN=*), INTENT(IN) :: output_prefix, edge
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_cell), INTENT(IN) :: cell
      LOGICAL, INTENT(IN) :: film
      TYPE(t_nococonv), INTENT(IN) :: nococonv

      TYPE(t_periodic_neighbor), ALLOCATABLE :: neighbors(:)
      TYPE(t_coordination_diagnostics) :: coordination
      TYPE(t_structural_frame_diagnostics) :: frame_diagnostics
      TYPE(t_rixs_state_site) :: site_record
      INTEGER, ALLOCATABLE :: ligand_atoms(:)
      REAL :: normalized_lattice(3,3), reference_frame(3,3), displacements(3,rixs_state_n_ligands)
      REAL :: singular_values(3), condition_number, norm_value
      COMPLEX :: orbital_transform(-local_d_l:local_d_l,-local_d_l:local_d_l)
      COMPLEX :: spin_transform(2,2), combined_transform(10,10)
      CHARACTER(LEN=512) :: filename
      CHARACTER(LEN=256) :: message
      INTEGER :: axis, atom, itype, iatom_l, ligand, site_index
      LOGICAL :: success, reflection_corrected, added

#ifndef CPP_HDF
      CALL juDFT_error("RIXS writeStateCharacter requires an HDF5-enabled FLEUR build", &
                       calledby="m_rixs_state_character")
#endif
      filename = rixs_state_shard_filename(output_prefix,edge,global_rank)
      CALL context%initialize(.TRUE.,ligand_z,global_rank,filename)
      CALL rixs_state_select_atoms_by_z(atoms%nz,atoms%itype,ligand_z,ligand_atoms)
      IF (SIZE(ligand_atoms) == 0) THEN
         CALL juDFT_error("No physical ligand atoms match the requested RIXS state-character ligand Z", &
                          calledby="m_rixs_state_character")
      END IF

      normalized_lattice = cell%amat
      DO axis = 1, 3
         norm_value = SQRT(DOT_PRODUCT(normalized_lattice(:,axis),normalized_lattice(:,axis)))
         IF (norm_value <= 0.0) CALL juDFT_error("Cannot construct a crystallographic reference from the lattice", &
                                                 calledby="m_rixs_state_character")
         normalized_lattice(:,axis) = normalized_lattice(:,axis)/norm_value
      END DO
      CALL closest_proper_rotation(normalized_lattice,reference_frame,singular_values,condition_number, &
                                   reflection_corrected,success,message)
      IF (.NOT.success) CALL juDFT_error("Cannot construct RIXS state-character reference frame: "//TRIM(message), &
                                         calledby="m_rixs_state_character")

      DO atom = 1, atoms%nat
         itype = atoms%itype(atom)
         IF (atoms%nz(itype) /= absorber_z) CYCLE
         iatom_l = atom-atoms%firstAtom(itype)+1
         IF (rixs_state_physical_atom(atoms%firstAtom,atoms%neq,atoms%itype,itype,iatom_l) /= atom) THEN
            CALL juDFT_error("Unexpected absorber-site mapping in RIXS state-character preparation", &
                             calledby="m_rixs_state_character")
         END IF
         CALL resolve_fleur_periodic_neighbours(cell,atoms,film,atom,ligand_atoms,rixs_state_n_ligands, &
                                                neighbors,coordination,success,message)
         IF (.NOT.success) THEN
            CALL juDFT_error("Cannot resolve RIXS state-character ligand shell: "//TRIM(message), &
                             calledby="m_rixs_state_character")
         END IF
         DO ligand = 1, rixs_state_n_ligands
            displacements(:,ligand) = neighbors(ligand)%displacement
         END DO
         CALL construct_structural_frame_from_displacements(displacements,reference_frame, &
                                                            site_record%local_to_global,frame_diagnostics, &
                                                            success,message)
         IF (.NOT.success) THEN
            CALL juDFT_error("Cannot construct RIXS state-character structural frame: "//TRIM(message), &
                             calledby="m_rixs_state_character")
         END IF

         site_record%physical_atom = atom
         site_record%atom_type = itype
         site_record%iatom_l = iatom_l
         site_record%atomic_number = atoms%nz(itype)
         site_record%reference_frame = reference_frame
         site_record%native_mt_to_global = nococonv%umat(itype)
         CALL local_d_spin_native_to_structural(site_record%local_to_global,site_record%native_mt_to_global, &
                                                orbital_transform,spin_transform,combined_transform)
         site_record%orbital_global_to_local = orbital_transform
         site_record%spin_native_to_local = spin_transform
         site_record%combined_native_to_local = combined_transform
         DO ligand = 1, rixs_state_n_ligands
            site_record%ligand_atom_ids(ligand) = neighbors(ligand)%atom_id
            site_record%ligand_translations(:,ligand) = neighbors(ligand)%translation
            site_record%ligand_displacements(:,ligand) = neighbors(ligand)%displacement
            site_record%bond_distances(ligand) = neighbors(ligand)%distance
         END DO
         site_record%next_neighbor_distance = coordination%next_neighbor_distance
         site_record%shell_gap = coordination%shell_gap
         site_record%next_to_last_ratio = coordination%next_to_last_ratio
         site_record%opposite_pairs = frame_diagnostics%opposite_pairs
         site_record%opposite_pair_angles = frame_diagnostics%opposite_pair_angles_deg
         site_record%opposite_pair_costs = frame_diagnostics%opposite_pair_costs
         site_record%total_pair_cost = frame_diagnostics%total_pair_cost
         site_record%raw_axis_angles = frame_diagnostics%raw_mutual_angles_deg
         site_record%raw_axes = frame_diagnostics%raw_axes
         site_record%raw_to_orthonormal_angles = frame_diagnostics%raw_to_orthonormal_angles_deg
         site_record%raw_condition_number = frame_diagnostics%condition_number
         site_record%reference_alignment_score = frame_diagnostics%reference_alignment_score
         site_record%reference_alignment_gap = frame_diagnostics%reference_alignment_gap
         CALL context%add_site(site_record,site_index,added)
         IF (.NOT.added) CALL juDFT_error("Duplicate absorber site in RIXS state-character preparation", &
                                          calledby="m_rixs_state_character")
      END DO
      IF (context%n_sites == 0) CALL juDFT_error("RIXS state-character preparation found no absorber sites", &
                                                 calledby="m_rixs_state_character")
   END SUBROUTINE rixs_prepare_state_character_context

   SUBROUTINE rixs_characterize_site_bands(context, abc_spin, radfun, atoms, itype, iatom_l, ikpt, &
                                           k_vector, k_weight, eig_band, occ_band, valence_min, valence_max, &
                                           intermediate_min, intermediate_max)
      TYPE(t_rixs_state_character_context), INTENT(INOUT) :: context
      TYPE(t_abc), INTENT(IN) :: abc_spin(:)
      TYPE(t_radfun), INTENT(IN) :: radfun
      TYPE(t_atoms), INTENT(IN) :: atoms
      INTEGER, INTENT(IN) :: itype, iatom_l, ikpt
      REAL, INTENT(IN) :: k_vector(3), k_weight, eig_band(:), occ_band(:)
      INTEGER, INTENT(IN) :: valence_min, valence_max, intermediate_min, intermediate_max

      TYPE(t_rixs_state_record) :: record
      COMPLEX :: rho_native(-2:2,2,-2:2,2), rho_structural(-2:2,2,-2:2,2)
      COMPLEX :: rho_real(5,2,5,2), rho_t2g(3,2,3,2)
      COMPLEX :: physical_l_real(5,5,3), physical_l_t2g(3,3,3), effective_l_t2g(3,3,3)
      INTEGER, ALLOCATABLE :: roles(:)
      INTEGER :: atom, band, site_index
      LOGICAL :: added

      IF (.NOT.context%enabled) RETURN
      IF (SIZE(eig_band) /= SIZE(occ_band)) CALL juDFT_error("RIXS state-character eigenvalue/occupation size mismatch", &
                                                            calledby="m_rixs_state_character")
      atom = rixs_state_physical_atom(atoms%firstAtom,atoms%neq,atoms%itype,itype,iatom_l)
      site_index = context%find_site(atom)
      IF (site_index == 0) CALL juDFT_error("RIXS state-character site was not prepared", &
                                           calledby="m_rixs_state_character")
      ALLOCATE(roles(SIZE(eig_band)))
      CALL rixs_state_band_roles(SIZE(eig_band),valence_min,valence_max,intermediate_min,intermediate_max,roles)
      CALL local_d_real_angular_momentum(physical_l_real,physical_l_t2g,effective_l_t2g)

      DO band = 1, SIZE(eig_band)
         IF (roles(band) == 0) CYCLE
         IF (context%has_record(ikpt,band,atom)) CYCLE
         CALL local_d_spin_density_matrix(abc_spin,radfun,atoms,itype,band,iatom_l,rho_native)
         CALL transform_local_d_spin_density_with_matrix(rho_native, &
            context%sites(site_index)%combined_native_to_local,rho_structural)
         CALL local_d_transform_to_real(rho_structural,rho_real)
         CALL local_d_extract_t2g(rho_real,rho_t2g)

         record = t_rixs_state_record()
         record%ikpt = ikpt
         record%band = band
         record%physical_atom = atom
         record%atom_type = itype
         record%iatom_l = iatom_l
         record%role_mask = roles(band)
         record%k_vector = k_vector
         record%k_weight = k_weight
         record%band_energy = eig_band(band)
         record%occupation = occ_band(band)
         CALL local_d_character_weights(rho_real,record%orbital_weights,record%t2g_weight, &
                                        record%eg_weight,record%d_weight)
         CALL local_d_jeff_weights(rho_t2g,record%jeff_mj_weights,record%jeff_weights)
         CALL local_d_spin_expectation(rho_real,record%spin_expectation)
         CALL local_d_orbital_expectation(rho_real,physical_l_real,record%orbital_expectation)
         CALL flatten_native_density(rho_native,record%rho_native)
         CALL validate_state_record(record,rho_native)
         CALL context%add_record(record,added)
         IF (.NOT.added) CALL juDFT_error("Duplicate RIXS state-character record key", &
                                          calledby="m_rixs_state_character")
      END DO
   END SUBROUTINE rixs_characterize_site_bands

   SUBROUTINE rixs_finalize_state_character_context(context)
      TYPE(t_rixs_state_character_context), INTENT(INOUT) :: context

      IF (.NOT.context%enabled) RETURN
      CALL context%write_hdf()
      WRITE(*,'(a,i0,a,a)') "RIXS state character wrote ",context%n_records," records to ",TRIM(context%filename)
      CALL context%clear()
   END SUBROUTINE rixs_finalize_state_character_context

   SUBROUTINE flatten_native_density(rho_native,rho_flat)
      COMPLEX, INTENT(IN) :: rho_native(-2:2,2,-2:2,2)
      COMPLEX, INTENT(OUT) :: rho_flat(rixs_state_matrix_size,rixs_state_matrix_size)
      INTEGER :: m, mp, spin, spinp, row, column

      DO m = -2, 2
         DO spin = 1, 2
            row = (m+2)*2+spin
            DO mp = -2, 2
               DO spinp = 1, 2
                  column = (mp+2)*2+spinp
                  rho_flat(row,column) = rho_native(m,spin,mp,spinp)
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE flatten_native_density

   SUBROUTINE validate_state_record(record,rho_native)
      TYPE(t_rixs_state_record), INTENT(IN) :: record
      COMPLEX, INTENT(IN) :: rho_native(-2:2,2,-2:2,2)
      COMPLEX :: rho_flat(10,10)
      REAL :: scale, tolerance, hermiticity_error

      CALL flatten_native_density(rho_native,rho_flat)
      scale = MAX(1.0,ABS(record%d_weight))
      tolerance = 10000.0*EPSILON(scale)*scale
      hermiticity_error = MAXVAL(ABS(rho_flat-CONJG(TRANSPOSE(rho_flat))))
      IF (hermiticity_error > tolerance) CALL juDFT_error("Non-Hermitian RIXS state-character density matrix", &
                                                          calledby="m_rixs_state_character")
      IF (ABS(record%t2g_weight+record%eg_weight-record%d_weight) > tolerance) THEN
         CALL juDFT_error("RIXS state-character t2g/eg reconstruction failed",calledby="m_rixs_state_character")
      END IF
      IF (ABS(SUM(record%jeff_weights)-record%t2g_weight) > tolerance) THEN
         CALL juDFT_error("RIXS state-character j_eff reconstruction failed",calledby="m_rixs_state_character")
      END IF
   END SUBROUTINE validate_state_record

END MODULE m_rixs_state_character
