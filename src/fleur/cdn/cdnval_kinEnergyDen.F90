!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_cdnval_kinEnergyDen
   !! Module to compute the kinetic energy density (KED) directly from the
   !! wavefunctions using the LAPW muffin-tin expansion.
   !!
   !! This implements the "direct" approach:
   !!   τ(r) = (1/2) Σ_i f_i |∇ψ_i(r)|²
   !!
   !! The calculation reuses the same density matrix as the charge density
   !! (rhonmt), but converts to the KED using modified radial products that
   !! include the radial derivative of the basis functions and the angular
   !! gradient coupling factor.
   !!
   !! Reference: Doumont et al., Phys. Rev. B 106, 235159 (2022)
   !!            "Implementation of self-consistent MGGA functionals in
   !!             augmented plane wave based methods"

   USE m_juDFT
#ifdef CPP_MPI
   use mpi
#endif
   implicit none

CONTAINS

   SUBROUTINE cdnval_kinEnergyDen(eig_id, fmpi, kpts, jspin, noco, nococonv, input, &
                                   cell, atoms, enpara, stars, vacuum, sphhar, sym, &
                                   vTot, results, kinEnergyDen)
      !! Compute the muffin-tin kinetic energy density for valence electrons.
      !!
      !! This routine follows the same structure as cdnval: it loops over
      !! k-points and atom types, computes the matching coefficients (abc),
      !! builds the density matrix, and then converts to the kinetic energy
      !! density using to_kinetic_energy_density (instead of to_full_density
      !! for the charge density).
      !!
      !! The interstitial KED is not computed here — that requires a separate
      !! plane-wave gradient approach.

      USE m_types
      USE m_constants
      USE m_eig66_io
      USE m_genMTBasis
      USE m_types_dos
      USE m_types_mcd
      USE m_types_slab
      USE m_types_jDOS
      USE m_types_vacDOS
      USE m_types_orbcomp
      USE m_types_denmatrix
      USE m_types_radfun
      USE m_types_cdnval
      use m_types_abc
      USE m_pwden_kinEnergyDen
#ifdef CPP_MPI
      USE m_mpi_col_den
#endif

      IMPLICIT NONE

      TYPE(t_results), INTENT(IN)    :: results
      TYPE(t_mpi), INTENT(IN)        :: fmpi
      TYPE(t_enpara), INTENT(IN)     :: enpara
      TYPE(t_input), INTENT(IN)      :: input
      TYPE(t_vacuum), INTENT(IN)     :: vacuum
      TYPE(t_noco), INTENT(IN)       :: noco
      TYPE(t_nococonv), INTENT(IN)   :: nococonv
      TYPE(t_sym), INTENT(IN)        :: sym
      TYPE(t_stars), INTENT(IN)      :: stars
      TYPE(t_cell), INTENT(IN)       :: cell
      TYPE(t_kpts), INTENT(IN)       :: kpts
      TYPE(t_sphhar), INTENT(IN)     :: sphhar
      TYPE(t_atoms), INTENT(IN)      :: atoms
      TYPE(t_potden), INTENT(IN)     :: vTot
      TYPE(t_potden), INTENT(INOUT)  :: kinEnergyDen

      ! Scalar Arguments
      INTEGER, INTENT(IN) :: eig_id, jspin

      ! Local Scalars
      INTEGER :: ikpt, ikpt_i, jsp, ispin, ispinpr, ispin123
      INTEGER :: iErr, nbands, noccbd, iType
      INTEGER :: skip_t, skip_tt, nbasfcn, abc_itype
      INTEGER :: max_length_k_list, nk
      LOGICAL :: l_real

      ! Local Arrays
      REAL, ALLOCATABLE :: we(:), eig(:)
      REAL              :: bkpt(3)
      INTEGER, ALLOCATABLE :: ev_list(:)

      TYPE(t_lapw)           :: lapw
      TYPE(t_usdus)          :: usdus
      TYPE(t_mat)            :: zMat
      TYPE(t_radfun)         :: radfun(atoms%ntype)
      TYPE(t_abc), allocatable :: abc(:, :)
      TYPE(t_denmatrix), allocatable :: denmatrix(:, :, :)
      TYPE(t_cdnvalJob)      :: cdnvalJob

      CALL timestart("cdnval_kinEnergyDen")

      l_real = sym%invs .AND. (.NOT. noco%l_soc) .AND. (.NOT. noco%l_noco) .AND. atoms%n_hia == 0

      ! For simplicity, this routine only handles the diagonal spin case
      ! (no off-diagonal / l_mperp). Extension is straightforward.
      allocate(denmatrix(jspin:jspin, jspin:jspin, atoms%ntype))
      DO iType = 1, atoms%ntype
         call denmatrix(jspin, jspin, iType)%init(iType, atoms, input, sphhar)
      END DO
      allocate(abc(jspin:jspin, merge(1, atoms%ntype, atoms%n_v == 0)))

      ! Initialize
      CALL usdus%init(atoms, input%jspins)
      CALL cdnvalJob%init(fmpi, input, kpts, noco, results, jspin)

      skip_tt = dot_product(enpara%skiplo(:atoms%ntype, jspin), atoms%neq(:atoms%ntype))
      IF (noco%l_soc .OR. noco%l_noco) skip_tt = 2*skip_tt

      jsp = MERGE(1, jspin, noco%l_noco)

      max_length_k_list = size(cdnvalJob%k_list)
#ifdef CPP_MPI
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, max_length_k_list, 1, MPI_INTEGER, MPI_MAX, fmpi%mpi_comm, ierr)
#endif

      ! ====================================================================
      ! K-point loop: compute density matrix (same as cdnval)
      ! ====================================================================
      DO ikpt_i = 1, size(cdnvalJob%k_list)
         ikpt = cdnvalJob%k_list(ikpt_i)
         bkpt = kpts%bk(:, ikpt)

         CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, ikpt, cell, fmpi)
         skip_t = skip_tt
         ev_list = cdnvaljob%compact_ev_list(ikpt_i, .false.)
         noccbd = SIZE(ev_list)
         we = cdnvalJob%weights(ev_list, ikpt)
         eig = results%eig(ev_list, ikpt, jsp)

         IF (cdnvalJob%l_evp) THEN
            IF (minval(ev_list) > skip_tt) skip_t = 0
            IF (maxval(ev_list) <= skip_tt) skip_t = noccbd
            IF ((minval(ev_list) <= skip_tt) .AND. (maxval(ev_list) > skip_tt)) skip_t = mod(skip_tt, noccbd)
         END IF

         nbasfcn = MERGE(lapw%nv(1) + lapw%nv(2) + 2*atoms%nlotot, lapw%nv(1) + atoms%nlotot, noco%l_noco)
         CALL zMat%init(l_real, nbasfcn, noccbd)
         CALL read_eig(eig_id, ikpt, jsp, list=ev_list, neig=nbands, zmat=zMat)
#ifdef CPP_MPI
         CALL MPI_BARRIER(fmpi%mpi_comm, iErr)
#endif
         IF (noccbd .LE. 0) CYCLE

         ! Loop over atom types: compute abc coefficients and density matrix
         DO itype = 1, atoms%ntype
            abc_itype = min(itype, size(abc, 2))
            call radfun(itype)%generate_radial_functions(atoms, input, enpara, fmpi, vtot, iType)

            ispin = jspin
            call abc(ispin, abc_itype)%init(input, atoms, radfun(itype)%n_r, noccbd, itype)
            call abc(ispin, abc_itype)%calc_abc(input, atoms, sym, cell, lapw, noccbd, usdus, noco, nococonv, ispin, itype, zMat)

            ! Build the density matrix (same coefficients as for charge density)
            call denmatrix(ispin, ispin, itype)%rhonmt(atoms, sphhar, we, noccbd, itype, &
                                                        sym, .true., abc(ispin, abc_itype), abc(ispin, abc_itype))
         END DO ! itype

         ! ---------------------------------------------------------------
         ! Interstitial KED: compute via plane-wave gradient FFT
         ! ---------------------------------------------------------------
         CALL pwden_kinEnergyDen(stars, kpts, input, cell, atoms, sym, &
                                  ikpt, jspin, lapw, noccbd, we, kinEnergyDen, zMat)
      END DO ! k-point loop

#ifdef CPP_MPI
      ! Synchronize remaining k-points
      DO nk = size(cdnvalJob%k_list) + 1, max_length_k_list
         CALL MPI_BARRIER(fmpi%MPI_COMM, ierr)
      END DO
      ! Collect density matrix across MPI ranks
      DO itype = 1, atoms%ntype
         call denmatrix(jspin, jspin, itype)%mpi_collect(fmpi)
      END DO
      ! Collect interstitial KED star coefficients across MPI ranks
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, kinEnergyDen%pw(:, jspin), &
                         size(kinEnergyDen%pw, 1), MPI_DOUBLE_COMPLEX, MPI_SUM, fmpi%mpi_comm, iErr)
#endif

      ! ====================================================================
      ! Convert density matrix to kinetic energy density (rank 0 only)
      ! ====================================================================
      IF (fmpi%irank == 0) THEN
         CALL timestart("denmatrix_to_kinEnergyDen")
         DO itype = 1, atoms%ntype
            call radfun(itype)%generate_radial_functions(atoms, input, enpara, fmpi, vtot, iType)
            call denmatrix(jspin, jspin, itype)%to_kinetic_energy_density( &
                 jspin, jspin, itype, input, sphhar, atoms, noco, sym, radfun(itype), &
                 kinEnergyDen%mt)
         END DO
         CALL timestop("denmatrix_to_kinEnergyDen")
      END IF

      CALL timestop("cdnval_kinEnergyDen")

   END SUBROUTINE cdnval_kinEnergyDen

END MODULE m_cdnval_kinEnergyDen
