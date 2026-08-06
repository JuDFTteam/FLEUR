!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_matrix_element_factory
    !> Driver for the evaluation of matrix elements between first-variation
    !> eigenstates. It reads the eigenvectors, computes the basis matching (abc)
    !> coefficients and the radial functions, keeps them in module variables so
    !> that several different matrix elements can be evaluated without redoing
    !> this setup, and finally calls the calc_matrix_elements routine of the
    !> matrix-element object.
    USE m_types_mat
    USE m_types_abc
    USE m_types_radfun
    USE m_types_spinor_layout, ONLY: radial_slot, melem_stack_spinor
    USE m_types_usdus
    IMPLICIT NONE
    PRIVATE

    !Cached data that depends only on the potential (not on the k-point).
    !Invalidated by matrix_element_factory_reset only, i.e. once per SCF iteration.
    TYPE(t_radfun), ALLOCATABLE :: radfun_store(:)  !(ntype)
    TYPE(t_usdus)               :: usdus_store      !all types and spins
    LOGICAL                     :: radfun_valid = .FALSE.

    !Cached data for one k-point, keyed on (eig_id, ikpt, band selection)
    TYPE(t_mat), ALLOCATABLE :: zmat_store(:)   !(nrec) eigenvectors as read
    !Allocated only when the two records are the halves of one spinor: the same state
    !as a single matrix of 2N rows. A consumer that addresses a spin block by row offset
    !needs that shape; one that contracts record by record uses zmat_store.
    TYPE(t_mat), ALLOCATABLE :: zmat_spinor(:)  !(1)
    TYPE(t_abc), ALLOCATABLE :: abc_store(:,:)  !(2,ntype) matching coefficients
    INTEGER              :: cached_eig_id = -1, cached_ikpt = -1
    INTEGER              :: cached_nrec = -1   !record layout the cache was filled with
    INTEGER, ALLOCATABLE :: cached_list(:)      !band selection; unallocated=all bands
    LOGICAL              :: cache_valid = .FALSE.

    PUBLIC :: matrix_element_factory, matrix_element_factory_reset

CONTAINS

    SUBROUTINE matrix_element_factory_reset()
        !> Invalidate all cached data. Must be called whenever the eig-file
        !> content or the potential changes, i.e. once per SCF iteration.
        IF (ALLOCATED(radfun_store)) DEALLOCATE(radfun_store)
        CALL usdus_store%free()
        radfun_valid = .FALSE.
        CALL reset_k_cache()
    END SUBROUTINE matrix_element_factory_reset

    SUBROUTINE reset_k_cache()
        IF (ALLOCATED(zmat_store))  DEALLOCATE(zmat_store)
        IF (ALLOCATED(zmat_spinor)) DEALLOCATE(zmat_spinor)
        IF (ALLOCATED(abc_store))   DEALLOCATE(abc_store)
        IF (ALLOCATED(cached_list)) DEALLOCATE(cached_list)
        cached_eig_id = -1
        cached_ikpt   = -1
        cached_nrec   = -1
        cache_valid   = .FALSE.
    END SUBROUTINE reset_k_cache

    LOGICAL FUNCTION cache_hit(eig_id, ikpt, nrec, ev_list)
        INTEGER, INTENT(IN)           :: eig_id, ikpt, nrec
        INTEGER, OPTIONAL, INTENT(IN) :: ev_list(:)

        cache_hit = cache_valid .AND. eig_id == cached_eig_id .AND. ikpt == cached_ikpt &
                    .AND. nrec == cached_nrec
        IF (.NOT.cache_hit) RETURN
        IF (PRESENT(ev_list) .NEQV. ALLOCATED(cached_list)) THEN
            cache_hit = .FALSE.
        ELSE IF (PRESENT(ev_list)) THEN
            cache_hit = SIZE(ev_list) == SIZE(cached_list)
            IF (cache_hit) cache_hit = ALL(ev_list == cached_list)
        END IF
    END FUNCTION cache_hit

    SUBROUTINE matrix_element_factory(matel, eig_id, ikpt, input, atoms, sym, cell, &
                                      noco, nococonv, enpara, lapw, vtot, fmpi, ev_list, &
                                      l_both_spinors, kpts)
        USE m_types_matelements
        USE m_types_input
        USE m_types_atoms
        USE m_types_sym
        USE m_types_cell
        USE m_types_noco
        USE m_types_nococonv
        USE m_types_enpara
        USE m_types_lapw
        USE m_types_kpts
        USE m_types_potden
        USE m_types_mpi
        USE m_eig66_io, ONLY: read_eig
        USE m_judft, ONLY: judft_error

        CLASS(t_matelements), INTENT(INOUT) :: matel
        INTEGER,           INTENT(IN) :: eig_id, ikpt
        TYPE(t_input),     INTENT(IN) :: input
        TYPE(t_atoms),     INTENT(IN) :: atoms
        TYPE(t_sym),       INTENT(IN) :: sym
        TYPE(t_cell),      INTENT(IN) :: cell
        TYPE(t_noco),      INTENT(IN) :: noco
        TYPE(t_nococonv),  INTENT(IN) :: nococonv
        TYPE(t_enpara),    INTENT(IN) :: enpara
        TYPE(t_lapw),      INTENT(IN) :: lapw
        TYPE(t_potden),    INTENT(IN) :: vtot
        TYPE(t_mpi),       INTENT(IN) :: fmpi
        INTEGER, OPTIONAL, INTENT(IN) :: ev_list(:) !select these states, default: all
        !> Read records 1 and 2 as the two components of one spinor, which is how the
        !> eigenvectors are stored once spin-orbit coupling has been applied: two records
        !> while input%jspins is 1. Absent, input%jspins records are read, each an
        !> independent spin channel. Which of the two an eig file holds depends on the
        !> stage of the calculation rather than on any flag, so it cannot be decided here.
        LOGICAL, OPTIONAL, INTENT(IN) :: l_both_spinors
        !> The k-point mesh. Given, a k-point beyond the irreducible zone is served by rotating
        !> the states of its parent, which is where they are stored. Absent, only the irreducible
        !> points can be asked for.
        TYPE(t_kpts), OPTIONAL, INTENT(IN) :: kpts

        INTEGER :: num_bands, n, jsp, jsp_rad, nrec, neig_actual, ikpt_stored
        LOGICAL :: l_real_zmat
        INTEGER, ALLOCATABLE :: read_list(:)
        LOGICAL :: l_spinor_records
        INTEGER :: n_comp

        l_spinor_records = .FALSE.
        IF (PRESENT(l_both_spinors)) l_spinor_records = l_both_spinors
        IF (l_spinor_records .AND. (.NOT.noco%l_soc .OR. noco%l_noco)) CALL judft_error( &
            'matrix_element_factory: spinor records only exist for l_soc=T, l_noco=F', &
            calledby='matrix_element_factory')

        !One record per potential, with two exceptions. Non-collinearly the Hamiltonian was
        !2N x 2N and the whole spinor sits in a single record whatever jspins says. After a
        !second variation there are two records even for one potential.
        nrec = input%jspins
        IF (noco%l_noco)      nrec = 1
        IF (l_spinor_records) nrec = 2

        !How many spin components the state HAS, which is not how many records it is stored
        !in: non-collinearly both components live inside the single record, and a state with
        !no spin structure has one component that serves as both.
        n_comp = 1
        IF (noco%l_noco .OR. l_spinor_records .OR. input%jspins == 2) n_comp = 2

        num_bands = input%neig
        IF (PRESENT(ev_list)) THEN
            num_bands = SIZE(ev_list)
        ELSE
            ! Only as many states as are actually stored for this k-point are
            ! available; requesting more would read uninitialized eig storage
            ! (harmless zeros in serial mem/DA, but stale window memory under
            ! MPI-RMA). Clamp to the stored count, as the old alineso did.
            ikpt_stored = ikpt
            IF (PRESENT(kpts)) THEN
                IF (ikpt > kpts%nkpt) ikpt_stored = kpts%bkp(ikpt)
            END IF
            CALL read_eig(eig_id, ikpt_stored, 1, neig=neig_actual)
            num_bands = MIN(num_bands, neig_actual)
        END IF

        !Radial functions and energy derivatives: k-independent, so generated
        !only once per SCF iteration (i.e. after each factory reset)
        IF (.NOT.radfun_valid) THEN
            ALLOCATE(radfun_store(atoms%ntype))
            DO n = 1, atoms%ntype
                CALL radfun_store(n)%generate_radial_functions(atoms, input, enpara, fmpi, &
                                                               vtot, n, usdus_out=usdus_store)
            END DO
            radfun_valid = .TRUE.
        END IF

        !Eigenvectors and abc coefficients for this k-point and band selection
        IF (.NOT.cache_hit(eig_id, ikpt, nrec, ev_list)) THEN
            CALL reset_k_cache()

            !Read and select the state vectors
            !TODO: at the secvar_soc call site the eigenvectors are also read in
            !t_secvar%initialize (needed there for the back-transform); a getter
            !for zmat_store could avoid this duplication.
            ALLOCATE(zmat_store(nrec))
            ! In SOC/noncollinear mode the eigenvectors are stored as complex
            ! objects in eig66 and must be read into complex matrices.
            l_real_zmat = input%l_real .AND. (.NOT.(noco%l_soc .OR. noco%l_noco))
            IF (PRESENT(ev_list)) THEN
                read_list = ev_list
            ELSE
                read_list = [(n, n = 1, num_bands)]
            END IF
            DO jsp = 1, nrec
                CALL zmat_store(jsp)%init(l_real_zmat, lapw%nmat, num_bands)
                IF (PRESENT(kpts)) THEN
                    CALL read_eig(eig_id, ikpt, jsp, list=read_list, zmat=zmat_store(jsp), &
                                  kpts=kpts, input=input, noco=noco, nococonv=nococonv, &
                                  sym=sym, atoms=atoms, cell=cell)
                ELSE
                    CALL read_eig(eig_id, ikpt, jsp, list=read_list, zmat=zmat_store(jsp))
                END IF
            END DO

            !Two records that are the halves of one spinor are also kept stacked, so that
            !what is handed over is the whole state in one matrix however the eig file
            !stored it. Without this a consumer that reads the spin-down block at
            !nv(1) + nlotot addresses one row past the end of an N-row record, and the
            !values it finds there are plausible rather than obviously wrong.
            IF (l_spinor_records) THEN
                ALLOCATE(zmat_spinor(1))
                CALL melem_stack_spinor(zmat_store(1), zmat_store(2), zmat_spinor(1))
            END IF

            !Calculate the abc coefficients for all atom types and both spins
            ALLOCATE(abc_store(2, atoms%ntype))
            DO n = 1, atoms%ntype
                DO jsp = 1, n_comp
                    !A single radial set serves both components of a spinor.
                    jsp_rad = radial_slot(radfun_store, jsp)
                    CALL abc_store(jsp,n)%init(input, atoms, num_bands, n)
                    !One record per component when they were stored separately; the single
                    !record for both when it holds the whole spinor, where the spin index
                    !is what selects the block.
                    CALL abc_store(jsp,n)%calc_abc(input, atoms, sym, cell, lapw, num_bands, &
                                                   usdus_store, noco, nococonv, jsp_rad, n, &
                                                   zmat_store(MIN(jsp, nrec)))
                END DO
                !A state with no spin structure is its own second channel. Two real
                !components must never be equated.
                IF (n_comp == 1) abc_store(2,n) = abc_store(1,n)
            END DO

            cached_eig_id = eig_id
            cached_ikpt   = ikpt
            cached_nrec   = nrec
            IF (PRESENT(ev_list)) cached_list = ev_list
            cache_valid = .TRUE.
        END IF

        !An operator whose result carries Cartesian components keeps them in a store
        !that has no distributed counterpart, so the two cannot be combined. Say so here
        !rather than let init_mat abort with the sub-communicator already chosen.
        IF (fmpi%n_size > 1 .AND. (matel%n_alpha > 1 .OR. matel%n_beta > 1)) CALL judft_error( &
            'matrix_element_factory: an operator with Cartesian components cannot be distributed', &
            hint='run k-parallel (n_size = 1); eigenvector parallelism is not supported for it', &
            calledby='matrix_element_factory')

        !Provide the result matrix; its spinor structure (2x2 or single block)
        !follows from the spinor flags of the matrix-element object
#ifdef CPP_MPI
        IF (fmpi%n_size > 1) THEN
            !Eigenvalue parallelization: distribute the matrix blocks row-cyclically
            !(full rows, columns distributed cyclically) over the sub_comm group.
            CALL matel%init_mat(num_bands, mpi_subcomm=fmpi%sub_comm)
        ELSE
            CALL matel%init_mat(num_bands)
        END IF
#else
        CALL matel%init_mat(num_bands)
#endif

        !Compute the matrix elements. zmat is the state at this k-point in as few matrices
        !as it takes: one when it is a whole spinor, two when the records are independent
        !spin channels. The abc coefficients stay per record either way, since each is
        !contracted with its own component.
        IF (ALLOCATED(zmat_spinor)) THEN
            CALL matel%calc_matrix_elements(zmat_spinor, abc_store, radfun_store, usdus_store)
        ELSE
            CALL matel%calc_matrix_elements(zmat_store, abc_store, radfun_store, usdus_store)
        END IF

    END SUBROUTINE matrix_element_factory

END MODULE m_matrix_element_factory
