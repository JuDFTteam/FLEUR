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
    USE m_types_spinor_layout, ONLY: t_spinor_layout, radial_slot, melem_stack_spinor, &
                                     LAYOUT_SCALAR
    USE m_types_usdus
    USE m_judft, ONLY: judft_error
    IMPLICIT NONE
    PRIVATE

    !Cached data that depends only on the potential (not on the k-point).
    !Invalidated by matrix_element_factory_reset only, i.e. once per SCF iteration.
    TYPE(t_radfun), ALLOCATABLE, TARGET :: radfun_store(:)  !(ntype)
    TYPE(t_usdus),  TARGET              :: usdus_store      !all types and spins
    LOGICAL                     :: radfun_valid = .FALSE.

    !Everything one k-point contributes, keyed on (eig_id, ikpt, nrec, band selection).
    TYPE :: t_k_slot
        TYPE(t_mat), ALLOCATABLE :: zmat(:)         !(nrec) eigenvectors as read
        !Allocated only when the two records are the halves of one spinor: the same state
        !as a single matrix of 2N rows. A consumer that addresses a spin block by row offset
        !needs that shape; one that contracts record by record uses zmat.
        TYPE(t_mat), ALLOCATABLE :: zmat_spinor(:)  !(1)
        TYPE(t_abc), ALLOCATABLE :: abc(:,:)        !(2,ntype) matching coefficients
        INTEGER              :: eig_id = -1, ikpt = -1
        INTEGER              :: nrec = -1           !record layout it was filled with
        INTEGER, ALLOCATABLE :: list(:)             !band selection; unallocated = all bands
        LOGICAL              :: valid = .FALSE.
        INTEGER              :: stamp = 0           !last use, for the eviction order
    END TYPE t_k_slot

    !> Three, because that is what the quantities between two neighbours need: the k itself
    !> and the pair k+b1, k+b2. One would be enough for an operator at a single k, and the
    !> price of the other two is their share of the eigenvectors and coefficients -- a
    !> consumer that walks k in order fills all three and reuses none of them.
    INTEGER, PARAMETER :: N_KSLOT = 3
    TYPE(t_k_slot), TARGET :: kslot(N_KSLOT)
    INTEGER        :: use_clock = 0
    !> A caller that holds on to one k-point while asking for others anchors it here, and
    !> that slot is not overwritten while the anchor stands. Order of last use is not enough
    !> for that case: the anchor is not touched again during the run of neighbours, so it
    !> would be the oldest by the third of them and go.
    INTEGER        :: anchor_slot = 0

    PUBLIC :: matrix_element_factory, matrix_element_factory_reset, matrix_element_states
    PUBLIC :: matrix_element_release_anchor, matrix_element_radial

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
        INTEGER :: is
        DO is = 1, N_KSLOT
            CALL clear_slot(kslot(is))
        END DO
        use_clock   = 0
        anchor_slot = 0
    END SUBROUTINE reset_k_cache

    SUBROUTINE matrix_element_release_anchor()
        !> Let the anchored k-point be overwritten again. Cheap to call more than once.
        anchor_slot = 0
    END SUBROUTINE matrix_element_release_anchor

    SUBROUTINE clear_slot(sl)
        TYPE(t_k_slot), INTENT(INOUT) :: sl
        IF (ALLOCATED(sl%zmat))        DEALLOCATE(sl%zmat)
        IF (ALLOCATED(sl%zmat_spinor)) DEALLOCATE(sl%zmat_spinor)
        IF (ALLOCATED(sl%abc))         DEALLOCATE(sl%abc)
        IF (ALLOCATED(sl%list))        DEALLOCATE(sl%list)
        sl%eig_id = -1; sl%ikpt = -1; sl%nrec = -1
        sl%valid  = .FALSE.
        sl%stamp  = 0
    END SUBROUTINE clear_slot

    LOGICAL FUNCTION slot_matches(sl, eig_id, ikpt, nrec, ev_list)
        TYPE(t_k_slot), INTENT(IN)    :: sl
        INTEGER, INTENT(IN)           :: eig_id, ikpt, nrec
        INTEGER, OPTIONAL, INTENT(IN) :: ev_list(:)

        slot_matches = sl%valid .AND. eig_id == sl%eig_id .AND. ikpt == sl%ikpt &
                       .AND. nrec == sl%nrec
        IF (.NOT.slot_matches) RETURN
        IF (PRESENT(ev_list) .NEQV. ALLOCATED(sl%list)) THEN
            slot_matches = .FALSE.
        ELSE IF (PRESENT(ev_list)) THEN
            slot_matches = SIZE(ev_list) == SIZE(sl%list)
            IF (slot_matches) slot_matches = ALL(ev_list == sl%list)
        END IF
    END FUNCTION slot_matches

    !> The slot holding this k, or the one to overwrite if none does. l_hit says which
    !> of the two happened, so the caller knows whether it still has to fill it. The
    !> victim is the slot used longest ago, which keeps a k that is revisited on every
    !> iteration -- the k of a k/k+b pair -- while its neighbours come and go.
    INTEGER FUNCTION acquire_slot(eig_id, ikpt, nrec, l_hit, ev_list)
        INTEGER, INTENT(IN)           :: eig_id, ikpt, nrec
        LOGICAL, INTENT(OUT)          :: l_hit
        INTEGER, OPTIONAL, INTENT(IN) :: ev_list(:)

        INTEGER :: is, oldest

        DO is = 1, N_KSLOT
            IF (slot_matches(kslot(is), eig_id, ikpt, nrec, ev_list)) THEN
                l_hit = .TRUE.
                acquire_slot = is
                use_clock = use_clock + 1
                kslot(is)%stamp = use_clock
                RETURN
            END IF
        END DO

        l_hit = .FALSE.
        oldest = 0
        DO is = 1, N_KSLOT
            IF (is == anchor_slot) CYCLE          ! the anchor is off limits
            IF (.NOT.kslot(is)%valid) THEN
                oldest = is
                EXIT
            END IF
            IF (oldest == 0) THEN
                oldest = is
            ELSE IF (kslot(is)%stamp < kslot(oldest)%stamp) THEN
                oldest = is
            END IF
        END DO
        IF (oldest == 0) CALL judft_error( &
            'matrix_element_factory: every slot is anchored, so there is none to fill', &
            calledby='acquire_slot')
        CALL clear_slot(kslot(oldest))
        use_clock = use_clock + 1
        kslot(oldest)%stamp = use_clock
        acquire_slot = oldest
    END FUNCTION acquire_slot

    !> Radial functions and their muffin-tin integrals. They depend on the potential and not
    !> on the k-point, so they are generated once and kept until the next factory reset --
    !> i.e. once per SCF iteration.
    SUBROUTINE ensure_radial(atoms, input, enpara, fmpi, vtot)
        USE m_types_input
        USE m_types_atoms
        USE m_types_enpara
        USE m_types_potden
        USE m_types_mpi

        TYPE(t_atoms),  INTENT(IN) :: atoms
        TYPE(t_input),  INTENT(IN) :: input
        TYPE(t_enpara), INTENT(IN) :: enpara
        TYPE(t_mpi),    INTENT(IN) :: fmpi
        TYPE(t_potden), INTENT(IN) :: vtot

        INTEGER :: n

        IF (radfun_valid) RETURN
        ALLOCATE(radfun_store(atoms%ntype))
        DO n = 1, atoms%ntype
            CALL radfun_store(n)%generate_radial_functions(atoms, input, enpara, fmpi, &
                                                           vtot, n, usdus_out=usdus_store)
        END DO
        radfun_valid = .TRUE.
    END SUBROUTINE ensure_radial

    !> The radial functions themselves, for a caller that needs them for something the factory
    !> does not do -- the projections A_mn, the Gaunt overlaps, deciding which radial set a
    !> spin component reads.
    !>
    !> What comes back points into the cache and lives until matrix_element_factory_reset.
    SUBROUTINE matrix_element_radial(atoms, input, enpara, fmpi, vtot, radfun, usdus)
        USE m_types_input
        USE m_types_atoms
        USE m_types_enpara
        USE m_types_potden
        USE m_types_mpi

        TYPE(t_atoms),  INTENT(IN) :: atoms
        TYPE(t_input),  INTENT(IN) :: input
        TYPE(t_enpara), INTENT(IN) :: enpara
        TYPE(t_mpi),    INTENT(IN) :: fmpi
        TYPE(t_potden), INTENT(IN) :: vtot
        TYPE(t_radfun), POINTER, INTENT(OUT) :: radfun(:)   !> (ntype)
        TYPE(t_usdus),  POINTER, INTENT(OUT) :: usdus

        CALL ensure_radial(atoms, input, enpara, fmpi, vtot)
        radfun => radfun_store
        usdus  => usdus_store
    END SUBROUTINE matrix_element_radial

    !> Make sure the states of one k-point, their matching coefficients and the radial
    !> functions are in the cache, and say which slot holds them and how many bands they
    !> carry. Everything the two public entry points share lives here.
    SUBROUTINE ensure_slot(eig_id, ikpt, input, atoms, sym, cell, &
                           noco, nococonv, enpara, lapw, vtot, fmpi, is, num_bands, ev_list, &
                           l_both_spinors, kpts)
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
        INTEGER, INTENT(OUT) :: is         !> the slot the data ended up in
        INTEGER, INTENT(OUT) :: num_bands  !> states selected at this k

        INTEGER :: n, jsp, jsp_rad, nrec, neig_actual, ikpt_stored
        LOGICAL :: l_real_zmat, l_hit
        INTEGER, ALLOCATABLE :: read_list(:)
        LOGICAL :: l_spinor_records
        INTEGER :: n_comp
        TYPE(t_spinor_layout) :: layout

        l_spinor_records = .FALSE.
        IF (PRESENT(l_both_spinors)) l_spinor_records = l_both_spinors

        !> How many records the eigenvectors of this k-point occupy, and how many spin
        !> components the state HAS, which is not the same number: non-collinearly both
        !> components live inside a single record. t_spinor_layout owns that classification
        !> and validates it, so it is asked here instead of derived a second time.
        !>
        !> radfun is deliberately not passed: the radial functions are ensured further down,
        !> so layout%n_radial would fall back to input%jspins rather than to the extent of
        !> the array. Only nrec and layout are taken from here.
        CALL layout%init(input, noco, lapw, atoms, l_both_spinors=l_spinor_records)
        nrec   = layout%nrec
        n_comp = MERGE(2, 1, layout%layout /= LAYOUT_SCALAR)

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

        CALL ensure_radial(atoms, input, enpara, fmpi, vtot)

        !Eigenvectors and abc coefficients for this k-point and band selection. Several
        !k-points are held at once, so asking for a neighbour does not discard the k it is
        !a neighbour of.
        is = acquire_slot(eig_id, ikpt, nrec, l_hit, ev_list)
        IF (.NOT.l_hit) THEN
            !Read and select the state vectors
            !TODO: at the secvar_soc call site the eigenvectors are also read in
            !t_secvar%initialize (needed there for the back-transform); a getter
            !for zmat_store could avoid this duplication.
            ALLOCATE(kslot(is)%zmat(nrec))
            ! In SOC/noncollinear mode the eigenvectors are stored as complex
            ! objects in eig66 and must be read into complex matrices.
            l_real_zmat = input%l_real .AND. (.NOT.(noco%l_soc .OR. noco%l_noco))
            IF (PRESENT(ev_list)) THEN
                read_list = ev_list
            ELSE
                read_list = [(n, n = 1, num_bands)]
            END IF
            DO jsp = 1, nrec
                CALL kslot(is)%zmat(jsp)%init(l_real_zmat, lapw%nmat, num_bands)
                IF (PRESENT(kpts)) THEN
                    CALL read_eig(eig_id, ikpt, jsp, list=read_list, zmat=kslot(is)%zmat(jsp), &
                                  kpts=kpts, input=input, noco=noco, nococonv=nococonv, &
                                  sym=sym, atoms=atoms, cell=cell)
                ELSE
                    CALL read_eig(eig_id, ikpt, jsp, list=read_list, zmat=kslot(is)%zmat(jsp))
                END IF
            END DO

            !Two records that are the halves of one spinor are also kept stacked, so that
            !what is handed over is the whole state in one matrix however the eig file
            !stored it. Without this a consumer that reads the spin-down block at
            !nv(1) + nlotot addresses one row past the end of an N-row record, and the
            !values it finds there are plausible rather than obviously wrong.
            IF (l_spinor_records) THEN
                ALLOCATE(kslot(is)%zmat_spinor(1))
                CALL melem_stack_spinor(kslot(is)%zmat(1), kslot(is)%zmat(2), kslot(is)%zmat_spinor(1))
            END IF

            !Calculate the abc coefficients for all atom types and both spins
            ALLOCATE(kslot(is)%abc(2, atoms%ntype))
            DO n = 1, atoms%ntype
                DO jsp = 1, n_comp
                    !A single radial set serves both components of a spinor.
                    jsp_rad = radial_slot(radfun_store, jsp)
                    CALL kslot(is)%abc(jsp,n)%init(input, atoms, num_bands, n)
                    !One record per component when they were stored separately; the single
                    !record for both when it holds the whole spinor, where the spin index
                    !is what selects the block.
                    CALL kslot(is)%abc(jsp,n)%calc_abc(input, atoms, sym, cell, lapw, num_bands, &
                                                   usdus_store, noco, nococonv, jsp_rad, n, &
                                                   kslot(is)%zmat(MIN(jsp, nrec)))
                END DO
                !A state with no spin structure is its own second channel. Two real
                !components must never be equated.
                IF (n_comp == 1) kslot(is)%abc(2,n) = kslot(is)%abc(1,n)
            END DO

            kslot(is)%eig_id = eig_id
            kslot(is)%ikpt   = ikpt
            kslot(is)%nrec   = nrec
            IF (PRESENT(ev_list)) kslot(is)%list = ev_list
            kslot(is)%valid = .TRUE.
        END IF

    END SUBROUTINE ensure_slot

    !> Evaluate one operator at one k-point.
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
        INTEGER, OPTIONAL, INTENT(IN) :: ev_list(:)
        LOGICAL, OPTIONAL, INTENT(IN) :: l_both_spinors
        TYPE(t_kpts), OPTIONAL, INTENT(IN) :: kpts

        INTEGER :: is, num_bands

        CALL ensure_slot(eig_id, ikpt, input, atoms, sym, cell, noco, nococonv, enpara, &
                         lapw, vtot, fmpi, is, num_bands, ev_list, l_both_spinors, kpts)

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
        IF (ALLOCATED(kslot(is)%zmat_spinor)) THEN
            CALL matel%calc_matrix_elements(kslot(is)%zmat_spinor, kslot(is)%abc, &
                                            radfun_store, usdus_store)
        ELSE
            CALL matel%calc_matrix_elements(kslot(is)%zmat, kslot(is)%abc, &
                                            radfun_store, usdus_store)
        END IF

    END SUBROUTINE matrix_element_factory

    !> The states at one k-point and their coefficients, for a consumer that combines two
    !> k-points itself and so cannot be written as an operator at a single one.
    !>
    !> What comes back points into the cache and stays valid until a later request evicts
    !> that slot. Slots are reused in order of last use, so the data of a k-point asked for
    !> on every iteration outlives its neighbours; N_KSLOT of them are alive at any time.
    SUBROUTINE matrix_element_states(eig_id, ikpt, input, atoms, sym, cell, &
                                     noco, nococonv, enpara, lapw, vtot, fmpi, &
                                     zmat, abc, ev_list, l_both_spinors, kpts, l_anchor)
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
        !> The states one record per element, which is how they were read: a consumer that
        !> reaches a spin block by row offset needs them apart. The stacked 2N form the
        !> operators are given is not handed out here.
        TYPE(t_mat), POINTER, INTENT(OUT) :: zmat(:)    !> (nrec)
        TYPE(t_abc), POINTER, INTENT(OUT) :: abc(:,:)   !> (2,ntype)
        INTEGER, OPTIONAL, INTENT(IN) :: ev_list(:)
        LOGICAL, OPTIONAL, INTENT(IN) :: l_both_spinors
        TYPE(t_kpts), OPTIONAL, INTENT(IN) :: kpts
        !> Keep this k-point until the anchor is released. A caller that holds on to what
        !> comes back while asking for other k-points needs it: without an anchor its slot
        !> is the oldest by the third of them, and the pointers stop being its own.
        LOGICAL, OPTIONAL, INTENT(IN) :: l_anchor

        INTEGER :: is, num_bands

        CALL ensure_slot(eig_id, ikpt, input, atoms, sym, cell, noco, nococonv, enpara, &
                         lapw, vtot, fmpi, is, num_bands, ev_list, l_both_spinors, kpts)

        IF (PRESENT(l_anchor)) THEN
            IF (l_anchor) anchor_slot = is
        END IF

        zmat => kslot(is)%zmat
        abc  => kslot(is)%abc
    END SUBROUTINE matrix_element_states

END MODULE m_matrix_element_factory
