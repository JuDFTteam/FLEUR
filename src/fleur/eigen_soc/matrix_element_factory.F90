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
    USE m_types_usdus
    IMPLICIT NONE
    PRIVATE

    !Cached data that depends only on the potential (not on the k-point).
    !Invalidated by matrix_element_factory_reset only, i.e. once per SCF iteration.
    TYPE(t_radfun), ALLOCATABLE :: radfun_store(:)  !(ntype)
    TYPE(t_usdus)               :: usdus_store      !all types and spins
    LOGICAL                     :: radfun_valid = .FALSE.

    !Cached data for one k-point, keyed on (eig_id, ikpt, band selection)
    TYPE(t_mat), ALLOCATABLE :: zmat_store(:)   !(jspins) eigenvectors as read
    TYPE(t_abc), ALLOCATABLE :: abc_store(:,:)  !(2,ntype) matching coefficients
    INTEGER              :: cached_eig_id = -1, cached_ikpt = -1
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
        IF (ALLOCATED(zmat_store)) DEALLOCATE(zmat_store)
        IF (ALLOCATED(abc_store))  DEALLOCATE(abc_store)
        IF (ALLOCATED(cached_list)) DEALLOCATE(cached_list)
        cached_eig_id = -1
        cached_ikpt   = -1
        cache_valid   = .FALSE.
    END SUBROUTINE reset_k_cache

    LOGICAL FUNCTION cache_hit(eig_id, ikpt, ev_list)
        INTEGER, INTENT(IN)           :: eig_id, ikpt
        INTEGER, OPTIONAL, INTENT(IN) :: ev_list(:)

        cache_hit = cache_valid .AND. eig_id == cached_eig_id .AND. ikpt == cached_ikpt
        IF (.NOT.cache_hit) RETURN
        IF (PRESENT(ev_list) .NEQV. ALLOCATED(cached_list)) THEN
            cache_hit = .FALSE.
        ELSE IF (PRESENT(ev_list)) THEN
            cache_hit = SIZE(ev_list) == SIZE(cached_list)
            IF (cache_hit) cache_hit = ALL(ev_list == cached_list)
        END IF
    END FUNCTION cache_hit

    SUBROUTINE matrix_element_factory(matel, eig_id, ikpt, input, atoms, sym, cell, &
                                      noco, nococonv, enpara, lapw, vtot, fmpi, ev_list)
        USE m_types_matelements
        USE m_types_input
        USE m_types_atoms
        USE m_types_sym
        USE m_types_cell
        USE m_types_noco
        USE m_types_nococonv
        USE m_types_enpara
        USE m_types_lapw
        USE m_types_potden
        USE m_types_mpi
        USE m_eig66_io, ONLY: read_eig

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

        INTEGER :: num_bands, n, jsp

        num_bands = input%neig
        IF (PRESENT(ev_list)) num_bands = SIZE(ev_list)

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
        IF (.NOT.cache_hit(eig_id, ikpt, ev_list)) THEN
            CALL reset_k_cache()

            !Read and select the state vectors
            !TODO: at the secvar_soc call site the eigenvectors are also read in
            !t_secvar%initialize (needed there for the back-transform); a getter
            !for zmat_store could avoid this duplication.
            ALLOCATE(zmat_store(input%jspins))
            DO jsp = 1, input%jspins
                CALL zmat_store(jsp)%init(input%l_real, lapw%nmat, num_bands)
                IF (PRESENT(ev_list)) THEN
                    CALL read_eig(eig_id, ikpt, jsp, list=ev_list, zmat=zmat_store(jsp))
                ELSE
                    CALL read_eig(eig_id, ikpt, jsp, zmat=zmat_store(jsp))
                END IF
            END DO

            !Calculate the abc coefficients for all atom types and both spins
            ALLOCATE(abc_store(2, atoms%ntype))
            DO n = 1, atoms%ntype
                DO jsp = 1, input%jspins
                    CALL abc_store(jsp,n)%init(input, atoms, num_bands, n)
                    CALL abc_store(jsp,n)%calc_abc(input, atoms, sym, cell, lapw, num_bands, &
                                                   usdus_store, noco, nococonv, jsp, n, zmat_store(jsp))
                END DO
                IF (input%jspins == 1) abc_store(2,n) = abc_store(1,n)
            END DO

            cached_eig_id = eig_id
            cached_ikpt   = ikpt
            IF (PRESENT(ev_list)) cached_list = ev_list
            cache_valid = .TRUE.
        END IF

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

        !Compute the matrix elements
        CALL matel%calc_matrix_elements(zmat_store, abc_store, radfun_store, usdus_store)

    END SUBROUTINE matrix_element_factory

END MODULE m_matrix_element_factory
