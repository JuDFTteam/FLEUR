!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_secvar_soc
#ifdef CPP_MPI
    use mpi
#endif
    IMPLICIT NONE
CONTAINS

    SUBROUTINE secvar_soc_kpts(atoms, noco, nococonv, input, sym, cell, enpara, lapw, vtot, rsoc, fmpi, ikpt, eig_id, ne_out, eigval_out, &
                               l_both_spinors)
        USE m_types_secvar
        USE m_types_matelements_soc
        USE m_matrix_element_factory
        USE m_types

        TYPE(t_atoms),   INTENT(IN) :: atoms
        TYPE(t_noco),    INTENT(IN) :: noco
        TYPE(t_nococonv),INTENT(IN) :: nococonv
        TYPE(t_input),   INTENT(IN) :: input
        TYPE(t_sym),     INTENT(IN) :: sym
        TYPE(t_cell),    INTENT(IN) :: cell
        TYPE(t_enpara),  INTENT(IN) :: enpara
        TYPE(t_lapw),    INTENT(IN) :: lapw
        TYPE(t_potden),  INTENT(IN) :: vtot
        TYPE(t_rsoc),    INTENT(IN) :: rsoc
        TYPE(t_mpi),     INTENT(IN) :: fmpi
        INTEGER,         INTENT(IN) :: ikpt, eig_id
        INTEGER,         OPTIONAL, INTENT(OUT) :: ne_out
        REAL,            OPTIONAL, INTENT(OUT) :: eigval_out(:)
        LOGICAL,         OPTIONAL, INTENT(IN)  :: l_both_spinors

        TYPE(t_secvar)          :: secvar
        TYPE(t_matelements_soc) :: matelements
        INTEGER :: ncopy

        ! Initialize the second variation matrix
        CALL secvar%initialize(.TRUE., ikpt, eig_id, input, fmpi, lapw, atoms, l_both_spinors)

        write(oUnit,*) "Non-SOC eigenvalues for k-point ", ikpt
        write(oUnit,*) secvar%eig(1:secvar%ne_first, 1)
        if (input%jspins == 2) then
           write(oUnit,*) secvar%eig(1:secvar%ne_first, 2)
        end if

        ! Initialize the SOC matrix element evaluator and compute matrix elements.
        ! The factory reads the eigenvectors, provides the abc coefficients and
        ! radial functions (cached for reuse with other matrix elements) and
        ! allocates the result matrix.
        CALL matelements%init(atoms, noco, input, sym, cell, enpara, lapw, vtot, rsoc,  fmpi, nococonv)
        CALL matrix_element_factory(matelements, eig_id, ikpt, input, atoms, sym, cell, &
                                    noco, nococonv, enpara, lapw, vtot, fmpi)
        ! Hand the computed operator matrix over to the second-variation problem
        CALL MOVE_ALLOC(matelements%mat, secvar%mat)

        ! Add the first-variation eigenvalues on the diagonal
        CALL secvar%add_diagonal_elements()

        ! Diagonalize the second-variation matrix
        CALL secvar%diagonalize()

        ! Return eigenvalues to caller if requested
        IF (PRESENT(ne_out))     ne_out = secvar%ne_second
        IF (PRESENT(eigval_out)) THEN
            ncopy = MIN(secvar%ne_second, SIZE(eigval_out), SIZE(secvar%eigval))
            eigval_out(:ncopy) = secvar%eigval(:ncopy)
        END IF
        write(oUnit,*) "SOC eigenvalues for k-point ", ikpt
        write(oUnit,*) secvar%eigval(1:secvar%ne_second)

        ! Back-transform eigenvectors and write them to file
        CALL secvar%store_eigvec()

    END SUBROUTINE secvar_soc_kpts


    SUBROUTINE secvar_soc(eig_id, fmpi, nococonv, vTot, enpara, fi, results)
        use m_types
        use m_matrix_element_factory

        INTEGER,            INTENT(IN)    :: eig_id
        TYPE(t_mpi),        INTENT(INOUT) :: fmpi
        TYPE(t_nococonv),   INTENT(IN)    :: nococonv
        TYPE(t_potden),     INTENT(IN)    :: vTot
        TYPE(t_enpara),     INTENT(INOUT) :: enpara
        TYPE(t_fleurinput), INTENT(IN)    :: fi
        TYPE(t_results),    INTENT(INOUT) :: results

        TYPE(t_lapw)  :: lapw
        TYPE(t_rsoc)  :: rsoc
        TYPE(t_sym)   :: sym_l
        INTEGER :: nk, nk_i, neig_soc
        LOGICAL :: l_both_spinors
#ifdef CPP_MPI
        INTEGER :: ierr
#endif


        results%eig   = 1.0e300
        results%neig  = 0

        ! Wannier/wannierlib need BOTH spinor components of the SOC eigenvectors. For jspins=1
        ! only one first-variation spin channel exists, so store_eigvec would write a single
        ! record; fleur.F90 already opens the eig file with 2 spin records under exactly this
        ! condition, so ask for the second (spin-down) record to be filled as well.
        ! (This is the l_wann_store behaviour of the removed alineso/eigenso.)
        l_both_spinors = fi%input%l_wann .OR. fi%wannierlib%l_wannierize

        ! The SOC angular matrix elements (soangl) are set up in the global frame.
        ! Hence the basis matching coefficients must not be rotated into the local
        ! frame of the representative atom, so we disable the atom rotation here
        ! (consistent with eigenso, where sym%ngopr is forced to 1). Otherwise
        ! hsmt_ab/the LO setup rotate the abc coefficients into each atom's local
        ! frame, which is inconsistent with the global-frame soangl and yields
        ! wrong SOC matrix elements for symmetry-related (e.g. -k) k-points.
        sym_l = fi%sym
        sym_l%ngopr = 1

        ! Compute radial spin-orbit matrix elements
        call rsoc%init(fi%atoms)
        call rsoc%rad_matrix(fi%atoms, fi%noco, nococonv, fi%input, fmpi, enpara, vTot)
        call rsoc%angles(fi%atoms,fmpi,nococonv%theta,nococonv%phi)
        ! The potential changed since the last call, so all cached data
        ! (eigenvectors, abc coefficients, radial functions) is invalid
        CALL matrix_element_factory_reset()
        ! Loop over k-points assigned to this MPI task
        DO nk_i = 1, SIZE(fmpi%k_list)
            nk = fmpi%k_list(nk_i)
            CALL lapw%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, sym_l, nk, fi%cell, fmpi)
            IF (fmpi%n_rank == 0) THEN
               CALL secvar_soc_kpts(fi%atoms, fi%noco, nococonv, fi%input, sym_l, fi%cell, &
                               enpara, lapw, vTot, rsoc, fmpi, nk, eig_id, &
                               ne_out=results%neig(nk, 1), eigval_out=results%eig(:, nk, 1), &
                               l_both_spinors=l_both_spinors)
            ELSE
               CALL secvar_soc_kpts(fi%atoms, fi%noco, nococonv, fi%input, sym_l, fi%cell, &
                               enpara, lapw, vTot, rsoc, fmpi, nk, eig_id, &
                               l_both_spinors=l_both_spinors)
            END IF
        END DO
        ! Free the cached data
        CALL matrix_element_factory_reset()

#ifdef CPP_MPI
    neig_soc = SIZE(results%eig, 1)
        CALL MPI_ALLREDUCE(MPI_IN_PLACE, results%neig(:, 1), &
                           fi%kpts%nkpt, MPI_INTEGER, MPI_MAX, fmpi%mpi_comm, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, results%eig(:neig_soc,:,1), &
               neig_soc*fi%kpts%nkpt, MPI_DOUBLE_PRECISION, MPI_MIN, fmpi%mpi_comm, ierr)
#endif
        if (fi%input%jspins == 2) then
       neig_soc = SIZE(results%eig, 1)
       results%eig(:neig_soc,:,2) = results%eig(:neig_soc,:,1)
           results%neig(:, 2) = results%neig(:, 1)
        end if
    END SUBROUTINE secvar_soc

END MODULE m_secvar_soc
