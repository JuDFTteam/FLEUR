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

    SUBROUTINE secvar_soc_kpts(atoms, noco, nococonv, input, sym, cell, enpara, lapw, vtot, rsoc, fmpi, ikpt, eig_id, ne_out, eigval_out)
        USE m_types_secvar
        USE m_types_matelements_soc
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

        TYPE(t_secvar)          :: secvar
        TYPE(t_matelements_soc) :: matelements

        ! Initialize the second variation matrix
        CALL secvar%initialize(.TRUE., ikpt, eig_id, input, fmpi, lapw, atoms)

        write(oUnit,*) "Non-SOC eigenvalues for k-point ", ikpt
        write(oUnit,*) secvar%eig(1:secvar%ne_first, 1)
        if (input%jspins == 2) then
           write(oUnit,*) secvar%eig(1:secvar%ne_first, 2)
        end if

        ! Initialize the SOC matrix element evaluator and compute matrix elements
        CALL matelements%init(atoms, noco, input, sym, cell, enpara, lapw, vtot, rsoc,  fmpi, nococonv)
        CALL matelements%add_matrix_elements(secvar%zmat, secvar%mat)

        ! Add the first-variation eigenvalues on the diagonal
        CALL secvar%add_diagonal_elements()

        ! Diagonalize the second-variation matrix
        CALL secvar%diagonalize()

        ! Return eigenvalues to caller if requested
        IF (PRESENT(ne_out))     ne_out = secvar%ne_second
        IF (PRESENT(eigval_out)) eigval_out(:secvar%ne_second) = secvar%eigval
        write(oUnit,*) "SOC eigenvalues for k-point ", ikpt
        write(oUnit,*) secvar%eigval(1:secvar%ne_second)

        ! Back-transform eigenvectors and write them to file
        CALL secvar%store_eigvec()

    END SUBROUTINE secvar_soc_kpts


    SUBROUTINE secvar_soc(eig_id, fmpi, nococonv, vTot, enpara, fi, results)
        use m_types 
        USE m_spnorb

        INTEGER,            INTENT(IN)    :: eig_id
        TYPE(t_mpi),        INTENT(INOUT) :: fmpi
        TYPE(t_nococonv),   INTENT(IN)    :: nococonv
        TYPE(t_potden),     INTENT(IN)    :: vTot
        TYPE(t_enpara),     INTENT(INOUT) :: enpara
        TYPE(t_fleurinput), INTENT(IN)    :: fi
        TYPE(t_results),    INTENT(INOUT) :: results

        TYPE(t_lapw)  :: lapw
        TYPE(t_rsoc)  :: rsoc
        INTEGER :: nk, nk_i
#ifdef CPP_MPI
        INTEGER :: ierr
#endif


        results%eig   = 1.0e300
        results%neig  = 0

        ! Compute radial spin-orbit matrix elements
        !CALL spnorb(fi%atoms, fi%noco, nococonv, fi%input, fmpi, enpara, vTot%mt, usdus, rsoc, .TRUE.)
        call rsoc%init(fi%atoms)
        call rsoc%rad_matrix(fi%atoms, fi%noco, nococonv, fi%input, fmpi, enpara, vTot)
        call rsoc%angles(fi%atoms,fmpi,nococonv%theta,nococonv%phi)
        ! Loop over k-points assigned to this MPI task
        DO nk_i = 1, SIZE(fmpi%k_list)
            nk = fmpi%k_list(nk_i)
            CALL lapw%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, nk, fi%cell, fmpi)
            CALL secvar_soc_kpts(fi%atoms, fi%noco, nococonv, fi%input, fi%sym, fi%cell, &
                            enpara, lapw, vTot, rsoc, fmpi, nk, eig_id, &
                            ne_out=results%neig(nk, 1), eigval_out=results%eig(:, nk, 1))
        END DO
        
#ifdef CPP_MPI
        CALL MPI_ALLREDUCE(MPI_IN_PLACE, results%neig(:, 1), &
                           fi%kpts%nkpt, MPI_INTEGER, MPI_SUM, fmpi%mpi_comm, ierr)
        CALL MPI_ALLREDUCE(MPI_IN_PLACE, results%eig(:2*fi%input%neig,:,1), &
                           2*fi%input%neig*fi%kpts%nkpt, MPI_DOUBLE_PRECISION, MPI_MIN, fmpi%mpi_comm, ierr)
#endif
        if (fi%input%jspins == 2) then
           results%eig(:2*fi%input%neig,:,2) = results%eig(:2*fi%input%neig,:,1)
           results%neig(:, 2) = results%neig(:, 1)
        end if
    END SUBROUTINE secvar_soc

END MODULE m_secvar_soc
