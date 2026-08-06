!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_matelements
    USE m_types_mat
    USE m_types_mpimat
    USE m_types_abc
    USE m_types_radfun
    USE m_types_usdus
    USE m_judft
    IMPLICIT NONE
    PRIVATE

    TYPE, ABSTRACT :: t_matelements
        LOGICAL :: spinorwavefcts = .FALSE. !> whether the matrix elements are for spinor wave functions
        LOGICAL :: spinoroperator = .FALSE. !> whether the matrix elements are for a spinor operator
        CLASS(t_mat), ALLOCATABLE :: mat(:,:) !> spinor matrix representation: 2x2 blocks if a spinor flag is set, a single block otherwise
        !> Result of an operator whose components carry Cartesian indices, which mat cannot
        !> hold: (band, band, spin, spin, component). The component runs over alpha, or over
        !> the pairs (alpha,beta) in row-major order for an operator with two of them.
        !>
        !> This storage and the distributed one are mutually exclusive, and init_mat refuses
        !> to combine them. Nothing needs both: the only operator that is distributed over the
        !> eigenvector sub-communicator is spin-orbit, which has no components, and the ones
        !> that have components run where that parallelism is ruled out.
        COMPLEX, ALLOCATABLE :: comp(:, :, :, :, :)
        INTEGER :: n_alpha = 1  !> Cartesian extents the operator declared; 1 means comp is unused
        INTEGER :: n_beta  = 1
    CONTAINS
        PROCEDURE :: init_mat    !> provide the result matrix with the dimensions implied by the spinor flags
        PROCEDURE :: add         !> accumulate the matrix elements of another object: this = this + other
        PROCEDURE :: full_matrix !> combine the spinor blocks into a single matrix
        PROCEDURE(calc_matrix_elements_iface), DEFERRED :: calc_matrix_elements
    END TYPE t_matelements

    ABSTRACT INTERFACE
        SUBROUTINE calc_matrix_elements_iface(this, zmat, abc, radfun, usdus)
            IMPORT :: t_matelements, t_mat, t_abc, t_radfun, t_usdus
            CLASS(t_matelements), INTENT(INOUT) :: this
            !> The state at this k-point in as few matrices as it takes: ONE when it is a
            !> whole spinor -- natively so, or stacked by the producer -- and TWO when the
            !> records are independent spin channels. SIZE(zmat) is therefore not jspins,
            !> and a consumer that addresses a spin block by row offset needs the first case.
            TYPE(t_mat),    INTENT(IN) :: zmat(:)
            TYPE(t_abc),    INTENT(IN) :: abc(:,:)  !> (2,ntype) matching coefficients
            TYPE(t_radfun), INTENT(IN) :: radfun(:) !> (ntype) radial functions and their MT integrals
            TYPE(t_usdus),  INTENT(IN) :: usdus     !> values/derivatives at the MT boundary, all types and spins
        END SUBROUTINE calc_matrix_elements_iface
    END INTERFACE

    PUBLIC :: t_matelements

CONTAINS

    SUBROUTINE init_mat(this, num_states, mpi_subcomm, l_real, n_alpha, n_beta)
        !> Provide the result matrix: a 2x2 block matrix in spin space if the
        !> wave functions or the operator are spinors, a single block otherwise.
        !> If mpi_subcomm is given, the blocks are t_mpimat distributed
        !> row-cyclically (full rows local, columns cyclic) over that
        !> communicator, as assumed by the calc_matrix_elements implementations.
        !> If the matrix is already allocated it is only cleared, so the
        !> allocation can be reused for several k-points.
        CLASS(t_matelements), INTENT(INOUT) :: this
        INTEGER, INTENT(IN)           :: num_states  !> global size of each block
        INTEGER, INTENT(IN), OPTIONAL :: n_alpha, n_beta !> Cartesian extents, default 1 each
        INTEGER, INTENT(IN), OPTIONAL :: mpi_subcomm !> distribute the blocks over this communicator
        LOGICAL, INTENT(IN), OPTIONAL :: l_real      !> default: complex matrix elements

        INTEGER :: nsp, i, j, na, nb
        LOGICAL :: l_real_local

        l_real_local = .FALSE.
        IF (PRESENT(l_real)) l_real_local = l_real

        na = 1; IF (PRESENT(n_alpha)) na = n_alpha
        nb = 1; IF (PRESENT(n_beta))  nb = n_beta
        IF ((na > 1 .OR. nb > 1) .AND. PRESENT(mpi_subcomm)) CALL judft_error( &
            "init_mat: a component-carrying operator cannot also be distributed", &
            hint="run it with one rank in the eigenvector sub-communicator", &
            calledby="init_mat")
        this%n_alpha = na; this%n_beta = nb

        !Components live in their own contiguous storage, cleared between k-points like mat.
        IF (na > 1 .OR. nb > 1) THEN
            IF (ALLOCATED(this%comp)) THEN
                this%comp = CMPLX(0.0, 0.0)
            ELSE
                nsp = MERGE(2, 1, this%spinorwavefcts.OR.this%spinoroperator)
                ALLOCATE(this%comp(num_states, num_states, nsp, nsp, na*nb), &
                         source=CMPLX(0.0, 0.0))
            END IF
        END IF

        nsp = MERGE(2, 1, this%spinorwavefcts.OR.this%spinoroperator)

        IF (ALLOCATED(this%mat)) THEN
            !Reuse the existing allocation (do not re-init a t_mpimat, that would
            !recreate its BLACS grid), just clear the data
            IF (SIZE(this%mat,1) /= nsp .OR. SIZE(this%mat,2) /= nsp) &
                CALL judft_bug("init_mat: the allocated matrix does not match the spinor flags")
            DO j = 1, nsp
                DO i = 1, nsp
                    CALL this%mat(i,j)%clear()
                END DO
            END DO
        ELSE
            IF (PRESENT(mpi_subcomm)) THEN
                ALLOCATE(t_mpimat::this%mat(nsp,nsp))
            ELSE
                ALLOCATE(t_mat::this%mat(nsp,nsp))
            END IF
            DO j = 1, nsp
                DO i = 1, nsp
                    !The distribution arguments are ignored by t_mat
                    CALL this%mat(i,j)%init(l_real_local, num_states, num_states, mpi_subcomm, MPIMAT_ROWCYCLIC)
                END DO
            END DO
        END IF
    END SUBROUTINE init_mat

    SUBROUTINE add(this, other)
        !> Accumulate the matrix elements of another object: this = this + other.
        !> Both objects must describe the same spinor structure and their blocks
        !> must match in type, size and distribution.
        CLASS(t_matelements), INTENT(INOUT) :: this
        CLASS(t_matelements), INTENT(IN)    :: other

        INTEGER :: i, j

        IF (this%spinorwavefcts.NEQV.other%spinorwavefcts .OR. this%spinoroperator.NEQV.other%spinoroperator) &
            CALL judft_error("Cannot add matrix elements with different spinor flags.")
        IF (.NOT.(ALLOCATED(this%mat).AND.ALLOCATED(other%mat))) &
            CALL judft_bug("add: the matrix elements have not been calculated yet")

        DO j = 1, SIZE(this%mat, 2)
            DO i = 1, SIZE(this%mat, 1)
                CALL this%mat(i,j)%add(other%mat(i,j))
            END DO
        END DO
    END SUBROUTINE add

    FUNCTION full_matrix(this) RESULT(fullmat)
        !> Combine the spinor blocks into a single matrix, i.e.
        !> ((M_uu,M_ud),(M_du,M_dd)) -> M of twice the size. Distributed blocks
        !> yield a t_mpimat on the same BLACS grid with the same distribution.
        CLASS(t_matelements), INTENT(IN) :: this
        CLASS(t_mat), ALLOCATABLE :: fullmat

        INTEGER :: i, j, nrows, ncols

        IF (.NOT.ALLOCATED(this%mat)) &
            CALL judft_bug("full_matrix: the matrix elements have not been calculated yet")

        SELECT TYPE (blk => this%mat(1,1))
        TYPE IS (t_mpimat)
            nrows = blk%global_size1
            ncols = blk%global_size2
            ALLOCATE(t_mpimat::fullmat)
            !Init from template: same grid and distribution, enlarged global size
            CALL fullmat%init(blk, SIZE(this%mat,1)*nrows, SIZE(this%mat,2)*ncols)
        CLASS DEFAULT
            nrows = blk%matsize1
            ncols = blk%matsize2
            ALLOCATE(fullmat, MOLD=blk)
            CALL fullmat%init(blk%l_real, SIZE(this%mat,1)*nrows, SIZE(this%mat,2)*ncols)
        END SELECT

        !Copy the blocks to their position (global indices); for t_mpimat the
        !copy also performs the required redistribution
        DO j = 1, SIZE(this%mat, 2)
            DO i = 1, SIZE(this%mat, 1)
                CALL fullmat%copy(this%mat(i,j), (i-1)*nrows+1, (j-1)*ncols+1)
            END DO
        END DO
    END FUNCTION full_matrix

END MODULE m_types_matelements
