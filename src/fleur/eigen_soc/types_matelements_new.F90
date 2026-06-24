!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_matelements
   USE m_types_mat
   IMPLICIT NONE
   PRIVATE

   TYPE, ABSTRACT :: t_matelements
      LOGICAL :: spinorwavefcts = .FALSE. !> whether the matrix elements are for spinor wave functions
      LOGICAL :: spinoroperator = .FALSE. !> whether the matrix elements are for spinor operators
      CLASS(t_mat), allocatable :: mat(:, :) !> 2x2 spinor matrix representation
   CONTAINS
      procedure  :: init   !init the matrix, set the spinor flags
      procedure  :: add    !add two matrix elements objects
      operator(+) :: add   !Make it simpler to add two matrix elements objects
      procedure  :: full_matrix !get the full matrix representation of the matrix elements
      PROCEDURE(calc_matrix_elements_iface), DEFERRED :: calc_matrix_elements
   END TYPE t_matelements

   ABSTRACT INTERFACE
      SUBROUTINE calc_matrix_elements_iface(this, zmat)
         IMPORT :: t_matelements, t_mat
         CLASS(t_matelements), INTENT(INOUT)         :: this
         TYPE(t_mat), INTENT(IN)            :: zmat(:)
      END SUBROUTINE calc_matrix_elements_iface
   END INTERFACE
   PUBLIC :: t_matelements
   CONTAINS
    SUBROUTINE init(this, spinorwf, spinorop,l_distributed)
        CLASS(t_matelements), INTENT(INOUT) :: this
        LOGICAL, INTENT(IN) :: spinorwf, spinorop,l_distributed
        integer:: size
        this%spinorwavefcts = spinorwf
        this%spinoroperator = spinorop
        size=merge(2,1,spinorwf .or. spinorop)
        if (l_distributed) then
           allocate(t_mpimat: this%mat(size, size))
        else
           allocate(t_mat:this%mat(size, size))
        end if
    END SUBROUTINE init
    SUBROUTINE add(this, other)
        CLASS(t_matelements), INTENT(INOUT) :: this
        CLASS(t_matelements), INTENT(IN)    :: other
        integer:: i,j
        if (this%spinorwavefcts .ne. other%spinorwavefcts .or. this%spinoroperator .ne. other%spinoroperator) &
           call judft_error("Cannot add matrix elements with different spinor flags.")
       
       DO i = 1, SIZE(this%mat, 1)
            DO j = 1, SIZE(this%mat, 2)
                this%mat(i, j) = this%mat(i, j) + other%mat(i, j)
            END DO
        END DO          
    END SUBROUTINE add
    FUNCTION full_matrix(this) RESULT(fullmat)
        CLASS(t_matelements), INTENT(IN) :: this
        CLASS(t_mat) :: fullmat(:)
        
        integer :: i, j, idx
        integer :: nrows, ncols
        nrows = SIZE(this%mat(1, 1)%data, 1)
        ncols = SIZE(this%mat(1, 1)%data, 2)
        allocate(fullmat(nrows * SIZE(this%mat, 1), ncols * SIZE(this%mat, 2)))
        fullmat = 0.0
        DO i = 1, SIZE(this%mat, 1)
            DO j = 1, SIZE(this%mat, 2)
                idx = (i - 1) * nrows + 1
                fullmat(idx:idx + nrows - 1, (j - 1) * ncols + 1:(j - 1) * ncols + ncols) = this%mat(i, j)%data
            END DO
        END DO
    END FUNCTION full_matrix
END MODULE m_types_matelements
