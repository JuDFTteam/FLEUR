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
        CLASS(t_mat), POINTER :: mat(:,:) !> 2x2 spinor matrix representation
    CONTAINS
        PROCEDURE(add_matrix_elements_iface), DEFERRED :: add_matrix_elements
    END TYPE t_matelements

    ABSTRACT INTERFACE
        SUBROUTINE add_matrix_elements_iface(this, zmat, mat_inout)
            IMPORT :: t_matelements, t_mat
            CLASS(t_matelements),           INTENT(INOUT)         :: this
            TYPE(t_mat),                    INTENT(IN)            :: zmat(:)
            CLASS(t_mat), OPTIONAL, TARGET, INTENT(INOUT)         :: mat_inout(:,:)
        END SUBROUTINE add_matrix_elements_iface
    END INTERFACE

    PUBLIC :: t_matelements

END MODULE m_types_matelements
