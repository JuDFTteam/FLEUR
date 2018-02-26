!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_eigen_redist_matrix
CONTAINS
  !> Collect Hamiltonian or overlap matrix to final form
  !!
  !! In the collinear case, this routine just copies mat(1,1) into the final matrix.
  !! If the matrices are distributed, the copy includes a redistribution into the block-cylic form needed by
  !! the diagonalization.
  !! In the non-collinear case, the 2x2 array of matrices is combined into the final matrix. Again a redistribution will happen in the parallel case
 
  
  SUBROUTINE eigen_redist_matrix(mpi,lapw,atoms,mat,mat_final)
    USE m_types
    USE m_types_mpimat
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)    :: mpi
    TYPE(t_lapw),INTENT(IN)   :: lapw
    TYPE(t_atoms),INTENT(IN)  :: atoms
    CLASS(t_mat),INTENT(INOUT):: mat(:,:)
    CLASS(t_mat),INTENT(INOUT):: mat_final

    INTEGER:: m

    !determine final matrix size and allocate the final matrix
    m=lapw%nv(1)+atoms%nlotot
    IF (SIZE(mat)>1) m=m+lapw%nv(2)+atoms%nlotot
    CALL mat_final%init(mat(1,1)%l_real,m,m,mpi%sub_comm,.TRUE.) !here the .true. creates a block-cyclic scalapack distribution
    
    !up-up component (or only component in collinear case)
    CALL mat_final%copy(mat(1,1),1,1)
    CALL mat(1,1)%free()
    IF (SIZE(mat)==1) RETURN

    !down-down component
    CALL mat_final%copy(mat(2,2),lapw%nv(1)+atoms%nlotot+1,lapw%nv(1)+atoms%nlotot+1)
    CALL mat(2,2)%free()

    !Now collect off-diagonal parts
    CALL mat(2,1)%add_transpose(mat(1,2))
    CALL mat_final%copy(mat(2,1),1,lapw%nv(1)+atoms%nlotot+1)
    CALL mat(1,2)%free()
    CALL mat(2,1)%free()
    
  END SUBROUTINE eigen_redist_matrix
END MODULE m_eigen_redist_matrix


