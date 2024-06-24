!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_dfpt_eigen_redist_matrix
CONTAINS
   SUBROUTINE dfpt_eigen_redist_matrix(fmpi,lapwq,lapw,atoms,mat,mat_final,mat_final_templ)
      !> Collect Hamiltonian or overlap matrix perturbation to final form
      !!
      !! In the collinear case, this routine just copies mat(1,1) into the final matrix.
      !! In the non-collinear case, the 2x2 array of matrices is combined into the final matrix.
      USE m_types
      USE m_types_mpimat

      IMPLICIT NONE

      TYPE(t_mpi),   INTENT(IN)    :: fmpi
      TYPE(t_lapw),  INTENT(IN)    :: lapwq, lapw
      TYPE(t_atoms), INTENT(IN)    :: atoms
      CLASS(t_mat),  INTENT(INOUT) :: mat(:,:)
      CLASS(t_mat),  INTENT(INOUT) :: mat_final

      CLASS(t_mat),  INTENT(IN), OPTIONAL :: mat_final_templ

      INTEGER :: mPr, m

      !  Determine final matrix size and allocate the final matrix
      m = lapw%nv(1) + atoms%nlotot
      IF (SIZE(mat)>1)  m = m + lapw%nv(2) + atoms%nlotot
      mPr = lapwq%nv(1) + atoms%nlotot
      IF (SIZE(mat)>1)  mPr = mPr + lapwq%nv(2) + atoms%nlotot

      IF (.NOT.PRESENT(mat_final_templ)) THEN
         CALL mat_final%init(mat(1,1)%l_real,mPr,m,fmpi%diag_sub_comm,.TRUE.) !here the .true. creates a block-cyclic scalapack distribution
      ELSE
         CALL mat_final%init(mat_final_templ)
      END IF

      !Collinear case --> only component
      IF (SIZE(mat)==1) THEN
         CALL mat_final%move(mat(1,1))
         CALL mat(1,1)%free()
         RETURN
      END IF

      !up-up
      CALL mat_final%copy(mat(1,1),1,1)
      CALL mat(1,1)%free()

      !down-down component
      CALL mat_final%copy(mat(2,2),lapwq%nv(1)+atoms%nlotot+1,lapw%nv(1)+atoms%nlotot+1)
      CALL mat(2,2)%free()

      !Now collect the off-diagonal parts
      !down-up
      CALL mat_final%copy(mat(2,1),lapwq%nv(1)+atoms%nlotot+1,1)
      CALL mat(2,1)%free()

      !up-down
      CALL mat_final%copy(mat(1,2),1,lapw%nv(1)+atoms%nlotot+1)
      CALL mat(1,2)%free()

   END SUBROUTINE dfpt_eigen_redist_matrix
END MODULE m_dfpt_eigen_redist_matrix
