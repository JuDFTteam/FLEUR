!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!<@brief
  !<This matrix_copy module is needed to copy a distributed into a non-distributed matrix or vice versa
  !<It calls the usual matrix-copy routines in other cases
MODULE m_matrix_copy
  IMPLICIT NONE
CONTAINS
  SUBROUTINE matrix_copy(mat_in,mat_out)
    USE m_types
    USE m_types_mpimat

    CLASS(t_mat),INTENT(IN)   ::mat_in
    CLASS(t_mat),INTENT(INOUT)::mat_out

    TYPE(t_mpimat)::tmp_mat

#ifdef CPP_MPI
    INCLUDE 'mpif.h'
#else
    INTEGER:: MPI_COMM_SELF
#endif
    
    SELECT TYPE(mat_in)
    CLASS is(t_mpimat)
       SELECT TYPE(mat_out)
       CLASS is(t_mpimat)
          CALL mat_out%copy(mat_in,1,1)
       CLASS is(t_mat)
          !Copy from t_mpimat to t_mat
          CALL tmp_mat%init(mat_in%l_real,mat_in%global_size1,mat_in%global_size2,MPI_COMM_SELF)
          CALL tmp_mat%copy(mat_in,1,1) !redistribute to single matrix
          IF (tmp_mat%l_real) THEN
             CALL move_ALLOC(tmp_mat%data_r,mat_out%data_r)
          ELSE
             CALL move_ALLOC(tmp_mat%data_r,mat_out%data_r)
          ENDIF
       END SELECT
    CLASS is (t_mat)
       SELECT TYPE(mat_out)
       CLASS is(t_mpimat)
          !Copy from non-distributed t_mat to t_mpimat
          CALL tmp_mat%init(mat_in%l_real,mat_in%matsize1,mat_in%matsize1,MPI_COMM_SELF)
          IF (tmp_mat%l_real) THEN
             tmp_mat%data_r=mat_in%data_r
          ELSE
             tmp_mat%data_c=mat_in%data_c
          ENDIF
          CALL mat_out%copy(tmp_mat,1,1)
       CLASS is(t_mat)
          CALL mat_out%copy(mat_in,1,1)
       END SELECT
    END SELECT
  END SUBROUTINE matrix_copy
END MODULE m_matrix_copy
