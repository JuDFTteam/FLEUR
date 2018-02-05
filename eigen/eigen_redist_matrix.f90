!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_eigen_redist_matrix
CONTAINS
  SUBROUTINE eigen_redist_matrix(mpi,lapw,atoms,mat,mat_final)
    USE m_types
    USE m_types_mpimat
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)    :: mpi
    TYPE(t_lapw),INTENT(IN)   :: lapw
    TYPE(t_atoms),INTENT(IN)  :: atoms
    CLASS(t_mat),INTENT(INOUT):: mat(:,:)
    CLASS(t_mat),INTENT(INOUT):: mat_final

    INTEGER:: i,j,m

    SELECT TYPE(mat_final)
    TYPE IS (t_mat)
       m=mat(1,1)%matsize1+mat(2,2)%matsize1
       CALL mat_final%alloc(.FALSE.,m,m)
       !up-up
       mat_final%data_c(:lapw%nv(1)+atoms%nlotot,:lapw%nv(1)+atoms%nlotot)=mat(1,1)%data_c
       !down-down
       mat_final%data_c(lapw%nv(1)+atoms%nlotot+1:,lapw%nv(1)+atoms%nlotot+1:)=mat(2,2)%data_c
       !off-diag
       DO i=1,mat(1,2)%matsize1 !First map U-part of smat&hmat(2,1) into smat(1,2)
          DO j=i+1,mat(1,2)%matsize2
             mat(2,1)%data_c(j,i)=CONJG(mat(1,2)%data_c(i,j))
          ENDDO
       ENDDO
       mat_final%data_c(:lapw%nv(1)+atoms%nlotot,lapw%nv(1)+atoms%nlotot+1:)=mat(2,1)%data_c
    TYPE IS (t_mpimat)
       STOP "ERROR"
!!$       
!!$    !UP/UP spin (or only spin in collinear case)
!!$    m1=mat(1,1)%matsize
!!$    CALL mat_final%create_blacsgrid(mpi%subcom,.FALSE.,m1,m2,1,desc)
!!$    CALL pdgemr2d(m1,m1,mat(1,1)%data_r,1,1,desc,mat_final%data_r,1,1,mat_final%blacs_desc,mat_final%blacs_ctext)
!!$
!!$    IF (nspins==1) RETURN !thats it rest only in noco case
!!$    
!!$    m2=mat(2,2)%matsize
!!$    !DOWN/DOWN
!!$    CALL mat_final%create_blacsgrid(mpi%subcom,.FALSE.,m2,m2,1,desc)
!!$    CALL pdgemr2d(m2,m2,mat(2,2)%data_r,1,1,desc,mat_final%data_r,m1+1,m1+1,mat_final%blacs_desc,mat_final%blacs_ctext)
!!$
!!$    !UP/DOWN
!!$    CALL mat_final%create_blacsgrid(mpi%subcom,.FALSE.,m2,m1,1,desc)
!!$    CALL mat_final%create_blacsgrid(mpi%subcom,.FALSE.,m1,m2,1,desc2)
!!$    
!!$    !first collect the two off-diagonal parts
!!$    CALL pdgeadd('t',m2,m1,1.0,mat(1,2)%data_r,1,1,desc2,1.0,mat(2,1)%data_c,1,1,desc)
!!$
!!$    !Now multiply the diagonal of the matrix by 1/2
!!$    ii=0
!!$    DO i=1,mpi%n_rank+1,MIN(m1,m2),mpi%n_size
!!$       ii=ii+1
!!$       mat(2,1)%data_r(i,ii)=mat(2,1)%data_r(i,ii)/2
!!$    ENDDO
!!$    !Communicate
!!$    CALL pdgemr2d(m2,m1,mat(2,1)%data_r,1,1,desc,mat_final%data_r,1,m1+1,mat_final%blacs_desc,mat_final%blacs_ctext)
    END SELECT
  END SUBROUTINE eigen_redist_matrix
END MODULE m_eigen_redist_matrix
