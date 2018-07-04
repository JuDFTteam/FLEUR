!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_symmetrize_matrix
  USE m_juDFT

CONTAINS
  SUBROUTINE symmetrize_matrix(mpi,noco,kpts,nk,hmat,smat)
    USE m_types
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)     :: mpi
    TYPE(t_noco),INTENT(in)    :: noco
    TYPE(t_kpts),INTENT(in)    :: kpts
    INTEGER,INTENT(in)         :: nk
    CLASS(t_mat),INTENT(inout) :: hmat,smat

    REAL :: max_imag
    !Check if we could exploit a real matrix even without inversion symmetry
    realcomplex:IF (.NOT.noco%l_noco.AND..NOT.hmat%l_real) THEN
       IF (ALL(ABS(kpts%bk(:,nk))<1E-10)) THEN
          max_imag=MAXVAL(ABS(AIMAG(hmat%data_c)))
          IF (max_imag>1e-10) THEN
                     PRINT *,"Real matrix expected but imaginary part is:",max_imag
                     RETURN
          ENDIF
          
          IF (mpi%irank==0) THEN
             PRINT *,"Complex matrix made real"
             WRITE(6,*) "Complex matrix made real"
          END IF
             
          !We are using Gamma point, so matrix should be real
          IF (ALLOCATED(smat%data_r)) DEALLOCATE(smat%data_r)
          ALLOCATE(smat%data_r(SIZE(smat%data_c,1),SIZE(smat%data_c,2)))
          smat%data_r=smat%data_c;smat%l_real=.TRUE.
          DEALLOCATE(smat%data_c)

          IF (ALLOCATED(hmat%data_r)) DEALLOCATE(hmat%data_r)
          ALLOCATE(hmat%data_r(SIZE(hmat%data_c,1),SIZE(hmat%data_c,2)))
          hmat%data_r=hmat%data_c;hmat%l_real=.TRUE.
          DEALLOCATE(hmat%data_c)

       ENDIF
    ENDIF realcomplex

  END SUBROUTINE symmetrize_matrix
END MODULE m_symmetrize_matrix
