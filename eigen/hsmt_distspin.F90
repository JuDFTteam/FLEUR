!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hsmt_distspins
  IMPLICIT NONE
CONTAINS
  SUBROUTINE hsmt_distspins(chi,mat_tmp,mat)
#include"cpp_double.h"
    USE m_types
    COMPLEX,INTENT(in)        :: chi(2,2)
    TYPE(t_mat),INTENT(IN)    :: mat_tmp
!    CLASS(t_mat),INTENT(INOUT):: mat(:,:)
    TYPE(t_mat),INTENT(INOUT) ::mat(:,:)
    INTEGER:: iintsp,jintsp,i,j
#ifdef _OPENACC
    DO iintsp=1,2
       DO jintsp=1,2
       !$acc parallel loop collapse(2) copyin(chi) present(mat_tmp,mat_tmp%data_c) &
          !$acc present(mat(jintsp,iintsp),mat(jintsp,iintsp)%data_c)
          DO j=1,SIZE(mat_tmp%data_c,2)
             DO i=1,SIZE(mat_tmp%data_c,1)
                mat(jintsp,iintsp)%data_c(i,j)=chi(jintsp,iintsp)*mat_tmp%data_c(i,j)+mat(jintsp,iintsp)%data_c(i,j)
             ENDDO
          ENDDO
          !$acc end parallel loop
       ENDDO
    ENDDO
#else
    !$OMP PARALLEL DO PRIVATE(j,iintsp,jintsp) SHARED(mat_tmp,chi,mat) DEFAULT(NONE)
    DO j=1,mat_tmp%matsize2
      DO iintsp=1,2
          DO jintsp=1,2
             call CPP_BLAS_caxpy(mat_tmp%matsize1,chi(jintsp,iintsp),mat_tmp%data_c(:,j),1,mat(jintsp,iintsp)%data_c(:,j),1)
          ENDDO
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO
#endif
  END SUBROUTINE hsmt_distspins
END MODULE m_hsmt_distspins
