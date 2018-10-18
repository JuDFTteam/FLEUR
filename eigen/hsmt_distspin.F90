!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hsmt_distspins
  IMPLICIT NONE
CONTAINS
  SUBROUTINE hsmt_distspins(chi,mat_tmp,mat)
    USE m_types
    COMPLEX,INTENT(in)        :: chi(2,2)
    TYPE(t_mat),INTENT(IN)    :: mat_tmp
    CLASS(t_mat),INTENT(INOUT):: mat(2,2)

    INTEGER:: iintsp,jintsp

    DO iintsp=1,2
       DO jintsp=1,2
          mat(jintsp,iintsp)%data_c(:,:)=chi(iintsp,jintsp)*mat_tmp%data_c(:,:)+mat(jintsp,iintsp)%data_c(:,:)
       ENDDO
    ENDDO
  END SUBROUTINE hsmt_distspins
END MODULE m_hsmt_distspins
