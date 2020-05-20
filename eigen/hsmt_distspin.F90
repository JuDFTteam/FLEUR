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
!    CLASS(t_mat),INTENT(INOUT):: mat(:,:)
    TYPE(t_mat),INTENT(INOUT) ::mat(:,:)
    INTEGER:: iintsp,jintsp,i,j
    !$acc parallel copyin(chi) present(mat%data_c,mat_tmp)
    DO iintsp=1,2
       DO jintsp=1,2
	!$acc loop independent collapse(2)	
          DO i=1,size(mat_tmp%data_c,1);do j=1,size(mat_tmp%data_c,2)	
          mat(jintsp,iintsp)%data_c(i,j)=chi(jintsp,iintsp)*mat_tmp%data_c(i,j)+mat(jintsp,iintsp)%data_c(i,j)
          enddo;enddo
        !$acc end loop
       ENDDO
    ENDDO
    !$acc end parallel
  END SUBROUTINE hsmt_distspins
END MODULE m_hsmt_distspins
