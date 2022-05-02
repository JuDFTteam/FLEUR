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
    CLASS(t_mat),INTENT(IN)    :: mat_tmp
!    CLASS(t_mat),INTENT(INOUT):: mat(:,:)
    CLASS(t_mat),INTENT(INOUT) ::mat(:,:)
    INTEGER:: igSpinPr,igSpin
    DO igSpinPr=1,2
      DO igSpin=1,2
        call mat(igSpinPr,igSpin)%add(mat_tmp,chi(igSpinPr,igSpin))
      enddo
    ENDDO

  END SUBROUTINE hsmt_distspins
END MODULE m_hsmt_distspins
