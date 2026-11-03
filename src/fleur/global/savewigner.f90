!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_savewigner
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: d_wgn
      COMPLEX, ALLOCATABLE :: d_wgn(:,:,:,:)
      END MODULE m_savewigner
