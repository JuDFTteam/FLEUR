!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_vac_const
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: nvac_mpb

   ! vacua in the mixed basis; for nvac = 1 the second one is the mirror image
   INTEGER, PARAMETER :: NVAC_MPB = 2
END MODULE m_vac_const
