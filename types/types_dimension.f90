!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_dimension
   TYPE t_dimension
      INTEGER :: nvd
      INTEGER :: nv2d
      INTEGER :: neigd !to be removed!!!
      INTEGER :: lmd
      INTEGER :: nbasfcn
   END TYPE t_dimension
 END MODULE m_types_dimension
