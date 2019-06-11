!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_sliceplot
    TYPE t_sliceplot
     LOGICAL :: iplot=.false.
     LOGICAL :: slice=.false.
     LOGICAL :: plpot=.false.
     INTEGER :: kk=0
     INTEGER :: nnne=0
     REAL    :: e1s=0.
     REAL    :: e2s=0.
  END TYPE t_sliceplot
END MODULE m_types_sliceplot
