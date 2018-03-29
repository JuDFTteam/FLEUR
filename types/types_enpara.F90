!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_enpara
   TYPE t_enpara
     REAL, ALLOCATABLE :: el0(:,:,:)
     REAL, ALLOCATABLE :: evac0(:,:)
     REAL, ALLOCATABLE :: ello0(:,:,:)
     REAL, ALLOCATABLE :: enmix(:)
     INTEGER, ALLOCATABLE :: skiplo(:,:)
     LOGICAL, ALLOCATABLE :: lchange(:,:,:)
     LOGICAL, ALLOCATABLE :: lchg_v(:,:)
     LOGICAL, ALLOCATABLE :: llochg(:,:,:)
     REAL                 :: epara_min
  END TYPE t_enpara
END MODULE m_types_enpara
