!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_kpts

  
  TYPE t_kpts
     INTEGER :: specificationType
     !no
     INTEGER :: nkpt
     INTEGER :: ntet
     REAL    :: posScale
     LOGICAL :: l_gamma
     !(3,nkpt) k-vectors internal units
     REAL,ALLOCATABLE ::bk(:,:)
     !(nkpts) weights
     REAL,ALLOCATABLE ::wtkpt(:)
     INTEGER               ::  nkptf !<k-vectors in full BZ
     INTEGER               ::  nkpt3(3)
     REAL                  ::  kPointDensity(3) ! only used if k point set is defined as density
     REAL   ,ALLOCATABLE   ::  bkf(:,:)
     INTEGER,ALLOCATABLE   ::  bkp(:)
     INTEGER,ALLOCATABLE   ::  bksym(:)
     INTEGER                       :: numSpecialPoints
     CHARACTER(LEN=50),ALLOCATABLE :: specialPointNames(:)
     REAL   ,ALLOCATABLE           :: specialPoints(:,:)
     INTEGER,ALLOCATABLE           :: ntetra(:,:)
     REAL   ,ALLOCATABLE           :: voltet(:)
  ENDTYPE t_kpts

 
END MODULE m_types_kpts
