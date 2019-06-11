!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_hybrid
     TYPE t_hybrid
      LOGICAL               ::  l_hybrid = .false.
      LOGICAL               ::  l_subvxc = .false.
      LOGICAL               ::  l_calhf = .false.
      LOGICAL               ::  l_addhf = .false.
      INTEGER               ::  ewaldlambda =3
      INTEGER               ::  lexp =16
      INTEGER               ::  bands1 !Only read in
      INTEGER               ::  nbasp
      INTEGER               ::  maxlcutm1
      INTEGER               ::  maxindxm1
      INTEGER               ::  maxbasm1
      INTEGER               ::  maxindxp1
      INTEGER               ::  maxgptm
      INTEGER               ::  maxgptm1
      INTEGER               ::  maxindx
      INTEGER               ::  maxlmindx
      INTEGER               ::  gptmd
      INTEGER, ALLOCATABLE   ::  nindx(:, :)
      INTEGER, ALLOCATABLE   ::  select1(:, :)
      INTEGER, ALLOCATABLE   ::  lcutm1(:)
      INTEGER, ALLOCATABLE   ::  nindxm1(:, :)
      INTEGER, ALLOCATABLE   ::  gptm(:, :)
      INTEGER, ALLOCATABLE   ::  ngptm1(:)
      INTEGER, ALLOCATABLE   ::  pgptm1(:, :)
      INTEGER, ALLOCATABLE   ::  ngptm(:)
      INTEGER, ALLOCATABLE   ::  pgptm(:, :)
      INTEGER, ALLOCATABLE   ::  lcutwf(:)
      INTEGER, ALLOCATABLE   ::  map(:, :)
      INTEGER, ALLOCATABLE   ::  tvec(:, :, :)
      INTEGER, ALLOCATABLE ::  nbasm(:)
      REAL                  ::  gcutm1
      REAL                  ::  tolerance1 = 1e-4 !only read in
      REAL, ALLOCATABLE   ::  basm1(:, :, :, :)
      COMPLEX, ALLOCATABLE   ::  d_wgn2(:, :, :, :)
      INTEGER, ALLOCATABLE   ::  ne_eig(:), nbands(:), nobd(:)                   !alloc in eigen_HF_init
      REAL, ALLOCATABLE   ::  div_vv(:, :, :)
   END TYPE t_hybrid
 END MODULE m_types_hybrid
