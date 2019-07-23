MODULE m_types_hybrid
   IMPLICIT NONE

   TYPE t_hybrid
      LOGICAL                ::  l_hybrid = .false.
      LOGICAL                ::  l_subvxc = .false.
      LOGICAL                ::  l_calhf = .false.
      LOGICAL                ::  l_addhf = .false.
      INTEGER                ::  ewaldlambda
      INTEGER                ::  lexp = 0
      INTEGER                ::  bands1 !Only read in
      INTEGER                ::  nbasp
      INTEGER                ::  maxlcutm1
      INTEGER                ::  maxindxm1
      INTEGER                ::  maxbasm1
      INTEGER                ::  maxindxp1
      INTEGER                ::  maxgptm
      INTEGER                ::  maxgptm1
      INTEGER                ::  maxindx
      INTEGER                ::  maxlmindx
      INTEGER                ::  gptmd
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
      INTEGER, ALLOCATABLE   ::  nbasm(:)
      REAL                   ::  gcutm1
      REAL                   ::  tolerance1  !only read in
      REAL, ALLOCATABLE      ::  basm1(:, :, :, :)
      COMPLEX, ALLOCATABLE   ::  d_wgn2(:, :, :, :)
      INTEGER, ALLOCATABLE   ::  ne_eig(:), nbands(:), nobd(:)  !alloc in eigen_HF_init
      REAL, ALLOCATABLE      ::  div_vv(:, :, :)
   END TYPE t_hybrid
   
   TYPE t_prodtype
      INTEGER :: l1,l2,n1,n2
   END TYPE t_prodtype

   TYPE t_hybdat
      INTEGER                     :: lmaxcd, maxindxc
      INTEGER                     :: maxfac
      REAL,           ALLOCATABLE :: gridf(:,:)     !alloc in util.F
      INTEGER,        ALLOCATABLE :: nindxc(:,:)     !alloc in eigen_HF_init
      INTEGER,        ALLOCATABLE :: lmaxc(:)      !alloc in eigen_HF_init
      REAL,           ALLOCATABLE :: core1(:,:,:,:),          core2(:, :, :, :) !alloc in eigen_HF_init
      REAL,           ALLOCATABLE :: eig_c(:,:,:)  !alloc in eigen_HF_init
      INTEGER,        ALLOCATABLE :: kveclo_eig(:,:)     !alloc in eigen_HF_setup
      REAL,           ALLOCATABLE :: sfac(:), fac(:) !alloc in eigen_HF_init
      REAL,           ALLOCATABLE :: gauntarr(:,:,:,:,:,:) !alloc in eigen_HF_init
      REAL,           ALLOCATABLE :: bas1(:,:,:,:), bas2(:,:,:,:) !alloc in eigen_HF_init
      REAL,           ALLOCATABLE :: bas1_MT(:,:,:), drbas1_MT(:, :,:) !alloc in eigen_HF_init
      REAL,           ALLOCATABLE :: prodm(:,:,:,:)           !alloc in eigen_HF_setup
      TYPE(t_PRODTYPE), ALLOCATABLE :: prod(:,:,:)  !alloc in eigen_HF_setup
      INTEGER,        ALLOCATABLE :: pntgptd(:)    !alloc in eigen_HF_setup
      INTEGER,        ALLOCATABLE :: pntgpt(:,:,:,:)           !alloc in eigen_HF_setup
      INTEGER,        ALLOCATABLE :: nindxp1(:,:)
      REAL,           ALLOCATABLE :: stepfunc_r(:,:,:)
      COMPLEX,ALLOCATABLE ::  stepfunc_c(:,:,:)
   END TYPE t_hybdat

CONTAINS
END MODULE m_types_hybrid
