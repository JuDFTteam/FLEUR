!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------


MODULE m_types_setup
  !*************************************************************
  !     This module contains definitions for all kind of types
  !*************************************************************
  USE m_types_dimension
  USE m_types_cell
  USE m_types_noco
  USE m_types_sym
  use m_types_input
  USE m_types_vacuum
  use m_types_atoms
  USE m_types_stars
  use m_types_sphhar
  use m_types_sliceplot
  USE m_types_banddos
  USE m_types_wannier
  USE m_types_corespecInput

  !Not converted to fleur_setup type yet
  use m_types_oneD


  TYPE t_hybrid
     LOGICAL               ::  l_hybrid=.false.
     LOGICAL               ::  l_subvxc=.false.
     LOGICAL               ::  l_calhf=.false.
     LOGICAL               ::  l_addhf=.false.
     INTEGER               ::  ewaldlambda
     INTEGER               ::  lexp
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
     INTEGER,ALLOCATABLE   ::  nindx(:,:)
     INTEGER,ALLOCATABLE   ::  select1(:,:)
     INTEGER,ALLOCATABLE   ::  lcutm1(:)
     INTEGER,ALLOCATABLE   ::  nindxm1(:,:)
     INTEGER,ALLOCATABLE   ::  gptm(:,:)
     INTEGER,ALLOCATABLE   ::  ngptm1(:)
     INTEGER,ALLOCATABLE   ::  pgptm1(:,:)
     INTEGER,ALLOCATABLE   ::  ngptm (:)
     INTEGER,ALLOCATABLE   ::  pgptm (:,:)
     INTEGER,ALLOCATABLE   ::  lcutwf(:)
     INTEGER,ALLOCATABLE   ::  map(:,:)
     INTEGER,ALLOCATABLE   ::  tvec(:,:,:)
     INTEGER , ALLOCATABLE ::  nbasm(:)                                     
     REAL                  ::  gcutm1
     REAL                  ::  tolerance1  !only read in
     REAL   ,ALLOCATABLE   ::  basm1(:,:,:,:)
     COMPLEX,ALLOCATABLE   ::  d_wgn2(:,:,:,:)
     INTEGER,ALLOCATABLE   ::  ne_eig(:),nbands(:),nobd(:)                   !alloc in eigen_HF_init
     REAL   ,ALLOCATABLE   ::  div_vv(:,:,:)
  END TYPE t_hybrid


 
  TYPE t_obsolete
     INTEGER:: lepr !floating energy parameters...
     INTEGER:: ndvgrd !remove
     REAL   :: chng   !remove
     LOGICAL :: lwb   !remove
  END TYPE t_obsolete



 
END MODULE m_types_setup
