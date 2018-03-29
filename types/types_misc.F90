!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_misc
  !*************************************************************
  !     This module contains definitions for all kind of types
  !*************************************************************
  !
  ! Types for orbital moment calculation:
  !
  !
  ! Types for spin-off-diagonal charge density:
  !
  TYPE t_mt21                          ! 'normal' contributions
     SEQUENCE
     REAL ::  uun,udn,dun,ddn           ! normes of radial overlaps
     COMPLEX :: uu,ud,du,dd             ! values
  END TYPE t_mt21

  TYPE t_lo21                          ! ocal orbitals & (u,d)
     SEQUENCE
     REAL ::  uulon,dulon,uloun,ulodn   ! normes of radial overlaps
     COMPLEX :: uulo,dulo,ulou,ulod     ! values
  END TYPE t_lo21

  TYPE t_orb                           ! 'normal' contributions
     SEQUENCE
     REAL :: uu,dd                      ! z   component
     COMPLEX :: uup,uum,ddp,ddm         ! +/- component
  END TYPE t_orb

  TYPE t_orbl                          ! local orbitals & (u,d)
     SEQUENCE
     REAL :: uulo,dulo
     COMPLEX :: uulop,uulom,dulop,dulom
  END TYPE t_orbl

  TYPE t_orblo                         ! lo,lo' contributions
     SEQUENCE
     REAL :: z
     COMPLEX :: p,m
  END TYPE t_orblo
  
  !
  ! Type for the HF total energy
  !
  TYPE t_energy_hf
     REAL :: valence
     REAL :: core
  END TYPE t_energy_hf

  
  TYPE prodtype
     INTEGER :: l1,l2,n1,n2
  END TYPE prodtype
  TYPE t_hybdat
     INTEGER              :: lmaxcd,maxindxc
     REAL,  ALLOCATABLE   ::  gridf(:,:)                                    !alloc in util.F
     INTEGER , ALLOCATABLE::  nindxc(:,:)                                   !alloc in eigen_HF_init
     INTEGER,ALLOCATABLE  :: lmaxc(:)                                       !alloc in eigen_HF_init
     REAL,    ALLOCATABLE ::  core1(:,:,:,:),core2(:,:,:,:)                 !alloc in eigen_HF_init
     REAL,    ALLOCATABLE ::  eig_c(:,:,:)                                  !alloc in eigen_HF_init
     INTEGER , ALLOCATABLE::  kveclo_eig(:,:)                               !alloc in eigen_HF_setup
     INTEGER              ::  maxfac
     REAL,    ALLOCATABLE ::  sfac(:),fac(:)                                !alloc in eigen_HF_init
     REAL,    ALLOCATABLE ::  gauntarr(:,:,:,:,:,:)                         !alloc in eigen_HF_init
     REAL,    ALLOCATABLE ::  bas1(:,:,:,:),bas2(:,:,:,:)                   !alloc in eigen_HF_init
     REAL ,   ALLOCATABLE ::  bas1_MT(:,:,:),drbas1_MT(:,:,:)               !alloc in eigen_HF_init
     REAL, ALLOCATABLE    ::  prodm(:,:,:,:)                                !alloc in eigen_HF_setup
     TYPE(PRODTYPE),ALLOCATABLE :: prod(:,:,:)                              !alloc in eigen_HF_setup
     INTEGER, ALLOCATABLE :: pntgptd(:)                                     !alloc in eigen_HF_setup
     INTEGER, ALLOCATABLE :: pntgpt(:,:,:,:)                                !alloc in eigen_HF_setup
     INTEGER,ALLOCATABLE   ::  nindxp1(:,:)
  END TYPE t_hybdat


  TYPE t_results
     REAL, ALLOCATABLE :: force(:,:,:)   !< Forces calculated on all atoms (for each spin)
     REAL, ALLOCATABLE :: force_old(:,:) !< Forces on all atoms from last iteration
     REAL              :: ef        !<Fermie energy
     REAL              :: seigc     !<sum of the core eigenvalues
     REAL              :: seigsc    !<weighted sum of the semi-core eigenvalues
     REAL              :: seigv     !<weighted sum of the occupied valence eigenvalues
     REAL              :: seigscv   !<sum of seigv and seigsc
     REAL              :: ts        !<entropy contribution to the free energy
     REAL              :: te_vcoul  !<charge density-coulomb potential integral
     REAL              :: te_veff   !<charge density-effective potential integral
     REAL              :: te_exc    !<charge density-ex-corr.energy density integral
     REAL              :: e_ldau    !<total energy contribution of LDA+U
     REAL              :: tote
     REAL              :: last_distance
     REAL              :: bandgap
     TYPE(t_energy_hf) ::  te_hfex
     REAL              ::  te_hfex_loc(2)
     REAL, ALLOCATABLE :: w_iks(:,:,:)
  END TYPE t_results


 
  TYPE t_zMat
     LOGICAL              :: l_real
     INTEGER              :: nbasfcn
     INTEGER              :: nbands
     REAL,    ALLOCATABLE :: z_r(:,:) ! z_r(nbasfcn,nbands)
     COMPLEX, ALLOCATABLE :: z_c(:,:) ! z_c(nbasfcn,nbands)
  END TYPE t_zMat

  TYPE t_hamOvlp
     LOGICAL              :: l_real
     INTEGER              :: matsize
     REAL,    ALLOCATABLE :: a_r(:), b_r(:)
     COMPLEX, ALLOCATABLE :: a_c(:), b_c(:)
  END TYPE t_hamOvlp
 
 
END MODULE m_types_misc
