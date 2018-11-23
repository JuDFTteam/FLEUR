!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------


MODULE m_types_setup
  !*************************************************************
  !     This module contains definitions for all kind of types
  !*************************************************************
  USE m_types_cell
  USE m_types_noco
  USE m_types_sym
  use m_types_input
  USE m_types_vacuum
  use m_types_atoms
  USE m_types_stars
  use m_types_sphhar

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


  TYPE t_dimension
     INTEGER :: nspd
     INTEGER :: nvd
     INTEGER :: nv2d
     INTEGER :: neigd
     INTEGER :: neigd2
     INTEGER :: ncvd
     INTEGER :: nstd
     INTEGER :: msh
     INTEGER :: lmd
     INTEGER :: lmplmd
     INTEGER :: nbasfcn
  END TYPE t_dimension


  TYPE t_sliceplot
     LOGICAL :: iplot
     LOGICAL :: slice
     LOGICAL :: plpot
     LOGICAL :: pallst
     INTEGER :: kk
     INTEGER :: nnne
     REAL    :: e1s
     REAL    :: e2s
  END TYPE t_sliceplot

  TYPE t_banddos
     LOGICAL :: dos
     LOGICAL :: band
     LOGICAL :: l_mcd
     LOGICAL :: l_orb
     LOGICAL :: vacdos
     INTEGER :: ndir
     INTEGER :: orbCompAtom
     REAL    :: e1_dos
     REAL    :: e2_dos
     REAL    :: sig_dos
     REAL    :: e_mcd_lo
     REAL    :: e_mcd_up
     LOGICAL :: unfoldband
     INTEGER :: s_cell_x
     INTEGER :: s_cell_y
     INTEGER :: s_cell_z     
  END TYPE t_banddos

  TYPE t_obsolete
     INTEGER:: lepr !floating energy parameters...
     INTEGER:: ndvgrd !remove
     REAL   :: chng   !remove
     LOGICAL :: lwb   !remove
  END TYPE t_obsolete


 ! type for the input to the calculation of the core spectrum (EELS)
  TYPE t_coreSpecInput
     integer :: verb  ! output verbosity
     integer :: atomType  ! atomic type used for calculation of core spectra
     character(LEN=1) :: edge  ! edge character (K,L,M,N,O,P)
     integer :: edgeidx(11)  ! l-j edges
     integer :: lx  ! maximum lmax considered in spectra calculation
     real :: ek0  ! kinetic energy of incoming electrons
     real :: emn  ! energy spectrum lower bound
     real :: emx  ! energy spectrum upper bound
     real :: ein  ! energy spectrum increment
  END TYPE t_coreSpecInput


  !
  ! type for wannier-functions
  !
  TYPE t_wann
     INTEGER :: wan90version
     INTEGER :: oc_num_orbs
     INTEGER,ALLOCATABLE :: oc_orbs(:)
     LOGICAL :: l_unformatted
     LOGICAL :: l_oc_f
     LOGICAL :: l_ndegen
     LOGICAL :: l_orbitalmom
     LOGICAL :: l_orbcomp
     LOGICAL :: l_orbcomprs
     LOGICAL :: l_denmat
     LOGICAL :: l_perturbrs
     LOGICAL :: l_perturb
     LOGICAL :: l_nedrho
     LOGICAL :: l_anglmomrs
     LOGICAL :: l_anglmom
     LOGICAL :: l_spindisp
     LOGICAL :: l_spindisprs
     LOGICAL :: l_socspicom
     LOGICAL :: l_socspicomrs
     LOGICAL :: l_offdiposoprs
     LOGICAL :: l_offdiposop
     LOGICAL :: l_torque
     LOGICAL :: l_torquers
     LOGICAL :: l_atomlist
     INTEGER :: atomlist_num
     INTEGER,ALLOCATABLE :: atomlist(:)
     LOGICAL :: l_berry
     LOGICAL :: l_perpmagrs
     LOGICAL :: l_perpmag
     LOGICAL :: l_perpmagat
     LOGICAL :: l_perpmagatrs
     LOGICAL :: l_socmatrs
     LOGICAL :: l_socmat
     LOGICAL :: l_soctomom
     LOGICAL :: l_kptsreduc2
     LOGICAL :: l_nablapaulirs
     LOGICAL :: l_nablars
     LOGICAL :: l_surfcurr
     LOGICAL :: l_updown
     LOGICAL :: l_ahe
     LOGICAL :: l_she
     LOGICAL :: l_rmat
     LOGICAL :: l_nabla
     LOGICAL :: l_socodi
     LOGICAL :: l_pauli
     LOGICAL :: l_pauliat
     LOGICAL :: l_potmat
     LOGICAL :: l_projgen
     LOGICAL :: l_plot_symm
     LOGICAL :: l_socmmn0
     LOGICAL :: l_bzsym
     LOGICAL :: l_hopping
     LOGICAL :: l_kptsreduc
     LOGICAL :: l_prepwan90
     LOGICAL :: l_plot_umdat
     LOGICAL :: l_wann_plot
     LOGICAL :: l_bynumber
     LOGICAL :: l_stopopt
     LOGICAL :: l_matrixmmn
     LOGICAL :: l_matrixamn
     LOGICAL :: l_projmethod
     LOGICAL :: l_wannierize
     LOGICAL :: l_plotw90
     LOGICAL :: l_byindex
     LOGICAL :: l_byenergy
     LOGICAL :: l_proj_plot
     LOGICAL :: l_bestproj
     LOGICAL :: l_ikptstart
     LOGICAL :: l_lapw
     LOGICAL :: l_plot_lapw
     LOGICAL :: l_fermi
     LOGICAL :: l_dipole
     LOGICAL :: l_dipole2
     LOGICAL :: l_dipole3
     LOGICAL :: l_mmn0
     LOGICAL :: l_mmn0at
     LOGICAL :: l_manyfiles
     LOGICAL :: l_collectmanyfiles
     LOGICAL :: l_ldauwan
     LOGICAL :: l_lapw_kpts
     LOGICAL :: l_lapw_gfleur
     LOGICAL :: l_kpointgen
     LOGICAL :: l_w90kpointgen
     LOGICAL :: l_finishnocoplot
     LOGICAL :: l_finishgwf
     LOGICAL :: l_skipkov
     LOGICAL :: l_matrixuHu
     LOGICAL :: l_matrixuHu_dmi
     INTEGER :: ikptstart
     INTEGER :: band_min(1:2)
     INTEGER :: band_max(1:2)
     INTEGER :: gfthick
     INTEGER :: gfcut
     INTEGER :: unigrid(6)
     INTEGER :: mhp(3)
     !---> gwf
     LOGICAL :: l_ms
     LOGICAL :: l_sgwf
     LOGICAL :: l_socgwf
     LOGICAL :: l_gwf
     LOGICAL :: l_bs_comf
     LOGICAL :: l_exist
     LOGICAL :: l_opened
     LOGICAL :: l_cleverskip
     LOGICAL :: l_dim(3)
     REAL    :: scale_param
     REAL    :: aux_latt_const
     REAL    :: hdwf_t1
     REAL    :: hdwf_t2
     INTEGER :: nparampts
     CHARACTER(len=20) :: fn_eig
     CHARACTER(len=20) :: param_file
     REAL,ALLOCATABLE :: param_vec(:,:)
     REAL,ALLOCATABLE :: param_alpha(:,:)
     CHARACTER(LEN=20), ALLOCATABLE :: jobList(:)
     !---> gwf

  END TYPE t_wann

 
END MODULE m_types_setup
