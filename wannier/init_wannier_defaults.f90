!--------------------------------------------------------------------------------
! Copyright (c) 2017 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_init_wannier_defaults

CONTAINS

SUBROUTINE initWannierDefaults(wann)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!
   !!!  This subroutine sets most of the attributes of the t_wann
   !!!  type to standard values.
   !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   USE m_types_setup

   IMPLICIT NONE

   TYPE(t_wann), INTENT(INOUT) :: wann

   wann%wan90version = 2 ! Set the standard to Wannier90-1.2

   wann%oc_num_orbs = 0
!        integer,allocatable :: oc_orbs(:)

   wann%l_unformatted = .FALSE.
   wann%l_oc_f = .FALSE.
   wann%l_ndegen = .FALSE.
   wann%l_orbitalmom = .FALSE.
   wann%l_orbcomp = .FALSE.
   wann%l_orbcomprs = .FALSE.
   wann%l_denmat = .FALSE.
   wann%l_perturbrs = .FALSE.
   wann%l_perturb = .FALSE.
   wann%l_nedrho = .FALSE.
   wann%l_anglmomrs = .FALSE.
   wann%l_anglmom = .FALSE.
   wann%l_spindisp = .FALSE.
   wann%l_spindisprs = .FALSE.
   wann%l_socspicom = .FALSE.
   wann%l_socspicomrs = .FALSE.
   wann%l_offdiposoprs = .FALSE.
   wann%l_offdiposop = .FALSE.
   wann%l_torque = .FALSE.
   wann%l_torquers = .FALSE.
   wann%l_atomlist = .FALSE.

   wann%atomlist_num = 0 ! has to be initialize to atoms%nat or something smaller at some point
!        integer,allocatable :: atomlist(:)

   wann%l_berry = .FALSE.
   wann%l_perpmagrs = .FALSE.
   wann%l_perpmag = .FALSE.
   wann%l_perpmagat = .FALSE.
   wann%l_perpmagatrs = .FALSE.
   wann%l_socmatrs = .FALSE.
   wann%l_socmat = .FALSE.
   wann%l_soctomom = .FALSE.
   wann%l_kptsreduc2 = .FALSE.
   wann%l_nablapaulirs = .FALSE.
   wann%l_nablars = .FALSE.
   wann%l_surfcurr = .FALSE.
   wann%l_updown = .FALSE.
   wann%l_ahe = .FALSE.
   wann%l_she = .FALSE.
   wann%l_rmat = .FALSE.
   wann%l_nabla = .FALSE.
   wann%l_socodi = .FALSE.
   wann%l_pauli = .FALSE.
   wann%l_pauliat = .FALSE.
   wann%l_potmat = .FALSE.
   wann%l_projgen = .FALSE.
   wann%l_plot_symm = .FALSE.
   wann%l_socmmn0 = .FALSE.
   wann%l_bzsym = .FALSE.
   wann%l_hopping = .FALSE.
   wann%l_kptsreduc = .FALSE.
   wann%l_prepwan90 = .FALSE.
   wann%l_plot_umdat = .FALSE.
   wann%l_wann_plot = .FALSE.
   wann%l_bynumber = .FALSE.
   wann%l_stopopt = .FALSE.
   wann%l_matrixmmn = .FALSE.
   wann%l_matrixamn = .FALSE.
   wann%l_projmethod = .FALSE.
   wann%l_wannierize = .FALSE.
   wann%l_plotw90 = .FALSE.
   wann%l_byindex = .FALSE.
   wann%l_byenergy = .FALSE.
   wann%l_proj_plot = .FALSE.
   wann%l_bestproj = .FALSE.
   wann%l_ikptstart = .FALSE.
   wann%l_lapw = .FALSE.
   wann%l_plot_lapw = .FALSE.
   wann%l_fermi = .FALSE.
   wann%l_dipole = .FALSE.
   wann%l_dipole2 = .FALSE.
   wann%l_dipole3 = .FALSE.
   wann%l_mmn0 = .FALSE.
   wann%l_mmn0at = .FALSE.
   wann%l_manyfiles = .FALSE.
   wann%l_collectmanyfiles = .FALSE.
   wann%l_ldauwan = .FALSE.
   wann%l_lapw_kpts = .FALSE.
   wann%l_lapw_gfleur = .FALSE.
   wann%l_kpointgen = .FALSE.
   wann%l_w90kpointgen = .FALSE.
   wann%l_finishnocoplot = .FALSE.
   wann%l_finishgwf = .FALSE.
   wann%l_skipkov = .FALSE.
   wann%l_matrixuHu = .FALSE.
   wann%l_matrixuHu_dmi = .FALSE.

   wann%ikptstart = 1
   wann%band_min(1:2) = -1
   wann%band_max(1:2) = -1
   wann%gfthick = 0
   wann%gfcut = 0
   wann%unigrid(6) = 0
   wann%mhp(3) = 0

!---> gwf
   wann%l_ms = .FALSE.
   wann%l_sgwf = .FALSE.
   wann%l_socgwf = .FALSE.
   wann%l_gwf = .FALSE.
   wann%l_bs_comf = .FALSE.
   wann%l_exist = .FALSE.
   wann%l_opened = .FALSE.
   wann%l_cleverskip = .FALSE.
   wann%l_dim(3) = .FALSE.

   wann%scale_param = 0.0
   wann%aux_latt_const = 0.0
   wann%hdwf_t1 = 0.0
   wann%hdwf_t2 = 0.0
   wann%nparampts = 0
   wann%fn_eig = ''
   wann%param_file = ''

   wann%scale_param = 1.0
   wann%aux_latt_const = 8.0!5.5!5.45886450 !5.98136400 !8.0725882513951497 !5.4170 !1.0
   wann%param_file='qpts'
   wann%l_dim=.false.

END SUBROUTINE initWannierDefaults

END MODULE m_init_wannier_defaults
