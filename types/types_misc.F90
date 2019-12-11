!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_misc

   IMPLICIT NONE

   !*************************************************************
   !     This module contains definitions for all kind of types
   !*************************************************************

   ! Type for the HF total energy
   TYPE t_energy_hf
      REAL :: valence
      REAL :: core
   END TYPE t_energy_hf

   TYPE t_results
      REAL, ALLOCATABLE    :: force(:,:,:)   !< Forces calculated on all atoms (for each spin)
      REAL, ALLOCATABLE    :: force_old(:,:) !< Forces on all atoms from last iteration
      REAL                 :: ef        !<Fermie energy
      REAL                 :: seigc     !<sum of the core eigenvalues
      REAL                 :: seigsc    !<weighted sum of the semi-core eigenvalues
      REAL                 :: seigv     !<weighted sum of the occupied valence eigenvalues
      REAL                 :: seigscv   !<sum of seigv and seigsc
      REAL                 :: ts        !<entropy contribution to the free energy
      REAL                 :: te_vcoul  !<charge density-coulomb potential integral
      REAL                 :: te_veff   !<charge density-effective potential integral
      REAL                 :: te_exc    !<charge density-ex-corr.energy density integral
      REAL                 :: e_ldau    !<total energy contribution of LDA+U
      REAL                 :: tote
      REAL                 :: last_distance
      REAL                 :: bandgap
      COMPLEX, ALLOCATABLE    :: unfolding_weights(:,:,:) !weights for unfolding a supercell bandstructure
      TYPE(t_energy_hf)    ::  te_hfex
      REAL                 ::  te_hfex_loc(2)
      REAL, ALLOCATABLE    :: w_iks(:,:,:)
      REAL, ALLOCATABLE    :: w_iksRDMFT(:,:,:)
      REAL, ALLOCATABLE    :: eig(:,:,:)
      INTEGER, ALLOCATABLE :: neig(:,:) ! neig(nkpts,jspins) number of calculated eigenvalues for each k point, spin

   CONTAINS
      PROCEDURE,PASS :: init => results_init
   END TYPE t_results

   TYPE t_zMat
      LOGICAL              :: l_real
      INTEGER              :: nbasfcn
      INTEGER              :: nbands
      REAL,    ALLOCATABLE :: z_r(:,:) ! z_r(nbasfcn,nbands)
      COMPLEX, ALLOCATABLE :: z_c(:,:) ! z_c(nbasfcn,nbands)

   CONTAINS
      PROCEDURE,PASS :: init => zMat_init
   END TYPE t_zMat

   TYPE t_hamOvlp
      LOGICAL              :: l_real
      INTEGER              :: matsize
      REAL,    ALLOCATABLE :: a_r(:), b_r(:)
      COMPLEX, ALLOCATABLE :: a_c(:), b_c(:)
   END TYPE t_hamOvlp

CONTAINS

   SUBROUTINE zMat_init(thisZMat,l_real,nbasfcn,nbands)

      IMPLICIT NONE

   CLASS(t_zMat),      INTENT(INOUT) :: thisZMat
      LOGICAL,            INTENT(IN)    :: l_real
      INTEGER,            INTENT(IN)    :: nbasfcn,nbands

      thisZMat%l_real = l_real
      thisZMat%nbasfcn = nbasfcn
      thisZMat%nbands = nbands

      IF (ALLOCATED(thisZMat%z_r)) DEALLOCATE(thisZMat%z_r)
      IF (ALLOCATED(thisZMat%z_c)) DEALLOCATE(thisZMat%z_c)
      IF (l_real) THEN
         ALLOCATE(thisZMat%z_r(nbasfcn,nbands))
         thisZMat%z_r = 0.0
      ELSE
         ALLOCATE(thisZMat%z_c(nbasfcn,nbands))
         thisZMat%z_c = CMPLX(0.0,0.0)
      END IF

   END SUBROUTINE zMat_init

   SUBROUTINE results_init(thisResults,input,atoms,kpts,noco)

     USE m_types_atoms
     USE m_types_input
     USE m_types_noco
     USE m_types_dimension
     USE m_types_kpts
     USE m_types_lapw
      IMPLICIT NONE

      CLASS(t_results),      INTENT(INOUT) :: thisResults

      TYPE(t_input),         INTENT(IN)    :: input
      TYPE(t_atoms),         INTENT(IN)    :: atoms
      TYPE(t_kpts),          INTENT(IN)    :: kpts
      TYPE(t_noco),          INTENT(IN)    :: noco

      INTEGER                              :: neigd2

      thisResults%seigc           = 0.0
      thisResults%seigsc          = 0.0
      thisResults%seigv           = 0.0
      thisResults%seigscv         = 0.0
      thisResults%e_ldau          = 0.0
      thisResults%ts              = 0.0

      thisResults%te_vcoul        = 0.0
      thisResults%te_veff         = 0.0
      thisResults%te_exc          = 0.0
      thisResults%te_hfex%valence = 0.0
      thisResults%te_hfex%core    = 0.0
      thisResults%te_hfex_loc     = 0.0

      thisResults%tote            = 0.0
      thisResults%last_distance   = -1.0
      thisResults%bandgap         = 0.0
      thisResults%ef              = 0.0

      neigd2 = MIN(input%neig,lapw_dim_nbasfcn)
!   neigd2 = input%neig
      IF (noco%l_soc.AND.(.NOT.noco%l_noco)) neigd2 = 2*neigd2

      ALLOCATE (thisResults%force(3,atoms%ntype,input%jspins));thisResults%force=0.0
      ALLOCATE (thisResults%force_old(3,atoms%ntype));thisResults%force_old=0.0
      ALLOCATE (thisResults%w_iks(neigd2,kpts%nkpt,input%jspins))
      ALLOCATE (thisResults%neig(kpts%nkpt,input%jspins))
      ALLOCATE (thisResults%eig(neigd2,kpts%nkpt,input%jspins))
      ALLOCATE (thisResults%unfolding_weights(neigd2,kpts%nkpt,input%jspins))

      thisResults%force = 0.0
      thisResults%force_old = 0.0
      thisResults%w_iks = 0.0
      thisResults%neig = 0
      thisResults%eig = 0.0

      IF(input%l_rdmft) THEN
         ALLOCATE (thisResults%w_iksRDMFT(neigd2,kpts%nkpt,input%jspins))
         thisResults%w_iksRDMFT = 0.0
      END IF

   END SUBROUTINE results_init

END MODULE m_types_misc
