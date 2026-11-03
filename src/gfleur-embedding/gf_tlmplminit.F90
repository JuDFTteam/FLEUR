!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_tlmplminit
      use m_juDFT
      use m_gf_types
      use m_local_Hamiltonian
      use m_constants, only: POTDEN_TYPE_POTTOT
      IMPLICIT NONE
      private
      public :: gf_tlmplminit
      CONTAINS
      SUBROUTINE gf_tlmplminit(ld, fmpi)
!-----------------------------------------------
!   Build the radial functions and the local Hamiltonian (both in
!   tlmplm, incl. LDA+U terms and Cholesky decomposition) of one layer
!   from its current potential by calling the standard FLEUR local_ham.
!   Replaces the old fleur_tlmplm/fleur_usetup wrappers and the manual
!   allocation of the v25-era tuu/tdd/ind arrays.
!           (last modified: 05-04-08) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      TYPE(t_gflayer),INTENT(INOUT) :: ld
      TYPE(t_mpi),INTENT(IN)        :: fmpi
      !>
      !<-- Locals
      TYPE(t_potden) :: vx
      LOGICAL        :: l_error
      !>

      IF (.NOT.ALLOCATED(ld%hub1data%mag_mom))                          &
     &     CALL ld%hub1data%init(ld%fi%atoms, ld%fi%input, ld%fi%hub1inp,&
     &          fmpi, 0.0, 0.0, l_error)

      !a zeroed exchange potential - gfleur runs no hybrid functionals
      CALL vx%init(ld%stars, ld%fi%atoms, ld%sphhar, ld%fi%vacuum,       &
     &             ld%fi%noco, ld%fi%input%jspins, POTDEN_TYPE_POTTOT)

      CALL local_ham(ld%sphhar, ld%fi%atoms, ld%fi%sym, ld%fi%noco,      &
     &               ld%nococonv, ld%enpara, fmpi, ld%vTot, vx,          &
     &               ld%cdn_new, ld%fi%input, ld%fi%hub1inp,             &
     &               ld%hub1data, ld%tlmplm, 0.0)

      END SUBROUTINE gf_tlmplminit
      END
