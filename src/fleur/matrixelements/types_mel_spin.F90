!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Spin (Pauli) matrix element as a t_matrixelement:
!>
!>      S0_alpha,mn(k) = < psi_mk | sigma_alpha | psi_nk > ,  alpha = x,y,z  (ncomp=3)
!>
!>  Example/skeleton derived type proving the abstract interface: everything shared
!>  (eigenstate read, abc coefficients, gauge rotation, FT, interpolation, IO) is
!>  inherited from t_matrixelement; only the deferred per-k Bloch builder is bound,
!>  and it is a one-line delegation to the existing provider m_melem_spin
!>  (muffin-tin blocks from abc x radfun%integral + interstitial from the spinor zMat).
!>
!>  Spinor mode only: the on-site Pauli operator needs both spinor components in one
!>  eigenvector, so init errors out for the collinear jspins=2 case (there the spin
!>  channels are wannierised separately and sigma_z is trivial per channel) -- same
!>  behaviour as the current pipeline, which skips the collinear coarse spin slice.
!>
!>  Planned (documented only): a per-atom variant t_mel_spin_peratom with
!>  init_melem(ctx, 'spin_peratom', 3*atoms%nat, nsites=atoms%nat), delegating to
!>  melem_spin_peratom and RESHAPEing (nb,nb,3,nat) -> (nb,nb,3*nat).
MODULE m_types_mel_spin
   USE m_juDFT
   USE m_types_matrixelement
   USE m_types_lapw
   USE m_types_mat
   USE m_types_abc
   USE m_melem_spin, ONLY : melem_spin_bloch
   IMPLICIT NONE
   PRIVATE

   TYPE, EXTENDS(t_matrixelement) :: t_mel_spin
   CONTAINS
      PROCEDURE :: init       => mel_spin_init
      PROCEDURE :: calc_bloch => mel_spin_calc_bloch
   END TYPE t_mel_spin

   PUBLIC :: t_mel_spin

CONTAINS

   SUBROUTINE mel_spin_init(this, ctx)
      CLASS(t_mel_spin),   INTENT(INOUT) :: this
      TYPE(t_mel_context), INTENT(IN)    :: ctx

      IF (.NOT. ctx%l_spinors) CALL juDFT_error( &
         't_mel_spin needs spinor wavefunctions (noco or SOC); collinear spin is per-channel', &
         calledby='mel_spin_init')
      CALL this%init_melem(ctx, 'spin', 3)
   END SUBROUTINE mel_spin_init

   !> Deferred builder: delegate to the existing provider. abc carries both local
   !> spin components (spinor convention of the base k-loop); the sum-rule check is
   !> printed for the first few k-points, as in the current coarse pass.
   SUBROUTINE mel_spin_calc_bloch(this, ctx, ik, lapw, zMat, abc, o0_k)
      CLASS(t_mel_spin),   INTENT(INOUT) :: this
      TYPE(t_mel_context), INTENT(IN)    :: ctx
      INTEGER,             INTENT(IN)    :: ik
      TYPE(t_lapw),        INTENT(IN)    :: lapw
      TYPE(t_mat),         INTENT(IN)    :: zMat
      TYPE(t_abc),         INTENT(IN)    :: abc(:, :)     ! (ntype, 2)
      COMPLEX,             INTENT(INOUT) :: o0_k(:, :, :) ! (nb, nb, 3)

      CALL melem_spin_bloch(ctx%atoms, abc, ctx%radfun, ctx%nococonv, ctx%stars, &
                                 lapw, zMat, this%nb, ik, o0_k, ik <= 3)
   END SUBROUTINE mel_spin_calc_bloch

END MODULE m_types_mel_spin
