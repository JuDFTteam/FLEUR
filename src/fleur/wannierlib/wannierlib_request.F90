!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grunberg Institut, Forschungszentrum Julich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> What the input asked the matrix-element layer for, resolved and refused by name.
!>
!> The two operator lists are spelled the way the wannierization spells them, which is not
!> a fact about the layer being asked, so they are resolved against the exposure tables
!> here and what leaves is the catalogue entry each name needs built. A name no table
!> carries stops the run with the accepted names, rather than leaving an operator silently
!> absent from the output.
!>
!> The refusals are by name and not as a whole, so that a request which is partly possible
!> says which part is not and why.
MODULE m_wannierlib_request
   USE m_juDFT
   USE m_types_input
   USE m_types_wannierlib, ONLY: t_wannierlib_wannierize
   USE m_types_melem_request, ONLY: t_melem_request
   USE m_types_melem_optable, ONLY: WANNIERLIB_INTERP, WANNIERLIB_OPR, &
                                    melem_exposed_find, melem_exposed_names
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wannierlib_build_request
CONTAINS

   SUBROUTINE wannierlib_build_request(this, input, l_spinors, request)
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
      TYPE(t_input), INTENT(IN) :: input
      LOGICAL, INTENT(IN) :: l_spinors       !> the states carry both spin components
      TYPE(t_melem_request), INTENT(OUT) :: request

      !> The catalogue entry each requested name needs built, one per name.
      CHARACTER(LEN=20), ALLOCATABLE :: op_needs(:), op_r_needs(:)
      INTEGER :: iop, krow

      !> The two lists are resolved against the exposure tables HERE, because those tables
      !> say how the wannierization spells things and that is not a fact about the layer
      !> being asked. What it receives is the catalogue entry each name needs built. A name
      !> no table carries stops the run with the accepted names, rather than leaving an
      !> operator silently absent from the output.
      ALLOCATE (op_needs(SIZE(this%ops)), op_r_needs(SIZE(this%op_r_name)))
      DO iop = 1, SIZE(this%ops)
         krow = melem_exposed_find(this%ops(iop)%name, WANNIERLIB_INTERP)
         IF (krow == 0) CALL juDFT_error('wannierlib: "'//TRIM(this%ops(iop)%name)// &
            '" is not an operator that can be interpolated', &
            hint=melem_exposed_names(WANNIERLIB_INTERP), calledby='wannierlib_main')
         op_needs(iop) = WANNIERLIB_INTERP(krow)%operator
      END DO
      DO iop = 1, SIZE(this%op_r_name)
         krow = melem_exposed_find(this%op_r_name(iop), WANNIERLIB_OPR)
         IF (krow == 0) CALL juDFT_error('wannierlib: "'//TRIM(this%op_r_name(iop))// &
            '" has no real-space export', &
            hint=melem_exposed_names(WANNIERLIB_OPR), calledby='wannierlib_main')
         op_r_needs(iop) = WANNIERLIB_OPR(krow)%operator
      END DO

      CALL request%init(this%l_spin, this%l_orbmom, this%l_socop, this%l_operators_r, &
                        this%op_r_name, this%ops%name, &
                        op_r_needs, op_needs, this%l_ws_distance)
      !> On a film the three regions are all covered now: the pair overlap, the spin
      !> operator and the momentum each have their vacuum half, so nothing is refused by
      !> geometry any more. What each of the two needed, and why one without the other
      !> would have been worse than nothing:
      !>
      !>   F is the pair overlap contracted with the gauge, so the vacuum expansions were
      !>   all it ever needed. C needs those AND a momentum; and adding one without the
      !>   other would not have been a partial fix but a wrong answer, because
      !>   <u_{k+b1}|H_k|u_{k+b2}> = E.M - b1.P - (b1^2/2).M replaces H by the eigenvalue of
      !>   the bra and so holds only over the whole cell. An M reaching the vacuum next to a
      !>   P that did not would have been off by exactly b1.P_vac, with nothing to say so.
      IF (input%film) THEN
         !> The spin operator does reach the vacuum: t_matelements_spin adds the third
         !> region to each of its four blocks, one t_melem_vacabc per spin channel, in the
         !> same pairing the interstitial uses. The orbital moment and the spin-orbit
         !> operator never needed it, being muffin-tin quantities by construction.
      END IF


      !> Two channels are two eigenproblems, wannierised one after the other, so there is
      !> no single spin matrix over them to interpolate: within a channel sigma_z is +/-1 by
      !> orthonormality, and the transverse part lives in the cross-spin block, which needs
      !> BOTH gauges and therefore both wannierizations. It exists only as the combined 2N
      !> operator in real space.
      IF (input%jspins == 2 .AND. .NOT.l_spinors) THEN
         IF (request%needs_op('spin', interp_only=.TRUE.)) CALL judft_error( &
            "the spin operator cannot be interpolated when the two spin channels are "// &
            "wannierised separately", &
            hint="ask for it in <operatorsRs>: there both channels are combined into the "// &
                 "2N rspauli.1, which is the only form this operator has here", &
            calledby="wannierlib_main")
         IF (request%needs_op('spinOrbit', interp_only=.TRUE.)) CALL judft_error( &
            "the spin-orbit operator cannot be interpolated when the two spin channels are "// &
            "wannierised separately", &
            hint="ask for it in <operatorsRs>: there both channels are combined into the "// &
                 "2N rssocmat.1, which is the only form this operator has here", &
            calledby="wannierlib_main")
         IF (request%needs_op('exchangeTorque', interp_only=.TRUE.)) CALL judft_error( &
            "the exchange torque cannot be interpolated when the two spin channels "// &
            "are wannierised separately", &
            hint="ask for it in <operatorsRs>: collinearly this operator is the "// &
                 "up-down Bloch block alone, written to updown.perpmag", &
            calledby="wannierlib_main")
      END IF

   END SUBROUTINE wannierlib_build_request

END MODULE m_wannierlib_request
