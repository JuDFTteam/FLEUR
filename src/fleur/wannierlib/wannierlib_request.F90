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
                        this%op_r_name, this%ops%name, this%ops%total, &
                        op_r_needs, op_needs, this%l_ws_distance)
      !> On a film, refuse by name rather than as a whole. The pair overlap has its vacuum
      !> half now, so the wannierization itself works and with it everything built from mmn:
      !> hamiltonian, eigenstates, position, velocity, bmn. Three quantities still stop at
      !> the interstitial and would come out short without saying so.
      IF (input%film) THEN
         !> wannierlib_uiu and wannierlib_uhu ask melem_overlap_states for the pair overlap
         !> without the vacuum expansions, and C also needs a momentum that has no vacuum
         !> part at all.
         IF (request%has_op_r('fmn') .OR. request%has_op_r('cmn')) CALL juDFT_error( &
            'wannierlib: fmn and cmn have no vacuum contribution and cannot be asked for on a film', &
            hint='everything built from the neighbour overlaps does work on a film; these two '// &
                 'need the vacuum halves of the pair overlap and of the momentum', &
            calledby='wannierlib_main')
         !> The spin operator sums the muffin tins and the interstitial. In a film there is a
         !> third region carrying spin density and this layer does not reach it. The orbital
         !> moment and the spin-orbit operator are muffin-tin quantities by construction, not
         !> by omission, so they are unaffected.
         IF (request%has_op('spin') .OR. request%has_op_r('spin')) CALL juDFT_error( &
            'wannierlib: the spin operator has no vacuum contribution and cannot be asked for on a film', &
            hint='it sums the muffin tins and the interstitial only; orbital and spin_orbit are '// &
                 'muffin-tin quantities and do work on a film', &
            calledby='wannierlib_main')
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
            hint="ask for it in <operators_r>: there both channels are combined into the "// &
                 "2N rspauli.1, which is the only form this operator has here", &
            calledby="wannierlib_main")
         IF (request%needs_op('spin_orbit', interp_only=.TRUE.)) CALL judft_error( &
            "the spin-orbit operator cannot be interpolated when the two spin channels are "// &
            "wannierised separately", &
            hint="ask for it in <operators_r>: there both channels are combined into the "// &
                 "2N rssocmat.1, which is the only form this operator has here", &
            calledby="wannierlib_main")
      END IF

   END SUBROUTINE wannierlib_build_request

END MODULE m_wannierlib_request
