!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_make_sym
   USE m_judft
   IMPLICIT NONE
   PRIVATE
   PUBLIC make_sym
CONTAINS
   SUBROUTINE make_sym(sym, cell, atoms, noco, oneD, input, gfinp)
      !Generates missing symmetry info.
      !tau,mrot and nop have to be specified alread
      USE m_dwigner
      USE m_angles !Phase factors for spin-offdiagonal lda+u
      USE m_constants
      USE m_mapatom
      USE m_od_mapatom
      use m_ptsym
      USE m_types_sym
      USE m_types_cell
      USE m_types_atoms
      USE m_types_noco
      USE m_types_oneD
      use m_types_input
      USE m_types_gfinp
      use m_types_fleurinput_base, only: REAL_NOT_INITALIZED, CMPLX_NOT_INITALIZED
      TYPE(t_sym), INTENT(INOUT) :: sym
      TYPE(t_cell), INTENT(IN)   :: cell
      TYPE(t_atoms), INTENT(IN)  :: atoms
      TYPE(t_noco), INTENT(IN)   :: noco
      TYPE(t_oneD), INTENT(INOUT):: oneD
      TYPE(t_input), INTENT(IN)  :: input
      TYPE(t_gfinp), INTENT(IN)  :: gfinp

      integer :: nsymt
      integer, allocatable::nrot(:), locops(:, :)

      !Check for additional time-reversal symmetry
      IF (sym%invs .OR. noco%l_soc) THEN
         sym%nsym = sym%nop
      ELSE
         ! combine time reversal symmetry with the spatial symmetry opera
         ! thus the symmetry operations are doubled
         sym%nsym = 2*sym%nop
      END IF

      !Generated wigner symbols for LDA+U (includes DFT+HubbardI)
      IF (ALLOCATED(sym%d_wgn)) DEALLOCATE (sym%d_wgn)
      ALLOCATE (sym%d_wgn(-lmaxU_const:lmaxU_const, -lmaxU_const:lmaxU_const, lmaxU_const, sym%nop), &
               source=CMPLX_NOT_INITALIZED)
      IF (atoms%n_denmat + gfinp%n .GT. 0) THEN !replace with atoms%n_u+gfinp%n

         CALL d_wigner(sym%nop, sym%mrot, cell%bmat, lmaxU_const, sym%d_wgn, write=.TRUE.)
         !For spin-offdiagonal parts, we need additional phase factors
         IF (noco%l_mperp) THEN
            IF (ALLOCATED(sym%phase)) DEALLOCATE (sym%phase)
            ALLOCATE (sym%phase(sym%nop), source=REAL_NOT_INITALIZED)
            CALL angles(sym)
         ENDIF
      END IF

      !Atom specific symmetries

      allocate (locops(sym%nop, atoms%nat), nrot(atoms%nat))
      if(.not. allocated(sym%ntypsy)) allocate (sym%ntypsy(atoms%nat))
      call ptsym(atoms%ntype, atoms%nat, atoms%neq, atoms%taual, sym%nop, sym%mrot, sym%tau, atoms%lmax, &
                 nsymt, sym%ntypsy, nrot, locops)

      IF (.NOT. oneD%odd%d1) THEN
         CALL mapatom(sym, atoms, cell, input, noco,gfinp)
         allocate (oneD%ngopr1(atoms%nat))
         oneD%ngopr1 = sym%ngopr
      ELSE
         CALL juDFT_error("The oneD version is broken here. Compare call to mapatom with old version")
         CALL mapatom(sym, atoms, cell, input, noco,gfinp)
         !CALL od_mapatom(oneD,atoms,sym,cell)
      END IF

   END SUBROUTINE make_sym
END MODULE m_make_sym
