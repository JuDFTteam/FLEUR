!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_fleur_psqpw
      IMPLICIT NONE
!----------------------------------------------------------------
!Interface to the pseudo-charge generation of FLEUR.
!In the port to current FLEUR the old ~45-argument psqpw call becomes
!the standard psqpw on a t_potden; the convergence parameters ncv are
!set by convn during the layer initialization. The no_core option
!(used by the preconditioner) is realized by a local atoms copy with
!zatom=0 - no FLEUR-side flag needed.
!----------------------------------------------------------------
      PRIVATE
      PUBLIC fleur_psqpw
      CONTAINS
      !<-- S: fleur_psqpw(fmpi,ld,den,ispin,psq,no_core)

      SUBROUTINE fleur_psqpw(fmpi,ld,den,ispin,psq,no_core)
!-----------------------------------------------
!calls psqpw from FLEUR
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      USE m_gf_types
      USE m_psqpw
            IMPLICIT NONE
      !<-- Arguments
      TYPE(t_mpi),INTENT(IN)      :: fmpi
      TYPE(t_gflayer),INTENT(IN)  :: ld
      TYPE(t_potden),INTENT(IN)   :: den
      INTEGER,INTENT(IN)          :: ispin
                                                !pseudo-charge
      COMPLEX,INTENT(OUT)         :: psq(:)
      logical,intent(in),optional :: no_core
      !>
      !<-- Locals
      TYPE(t_atoms) :: atoms_loc
      LOGICAL       :: l_no_core
      !>
      l_no_core = .FALSE.
      IF (PRESENT(no_core)) l_no_core = no_core

      IF (l_no_core) THEN
         atoms_loc = ld%fi%atoms
         atoms_loc%zatom = 0.0
         CALL psqpw(fmpi,atoms_loc,ld%sphhar,ld%stars,ld%fi%vacuum,      &
     &              ld%fi%cell,ld%fi%input,ld%fi%sym,den,ispin,          &
     &              .FALSE.,den%potdenType,psq)
      ELSE
         CALL psqpw(fmpi,ld%fi%atoms,ld%sphhar,ld%stars,ld%fi%vacuum,    &
     &              ld%fi%cell,ld%fi%input,ld%fi%sym,den,ispin,          &
     &              .FALSE.,den%potdenType,psq)
      ENDIF
      END SUBROUTINE

      !>
      END
