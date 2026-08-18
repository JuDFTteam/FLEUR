!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_melem_optable
   !> Two tables with different owners. Keeping them apart is the point.
   !>
   !> MELEM_OPERATORS is the catalogue of what can be COMPUTED: Bloch matrices O(k) built
   !> from the eigenstates through the factory. Three entries, and they know nothing about
   !> Wannier functions. This is the table that would move with the operator layer if it
   !> ever became a service of its own.
   !>
   !> Its job today is to be what the exposure rows POINT AT. That is what turns "which
   !> coarse matrix does this name need" from a list to be remembered in four places into a
   !> lookup. Columns are added when something reads them, not in anticipation.
   !>
   !> WANNIERLIB_* say how the wannierization EXPOSES things -- which names its two input
   !> blocks accept, which output files each writes, and which catalogue entry (if any) it
   !> needs built for it. Most of what it exposes is not an operator at all: the interpolated
   !> Hamiltonian comes from the eigenvalues and the gauge, the position and the velocity
   !> from the neighbour overlaps, the currents from the algebra of the two, the eigenstates
   !> from a diagonalisation. None of those is a contraction over the states, so none of them
   !> belongs in the catalogue.
   !>
   !> ADDING AN OPERATOR: a row in MELEM_OPERATORS, plus a row in whichever exposure table
   !> should offer it. Then the schema enumeration, and a branch only if it is not served by
   !> the generic driver.

   IMPLICIT NONE
   PRIVATE

   ! ---------------------------------------------------------------- what can be computed
   TYPE, PUBLIC :: t_melem_operator
      CHARACTER(LEN=20) :: name = ''
      !> Cartesian components the coarse pass produces.
      INTEGER :: ncomp = 0
   END TYPE t_melem_operator

   TYPE(t_melem_operator), PARAMETER, PUBLIC :: MELEM_OPERATORS(3) = [ &
      t_melem_operator('spin',       3), &
      t_melem_operator('orbital',    3), &
      t_melem_operator('spin_orbit', 4)]

   ! ------------------------------------------------- how the wannierization exposes things
   TYPE, PUBLIC :: t_melem_exposed
      CHARACTER(LEN=20) :: name = ''
      !> Which catalogue entry has to be built for it; '' when nothing does.
      CHARACTER(LEN=20) :: operator = ''
      !> Components the generic interpolation driver writes; 0 = it has a driver of its own.
      INTEGER :: ncomp = 0
      !> Output basenames, without .dat and without the domain or channel suffix.
      CHARACTER(LEN=24) :: out1 = ''
      CHARACTER(LEN=24) :: out2 = ''
      !> Honours the total="T" attribute, i.e. it has a site-summed projection to choose.
      !> Only the two site-resolved ones do; soc and the currents have nothing to sum over.
      LOGICAL :: honours_total = .FALSE.
   END TYPE t_melem_exposed

   !> <interpolation>/<operator name="..">
   TYPE(t_melem_exposed), PARAMETER, PUBLIC :: WANNIERLIB_INTERP(6) = [ &
      t_melem_exposed('hamiltonian',    '',           0, 'bands_wann_interpol',     &
                      'bands_wann_interpol_ev', .FALSE.), &
      t_melem_exposed('spin',           'spin',       3, 'bands_wann_spin',         '', .TRUE.), &
      t_melem_exposed('orbital',        'orbital',    3, 'bands_wann_orbmom',       '', .TRUE.), &
      t_melem_exposed('soc',            'spin_orbit', 1, 'bands_wann_soc',          '', .FALSE.), &
      t_melem_exposed('velocity',       '',           0, 'bands_wann_velocity',     &
                      'bands_wann_berrycurv', .FALSE.), &
      t_melem_exposed('eigenstates',    '',           0, 'bands_wann_eigenstates',  '', .FALSE.)]

   !> <operators_r>/<operator name="..">. No output column: these write legacy names that
   !> carry a channel or a spin-block index, so the writer builds them itself.
   TYPE(t_melem_exposed), PARAMETER, PUBLIC :: WANNIERLIB_OPR(9) = [ &
      t_melem_exposed('hamiltonian', '',           0, '', '', .FALSE.), &
      t_melem_exposed('position',    '',           0, '', '', .FALSE.), &
      !> Second Berry connection, in the postw90 convention: Eq. (44) WYSV06 on the
      !> diagonal as well (no log branch) and A <- (A + A^dagger)/2 at each k. It is a
      !> separate name rather than a flag on `position` because both are needed at once:
      !> `position` writes the _r.dat that orbitrans reads and that carries the Wannier
      !> centres, while berry.F90 refuses that form for the orbital magnetisation
      !> ("transl_inv=T disabled for morb"). Writes WF*_rpw.dat.
      t_melem_exposed('position_pw90', '',         0, '', '', .FALSE.), &
      t_melem_exposed('bmn',         '',           0, '', '', .FALSE.), &
      t_melem_exposed('fmn',         '',           0, '', '', .FALSE.), &
      t_melem_exposed('cmn',         '',           0, '', '', .FALSE.), &
      t_melem_exposed('spin',        'spin',       3, '', '', .FALSE.), &
      t_melem_exposed('orbital',     'orbital',    3, '', '', .FALSE.), &
      t_melem_exposed('spin_orbit',  'spin_orbit', 4, '', '', .FALSE.)]

   PUBLIC :: melem_exposed_find, melem_operator_find, melem_exposed_names

CONTAINS

   !> Where a name sits in an exposure table, or 0 if that block does not accept it.
   PURE INTEGER FUNCTION melem_exposed_find(name, table) RESULT(idx)
      CHARACTER(LEN=*), INTENT(IN) :: name
      TYPE(t_melem_exposed), INTENT(IN) :: table(:)
      INTEGER :: i
      idx = 0
      DO i = 1, SIZE(table)
         IF (TRIM(name) == TRIM(table(i)%name)) THEN
            idx = i
            RETURN
         END IF
      END DO
   END FUNCTION melem_exposed_find

   !> Where an operator sits in the catalogue, or 0 for a name that is not one.
   PURE INTEGER FUNCTION melem_operator_find(name) RESULT(idx)
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER :: i
      idx = 0
      IF (LEN_TRIM(name) == 0) RETURN
      DO i = 1, SIZE(MELEM_OPERATORS)
         IF (TRIM(name) == TRIM(MELEM_OPERATORS(i)%name)) THEN
            idx = i
            RETURN
         END IF
      END DO
   END FUNCTION melem_operator_find

   !> The accepted names on one line, so an error can say what to write instead.
   PURE FUNCTION melem_exposed_names(table) RESULT(txt)
      TYPE(t_melem_exposed), INTENT(IN) :: table(:)
      CHARACTER(LEN=200) :: txt
      INTEGER :: i
      txt = 'accepted:'
      DO i = 1, SIZE(table)
         txt = TRIM(txt)//' '//TRIM(table(i)%name)
      END DO
   END FUNCTION melem_exposed_names

END MODULE m_types_melem_optable
