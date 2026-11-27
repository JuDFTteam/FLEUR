!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> The Monkhorst-Pack divisions of the wannierisation mesh, handed out only after
!> they have been checked against the k-point list they claim to describe.
MODULE m_wgauge_mesh
   USE m_juDFT
   USE m_constants, ONLY: oUnit
   USE m_types_kpts

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wgauge_mp_grid

CONTAINS

   !> WHY THIS IS NOT kpts%nkpt3 READ DIRECTLY. That field carries a warning of its own
   !> in types_kpts: "This variable is not supposed to be reliable. It is only meant as
   !> additional user information with the available data at IO time." It is filled from
   !> the @nx, @ny and @nz of whichever kPointList was selected, and a list that has none
   !> -- a path, or one written by hand -- leaves it at its default of zero.
   !>
   !> Six places need these divisions: Wannier90 takes them as mp_grid, the exported .win
   !> writes them, and the Wigner-Seitz construction builds the R mesh of every O(R) from
   !> them, including the minimum-distance replicas and the wig_vectors file. Handing them
   !> out from one place is what keeps the six from drifting, and is why this is a
   !> function rather than a check someone has to remember to call first.
   !>
   !> The test is exact, not a heuristic. Wannier90 is given num_kpts = kpts%nkptf and
   !> kpoints = kpts%bkf beside mp_grid, so the divisions have to multiply out to that
   !> count. Wannier90 makes the same test itself, in w90_readwrite_read_mp_grid, but only
   !> on the path that READS a .win; library mode sets the options directly and never
   !> reaches it, so without this there is no test anywhere.
   !>
   !> Every rank calls this with the same numbers -- nkpt3 and nkptf are both broadcast --
   !> so the verdict is the same everywhere and the abort is collective.
   FUNCTION wgauge_mp_grid(kpts) RESULT(mp_grid)
      TYPE(t_kpts), INTENT(IN) :: kpts
      INTEGER :: mp_grid(3)

      mp_grid = kpts%nkpt3

      !> Two failures, two messages. All-zero is "the selected list is not a mesh at all";
      !> a product that does not match is "it is a mesh, but not this one".
      IF (ANY(mp_grid < 1)) CALL juDFT_error( &
         'wannierlib: the selected k-point list carries no mesh divisions', &
         hint='wannierise on a kPointList that has nx, ny and nz -- a full-zone, '// &
              'Gamma-centred mesh. A path has no divisions, and the operators in R are '// &
              'built on the mesh those divisions describe.', &
         calledby='wgauge_mp_grid')

      IF (PRODUCT(mp_grid) /= kpts%nkptf) THEN
         WRITE (oUnit, '(/,a)')   'wannierlib: the mesh divisions are not those of the k-point list'
         WRITE (oUnit, '(a,3i6)') '   nx, ny, nz of the list   :', mp_grid
         WRITE (oUnit, '(a,i8)')  '   their product            :', PRODUCT(mp_grid)
         WRITE (oUnit, '(a,i8)')  '   k-points in the full zone:', kpts%nkptf
         CALL juDFT_error( &
            'wannierlib: nx*ny*nz is not the number of k-points in the full zone', &
            hint='the divisions and the list have to describe the same mesh; out.xml '// &
                 'carries both numbers', &
            calledby='wgauge_mp_grid')
      END IF
   END FUNCTION wgauge_mp_grid

END MODULE m_wgauge_mesh
