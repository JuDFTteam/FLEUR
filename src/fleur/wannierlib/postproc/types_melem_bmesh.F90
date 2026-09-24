!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Neighbour / b-shell description of the coarse k-mesh, in a form that carries no
!>  Wannier90 types.
!>
!>  The position (Berry connection) matrix element is the one operator that is not a
!>  single-k Bloch matrix: it is built from the neighbour overlaps M^(0,b)(k) and needs
!>  the finite-difference b-vectors and their shell weights,
!>
!>      A^(W)_alpha(k) = i sum_b w_b b_alpha ( M^(W,b)(k) - delta ) .
!>
!>  Those weights are produced by the Wannier90 kmesh setup. Taking them from there
!>  directly would make the whole matrixelements layer depend on the Wannier90 library,
!>  so the wannierization side fills this plain bundle instead and passes
!>  it in. Any other source of a b-mesh can be plugged in the same way.
MODULE m_types_melem_bmesh
   USE m_judft, ONLY: juDFT_error
   IMPLICIT NONE
   PRIVATE

   TYPE :: t_melem_bmesh
      !> Topology: known as soon as the mesh is, which is before anything is wannierised.
      INTEGER :: nntot = 0                   !< number of b-vectors per k-point
      INTEGER, ALLOCATABLE :: nnlist(:, :)   !< (nk, nntot) global k index of neighbour b of k
      INTEGER, ALLOCATABLE :: gkpb(:, :, :)  !< (3, nk, nntot) G that brings the neighbour back
      !> (3, nntot) the distinct b vectors, in internal coordinates. Deduplicated by value,
      !> so its second index is a storage slot and not a neighbour slot: the two orders
      !> coincide only by accident. Look an entry up by value, the way melem_mmkb_sph does,
      !> or ask shell_vector for the b of a neighbour slot.
      REAL,    ALLOCATABLE :: kdiff(:, :)
      !> Shell weights: these the wannierisation produces, so they arrive later. Both stay
      !> unallocated for as long as the neighbour overlaps are still being built, which is
      !> where shell_vector is the only way to a b vector.
      REAL,    ALLOCATABLE :: wb(:)          !< (nntot) b-shell weights
      !> (3, nntot, nk) the same b vectors CARTESIAN, which is what the weights above are
      !> paired with. Neither the coordinates nor the index order match kdiff or gkpb.
      REAL,    ALLOCATABLE :: bk(:, :, :)
      !> Reference Wannier centres, used only by the R=0 Berry-centre calibration check.
      !> Left unallocated when no reference is available -> the check is skipped.
      REAL,    ALLOCATABLE :: centres(:, :)  !< (3, num_wann)
   CONTAINS
      PROCEDURE :: free           => melem_bmesh_free
      PROCEDURE :: set_neighbours => melem_bmesh_set_neighbours
      PROCEDURE :: shell_vector   => melem_bmesh_shell_vector
      PROCEDURE :: pair_diffs     => melem_bmesh_pair_diffs
   END TYPE t_melem_bmesh

   PUBLIC :: t_melem_bmesh

CONTAINS

   SUBROUTINE melem_bmesh_free(this)
      CLASS(t_melem_bmesh), INTENT(INOUT) :: this

      IF (ALLOCATED(this%nnlist)) DEALLOCATE (this%nnlist)
      IF (ALLOCATED(this%gkpb)) DEALLOCATE (this%gkpb)
      IF (ALLOCATED(this%kdiff)) DEALLOCATE (this%kdiff)
      IF (ALLOCATED(this%wb)) DEALLOCATE (this%wb)
      IF (ALLOCATED(this%bk)) DEALLOCATE (this%bk)
      IF (ALLOCATED(this%centres)) DEALLOCATE (this%centres)
      this%nntot = 0
   END SUBROUTINE melem_bmesh_free

   !> The neighbour topology, from whoever knows the mesh. Separate from the weights because
   !> it is available earlier: the overlaps between neighbours are what the wannierisation is
   !> given, so their topology cannot wait for its output.
   SUBROUTINE melem_bmesh_set_neighbours(this, nntot, nnlist, gkpb, kdiff)
      CLASS(t_melem_bmesh), INTENT(INOUT) :: this
      INTEGER, INTENT(IN) :: nntot
      INTEGER, INTENT(IN) :: nnlist(:, :)   !> (nk, nntot)
      INTEGER, INTENT(IN) :: gkpb(:, :, :)  !> (3, nk, nntot)
      REAL,    INTENT(IN) :: kdiff(:, :)    !> (3, nntot)

      this%nntot  = nntot
      this%nnlist = nnlist
      this%gkpb   = gkpb
      this%kdiff  = kdiff
   END SUBROUTINE melem_bmesh_set_neighbours

   !> The b vector of one neighbour slot, in the internal coordinates kdiff is written in:
   !>
   !>     b(k, nn) = bkf(nnlist(k, nn)) + gkpb(:, k, nn) - bkf(k)
   !>
   !> Needs the topology and nothing else, so it answers from the moment the mesh is known
   !> and long before the shell weights and their cartesian bk exist. Takes the mesh points
   !> as a plain array to keep this bundle free of any type.
   PURE FUNCTION melem_bmesh_shell_vector(this, bkf, k, nn) RESULT(b)
      CLASS(t_melem_bmesh), INTENT(IN) :: this
      REAL,    INTENT(IN) :: bkf(:, :)   !> (3, nk) the mesh points
      INTEGER, INTENT(IN) :: k, nn
      REAL :: b(3)

      b = bkf(:, this%nnlist(k, nn)) + this%gkpb(:, k, nn) - bkf(:, k)
   END FUNCTION melem_bmesh_shell_vector

   !> The distinct b2 - b1 vectors, which are the ones the muffin-tin half of a pair overlap
   !> needs a radial table for. Deduplicated by value and in internal coordinates, so the
   !> result goes straight into melem_ujugaunt and is matched by melem_mmkb_sph the way the
   !> neighbour table is.
   !>
   !> Swept over every k rather than assumed uniform: on a uniform mesh every k gives the same
   !> set and the sweep is redundant, and where it is not the table has to cover all of them
   !> or a pair overlap stops the run.
   !>
   !> Lives on the mesh because that is all it needs. uHu and uIu each carried an identical
   !> private copy, and which difference vectors exist is a property of the topology rather
   !> than of either consumer.
   SUBROUTINE melem_bmesh_pair_diffs(this, bkf, kdiff_pair, npair)
      CLASS(t_melem_bmesh), INTENT(IN) :: this
      REAL, INTENT(IN) :: bkf(:, :)                    !> (3, nk) the mesh points
      REAL, ALLOCATABLE, INTENT(OUT) :: kdiff_pair(:, :)
      INTEGER, INTENT(OUT) :: npair

      REAL :: d(3)
      INTEGER :: k, b1, b2, i
      LOGICAL :: seen

      ALLOCATE (kdiff_pair(3, this%nntot**2))
      kdiff_pair = 0.0
      npair = 0

      DO k = 1, SIZE(this%nnlist, 1)
         DO b1 = 1, this%nntot
            DO b2 = 1, this%nntot
               d = this%shell_vector(bkf, k, b2) - this%shell_vector(bkf, k, b1)
               seen = .FALSE.
               DO i = 1, npair
                  IF (ALL(ABS(kdiff_pair(:, i) - d) <= 1.0e-4)) THEN
                     seen = .TRUE.
                     EXIT
                  END IF
               END DO
               IF (seen) CYCLE
               IF (npair == SIZE(kdiff_pair, 2)) CALL juDFT_error( &
                  'wannierlib: more distinct b2-b1 vectors than pairs of neighbours', &
                  calledby='melem_bmesh_pair_diffs')
               npair = npair + 1
               kdiff_pair(:, npair) = d
            END DO
         END DO
      END DO
   END SUBROUTINE melem_bmesh_pair_diffs

END MODULE m_types_melem_bmesh
