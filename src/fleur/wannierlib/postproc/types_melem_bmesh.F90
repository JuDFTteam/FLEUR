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
!>  Those weights are produced by the Wannier90 kmesh setup and used to live in
!>  m_wannierlib_w90_adapter's saved lib_common_type. Reaching into that from here would
!>  make the whole matrixelements layer depend on the Wannier90 library, so the
!>  wannierization side fills this plain bundle instead (wannierlib_get_bmesh) and passes
!>  it in. Any other source of a b-mesh can be plugged in the same way.
MODULE m_types_melem_bmesh
   IMPLICIT NONE
   PRIVATE

   TYPE :: t_melem_bmesh
      !> Topology: known as soon as the mesh is, which is before anything is wannierised.
      INTEGER :: nntot = 0                   !< number of b-vectors per k-point
      INTEGER, ALLOCATABLE :: nnlist(:, :)   !< (nk, nntot) global k index of neighbour b of k
      INTEGER, ALLOCATABLE :: gkpb(:, :, :)  !< (3, nk, nntot) G that brings the neighbour back
      REAL,    ALLOCATABLE :: kdiff(:, :)    !< (3, nntot) the distinct b vectors
      !> Shell weights: these the wannierisation produces, so they arrive later.
      REAL,    ALLOCATABLE :: wb(:)          !< (nntot) b-shell weights
      REAL,    ALLOCATABLE :: bk(:, :, :)    !< (3, nntot, nk) cartesian b vectors
      !> Reference Wannier centres, used only by the R=0 Berry-centre calibration check.
      !> Left unallocated when no reference is available -> the check is skipped.
      REAL,    ALLOCATABLE :: centres(:, :)  !< (3, num_wann)
   CONTAINS
      PROCEDURE :: free           => melem_bmesh_free
      PROCEDURE :: set_neighbours => melem_bmesh_set_neighbours
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

END MODULE m_types_melem_bmesh
