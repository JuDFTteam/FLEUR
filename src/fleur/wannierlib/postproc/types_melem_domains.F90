!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_melem_domains
   !> Where an interpolation is to be evaluated: one entry per <domain>, each a k-set and
   !> the suffix its output files carry.
   !>
   !> There is nothing else to distinguish them. Whether a k-set traces a line, covers a
   !> plane or fills the zone is a property of the list the user named, not something this
   !> code needs to know: every domain is interpolated the same way and differs only in
   !> where its output lands.
   !>
   !> IMPORTANT: the k-sets and suffixes exist ONLY on rank 0, because only rank 0 writes
   !> the k-set file and renames the outputs. n is broadcast and bounds the domain loop on
   !> every rank; the arrays are not, so they are unallocated off rank 0. Do not read or
   !> validate them there -- a check on nkpt in this init aborted every non-root rank.

   USE m_judft
   USE m_types_kpts

   IMPLICIT NONE
   PRIVATE

   TYPE t_melem_domains
      INTEGER :: n = 0                              !> number of domains (broadcast)
      TYPE(t_kpts),      ALLOCATABLE :: kset(:)     !> (n) rank 0 only
      CHARACTER(LEN=64), ALLOCATABLE :: suffix(:)   !> (n) rank 0 only
   CONTAINS
      PROCEDURE :: init => melem_domains_init
   END TYPE t_melem_domains

   PUBLIC :: t_melem_domains

CONTAINS

   !> kset and suffix are the input type's arrays, which only rank 0 has. They are copied
   !> when present and left alone otherwise, so this runs on every rank without asking who
   !> it is: n alone decides how many times the domain loop turns.
   SUBROUTINE melem_domains_init(this, n, kset, suffix)
      CLASS(t_melem_domains), INTENT(OUT) :: this
      INTEGER,                INTENT(IN)  :: n
      TYPE(t_kpts),      ALLOCATABLE, INTENT(IN) :: kset(:)
      CHARACTER(LEN=64), ALLOCATABLE, INTENT(IN) :: suffix(:)

      this%n = n
      IF (ALLOCATED(kset))   this%kset   = kset
      IF (ALLOCATED(suffix)) this%suffix = suffix
   END SUBROUTINE melem_domains_init

END MODULE m_types_melem_domains
