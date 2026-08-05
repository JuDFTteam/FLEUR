!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_melem_domains
   !> Where an interpolation is to be evaluated: which output domains were asked for and
   !> the k-set of each one.
   !>
   !> The three domains are not interchangeable. A plane and a grid carry their k-points
   !> with them, so an empty k-set there means nothing was resolved. A path may instead
   !> name a file, and its default name is the one the interpolation reads anyway, so a
   !> path with neither k-set nor a name of its own is a request to use what is already
   !> on disk rather than an incomplete one.

   USE m_judft
   USE m_types_kpts

   IMPLICIT NONE
   PRIVATE

   TYPE t_melem_domains
      LOGICAL :: l_path  = .FALSE.   !> 1D path
      LOGICAL :: l_plane = .FALSE.   !> 2D plane
      LOGICAL :: l_grid  = .FALSE.   !> 3D uniform mesh
      CHARACTER(LEN=64) :: path_file = 'kpts_interpol'  !> k-list file of the path domain
      TYPE(t_kpts) :: path_kset      !> empty unless the path was given as a named list
      TYPE(t_kpts) :: plane_kset
      TYPE(t_kpts) :: grid_kset
   CONTAINS
      PROCEDURE :: init => melem_domains_init
   END TYPE t_melem_domains

   PUBLIC :: t_melem_domains

CONTAINS

   SUBROUTINE melem_domains_init(this, l_path, l_plane, l_grid, path_file, &
                                 path_kset, plane_kset, grid_kset)
      CLASS(t_melem_domains), INTENT(OUT) :: this
      LOGICAL,                INTENT(IN)  :: l_path, l_plane, l_grid
      CHARACTER(LEN=*),       INTENT(IN)  :: path_file
      TYPE(t_kpts),           INTENT(IN)  :: path_kset, plane_kset, grid_kset

      this%l_path  = l_path
      this%l_plane = l_plane
      this%l_grid  = l_grid
      this%path_file  = path_file
      this%path_kset  = path_kset
      this%plane_kset = plane_kset
      this%grid_kset  = grid_kset

      !> An empty name would send the path domain to read the current directory itself.
      IF (l_path .AND. LEN_TRIM(path_file) == 0) &
         CALL judft_error("t_melem_domains: the path domain has neither k-set nor file name", &
                          calledby="melem_domains_init")
   END SUBROUTINE melem_domains_init

END MODULE m_types_melem_domains
