!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_local_coordination_core
   USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
   IMPLICIT NONE
   PRIVATE

   INTEGER, PARAMETER :: default_maximum_image_range = 32

   TYPE, PUBLIC :: t_periodic_neighbor
      ! atom_id identifies a physical atom in the supplied reference cell.
      ! translation is the lattice-vector shift applied to that atom.
      INTEGER :: atom_id = 0
      INTEGER :: translation(3) = 0
      REAL :: displacement(3) = 0.0
      REAL :: distance = 0.0
   END TYPE t_periodic_neighbor

   TYPE, PUBLIC :: t_coordination_diagnostics
      INTEGER :: search_range = 0
      INTEGER :: enumerated_images = 0
      LOGICAL :: has_next_neighbor = .FALSE.
      REAL :: next_neighbor_distance = HUGE(1.0)
      REAL :: shell_gap = HUGE(1.0)
      REAL :: next_to_last_ratio = HUGE(1.0)
      REAL :: minimum_unseen_distance_bound = HUGE(1.0)
   END TYPE t_coordination_diagnostics

   PUBLIC :: resolve_periodic_neighbours_cartesian

CONTAINS

   SUBROUTINE resolve_periodic_neighbours_cartesian(lattice, periodic, center, candidate_positions, candidate_ids, &
                                                    nneighbors, neighbors, diagnostics, success, message, &
                                                    maximum_image_range)
      ! Resolve distinct (physical atom, lattice translation) images by their
      ! Cartesian distance from center. lattice(:,i) is lattice vector i.
      !
      ! The image search is adaptive. For a periodic basis B and any integer
      ! translation outside [-r,r]^p,
      !
      !   |B n + delta| >= sigma_min(B)*(r+1) - max_candidate|delta|.
      !
      ! Enumeration stops only when this lower bound exceeds the (N+1)-th
      ! distance. Thus the returned N images are complete even for skew cells;
      ! no component-wise fractional-coordinate wrapping is used.
      REAL, INTENT(IN) :: lattice(3,3), center(3), candidate_positions(:,:)
      LOGICAL, INTENT(IN) :: periodic(3)
      INTEGER, INTENT(IN) :: candidate_ids(:), nneighbors
      TYPE(t_periodic_neighbor), ALLOCATABLE, INTENT(OUT) :: neighbors(:)
      TYPE(t_coordination_diagnostics), INTENT(OUT) :: diagnostics
      LOGICAL, INTENT(OUT) :: success
      CHARACTER(LEN=*), INTENT(OUT) :: message
      INTEGER, OPTIONAL, INTENT(IN) :: maximum_image_range

      TYPE(t_periodic_neighbor), ALLOCATABLE :: images(:)
      REAL :: sigma_minimum, maximum_offset, lower_bound, zero_tolerance, scale
      INTEGER :: image_range, range_limit, nimages, candidate
      LOGICAL :: basis_success, enumeration_success
      CHARACTER(LEN=256) :: basis_message, enumeration_message

      diagnostics = t_coordination_diagnostics()
      success = .FALSE.
      message = ""
      ALLOCATE(neighbors(0))

      IF (SIZE(candidate_positions,1) /= 3 .OR. SIZE(candidate_positions,2) /= SIZE(candidate_ids)) THEN
         message = "candidate_positions must have shape (3,size(candidate_ids))"
         RETURN
      END IF
      IF (nneighbors <= 0) THEN
         message = "nneighbors must be positive"
         RETURN
      END IF
      IF (SIZE(candidate_ids) == 0) THEN
         message = "candidate set is empty"
         RETURN
      END IF
      IF (.NOT.ALL(ieee_is_finite(lattice)) .OR. .NOT.ALL(ieee_is_finite(center)) .OR. &
          .NOT.ALL(ieee_is_finite(candidate_positions))) THEN
         message = "lattice and Cartesian coordinates must be finite"
         RETURN
      END IF
      IF (ANY(candidate_ids <= 0)) THEN
         message = "candidate physical-atom identifiers must be positive"
         RETURN
      END IF
      DO candidate = 1, SIZE(candidate_ids)-1
         IF (ANY(candidate_ids(candidate+1:) == candidate_ids(candidate))) THEN
            message = "candidate physical-atom identifiers must be unique"
            RETURN
         END IF
      END DO

      range_limit = default_maximum_image_range
      IF (PRESENT(maximum_image_range)) range_limit = maximum_image_range
      IF (range_limit < 0) THEN
         message = "maximum_image_range must be nonnegative"
         RETURN
      END IF

      CALL periodic_basis_smallest_singular_value(lattice, periodic, sigma_minimum, basis_success, basis_message)
      IF (.NOT.basis_success) THEN
         message = TRIM(basis_message)
         RETURN
      END IF

      maximum_offset = 0.0
      DO candidate = 1, SIZE(candidate_ids)
         maximum_offset = MAX(maximum_offset, vector_norm(candidate_positions(:,candidate)-center))
      END DO
      scale = MAX(1.0, MAX(MAXVAL(ABS(lattice)), maximum_offset))
      zero_tolerance = 1000.0*EPSILON(scale)*scale

      DO image_range = 0, range_limit
         CALL enumerate_images(lattice, periodic, center, candidate_positions, candidate_ids, image_range, &
                               zero_tolerance, images, nimages, enumeration_success, enumeration_message)
         IF (.NOT.enumeration_success) THEN
            message = TRIM(enumeration_message)
            RETURN
         END IF
         IF (nimages < nneighbors) THEN
            IF (.NOT.ANY(periodic)) THEN
               message = "candidate set contains fewer nonzero images than requested neighbours"
               RETURN
            END IF
            CYCLE
         END IF
         CALL sort_neighbors(images, nimages)

         IF (.NOT.ANY(periodic)) THEN
            lower_bound = HUGE(1.0)
            EXIT
         END IF
         IF (nimages < nneighbors+1) CYCLE
         lower_bound = sigma_minimum*REAL(image_range+1) - maximum_offset
         IF (lower_bound > images(nneighbors+1)%distance + &
             1000.0*EPSILON(MAX(1.0,images(nneighbors+1)%distance))*MAX(1.0,images(nneighbors+1)%distance)) EXIT
      END DO

      IF (image_range > range_limit) THEN
         WRITE(message,'(a,i0,a)') "Periodic image search did not prove completeness within range ", range_limit, &
                                  "; increase maximum_image_range or inspect the lattice"
         RETURN
      END IF

      DEALLOCATE(neighbors)
      ALLOCATE(neighbors(nneighbors))
      neighbors = images(:nneighbors)
      diagnostics%search_range = image_range
      diagnostics%enumerated_images = nimages
      diagnostics%minimum_unseen_distance_bound = lower_bound
      IF (nimages >= nneighbors+1) THEN
         diagnostics%has_next_neighbor = .TRUE.
         diagnostics%next_neighbor_distance = images(nneighbors+1)%distance
         diagnostics%shell_gap = images(nneighbors+1)%distance-images(nneighbors)%distance
         IF (images(nneighbors)%distance > zero_tolerance) THEN
            diagnostics%next_to_last_ratio = images(nneighbors+1)%distance/images(nneighbors)%distance
         END IF
      END IF
      success = .TRUE.
      message = "OK"
   END SUBROUTINE resolve_periodic_neighbours_cartesian

   SUBROUTINE periodic_basis_smallest_singular_value(lattice, periodic, sigma_minimum, success, message)
      REAL, INTENT(IN) :: lattice(3,3)
      LOGICAL, INTENT(IN) :: periodic(3)
      REAL, INTENT(OUT) :: sigma_minimum
      LOGICAL, INTENT(OUT) :: success
      CHARACTER(LEN=*), INTENT(OUT) :: message

      REAL :: basis(3,3), singular_values(3), u_dummy(1,1), vt_dummy(1,1), work(64), scale
      INTEGER :: direction, column, ndim, info

      basis = 0.0
      singular_values = 0.0
      ndim = COUNT(periodic)
      sigma_minimum = HUGE(1.0)
      success = .FALSE.
      message = ""
      IF (ndim == 0) THEN
         success = .TRUE.
         message = "OK"
         RETURN
      END IF

      column = 0
      DO direction = 1, 3
         IF (.NOT.periodic(direction)) CYCLE
         column = column + 1
         basis(:,column) = lattice(:,direction)
      END DO
      scale = MAXVAL(ABS(basis(:,:ndim)))
      IF (scale <= 0.0) THEN
         message = "periodic lattice basis contains a zero direction"
         RETURN
      END IF
      CALL DGESVD('N','N',3,ndim,basis,3,singular_values,u_dummy,1,vt_dummy,1,work,SIZE(work),info)
      IF (info /= 0) THEN
         WRITE(message,'(a,i0)') "SVD of periodic lattice basis failed with INFO=", info
         RETURN
      END IF
      sigma_minimum = singular_values(ndim)
      IF (sigma_minimum <= 1000.0*EPSILON(scale)*scale) THEN
         message = "periodic lattice vectors are linearly dependent or ill conditioned"
         RETURN
      END IF
      success = .TRUE.
      message = "OK"
   END SUBROUTINE periodic_basis_smallest_singular_value

   SUBROUTINE enumerate_images(lattice, periodic, center, candidate_positions, candidate_ids, image_range, &
                               zero_tolerance, images, nimages, success, message)
      REAL, INTENT(IN) :: lattice(3,3), center(3), candidate_positions(:,:), zero_tolerance
      LOGICAL, INTENT(IN) :: periodic(3)
      INTEGER, INTENT(IN) :: candidate_ids(:), image_range
      TYPE(t_periodic_neighbor), ALLOCATABLE, INTENT(OUT) :: images(:)
      INTEGER, INTENT(OUT) :: nimages
      LOGICAL, INTENT(OUT) :: success
      CHARACTER(LEN=*), INTENT(OUT) :: message

      INTEGER(KIND=8) :: maximum_images_8
      INTEGER :: lower(3), upper(3), translation(3), candidate, maximum_images, ix, iy, iz
      REAL :: displacement(3), distance

      WHERE (periodic)
         lower = -image_range
         upper = image_range
      ELSEWHERE
         lower = 0
         upper = 0
      END WHERE
      maximum_images_8 = INT(SIZE(candidate_ids),8)
      maximum_images_8 = maximum_images_8*INT(upper(1)-lower(1)+1,8) &
                       *INT(upper(2)-lower(2)+1,8)*INT(upper(3)-lower(3)+1,8)
      success = .FALSE.
      message = ""
      IF (maximum_images_8 > HUGE(maximum_images)) THEN
         ALLOCATE(images(0))
         nimages = 0
         message = "periodic image enumeration exceeds integer capacity"
         RETURN
      END IF
      maximum_images = INT(maximum_images_8)
      ALLOCATE(images(maximum_images))
      nimages = 0
      DO ix = lower(1), upper(1)
         DO iy = lower(2), upper(2)
            DO iz = lower(3), upper(3)
               translation = [ix,iy,iz]
               DO candidate = 1, SIZE(candidate_ids)
                  displacement = candidate_positions(:,candidate)-center &
                               + MATMUL(lattice,REAL(translation))
                  distance = vector_norm(displacement)
                  ! Suppress the central atom itself if it occurs in the
                  ! candidate set; nonzero periodic images remain legitimate.
                  IF (distance <= zero_tolerance) CYCLE
                  nimages = nimages + 1
                  images(nimages)%atom_id = candidate_ids(candidate)
                  images(nimages)%translation = translation
                  images(nimages)%displacement = displacement
                  images(nimages)%distance = distance
               END DO
            END DO
         END DO
      END DO
      success = .TRUE.
      message = "OK"
   END SUBROUTINE enumerate_images

   SUBROUTINE sort_neighbors(images, nimages)
      TYPE(t_periodic_neighbor), INTENT(INOUT) :: images(:)
      INTEGER, INTENT(IN) :: nimages

      IF (nimages > 1) CALL quicksort_neighbors(images, 1, nimages)
   END SUBROUTINE sort_neighbors

   RECURSIVE SUBROUTINE quicksort_neighbors(images, left, right)
      TYPE(t_periodic_neighbor), INTENT(INOUT) :: images(:)
      INTEGER, INTENT(IN) :: left, right
      TYPE(t_periodic_neighbor) :: pivot, temporary
      INTEGER :: i, j

      i = left
      j = right
      pivot = images((left+right)/2)
      DO
         DO WHILE (neighbor_precedes(images(i),pivot))
            i = i+1
         END DO
         DO WHILE (neighbor_precedes(pivot,images(j)))
            j = j-1
         END DO
         IF (i <= j) THEN
            temporary = images(i)
            images(i) = images(j)
            images(j) = temporary
            i = i+1
            j = j-1
         END IF
         IF (i > j) EXIT
      END DO
      IF (left < j) CALL quicksort_neighbors(images,left,j)
      IF (i < right) CALL quicksort_neighbors(images,i,right)
   END SUBROUTINE quicksort_neighbors

   LOGICAL FUNCTION neighbor_precedes(first, second)
      TYPE(t_periodic_neighbor), INTENT(IN) :: first, second
      INTEGER :: direction

      IF (first%distance < second%distance) THEN
         neighbor_precedes = .TRUE.
         RETURN
      ELSE IF (first%distance > second%distance) THEN
         neighbor_precedes = .FALSE.
         RETURN
      END IF
      IF (first%atom_id < second%atom_id) THEN
         neighbor_precedes = .TRUE.
         RETURN
      ELSE IF (first%atom_id > second%atom_id) THEN
         neighbor_precedes = .FALSE.
         RETURN
      END IF
      DO direction = 1, 3
         IF (first%translation(direction) < second%translation(direction)) THEN
            neighbor_precedes = .TRUE.
            RETURN
         ELSE IF (first%translation(direction) > second%translation(direction)) THEN
            neighbor_precedes = .FALSE.
            RETURN
         END IF
      END DO
      neighbor_precedes = .FALSE.
   END FUNCTION neighbor_precedes

   REAL FUNCTION vector_norm(vector)
      REAL, INTENT(IN) :: vector(3)
      vector_norm = SQRT(MAX(0.0,DOT_PRODUCT(vector,vector)))
   END FUNCTION vector_norm

END MODULE m_local_coordination_core
