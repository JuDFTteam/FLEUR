!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_kpts
   INTEGER, PARAMETER:: kpts_by_number = 1
   INTEGER, PARAMETER:: kpts_by_mesh = 2
   INTEGER, PARAMETER:: kpts_by_list = 3

   TYPE t_kpts
      INTEGER :: specificationType

      !no
      INTEGER               :: nkpt
      INTEGER               :: ntet
      INTEGER               :: ntria
      REAL                  :: posScale
      LOGICAL               :: l_gamma
      !(3,nkpt) k-vectors internal units
      REAL, ALLOCATABLE      :: bk(:, :)
      !(nkpts) weights
      REAL, ALLOCATABLE      :: wtkpt(:)
      INTEGER               :: nkptf = 0!<k-vectors in full BZ
      INTEGER               :: nkpt3(3)
      REAL                  :: kPointDensity(3) ! only used if k point set is defined as density
      REAL, ALLOCATABLE   :: bkf(:, :)
      INTEGER, ALLOCATABLE   :: bkp(:)
      INTEGER, ALLOCATABLE   :: bksym(:)
      INTEGER                       :: numSpecialPoints
      INTEGER, ALLOCATABLE          :: specialPointIndices(:)
      CHARACTER(LEN=50), ALLOCATABLE :: specialPointNames(:)
      REAL, ALLOCATABLE           :: specialPoints(:, :)
      INTEGER, ALLOCATABLE           :: ntetra(:, :)
      INTEGER, ALLOCATABLE           :: itria(:,:)
      REAL, ALLOCATABLE           :: voltria(:)
      REAL, ALLOCATABLE           :: voltet(:)
      REAL, ALLOCATABLE           :: sc_list(:, :) !list for all information about folding of bandstructure (need for unfoldBandKPTS)((k(x,y,z),K(x,y,z),m(g1,g2,g3)),(nkpt),k_original(x,y,z))
   CONTAINS
      procedure :: get_nk => kpts_get_nk
      procedure :: to_first_bz => kpts_to_first_bz
      procedure :: is_kpt => kpts_is_kpt
      procedure :: write_kpts
      generic   :: write(unformatted) => write_kpts
      procedure :: read_kpts
      generic   :: read(unformatted) => read_kpts
   ENDTYPE t_kpts
contains
   subroutine read_kpts(kpts, unit, iostat, iomsg)
      implicit NONE
      class(t_kpts), intent(inout)   :: kpts
      integer, intent(in)         :: unit         ! Internal unit to write to.
      integer, intent(out)        :: iostat      ! non zero on error, etc.
      character(*), intent(inout) :: iomsg  ! define if iostat non zero.

      integer :: shape_1d(1), shape_2d(2), shape_3d(3)

      read(unit, iostat=iostat, iomsg=iomsg) kpts%nkpt, kpts%ntet, kpts%posScale, kpts%l_gamma

      read(unit, iostat=iostat, iomsg=iomsg) shape_2d
      allocate(kpts%bk(shape_2d(1), shape_2d(2)))
      read(unit, iostat=iostat, iomsg=iomsg) kpts%bk

      read(unit, iostat=iostat, iomsg=iomsg) shape_1d
      allocate(kpts%wtkpt(shape_1d(1)))
      read(unit, iostat=iostat, iomsg=iomsg) kpts%wtkpt

      read(unit, iostat=iostat, iomsg=iomsg) kpts%nkptf, kpts%nkpt3, kpts%kPointDensity

      read(unit, iostat=iostat, iomsg=iomsg) shape_2d
      allocate(kpts%bkf(shape_2d(1), shape_2d(2)))
      read(unit, iostat=iostat, iomsg=iomsg) kpts%bkf

      read(unit, iostat=iostat, iomsg=iomsg) shape_1d
      allocate(kpts%bkp(shape_1d(1)))
      read(unit, iostat=iostat, iomsg=iomsg) kpts%bkp

      read(unit, iostat=iostat, iomsg=iomsg) shape_1d
      allocate(kpts%bksym(shape_1d(1)))
      read(unit, iostat=iostat, iomsg=iomsg) kpts%bksym

      read(unit, iostat=iostat, iomsg=iomsg) kpts%numSpecialPoints

      read(unit, iostat=iostat, iomsg=iomsg) shape_1d
      allocate(kpts%specialPointIndices(shape_1d(1)))
      read(unit, iostat=iostat, iomsg=iomsg) kpts%specialPointIndices

      read(unit, iostat=iostat, iomsg=iomsg) shape_1d
      allocate(kpts%specialPointNames(shape_1d(1)))
      read(unit, iostat=iostat, iomsg=iomsg) kpts%specialPointNames

      read(unit, iostat=iostat, iomsg=iomsg) shape_2d
      allocate(kpts%specialPoints(shape_2d(1), shape_2d(2)))
      read(unit, iostat=iostat, iomsg=iomsg) kpts%specialPoints

      read(unit, iostat=iostat, iomsg=iomsg) shape_2d
      allocate(kpts%ntetra(shape_2d(1), shape_2d(2)))
      read(unit, iostat=iostat, iomsg=iomsg) kpts%ntetra

      read(unit, iostat=iostat, iomsg=iomsg) shape_1d
      allocate(kpts%voltet(shape_2d(1)))
      read(unit, iostat=iostat, iomsg=iomsg) kpts%voltet

      read(unit, iostat=iostat, iomsg=iomsg) shape_2d
      allocate(kpts%sc_list(shape_2d(1), shape_2d(2)))
      read(unit, iostat=iostat, iomsg=iomsg) kpts%sc_list
   end subroutine read_kpts

   subroutine write_kpts(kpts, unit, iostat, iomsg)
      implicit NONE
      class(t_kpts), intent(in)   :: kpts
      integer, intent(in)         :: unit         ! Internal unit to write to.
      integer, intent(out)        :: iostat      ! non zero on error, etc.
      character(*), intent(inout) :: iomsg  ! define if iostat non zero.

      write(unit, iostat=iostat, iomsg=iomsg) kpts%nkpt, kpts%ntet, kpts%posScale, kpts%l_gamma

      write(unit, iostat=iostat, iomsg=iomsg) shape(kpts%bk)
      write(unit, iostat=iostat, iomsg=iomsg) kpts%bk

      write(unit, iostat=iostat, iomsg=iomsg) shape(kpts%wtkpt)
      write(unit, iostat=iostat, iomsg=iomsg) kpts%wtkpt

      write(unit, iostat=iostat, iomsg=iomsg) kpts%nkptf, kpts%nkpt3, kpts%kPointDensity

      write(unit, iostat=iostat, iomsg=iomsg) shape(kpts%bkf)
      write(unit, iostat=iostat, iomsg=iomsg) kpts%bkf

      write(unit, iostat=iostat, iomsg=iomsg) shape(kpts%bkp)
      write(unit, iostat=iostat, iomsg=iomsg) kpts%bkp

      write(unit, iostat=iostat, iomsg=iomsg) shape(kpts%bksym)
      write(unit, iostat=iostat, iomsg=iomsg) kpts%bksym

      write(unit, iostat=iostat, iomsg=iomsg) kpts%numSpecialPoints

      write(unit, iostat=iostat, iomsg=iomsg) shape(kpts%specialPointIndices)
      write(unit, iostat=iostat, iomsg=iomsg) kpts%specialPointIndices

      write(unit, iostat=iostat, iomsg=iomsg) shape(kpts%specialPointNames)
      write(unit, iostat=iostat, iomsg=iomsg) kpts%specialPointNames

      write(unit, iostat=iostat, iomsg=iomsg) shape(kpts%specialPoints)
      write(unit, iostat=iostat, iomsg=iomsg) kpts%specialPoints

      write(unit, iostat=iostat, iomsg=iomsg) shape(kpts%ntetra)
      write(unit, iostat=iostat, iomsg=iomsg) kpts%ntetra

      write(unit, iostat=iostat, iomsg=iomsg) shape(kpts%voltet)
      write(unit, iostat=iostat, iomsg=iomsg) kpts%voltet

      write(unit, iostat=iostat, iomsg=iomsg) shape(kpts%sc_list)
      write(unit, iostat=iostat, iomsg=iomsg) kpts%sc_list
   end subroutine write_kpts

   function kpts_get_nk(kpts, kpoint) result(ret_idx)
      ! get the index of a kpoint
      implicit NONE
      class(t_kpts), intent(in)    :: kpts
      real, intent(in)            :: kpoint(3)
      integer                     :: idx, ret_idx

      DO idx = 1, kpts%nkptf
         IF (all(abs(kpoint - kpts%bkf(:,idx)) < 1E-06)) THEN
            ret_idx = idx
            return
         END IF
      END DO
      ret_idx = 0
   end function kpts_get_nk

   function kpts_to_first_bz(kpts, kpoint) result(out_point)
      implicit NONE
      class(t_kpts), intent(in)  :: kpts
      real, intent(in)           :: kpoint(3)
      real                       :: out_point(3)

      out_point = kpoint - floor(kpoint)
   end function kpts_to_first_bz

   function kpts_is_kpt(kpts, kpoint) result(is_kpt)
      implicit none
      class(t_kpts), intent(in)  :: kpts
      real, intent(in)           :: kpoint(3)
      logical                    :: is_kpt

      is_kpt = kpts%get_nk(kpoint) > 0
   end function kpts_is_kpt
END MODULE m_types_kpts
