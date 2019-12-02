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
      REAL, ALLOCATABLE           :: voltet(:)
      REAL, ALLOCATABLE           :: sc_list(:, :) !list for all information about folding of bandstructure (need for unfoldBandKPTS)((k(x,y,z),K(x,y,z),m(g1,g2,g3)),(nkpt),k_original(x,y,z))
   CONTAINS
      procedure :: get_nk => kpts_get_nk
      procedure :: to_first_bz => kpts_to_first_bz
      procedure :: is_kpt => kpts_is_kpt
      procedure :: write_kpts
      generic   :: write(unformatted) => write_kpts
   ENDTYPE t_kpts
contains
   subroutine write_kpts(dtv, unit, iostat, iomsg)
      implicit NONE
      class(t_kpts), intent(in)   :: dtv
      integer, intent(in)         :: unit         ! Internal unit to write to.
      integer, intent(out)        :: iostat      ! non zero on error, etc.
      character(*), intent(inout) :: iomsg  ! define if iostat non zero.

      write(unit, iostat=iostat, iomsg=iomsg) dtv%nkpt, dtv%ntet, dtv%posScale, dtv%l_gamma

      write(unit, iostat=iostat, iomsg=iomsg) shape(dtv%bk)
      write(unit, iostat=iostat, iomsg=iomsg) dtv%bk

      write(unit, iostat=iostat, iomsg=iomsg) shape(dtv%wtkpt)
      write(unit, iostat=iostat, iomsg=iomsg) dtv%wtkpt

      write(unit, iostat=iostat, iomsg=iomsg) dtv%nkptf, dtv%nkpt3, dtv%kPointDensity

      write(unit, iostat=iostat, iomsg=iomsg) shape(dtv%bkf)
      write(unit, iostat=iostat, iomsg=iomsg) dtv%bkf

      write(unit, iostat=iostat, iomsg=iomsg) shape(dtv%bkp)
      write(unit, iostat=iostat, iomsg=iomsg) dtv%bkp

      write(unit, iostat=iostat, iomsg=iomsg) shape(dtv%bksym)
      write(unit, iostat=iostat, iomsg=iomsg) dtv%bksym

      write(unit, iostat=iostat, iomsg=iomsg) dtv%numSpecialPoints

      write(unit, iostat=iostat, iomsg=iomsg) shape(dtv%specialPointIndices)
      write(unit, iostat=iostat, iomsg=iomsg) dtv%specialPointIndices

      write(unit, iostat=iostat, iomsg=iomsg) shape(dtv%specialPointNames)
      write(unit, iostat=iostat, iomsg=iomsg) dtv%specialPointNames

      write(unit, iostat=iostat, iomsg=iomsg) shape(dtv%specialPoints)
      write(unit, iostat=iostat, iomsg=iomsg) dtv%specialPoints

      write(unit, iostat=iostat, iomsg=iomsg) shape(dtv%ntetra)
      write(unit, iostat=iostat, iomsg=iomsg) dtv%ntetra

      write(unit, iostat=iostat, iomsg=iomsg) shape(dtv%voltet)
      write(unit, iostat=iostat, iomsg=iomsg) dtv%voltet

      write(unit, iostat=iostat, iomsg=iomsg) shape(dtv%sc_list)
      write(unit, iostat=iostat, iomsg=iomsg) dtv%sc_list
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
