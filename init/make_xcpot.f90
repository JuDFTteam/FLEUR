!--------------------------------------------------------------------------------
! Copyright (c) 2017 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_make_xcpot
   use m_juDFT
   implicit none
   private
   public make_xcpot

contains
   subroutine make_xcpot(fmpi,xcpot, atoms, input)
      use m_types_xcpot
      use m_types_atoms
      use m_types_input
      USE m_types_xcpot_libxc
      USE m_types_xcpot_inbuild
      USE m_types_xcpot_inbuild_nofunction
      USE m_types_mpi

      TYPE(t_mpi),INTENT(IN)        :: fmpi
      TYPE(t_input), INTENT(IN)     :: input
      TYPE(t_atoms), INTENT(IN)     :: atoms
      CLASS(t_xcpot), ALLOCATABLE, INTENT(INOUT) :: xcpot

      INTEGER              :: func_vxc_id_c, func_vxc_id_x, func_exc_id_c, func_exc_id_x
      REAL                 :: gmaxxc
      LOGICAL              :: l_libxc
      LOGICAL              :: l_inbuild
      CHARACTER(len=10)    :: inbuild_name
      LOGICAL              :: l_relativistic

      !Finish setup of xcpot
      xcpot%l_libxc = (xcpot%inbuild_name == "LibXC")
      IF (xcpot%l_libxc) THEN
         write (*,*) "func_ids", xcpot%func_vxc_id_c, xcpot%func_vxc_id_x, xcpot%func_exc_id_c, xcpot%func_exc_id_x

         func_vxc_id_c  = xcpot%func_vxc_id_c
         func_vxc_id_x  = xcpot%func_vxc_id_x
         func_exc_id_c  = xcpot%func_exc_id_c
         func_exc_id_x  = xcpot%func_exc_id_x
         gmaxxc         = xcpot%gmaxxc
         l_libxc        = .TRUE.
         l_inbuild      = .FALSE.
         inbuild_name   = xcpot%inbuild_name
         l_relativistic = xcpot%l_relativistic

         DEALLOCATE (xcpot)
         ALLOCATE (t_xcpot_libxc::xcpot)
         SELECT TYPE (xcpot)
         CLASS is (t_xcpot_libxc)!just allocated like this
            CALL xcpot%init(func_vxc_id_x, func_vxc_id_c, func_exc_id_x, func_exc_id_c, input%jspins)
         END SELECT
         xcpot%gmaxxc         = gmaxxc
         xcpot%l_libxc        = l_libxc
         xcpot%l_inbuild      = l_inbuild
         xcpot%inbuild_name   = inbuild_name
         xcpot%l_relativistic = l_relativistic
      ELSE
        SELECT TYPE (xcpot)
        CLASS is (t_xcpot_inbuild_nf)
          Call Xcpot%Init(Atoms%Ntype)
          CLASS DEFAULT
          CALL judft_error("Error in setup xcpot")
        END SELECT
      END IF

   end subroutine make_xcpot
end MODULE m_make_xcpot
