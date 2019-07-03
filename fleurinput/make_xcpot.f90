!--------------------------------------------------------------------------------
! Copyright (c) 2017 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_make_xcpot
  implicit none
  use m_juDFT
  private
  public make_xcpot
 
contains
  subroutine make_xcpot(xcpot,atoms,input)
    use m_types
    USE m_types_forcetheo_extended
    USE m_types_xcpot_libxc
    USE m_types_xcpot_inbuild
    USE m_types_xcpot_inbuild_nofunction

    TYPE(t_input),    INTENT(IN)     :: input
    TYPE(t_atoms),    INTENT(IN)     :: atoms
    CLASS(t_xcpot),ALLOCATABLE,INTENT(INOUT) :: xcpot


    INTEGER              ::func_vxc_id_c,func_vxc_id_x,func_exc_id_c,func_exc_id_x
 
    !Finish setup of xcpot
    IF (xcpot%l_libxc) THEN
       func_vxc_id_c=xcpot%func_vxc_id_c
       func_vxc_id_x=xcpot%func_vxc_id_x
       func_exc_id_c=xcpot%func_exc_id_c
       func_exc_id_x=xcpot%func_exc_id_x
       DEALLOCATE(xcpot)
       ALLOCATE(t_xcpot_libxc::xcpot)
       SELECT TYPE(xcpot)
       CLASS is (t_xcpot_libxc)!just allocated like this
          CALL xcpot%init(func_vxc_id_x,func_vxc_id_c,func_exc_id_x,func_exc_id_c,input%jspins)
       END SELECT
    ELSE
       SELECT TYPE(xcpot)
       CLASS is (t_xcpot_inbuild_nf)
          CALL xcpot%init(atoms%ntype)
       CLASS DEFAULT
          CALL judft_error("Error in setup xcpot")
       END SELECT
    END IF

end subroutine make_xcpot
end MODULE m_make_xcpot
