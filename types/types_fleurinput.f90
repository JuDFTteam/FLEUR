!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_fleurinput
  implicit none
  private

  type,abstract:: t_fleurinput
   contains
     procedure :: read_xml
     procedure :: mpi_bc
  end type t_fleurinput
  public t_fleurinput
  
contains
  subroutine read_xml(fleur_base,xml)
    class(t_fleur_base),intent(out):: fleur_base
    type(t_xml),INTENT(IN)         :: xml
  end subroutine read_xml

  subroutine mpi_bc(fleur_base,mpi)
    class(t_fleur_base),intent(INOUT):: fleur_base
    type(t_mpit),INTENT(IN)          :: mpi
  end subroutine mpi_bc
end MODULE m_types_fleurinput

