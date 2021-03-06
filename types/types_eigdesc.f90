!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_eigdesc
  USE m_juDFT
  IMPLICIT NONE
  PRIVATE
  PUBLIC:: t_eigdesc,t_eigdesc_list

  TYPE:: t_eigdesc
    !each eigenvalue might be described by weights
    CHARACTER(len=20),ALLOCATABLE:: weight_names(:)!This must be allocated in init of derived type
  CONTAINS
    procedure :: get_weight_name
    procedure :: get_weight !should be overwritten in derived type
    procedure :: get_num_weights
    !procedure :: write_hdf5 TODO!!
    !procedure :: read_hdf5
    procedure :: write
    procedure :: read
  END TYPE

  TYPE:: t_eigdesc_list
    CLASS(t_eigdesc),POINTER :: p
  END TYPE

CONTAINS

  subroutine write(this,id)
    class(t_eigdesc),INTENT(IN):: this
    INTEGER,INTENT(IN)         :: id

    INTEGER:: n
    real,allocatable :: w(:,:,:)

    if (this%get_num_weights()<1) return
    w=this%get_weight(1)
    write(id) this%get_num_weights(),shape(w)
    write(id) this%get_weight_name(id),w
    do n=2,this%get_num_weights()
      w=this%get_weight(1)
      write(id) this%get_weight_name(id),w
    enddo
  end subroutine

  subroutine read(this,id)
    class(t_eigdesc),INTENT(IN):: this
    INTEGER,INTENT(IN)         :: id
  end subroutine

  function get_weight(this,id)
    class(t_eigdesc),intent(in):: this
    INTEGER,intent(in)         :: id
    real,allocatable:: get_weight(:,:,:)
  end function

  integer function get_num_weights(this)
    class(t_eigdesc),intent(in):: this
    get_num_weights=0
    if (allocated(this%weight_names)) get_num_weights=size(this%weight_names)
  end function

  character(len=20) function get_weight_name(this,id)
    class(t_eigdesc),intent(in):: this
    INTEGER,intent(in)         :: id
    if (.not.allocated(this%weight_names)) call judft_error("No weight names in t_eigdesc")
    if (id>size(this%weight_names)) call judft_error("Not enough weight names in t_eigdesc")
    get_weight_name=this%weight_names(id)
  end function

END MODULE m_types_eigdesc
