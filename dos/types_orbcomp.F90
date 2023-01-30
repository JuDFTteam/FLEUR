!--------------------------------------------------------------------------------
! Copyright (c) 2020 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_orbcomp
  use m_judft
  use m_types_eigdos
  implicit none
  PRIVATE
  public t_orbcomp
   TYPE,extends(t_eigdos):: t_orbcomp

      REAL, ALLOCATABLE    :: comp(:,:,:,:,:)
      REAL, ALLOCATABLE    :: qmtp(:,:,:,:)
      INTEGER,ALLOCATABLE  :: n_dos_to_na(:)
      CONTAINS
         PROCEDURE,PASS :: init => orbcomp_init
         PROCEDURE      :: get_num_weights
         PROCEDURE      :: get_weight_eig
         PROCEDURE      :: get_weight_name
         PROCEDURE      :: sym_weights
   END TYPE t_orbcomp
CONTAINS

  subroutine sym_weights(this)
    class(t_orbcomp),intent(inout):: this
    integer:: i,j
    return !This is done later in get_weights
    end subroutine


  integer function get_num_weights(this)
    class(t_orbcomp),intent(in):: this
    get_num_weights=23*size(this%comp,3)
  END function

  character(len=20) function get_weight_name(this,id)
    class(t_orbcomp),intent(in):: this
    INTEGER,intent(in)         :: id

    INTEGER :: ind,na,nc
    ind=0
    DO na=1,size(this%comp,3)
      DO nc=1,23
        ind=ind+1
        if (ind==id) THEN
          write(get_weight_name,"(a,i0,a,i0)") "ORB:",this%n_dos_to_na(na),",ind:",nc
          RETURN
        ELSE IF(ind>id) then
          CALL judft_error("Types_mcd: data not found")
        ENDIF
      ENDDO
    ENDDO
  end function

  function get_weight_eig(this,id)
    class(t_orbcomp),intent(in):: this
    INTEGER,intent(in)      :: id
    real,allocatable:: get_weight_eig(:,:,:)

    integer :: i,ind,na

    ind = 0
    DO na=1,size(this%comp,3)
      DO i= 1, 23
        ind = ind+1
        if (ind==id) THEN
            get_weight_eig=this%comp(:,i,na,:,:)*this%qmtp(:,na,:,:)/10000.
            call this%sym_weights_eigdos(get_weight_eig)
            return
        ENDIF
      ENDDO
    ENDDO
  end function


SUBROUTINE orbcomp_init(thisOrbcomp,input,banddos,atoms,kpts,eig)

   USE m_types_setup
   USE m_types_kpts

   IMPLICIT NONE

   CLASS(t_orbcomp),      INTENT(INOUT) :: thisOrbcomp
   TYPE(t_input),         INTENT(IN)    :: input
   TYPE(t_banddos),       INTENT(IN)    :: banddos

   TYPE(t_atoms),         INTENT(IN)    :: atoms
   TYPE(t_kpts),          INTENT(IN)    :: kpts
   REAL,INTENT(IN)                      :: eig(:,:,:)
   thisOrbcomp%n_dos_to_na=banddos%dos_atomlist
   IF ((banddos%l_orb).AND.banddos%dos) THEN
      ALLOCATE(thisOrbcomp%comp(input%neig,23,size(banddos%dos_atomlist),kpts%nkpt,input%jspins))
      ALLOCATE(thisOrbcomp%qmtp(input%neig,size(banddos%dos_atomlist),kpts%nkpt,input%jspins))
      thisOrbcomp%eig=eig
   ELSE
      ALLOCATE(thisOrbcomp%dos(0,0,0))
      ALLOCATE(thisOrbcomp%comp(1,1,0,1,input%jspins))
      ALLOCATE(thisOrbcomp%qmtp(1,0,1,input%jspins))
   END IF

   thisOrbcomp%comp = 0.0
   thisOrbcomp%qmtp = 0.0
   thisOrbcomp%name_of_dos="Orbcomp"
END SUBROUTINE orbcomp_init

end module m_types_orbcomp
