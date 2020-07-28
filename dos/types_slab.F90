!--------------------------------------------------------------------------------
! Copyright (c) 2020 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_slab
  use m_judft
  use m_types_eigdos
  implicit none
  PRIVATE
  public t_slab
   TYPE,extends(t_eigdos):: t_slab
      INTEGER              :: nsld, nsl

      INTEGER, ALLOCATABLE :: nmtsl(:,:)
      INTEGER, ALLOCATABLE :: nslat(:,:)
      REAL,    ALLOCATABLE :: zsl(:,:)
      REAL,    ALLOCATABLE :: volsl(:)
      REAL,    ALLOCATABLE :: volintsl(:)
      REAL,    ALLOCATABLE :: qintsl(:,:,:,:)
      REAL,    ALLOCATABLE :: qmtsl(:,:,:,:)

      CONTAINS
         PROCEDURE,PASS :: init => slab_init
         PROCEDURE      :: get_num_weights
         PROCEDURE      :: get_weight_eig
         PROCEDURE      :: get_weight_name
   END TYPE t_slab
CONTAINS
  integer function get_num_weights(this)
    class(t_slab),intent(in):: this
    get_num_weights=2*this%nsl
  END function

  character(len=20) function get_weight_name(this,id)
    class(t_slab),intent(in):: this
    INTEGER,intent(in)         :: id

    INTEGER :: ind,n
    ind=0
    DO n=1,this%nsl
      ind=ind+1
      if (ind==id) write(get_weight_name,"(a,i0)") "SLAB(INT):",n
      ind=ind+1
      if (ind==id) write(get_weight_name,"(a,i0)") "SLAB(MT):",n
      IF(ind>id) return
    ENDDO
  end function

  function get_weight_eig(this,id)
    class(t_slab),intent(in):: this
    INTEGER,intent(in)      :: id
    real,allocatable:: get_weight_eig(:,:,:)


    INTEGER :: ind,n
    ind=0
    DO n=1,this%nsl
      ind=ind+1
      if (ind==id) get_weight_eig=this%qintsl(n,:,:,:)
      ind=ind+1
      if (ind==id) get_weight_eig=this%qmtsl(n,:,:,:)
    ENDDO
  end function



  SUBROUTINE slab_init(thisSlab,banddos,atoms,cell,input,kpts)

   USE m_types_setup
   USE m_types_kpts
   USE m_slabdim
   USE m_slabgeom

   IMPLICIT NONE

   CLASS(t_slab),      INTENT(INOUT) :: thisSlab
   TYPE(t_banddos),    INTENT(IN)    :: banddos

   TYPE(t_atoms),      INTENT(IN)    :: atoms
   TYPE(t_cell),       INTENT(IN)    :: cell
   TYPE(t_input),      INTENT(IN)    :: input
   TYPE(t_kpts),       INTENT(IN)    :: kpts

   INTEGER :: nsld
   thisSlab%name_of_dos="SLAB"
   nsld=1
   IF (banddos%l_slab.AND.banddos%dos) THEN
      CALL slab_dim(atoms, nsld)
      ALLOCATE (thisSlab%nmtsl(atoms%ntype,nsld))
      ALLOCATE (thisSlab%nslat(atoms%nat,nsld))
      ALLOCATE (thisSlab%zsl(2,nsld))
      ALLOCATE (thisSlab%volsl(nsld))
      ALLOCATE (thisSlab%volintsl(nsld))
      ALLOCATE (thisSlab%qintsl(nsld,input%neig,kpts%nkpt,input%jspins))
      ALLOCATE (thisSlab%qmtsl(nsld,input%neig,kpts%nkpt,input%jspins))
      CALL slabgeom(atoms,cell,nsld,thisSlab%nsl,thisSlab%zsl,thisSlab%nmtsl,&
                    thisSlab%nslat,thisSlab%volsl,thisSlab%volintsl)
   ELSE
     allocate(thisSlab%dos(0,0,0))
      ALLOCATE (thisSlab%nmtsl(1,1))
      ALLOCATE (thisSlab%nslat(1,1))
      ALLOCATE (thisSlab%zsl(1,1))
      ALLOCATE (thisSlab%volsl(1))
      ALLOCATE (thisSlab%volintsl(1))
      ALLOCATE (thisSlab%qintsl(1,1,1,input%jspins))
      ALLOCATE (thisSlab%qmtsl(1,1,1,input%jspins))
   END IF
   thisSlab%nsld = nsld

   thisSlab%nmtsl = 0
   thisSlab%nslat = 0
   thisSlab%zsl = 0.0
   thisSlab%volsl = 0.0
   thisSlab%volintsl = 0.0
   thisSlab%qintsl = 0.0
   thisSlab%qmtsl = 0.0

END SUBROUTINE slab_init

end module m_types_slab
