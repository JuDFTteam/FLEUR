!--------------------------------------------------------------------------------
! Copyright (c) 2020 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_jdos
   use m_judft
   use m_types_eigdos
   implicit none
   PRIVATE
   public t_jdos
   TYPE,extends(t_eigdos):: t_jDOS

      REAL, ALLOCATABLE    :: comp(:,:,:,:,:)  !decomposition in percent
      REAL, ALLOCATABLE    :: qmtp(:,:,:)      !How much of the state is in the muffin-tin sphere
      REAL, ALLOCATABLE    :: occ(:,:,:)       !Occupation of the j-states
      INTEGER,ALLOCATABLE  :: n_dos_to_na(:)

      CONTAINS
         PROCEDURE,PASS :: init => jDOS_init
         PROCEDURE      :: get_weight_eig
         PROCEDURE      :: get_num_weights
         PROCEDURE      :: get_weight_name
         PROCEDURE      :: get_spins
   END TYPE t_jDOS
   CONTAINS

   pure integer function get_spins(this)
      CLASS(t_jdos),INTENT(IN)::this
      get_spins=1
   END function

   function get_weight_eig(this,id)
      class(t_jdos),intent(in):: this
      INTEGER,intent(in)      :: id
      real,allocatable:: get_weight_eig(:,:,:)

      integer :: i,l,jj,na

      ALLOCATE(get_weight_eig(size(this%comp,1),size(this%comp,5),1))

      i = 0
      DO na=1,size(this%comp,4)
         DO l= 0, 3
            DO jj = 1, MERGE(1,2,l==0)
               i = i+1
               if (i==id) get_weight_eig(:,:,1)=this%comp(:,l,jj,na,:)*this%qmtp(:,na,:)/10000.
               if (i>id) RETURN
            ENDDO
         ENDDO
      ENDDO
   end function

   integer function get_num_weights(this)
      class(t_jdos),intent(in):: this
      get_num_weights = 7*size(this%comp,4)
   end function


   character(len=20) function get_weight_name(this,id)
      class(t_jdos),intent(in):: this
      INTEGER,intent(in)         :: id
      integer :: i,l,jj,na
      character :: spdfg(0:4)=["s","p","d","f","g"]
      character(len=3) :: jname

      i = 0
      DO na=1,size(this%comp,4)
         DO l= 0, 3
            DO jj = -1, MERGE(-1,1,l==0), 2
               i = i+1
               WRITE(jname,'(i1,a,i1)') INT(2*l+jj),'-',2
               if (i==id) THEN
                  IF(l.EQ.0) write(get_weight_name,"(a,i0,a)") "jDOS:",this%n_dos_to_na(na),spdfg(l)
                  IF(l.NE.0) write(get_weight_name,"(a,i0,a,a)") "jDOS:",this%n_dos_to_na(na),spdfg(l),jname
               endif
               if (i>id) RETURN
            ENDDO
         ENDDO
      ENDDO

   end function


  SUBROUTINE jDOS_init(thisjDOS,input,banddos,atoms,kpts,eig)

      USE m_types_setup
      USE m_types_kpts

      IMPLICIT NONE

      CLASS(t_jDOS),         INTENT(INOUT) :: thisjDOS
      TYPE(t_input),         INTENT(IN)    :: input
      TYPE(t_banddos),       INTENT(IN)    :: banddos

      TYPE(t_atoms),         INTENT(IN)    :: atoms
      TYPE(t_kpts),          INTENT(IN)    :: kpts
      REAL,INTENT(IN)                      :: eig(:,:,:)

      thisjDOS%n_dos_to_na=banddos%dos_atomlist
      IF (banddos%l_jdos.AND.banddos%dos) THEN
         ALLOCATE(thisjDOS%comp(input%neig,0:3,2,size(banddos%dos_atomlist),kpts%nkpt),source = 0.0)
         ALLOCATE(thisjDOS%qmtp(input%neig,size(banddos%dos_atomlist),kpts%nkpt),source = 0.0)
         ALLOCATE(thisjDOS%occ(0:3,2,size(banddos%dos_atomlist)),source=0.0)
         thisjDOS%eig = eig
      ELSE
         ALLOCATE(thisjDOS%dos(0,0,0))
         ALLOCATE(thisjDOS%comp(1,1,1,0,1),source = 0.0)
         ALLOCATE(thisjDOS%qmtp(1,0,1),source = 0.0)
         ALLOCATE(thisjDOS%occ(1,1,0),source=0.0)
      END IF

      thisjDOS%name_of_dos="jDOS"

   END SUBROUTINE jDOS_init
end module m_types_jDOS
