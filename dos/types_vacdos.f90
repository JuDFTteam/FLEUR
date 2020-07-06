!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_vacdos
  USE m_juDFT
  USE m_types_eigdos
  IMPLICIT NONE
  PRIVATE
  PUBLIC:: t_vacdos
  TYPE,extends(t_eigdos):: t_vacdos
     REAL,    ALLOCATABLE :: qvac(:,:,:,:)
     REAL,    ALLOCATABLE :: qvlay(:,:,:,:,:)
     COMPLEX, ALLOCATABLE :: qstars(:,:,:,:,:,:)
     CHARACTER(len=20),ALLOCATABLE:: weight_names(:)!This must be allocated in init of derived type

   CONTAINS
     PROCEDURE,PASS :: init => dos_init
     PROCEDURE      :: get_weight_eig
     PROCEDURE      :: get_num_spins
     PROCEDURE      :: get_num_weights
     PROCEDURE      :: get_weight_name
  END TYPE t_vacdos

CONTAINS

  integer function get_num_weights(this)
    class(t_vacdos),intent(in):: this
    get_num_weights=0
    if (allocated(this%weight_names)) get_num_weights=size(this%weight_names)
  end function

  character(len=20) function get_weight_name(this,id)
    class(t_vacdos),intent(in):: this
    INTEGER,intent(in)         :: id
    if (.not.allocated(this%weight_names)) call judft_error("No weight names in t_eigdos")
    if (id>size(this%weight_names)) call judft_error("Not enough weight names in t_eigdos")
    get_weight_name=this%weight_names(id)
  end function


  integer function get_num_spins(this)
    class(t_vacdos),intent(in):: this
    get_num_spins= size(this%qvac,4)
  end function

  function get_weight_eig(this,id)
    class(t_vacdos),intent(in):: this
    INTEGER,intent(in)     :: id
    real,allocatable:: get_weight_eig(:,:,:)

    INTEGER :: ind,l,ntype,i
    allocate(get_weight_eig(size(this%qvac,1),size(this%qvac,3),size(this%qvac,4)))

    ind=0
    do i=1,2
      ind=ind+1
      if (ind==id) get_weight_eig=this%qvac(:,i,:,:)
    end do
    do i=1,size(this%qvlay,2)
      ind=ind+1
      if (ind==id) get_weight_eig=this%qvlay(:,i,1,:,:)
      ind=ind+1
      if (ind==id) get_weight_eig=this%qvlay(:,i,2,:,:)
    end do
    DO l=1,size(this%qstars,3)
      do i=1,size(this%qstars,1)
        ind=ind+1
        if (ind==id) get_weight_eig=real(this%qstars(i,:,l,1,:,:))
        ind=ind+1
        if (ind==id) get_weight_eig=aimag(this%qstars(i,:,l,1,:,:))
        ind=ind+1
        if (ind==id) get_weight_eig=real(this%qstars(i,:,l,2,:,:))
        ind=ind+1
        if (ind==id) get_weight_eig=aimag(this%qstars(i,:,l,2,:,:))
      end do
    end do
  end function

SUBROUTINE dos_init(thisDOS,input,atoms,kpts,banddos,eig)
  USE m_types_input
  USE m_types_atoms
  USE m_types_banddos
  USE m_types_kpts
  IMPLICIT NONE
  CLASS(t_vacdos),           INTENT(INOUT) :: thisDOS
  TYPE(t_input),          INTENT(IN)    :: input
  TYPE(t_atoms),          INTENT(IN)    :: atoms
  TYPE(t_kpts),           INTENT(IN)    :: kpts
  TYPE(t_banddos),         INTENT(IN)    :: banddos
  real,intent(in)                       :: eig(:,:,:)

  INTEGER :: ntype,l,i,ind
  character :: spdfg(0:4)=["s","p","d","f","g"]
  thisDOS%name_of_dos="Vacuum"
  thisDOS%eig=eig
  ALLOCATE(thisDOS%qvac(input%neig,2,kpts%nkpt,input%jspins))
  ALLOCATE(thisDOS%qvlay(input%neig,banddos%layers,2,kpts%nkpt,input%jspins))
  ALLOCATE(thisDOS%qstars(banddos%nstars,input%neig,banddos%layers,2,kpts%nkpt,input%jspins))

  thisDOS%qvac = 0.0
  thisDOS%qvlay = 0.0
  thisDOS%qstars = CMPLX(0.0,0.0)

  if (.not.banddos%vacdos) THEN
    allocate(thisDOS%weight_names(0))
    RETURN
  endif
  allocate(thisDOS%weight_names(2+banddos%layers*(4*banddos%nstars+2)))
  ind=1
  thisDOS%weight_names(ind)="VAC1"
  ind=ind+1
  thisDOS%weight_names(ind)="VAC2"
  do i=1,banddos%layers
    ind=ind+1
    write(thisDOS%weight_names(ind),"(a,i0)") "LAYER1-",i
    ind=ind+1
    write(thisDOS%weight_names(ind),"(a,i0)") "LAYER2-",i
  end do
  DO l=1,banddos%layers
    do i=1,banddos%nstars
      ind=ind+1
      write(thisDOS%weight_names(ind),"(a,i0,a,i0)") "R(gVAC1)-",l,"-",i
      ind=ind+1
      write(thisDOS%weight_names(ind),"(a,i0,a,i0)") "I(gVAC1)-",l,"-",i
      ind=ind+1
      write(thisDOS%weight_names(ind),"(a,i0,a,i0)") "R(gVAC2)-",l,"-",i
      ind=ind+1
      write(thisDOS%weight_names(ind),"(a,i0,a,i0)") "I(gVAC2)-",l,"-",i
    end do
  end do


END SUBROUTINE dos_init



END MODULE m_types_vacdos
