!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_dos
  USE m_juDFT
  USE m_types_eigdos
  IMPLICIT NONE
  PRIVATE
  PUBLIC:: t_dos
  TYPE,extends(t_eigdos):: t_dos
     INTEGER,ALLOCATABLE :: neq(:)
     INTEGER, ALLOCATABLE :: jsym(:,:,:)
     INTEGER, ALLOCATABLE :: ksym(:,:,:)
     REAL,    ALLOCATABLE :: qis(:,:,:)
     REAL,    ALLOCATABLE :: qal(:,:,:,:,:)
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
  END TYPE t_dos

CONTAINS

  integer function get_num_weights(this)
    class(t_dos),intent(in):: this
    get_num_weights=0
    if (allocated(this%weight_names)) get_num_weights=size(this%weight_names)
  end function

  character(len=20) function get_weight_name(this,id)
    class(t_dos),intent(in):: this
    INTEGER,intent(in)         :: id
    if (.not.allocated(this%weight_names)) call judft_error("No weight names in t_eigdos")
    if (id>size(this%weight_names)) call judft_error("Not enough weight names in t_eigdos")
    get_weight_name=this%weight_names(id)
  end function




  integer function get_num_spins(this)
    class(t_dos),intent(in):: this
    get_num_spins= size(this%qis,3)
  end function

  function get_weight_eig(this,id)
    class(t_dos),intent(in):: this
    INTEGER,intent(in)     :: id
    real,allocatable:: get_weight_eig(:,:,:)

    INTEGER :: ind,l,ntype,i
    allocate(get_weight_eig,mold=this%qis)

    if (id==1) get_weight_eig=1.0
    if (id==2) THEN
      get_weight_eig=this%qis
      if (all(get_weight_eig==0.0))  then
        !No INT dos calculated so far...
        get_weight_eig=1.0
        DO ntype=1,size(this%qal,2)
          DO l=0,3
            get_weight_eig=get_weight_eig-this%qal(l,ntype,:,:,:)*this%neq(ntype)
          ENDDO
        ENDDO
      endif
    endif
    ind=2
    DO ntype=1,size(this%qal,2)
      DO l=0,3
        ind=ind+1
        if (ind==id) get_weight_eig=this%qal(l,ntype,:,:,:)
      ENDDO
    ENDDO
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

SUBROUTINE dos_init(thisDOS,input,atoms,kpts,vacuum,eig)
  USE m_types_input
  USE m_types_atoms
  USE m_types_vacuum
  USE m_types_kpts
  IMPLICIT NONE
  CLASS(t_dos),           INTENT(INOUT) :: thisDOS
  TYPE(t_input),          INTENT(IN)    :: input
  TYPE(t_atoms),          INTENT(IN)    :: atoms
  TYPE(t_kpts),           INTENT(IN)    :: kpts
  TYPE(t_vacuum),         INTENT(IN)    :: vacuum
  real,intent(in)                       :: eig(:,:,:)

  INTEGER :: ntype,l,i,ind
  character :: spdfg(0:4)=["s","p","d","f","g"]
  thisDOS%name_of_dos="Local"
  thisDOS%neq=atoms%neq
  thisDOS%eig=eig
  ALLOCATE(thisDOS%jsym(input%neig,kpts%nkpt,input%jspins))
  ALLOCATE(thisDOS%ksym(input%neig,kpts%nkpt,input%jspins))
  ALLOCATE(thisDOS%qis(input%neig,kpts%nkpt,input%jspins))
  ALLOCATE(thisDOS%qal(0:3,atoms%ntype,input%neig,kpts%nkpt,input%jspins))
  ALLOCATE(thisDOS%qvac(input%neig,2,kpts%nkpt,input%jspins))
  ALLOCATE(thisDOS%qvlay(input%neig,vacuum%layerd,2,kpts%nkpt,input%jspins))
  ALLOCATE(thisDOS%qstars(vacuum%nstars,input%neig,vacuum%layerd,2,kpts%nkpt,input%jspins))

  thisDOS%jsym = 0
  thisDOS%ksym = 0
  thisDOS%qis = 0.0
  thisDOS%qal = 0.0
  thisDOS%qvac = 0.0
  thisDOS%qvlay = 0.0
  thisDOS%qstars = CMPLX(0.0,0.0)

  allocate(thisDOS%weight_names(4+4*atoms%ntype+vacuum%layerd*(vacuum%nstars+1)))
  thisDOS%weight_names(1)="Total"
  thisDOS%weight_names(2)="INT"
  ind=2
  DO ntype=1,atoms%ntype
    DO l=0,3
      ind=ind+1
      write(thisDOS%weight_names(ind),"(a,i0,a)") "MT:",ntype,spdfg(l)
    ENDDO
  ENDDO
  ind=ind+1
  thisDOS%weight_names(ind)="VAC1"
  ind=ind+1
  thisDOS%weight_names(ind)="VAC2"
  do i=1,vacuum%layerd
    ind=ind+1
    write(thisDOS%weight_names(ind),"(a,i0)") "LAYER1-",i
    ind=ind+1
    write(thisDOS%weight_names(ind),"(a,i0)") "LAYER2-",i
  end do
  DO l=1,vacuum%layerd
    do i=1,vacuum%nstars
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



END MODULE m_types_dos
