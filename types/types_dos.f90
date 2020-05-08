!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_dos
  USE m_juDFT
  USE m_types_eigdesc
  IMPLICIT NONE
  PRIVATE
  PUBLIC:: t_dos
  TYPE,extends(t_eigdesc):: t_dos
     INTEGER, ALLOCATABLE :: jsym(:,:,:)
     INTEGER, ALLOCATABLE :: ksym(:,:,:)
     REAL,    ALLOCATABLE :: qis(:,:,:)
     REAL,    ALLOCATABLE :: qal(:,:,:,:,:)
     REAL,    ALLOCATABLE :: qvac(:,:,:,:)
     REAL,    ALLOCATABLE :: qvlay(:,:,:,:,:)
     COMPLEX, ALLOCATABLE :: qstars(:,:,:,:,:,:)
   CONTAINS
     PROCEDURE,PASS :: init => dos_init
     PROCEDURE      :: get_weight
  END TYPE t_dos

CONTAINS
  function get_weight(this,id)
    class(t_dos),intent(in):: this
    INTEGER,intent(in)     :: id
    real,allocatable:: get_weight(:,:,:)

    INTEGER :: ind,l,ntype,i
    allocate(get_weight,mold=this%qis)

    if (id==1) get_weight=this%qis
    ind=1
    DO ntype=1,size(this%qal,2)
      DO l=0,3
        ind=ind+1
        if (ind==id) get_weight=this%qal(l,ntype,:,:,:)
      ENDDO
    ENDDO
    do i=1,2
      ind=ind+1
      if (ind==id) get_weight=this%qvac(:,i,:,:)
    end do
    do i=1,size(this%qvlay,2)
      ind=ind+1
      if (ind==id) get_weight=this%qvlay(:,i,1,:,:)
      ind=ind+1
      if (ind==id) get_weight=this%qvlay(:,i,2,:,:)
    end do
    DO l=1,size(this%qstars,3)
      do i=1,size(this%qstars,1)
        ind=ind+1
        if (ind==id) get_weight=real(this%qstars(i,:,l,1,:,:))
        ind=ind+1
        if (ind==id) get_weight=aimag(this%qstars(i,:,l,1,:,:))
        ind=ind+1
        if (ind==id) get_weight=real(this%qstars(i,:,l,2,:,:))
        ind=ind+1
        if (ind==id) get_weight=aimag(this%qstars(i,:,l,2,:,:))
      end do
    end do
  end function

SUBROUTINE dos_init(thisDOS,input,atoms,kpts,vacuum)
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

  INTEGER :: ntype,l,i,ind
  character :: spdfg(0:4)=["s","p","d","f","g"]
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
  ind=1

  allocate(thisDOS%weight_names(3+3*atoms%ntype+vacuum%layerd*(vacuum%nstars+1)))
  thisDOS%weight_names(ind)="INT"
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
    write(thisDOS%weight_names(ind),"(a,i0)") "LAYER1:",i
    ind=ind+1
    write(thisDOS%weight_names(ind),"(a,i0)") "LAYER2:",i
  end do
  DO l=1,vacuum%layerd
    do i=1,vacuum%nstars
      ind=ind+1
      write(thisDOS%weight_names(ind),"(a,i0,a,i0)") "R(gVAC1):",l,"-",i
      ind=ind+1
      write(thisDOS%weight_names(ind),"(a,i0,a,i0)") "I(gVAC1):",l,"-",i
      ind=ind+1
      write(thisDOS%weight_names(ind),"(a,i0,a,i0)") "R(gVAC2):",l,"-",i
      ind=ind+1
      write(thisDOS%weight_names(ind),"(a,i0,a,i0)") "I(gVAC2):",l,"-",i
    end do
  end do


END SUBROUTINE dos_init

END MODULE m_types_dos
