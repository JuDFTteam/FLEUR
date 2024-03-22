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
     REAL,    ALLOCATABLE :: qis(:,:,:)
     REAL,    ALLOCATABLE :: qal(:,:,:,:,:)
     REAL,    ALLOCATABLE :: qTot(:,:,:)
     CHARACTER(len=20),ALLOCATABLE:: weight_names(:)!This must be allocated in init of derived type

   CONTAINS
     PROCEDURE,PASS :: init => dos_init
     PROCEDURE      :: get_weight_eig
     PROCEDURE      :: get_num_spins
     PROCEDURE      :: get_num_weights
     PROCEDURE      :: get_weight_name
     PROCEDURE      :: sym_weights
  END TYPE t_dos

CONTAINS

subroutine sym_weights(this)
  class(t_dos),intent(inout):: this

  integer:: i,j

  DO i=0,size(this%qal,1)-1
    DO j=1,size(this%qal,2)
      call this%sym_weights_eigdos(this%qal(i,j,:,:,:))
    enddo
  ENDDO  
  call this%sym_weights_eigdos(this%qis(:,:,:))
  call this%sym_weights_eigdos(this%qtot(:,:,:))
end subroutine

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

    if (id==1) THEN
       get_weight_eig=this%qTot
       if (all(this%qis==0.0))  then
          get_weight_eig= 1.0
       END IF
    END IF
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
    if (id==3) get_weight_eig=1.*this%jsym
    ind=3
    DO ntype=1,size(this%qal,2)
      DO l=0,3
        ind=ind+1
        if (ind==id) get_weight_eig=this%qal(l,ntype,:,:,:)
      ENDDO
    ENDDO
  end function

SUBROUTINE dos_init(thisDOS,input,atoms,kpts,banddos,eig)
  USE m_types_input
  USE m_types_atoms
  USE m_types_banddos
  USE m_types_kpts
  IMPLICIT NONE
  CLASS(t_dos),           INTENT(INOUT) :: thisDOS
  TYPE(t_input),          INTENT(IN)    :: input
  TYPE(t_atoms),          INTENT(IN)    :: atoms
  TYPE(t_kpts),           INTENT(IN)    :: kpts
  TYPE(t_banddos),         INTENT(IN)    :: banddos
  real,intent(in)                       :: eig(:,:,:)

  INTEGER :: ntype,l,i,ind
  character :: spdfg(0:4)=["s","p","d","f","g"]
  thisDOS%name_of_dos="Local"
  thisDOS%neq=atoms%neq(banddos%dos_typelist)
  thisDOS%eig=eig
  ALLOCATE(thisDOS%jsym(input%neig,kpts%nkpt,input%jspins))
  ALLOCATE(thisDOS%qis(input%neig,kpts%nkpt,input%jspins))
  ALLOCATE(thisDOS%qal(0:3,size(banddos%dos_typelist),input%neig,kpts%nkpt,input%jspins))
  ALLOCATE(thisDOS%qTot(input%neig,kpts%nkpt,input%jspins))

  thisDOS%jsym = 0
  thisDOS%qis = 0.0
  thisDOS%qal = 0.0
  thisDOS%qTot = 0.0

  allocate(thisDOS%weight_names(3+4*size(banddos%dos_typelist)))
  thisDOS%weight_names(1)="Total"
  thisDOS%weight_names(2)="INT"
  thisDOS%weight_names(3)="Sym"
  ind=3
  DO ntype=1,size(banddos%dos_typelist)
    DO l=0,3
      ind=ind+1
      write(thisDOS%weight_names(ind),"(a,i0,a)") "MT:",banddos%dos_typelist(ntype),spdfg(l)
    ENDDO
  ENDDO


END SUBROUTINE dos_init



END MODULE m_types_dos
