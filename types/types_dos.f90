!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_dos
  USE m_juDFT
  IMPLICIT NONE
  PRIVATE
  PUBLIC:: t_dos
  TYPE:: t_dos
     INTEGER, ALLOCATABLE :: jsym(:,:,:)
     INTEGER, ALLOCATABLE :: ksym(:,:,:)
     REAL,    ALLOCATABLE :: qis(:,:,:)
     REAL,    ALLOCATABLE :: qal(:,:,:,:,:)
     REAL,    ALLOCATABLE :: qvac(:,:,:,:)
     REAL,    ALLOCATABLE :: qvlay(:,:,:,:,:)
     COMPLEX, ALLOCATABLE :: qstars(:,:,:,:,:,:)
   CONTAINS
     PROCEDURE,PASS :: init => dos_init
  END TYPE t_dos

  TYPE,ABSTRACT:: t_eigdesc
    !each eigenvalue might be described by weights
    CHARACTER(len=20),ALLOCATABLE:: weight_names(:)
  CONTAINS
    procedure :: get_weight_name
    procedure,DEFERRED :: get_weight
    procedure :: get_num_weights
    procedure :: write_hdf5
    procedure :: read_hdf5
    procedure :: write
    procedure :: read
  END TYPE


CONTAINS

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

END SUBROUTINE dos_init

END MODULE m_types_dos
