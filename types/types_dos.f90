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

CONTAINS

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
