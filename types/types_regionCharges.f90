!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_regionCharges

IMPLICIT NONE

PRIVATE

   TYPE t_regionCharges

      REAL,    ALLOCATABLE :: qis(:,:,:)

      REAL,    ALLOCATABLE :: qal(:,:,:,:)
      REAL,    ALLOCATABLE :: sqal(:,:,:)
      REAL,    ALLOCATABLE :: ener(:,:,:)

      REAL,    ALLOCATABLE :: sqlo(:,:,:)
      REAL,    ALLOCATABLE :: enerlo(:,:,:)

      REAL,    ALLOCATABLE :: qvac(:,:,:,:)
      REAL,    ALLOCATABLE :: svac(:,:)
      REAL,    ALLOCATABLE :: pvac(:,:)
      REAL,    ALLOCATABLE :: qvlay(:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: qstars(:,:,:,:)

      CONTAINS
         PROCEDURE,PASS :: init => regionCharges_init
   END TYPE t_regionCharges

PUBLIC t_regionCharges

CONTAINS

SUBROUTINE regionCharges_init(thisRegCharges,input,atoms,dimension,kpts,vacuum)

   USE m_types_setup
   USE m_types_kpts

   IMPLICIT NONE

   CLASS(t_regionCharges), INTENT(INOUT) :: thisRegCharges
   TYPE(t_input),          INTENT(IN)    :: input
   TYPE(t_atoms),          INTENT(IN)    :: atoms
   TYPE(t_dimension),      INTENT(IN)    :: dimension
   TYPE(t_kpts),           INTENT(IN)    :: kpts
   TYPE(t_vacuum),         INTENT(IN)    :: vacuum

   ALLOCATE(thisRegCharges%qis(dimension%neigd,kpts%nkpt,input%jspins))

   ALLOCATE(thisRegCharges%qal(0:3,atoms%ntype,dimension%neigd,input%jspins))
   ALLOCATE(thisRegCharges%sqal(0:3,atoms%ntype,input%jspins))
   ALLOCATE(thisRegCharges%ener(0:3,atoms%ntype,input%jspins))

   ALLOCATE(thisRegCharges%sqlo(atoms%nlod,atoms%ntype,input%jspins))
   ALLOCATE(thisRegCharges%enerlo(atoms%nlod,atoms%ntype,input%jspins))

   ALLOCATE(thisRegCharges%qvac(dimension%neigd,2,kpts%nkpt,input%jspins))
   ALLOCATE(thisRegCharges%svac(2,input%jspins))
   ALLOCATE(thisRegCharges%pvac(2,input%jspins))
   ALLOCATE(thisRegCharges%qvlay(dimension%neigd,vacuum%layerd,2,kpts%nkpt,input%jspins))
   ALLOCATE(thisRegCharges%qstars(vacuum%nstars,dimension%neigd,vacuum%layerd,2))

   thisRegCharges%qis = 0.0

   thisRegCharges%qal = 0.0
   thisRegCharges%sqal = 0.0
   thisRegCharges%ener = 0.0

   thisRegCharges%sqlo = 0.0
   thisRegCharges%enerlo = 0.0

   thisRegCharges%qvac = 0.0
   thisRegCharges%svac = 0.0
   thisRegCharges%pvac = 0.0
   thisRegCharges%qvlay = 0.0
   thisRegCharges%qstars = CMPLX(0.0,0.0)

END SUBROUTINE regionCharges_init

END MODULE m_types_regionCharges
