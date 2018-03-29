!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_cdnval

IMPLICIT NONE

PRIVATE

   TYPE t_orb
      REAL, ALLOCATABLE    :: uu(:,:,:,:)
      REAL, ALLOCATABLE    :: dd(:,:,:,:)
      COMPLEX, ALLOCATABLE :: uup(:,:,:,:)
      COMPLEX, ALLOCATABLE :: uum(:,:,:,:)
      COMPLEX, ALLOCATABLE :: ddp(:,:,:,:)
      COMPLEX, ALLOCATABLE :: ddm(:,:,:,:)

      REAL, ALLOCATABLE    :: uulo(:,:,:,:)
      REAL, ALLOCATABLE    :: dulo(:,:,:,:)
      COMPLEX, ALLOCATABLE :: uulop(:,:,:,:)
      COMPLEX, ALLOCATABLE :: uulom(:,:,:,:)
      COMPLEX, ALLOCATABLE :: dulop(:,:,:,:)
      COMPLEX, ALLOCATABLE :: dulom(:,:,:,:)

      REAL, ALLOCATABLE    :: z(:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: p(:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: m(:,:,:,:,:)

      CONTAINS
         PROCEDURE,PASS :: init => orb_init
   END TYPE t_orb

PUBLIC t_orb

CONTAINS

SUBROUTINE orb_init(thisOrb, atoms, noco, jsp_start, jsp_end)

   USE m_types_setup

   IMPLICIT NONE

   CLASS(t_orb), INTENT(INOUT)    :: thisOrb
   TYPE(t_atoms), INTENT(IN)      :: atoms
   TYPE(t_noco), INTENT(IN)       :: noco
   INTEGER, INTENT(IN)            :: jsp_start
   INTEGER, INTENT(IN)            :: jsp_end

   INTEGER                        :: dim1, dim2, dim3

   IF(ALLOCATED(thisOrb%uu)) DEALLOCATE(thisOrb%uu)
   IF(ALLOCATED(thisOrb%dd)) DEALLOCATE(thisOrb%dd)
   IF(ALLOCATED(thisOrb%uup)) DEALLOCATE(thisOrb%uup)
   IF(ALLOCATED(thisOrb%uum)) DEALLOCATE(thisOrb%uum)
   IF(ALLOCATED(thisOrb%ddp)) DEALLOCATE(thisOrb%ddp)
   IF(ALLOCATED(thisOrb%ddm)) DEALLOCATE(thisOrb%ddm)

   IF(ALLOCATED(thisOrb%uulo)) DEALLOCATE(thisOrb%uulo)
   IF(ALLOCATED(thisOrb%dulo)) DEALLOCATE(thisOrb%dulo)
   IF(ALLOCATED(thisOrb%uulop)) DEALLOCATE(thisOrb%uulop)
   IF(ALLOCATED(thisOrb%uulom)) DEALLOCATE(thisOrb%uulom)
   IF(ALLOCATED(thisOrb%dulop)) DEALLOCATE(thisOrb%dulop)
   IF(ALLOCATED(thisOrb%dulom)) DEALLOCATE(thisOrb%dulom)

   IF(ALLOCATED(thisOrb%z)) DEALLOCATE(thisOrb%z)
   IF(ALLOCATED(thisOrb%p)) DEALLOCATE(thisOrb%p)
   IF(ALLOCATED(thisOrb%m)) DEALLOCATE(thisOrb%m)

   dim1 = 0
   dim2 = 1
   dim3 = 1
   IF (noco%l_soc) THEN
      dim1 = atoms%lmaxd
      dim2 = atoms%ntype
      dim3 = atoms%nlod
   END IF

   ALLOCATE(thisOrb%uu(0:dim1,-atoms%lmaxd:atoms%lmaxd,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%dd(0:dim1,-atoms%lmaxd:atoms%lmaxd,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%uup(0:dim1,-atoms%lmaxd:atoms%lmaxd,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%uum(0:dim1,-atoms%lmaxd:atoms%lmaxd,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%ddp(0:dim1,-atoms%lmaxd:atoms%lmaxd,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%ddm(0:dim1,-atoms%lmaxd:atoms%lmaxd,dim2,jsp_start:jsp_end))

   ALLOCATE(thisOrb%uulo(dim3,-atoms%llod:atoms%llod,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%dulo(dim3,-atoms%llod:atoms%llod,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%uulop(dim3,-atoms%llod:atoms%llod,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%uulom(dim3,-atoms%llod:atoms%llod,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%dulop(dim3,-atoms%llod:atoms%llod,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%dulom(dim3,-atoms%llod:atoms%llod,dim2,jsp_start:jsp_end))

   ALLOCATE(thisOrb%z(dim3,dim3,-atoms%llod:atoms%llod,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%p(dim3,dim3,-atoms%llod:atoms%llod,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%m(dim3,dim3,-atoms%llod:atoms%llod,dim2,jsp_start:jsp_end))

   thisOrb%uu = 0.0
   thisOrb%dd = 0.0
   thisOrb%uup = CMPLX(0.0,0.0)
   thisOrb%uum = CMPLX(0.0,0.0)
   thisOrb%ddp = CMPLX(0.0,0.0)
   thisOrb%ddm = CMPLX(0.0,0.0)

   thisOrb%uulo = 0.0
   thisOrb%dulo = 0.0
   thisOrb%uulop = CMPLX(0.0,0.0)
   thisOrb%uulom = CMPLX(0.0,0.0)
   thisOrb%dulop = CMPLX(0.0,0.0)
   thisOrb%dulom = CMPLX(0.0,0.0)

   thisOrb%z = 0.0
   thisOrb%p = CMPLX(0.0,0.0)
   thisOrb%m = CMPLX(0.0,0.0)

END SUBROUTINE orb_init

END MODULE m_types_cdnval
