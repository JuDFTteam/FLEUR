!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_usdus
  TYPE t_usdus
     REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: us
     REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: dus
     REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: uds
     REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: duds !(0:lmaxd,ntype,jspd)
     REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: ddn  !(0:lmaxd,ntype,jspd)
     REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: ulos
     REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: dulos
     REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: uulon
     REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: dulon     ! (nlod,ntype,jspd)
     REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: uloulopn  ! (nlod,nlod,ntypd,jspd)
     REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: uuilon
     REAL,ALLOCATABLE,DIMENSION(:,:,:)   :: duilon    ! (nlod,ntype,jspd)
     REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: ulouilopn ! (nlod,nlod,ntypd,jspd)
   CONTAINS
     PROCEDURE :: init => usdus_init
  END TYPE t_usdus


 
CONTAINS
  SUBROUTINE usdus_init(ud,atoms,jsp)
    USE m_judft
    USE m_types_setup
    IMPLICIT NONE
    CLASS(t_usdus)           :: ud
    TYPE(t_atoms),INTENT(IN) :: atoms
    INTEGER,INTENT(IN)       :: jsp

    INTEGER :: err(13)
    ALLOCATE ( ud%uloulopn(atoms%nlod,atoms%nlod,atoms%ntype,jsp),stat=err(1) )
    ALLOCATE ( ud%ddn(0:atoms%lmaxd,atoms%ntype,jsp),stat=err(2) )
    ALLOCATE ( ud%us(0:atoms%lmaxd,atoms%ntype,jsp),stat=err(3))
    ALLOCATE ( ud%uds(0:atoms%lmaxd,atoms%ntype,jsp),stat=err(4) )
    ALLOCATE ( ud%dus(0:atoms%lmaxd,atoms%ntype,jsp),stat=err(5))
    ALLOCATE ( ud%duds(0:atoms%lmaxd,atoms%ntype,jsp),stat=err(6))
    ALLOCATE ( ud%ulos(atoms%nlod,atoms%ntype,jsp ),stat=err(7))
    ALLOCATE (ud%dulos(atoms%nlod,atoms%ntype,jsp ),stat=err(8) )
    ALLOCATE (ud%uulon(atoms%nlod,atoms%ntype,jsp ),stat=err(9))
    ALLOCATE (ud%dulon(atoms%nlod,atoms%ntype,jsp) ,stat=err(10))
    ALLOCATE (ud%uuilon(atoms%nlod,atoms%ntype,jsp),stat=err(11))
    ALLOCATE (ud%duilon(atoms%nlod,atoms%ntype,jsp),stat=err(12))
    ALLOCATE (ud%ulouilopn(atoms%nlod,atoms%nlod,atoms%ntype,jsp),stat=err(13))

    IF (ANY(err>0)) CALL judft_error("Not enough memory allocating usdus datatype")

  END SUBROUTINE usdus_init
  
 
END MODULE m_types_usdus
