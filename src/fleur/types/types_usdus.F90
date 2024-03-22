!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_usdus
   TYPE t_usdus
      REAL, ALLOCATABLE, DIMENSION(:, :, :)    :: us
      REAL, ALLOCATABLE, DIMENSION(:, :, :)    :: dus
      REAL, ALLOCATABLE, DIMENSION(:, :, :)    :: uds
      REAL, ALLOCATABLE, DIMENSION(:, :, :)    :: duds !(0:lmaxd,ntype,jspd)
      REAL, ALLOCATABLE, DIMENSION(:, :, :)    :: ddn  !(0:lmaxd,ntype,jspd)
      REAL, ALLOCATABLE, DIMENSION(:, :, :)   :: ulos
      REAL, ALLOCATABLE, DIMENSION(:, :, :)   :: dulos
      REAL, ALLOCATABLE, DIMENSION(:, :, :)   :: uulon
      REAL, ALLOCATABLE, DIMENSION(:, :, :)   :: dulon     ! (nlod,ntype,jspd)
      REAL, ALLOCATABLE, DIMENSION(:, :, :, :) :: uloulopn  ! (nlod,nlod,ntypd,jspd)
      REAL, ALLOCATABLE, DIMENSION(:, :, :)   :: uuilon
      REAL, ALLOCATABLE, DIMENSION(:, :, :)   :: duilon    ! (nlod,ntype,jspd)
      REAL, ALLOCATABLE, DIMENSION(:, :, :, :) :: ulouilopn ! (nlod,nlod,ntypd,jspd)
   CONTAINS
      PROCEDURE :: init => usdus_init
      PROCEDURE :: free => usdus_free
   END TYPE t_usdus

CONTAINS
   SUBROUTINE usdus_init(ud, atoms, jsp)
      USE m_judft
      USE m_types_setup
      IMPLICIT NONE
      CLASS(t_usdus)           :: ud
      TYPE(t_atoms), INTENT(IN) :: atoms
      INTEGER, INTENT(IN)       :: jsp

      INTEGER :: err(13)

      err = 0
      if(.not. allocated(ud%uloulopn)) &
         allocate(ud%uloulopn(atoms%nlod, atoms%nlod, atoms%ntype, jsp), stat=err(1))
      if(.not. allocated(ud%ddn)) &
         allocate(ud%ddn(0:atoms%lmaxd, atoms%ntype, jsp), stat=err(2))
      if(.not. allocated(ud%us)) &
         allocate(ud%us(0:atoms%lmaxd, atoms%ntype, jsp), stat=err(3))
      if(.not. allocated(ud%uds)) &
         allocate(ud%uds(0:atoms%lmaxd, atoms%ntype, jsp), stat=err(4))
      if(.not. allocated(ud%dus)) &
         allocate(ud%dus(0:atoms%lmaxd, atoms%ntype, jsp), stat=err(5))
      if(.not. allocated(ud%duds)) &
         allocate(ud%duds(0:atoms%lmaxd, atoms%ntype, jsp), stat=err(6))
      if(.not. allocated(ud%ulos)) &
         allocate(ud%ulos(atoms%nlod, atoms%ntype, jsp), stat=err(7))
      if(.not. allocated(ud%dulos)) &
         allocate(ud%dulos(atoms%nlod, atoms%ntype, jsp), stat=err(8))
      if(.not. allocated(ud%uulon)) &
         allocate(ud%uulon(atoms%nlod, atoms%ntype, jsp), stat=err(9))
      if(.not. allocated(ud%dulon)) &
         allocate(ud%dulon(atoms%nlod, atoms%ntype, jsp), stat=err(10))
      if(.not. allocated(ud%uuilon)) &
         allocate(ud%uuilon(atoms%nlod, atoms%ntype, jsp), stat=err(11))
      if(.not. allocated(ud%duilon)) &
         allocate(ud%duilon(atoms%nlod, atoms%ntype, jsp), stat=err(12))
      if(.not. allocated(ud%ulouilopn)) &
         allocate(ud%ulouilopn(atoms%nlod, atoms%nlod, atoms%ntype, jsp), stat=err(13))

      !write (*,*) "err array", err
      IF(ANY(err > 0)) CALL judft_error("Not enough memory allocating usdus datatype")

      ud%uloulopn  = 0; ud%ddn       = 0; ud%us        = 0
      ud%uds       = 0; ud%dus       = 0; ud%duds      = 0
      ud%ulos      = 0; ud%dulos     = 0; ud%uulon     = 0
      ud%dulon     = 0; ud%uuilon    = 0; ud%duilon    = 0
      ud%ulouilopn = 0
   END SUBROUTINE usdus_init

   SUBROUTINE usdus_free(ud)
      IMPLICIT NONE
      CLASS(t_usdus)           :: ud

      if(allocated(ud%uloulopn)) deallocate(ud%uloulopn)
      if(allocated(ud%ddn)) deallocate(ud%ddn)
      if(allocated(ud%us)) deallocate(ud%us)
      if(allocated(ud%uds)) deallocate(ud%uds)
      if(allocated(ud%dus)) deallocate(ud%dus)
      if(allocated(ud%duds)) deallocate(ud%duds)
      if(allocated(ud%ulos)) deallocate(ud%ulos)
      if(allocated(ud%dulos)) deallocate(ud%dulos)
      if(allocated(ud%uulon)) deallocate(ud%uulon)
      if(allocated(ud%dulon)) deallocate(ud%dulon)
      if(allocated(ud%uuilon)) deallocate(ud%uuilon)
      if(allocated(ud%duilon)) deallocate(ud%duilon)
      if(allocated(ud%ulouilopn)) deallocate(ud%ulouilopn)

   END SUBROUTINE usdus_free

END MODULE m_types_usdus
