!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_hub1data

   USE m_constants
   USE m_types_setup

   IMPLICIT NONE

   PRIVATE

   TYPE t_hub1data
      !Contains the information for the hubbard 1 solver,
      !which is calculated on the fly (fixed parameters are found in types_hub1inp)
      INTEGER           :: iter=0
      LOGICAL           :: l_runthisiter=.FALSE.   !switch which determines wether Hubbard 1 will be run in the current iteration
      LOGICAL           :: l_performSpinavg = .TRUE.


      REAL, ALLOCATABLE :: mag_mom(:,:)    !magnetic moment (for exchange splitting)
      REAL, ALLOCATABLE :: xi(:)           !Spin-orbit coupling parameter
      REAL, ALLOCATABLE :: ccfmat(:,:,:)   !crystal field splitting matrix

      REAL, ALLOCATABLE :: cdn_spherical(:,:,:)

      CONTAINS

      PROCEDURE, PASS :: init => hub1data_init

   END TYPE t_hub1data

   PUBLIC t_hub1data

   CONTAINS

   SUBROUTINE hub1data_init(this,atoms,hub1inp)

      CLASS(t_hub1data),   INTENT(INOUT) :: this
      TYPE(t_atoms),       INTENT(IN)    :: atoms
      TYPE(t_hub1inp),     INTENT(IN)    :: hub1inp

      INTEGER :: i_hia

      this%l_performSpinavg = .NOT.hub1inp%l_dftSpinpol

      ALLOCATE (this%mag_mom(4*atoms%ntype,lmaxU_const-1),source=0.0)
      ALLOCATE (this%xi(4*atoms%ntype),source=0.0)
      ALLOCATE (this%ccfmat(4*atoms%ntype,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const),source=0.0)
      ALLOCATE (this%cdn_spherical(atoms%jmtd,0:lmaxU_const,atoms%ntype),source=0.0)

      DO i_hia = 1, atoms%n_hia
         IF(hub1inp%l_soc_given(i_hia)) THEN
            this%xi(i_hia) = hub1inp%xi_par(i_hia)
         ENDIF
      ENDDO

   END SUBROUTINE hub1data_init

END MODULE m_types_hub1data