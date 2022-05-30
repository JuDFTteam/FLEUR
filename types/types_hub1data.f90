!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_hub1data

   USE m_constants
   USE m_juDFT

   IMPLICIT NONE

   PRIVATE

   TYPE t_hub1data
      !Contains the information for the hubbard 1 solver,
      !which is calculated on the fly (fixed parameters are found in types_hub1inp)
      INTEGER           :: iter=0
      INTEGER           :: overallIteration=0
      LOGICAL           :: l_runthisiter=.FALSE.   !switch which determines wether Hubbard 1 will be run in the current iteration
      LOGICAL           :: l_performSpinavg = .TRUE.


      REAL, ALLOCATABLE :: mag_mom(:,:)    !magnetic moment (for exchange splitting)
      REAL, ALLOCATABLE :: xi(:)           !Spin-orbit coupling parameter
      REAL, ALLOCATABLE :: ccfmat(:,:,:)   !crystal field splitting matrix

      REAL, ALLOCATABLE :: cdn_atomic(:,:,:,:) !atomic contribution to the charge density
                                               !is used to calculate CF coefficients in writeCFoutput

      CONTAINS

      PROCEDURE, PASS :: init   => hub1data_init
      PROCEDURE, PASS :: mpi_bc => hub1data_mpi_bc

   END TYPE t_hub1data

   PUBLIC t_hub1data

   CONTAINS

   SUBROUTINE hub1data_init(this,atoms,input,hub1inp,fmpi,mmpmatDistancePrev,occDistancePrev,l_error)

      USE m_types_mpi
      USE m_types_atoms
      USE m_types_input
      USE m_types_hub1inp
      USE m_gaunt

      CLASS(t_hub1data),   INTENT(INOUT) :: this
      TYPE(t_atoms),       INTENT(IN)    :: atoms
      TYPE(t_input),       INTENT(IN)    :: input
      TYPE(t_hub1inp),     INTENT(IN)    :: hub1inp
      TYPE(t_mpi),         INTENT(IN)    :: fmpi
      REAL,                INTENT(IN)    :: mmpmatDistancePrev,occDistancePrev
      LOGICAL,             INTENT(IN)    :: l_error

      INTEGER :: i_hia, l, m, mp, lcoeff, mcoeff
      REAL :: gaunt_coef


      this%l_performSpinavg = .FALSE.
      this%iter = 0
      this%l_runthisiter = .FALSE.
      IF(atoms%n_hia>0) THEN
         IF(fmpi%irank == 0) THEN
            this%l_performSpinavg = .NOT.hub1inp%l_dftSpinpol

            IF(.NOT.l_error.AND..NOT.hub1inp%l_forceHIAiteration) THEN
               IF(hub1inp%l_correctEtot.AND..NOT.hub1inp%l_dftSpinpol.AND.&
                  mmpmatDistancePrev<hub1inp%minmatDistance.AND.&
                  occDistancePrev<hub1inp%minoccDistance) THEN
                  !If we read converged distances it is assumed that the correction should kick in
                  WRITE(*,*) "Previous density matrix was converged"
                  WRITE(*,*) "Switching off spin averaging"
                  this%l_performSpinavg = .FALSE.
               ENDIF
            ELSE
               IF(hub1inp%l_correctEtot) THEN
                  IF(l_error) THEN
                     WRITE(*,*) "No previous density matrix distances found"
                     WRITE(*,*) "setting spin averaging according to dftSpinpol"
                  ELSE IF(hub1inp%l_forceHIAiteration) THEN
                     WRITE(*,*) "Previous density matrix distances are ignored"
                     WRITE(*,*) "setting spin averaging according to dftSpinpol"
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF


      ALLOCATE (this%mag_mom(MAX(1,atoms%n_hia),lmaxU_const-1),source=0.0)
      ALLOCATE (this%xi(MAX(1,atoms%n_hia)),source=0.0)
      DO i_hia = 1, atoms%n_hia
         IF(hub1inp%l_soc_given(i_hia)) THEN
            this%xi(i_hia) = hub1inp%xi_par(i_hia)
         ENDIF
      ENDDO

      ALLOCATE (this%cdn_atomic(atoms%jmtd,0:lmaxU_const,atoms%ntype,input%jspins),source=0.0)

      ALLOCATE (this%ccfmat(MAX(1,atoms%n_hia),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const),source=0.0)
      IF(ANY(ABS(hub1inp%ccf(:)).GT.1e-12)) THEN
         DO i_hia = 1, atoms%n_hia
            l = atoms%lda_u(atoms%n_u+i_hia)%l
            DO lcoeff = 0, 6
               DO mcoeff = -lcoeff,lcoeff
                  IF(ABS(hub1inp%cfCoeffs(i_hia,lcoeff,mcoeff)).LT.1e-12) CYCLE
                  DO m = -l,l
                     DO mp = -l,l
                        !TODO: The (-1)**mp term comes from the fact that gaunt1
                        !should calculate the integral of three spherical harmonics
                        !with the first one complex conjugated. However, the complex
                        !conjugation is only reflected in the selection rules at the top
                        !and the prefactor is omitted
                        gaunt_coef = (-1)**mp * gaunt1(l,lcoeff,l,mp,mcoeff,m,6)
                        IF(ABS(gaunt_coef).LT.1e-12) CYCLE
                        this%ccfmat(i_hia,m,mp) = &
                           this%ccfmat(i_hia,m,mp) + sfp_const/sqrt(2*lcoeff+1.0) * gaunt_coef &
                                                    * hub1inp%cfCoeffs(i_hia,lcoeff,mcoeff) &
                                                    * boltzmann_const
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF

   END SUBROUTINE hub1data_init

   SUBROUTINE hub1data_mpi_bc(this, mpi_comm, irank)
      USE m_mpi_bc_tool
      CLASS(t_hub1data), INTENT(INOUT)::this
      INTEGER, INTENT(IN):: mpi_comm
      INTEGER, INTENT(IN), OPTIONAL::irank
      INTEGER ::rank,myrank,n,ierr
      IF (PRESENT(irank)) THEN
         rank = irank
      ELSE
         rank = 0
      END IF
      CALL mpi_bc(this%iter,rank,mpi_comm)
      CALL mpi_bc(this%overallIteration,rank,mpi_comm)
      CALL mpi_bc(this%l_runthisiter,rank,mpi_comm)
      CALL mpi_bc(this%l_performSpinavg,rank,mpi_comm)
      CALL mpi_bc(this%mag_mom,rank,mpi_comm)
      CALL mpi_bc(this%xi,rank,mpi_comm)
      CALL mpi_bc(this%ccfmat,rank,mpi_comm)
      CALL mpi_bc(this%cdn_atomic,rank,mpi_comm)

   END SUBROUTINE hub1data_mpi_bc

END MODULE m_types_hub1data