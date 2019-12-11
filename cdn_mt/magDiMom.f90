! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_magDiMom

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine calculates intraatomic magnetic dipole moments.
!
!                                           GM'2018
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE magDiMom(input,atoms,sphhar,noco,l_fmpl2,rho,magDipoles,elecDipoles)

   USE m_constants
   USE m_types
   USE m_juDFT
   USE m_rotdenmat
   USE m_lattHarmsSphHarmsConv
   USE m_gaunt
   USE m_intgr

   IMPLICIT NONE

   TYPE(t_input),         INTENT(IN)    :: input
   TYPE(t_atoms),         INTENT(IN)    :: atoms
   TYPE(t_sphhar),        INTENT(IN)    :: sphhar
   TYPE(t_noco),          INTENT(IN)    :: noco
   REAL,                  INTENT(IN)    :: rho(:,0:,:,:) ! if l_fmpl last dimension is 4, otherwise 2.

   LOGICAL,               INTENT(IN)    :: l_fmpl2
   REAL,                  INTENT(INOUT) :: magDipoles(:,:)
   REAL,                  INTENT(INOUT) :: elecDipoles(:,:)


   REAL,    ALLOCATABLE :: inRho(:,:,:,:)
   COMPLEX, ALLOCATABLE :: rhoSphHarms(:,:,:,:), rhoTempSphHarms(:,:,:,:)
   COMPLEX, ALLOCATABLE :: rhoSphHarmsR(:,:,:)

   INTEGER :: iType, ilh, l, m, lm, lp, mp, lmp, i
   REAL    :: theta, phi, cdn11, cdn22
   REAL    :: magDipole(3), myCharge, elecDipole(3)
   REAL    :: constA
   COMPLEX :: gauntA, gauntB, gauntC
   COMPLEX :: cdn21

   IF(input%jspins.EQ.1) RETURN
   IF(.NOT.noco%l_noco) RETURN

   !---> calculate the charge and magnetization density in the muffin tins
   ALLOCATE(inRho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,4))
   DO iType = 1,atoms%ntype
      IF (.NOT.l_fmpl2) THEN
         theta = noco%beta(iType)
         phi   = noco%alph(iType)
         inRho(:,:,iType,1) = rho(:,:,iType,1) + rho(:,:,iType,2)
         inRho(:,:,iType,2) = rho(:,:,iType,1) - rho(:,:,iType,2)
         inRho(:,:,iType,3) = inRho(:,:,iType,2) * SIN(phi)*SIN(theta)
         inRho(:,:,iType,4) = inRho(:,:,iType,2) * COS(theta)
         inRho(:,:,iType,2) = inRho(:,:,iType,2) * COS(phi)*SIN(theta)
      ELSE
         DO ilh = 0,sphhar%nlh(atoms%ntypsy(iType))
            DO i = 1,atoms%jri(iType)
               
               cdn11 = rho(i,ilh,iType,1)
               cdn22 = rho(i,ilh,iType,2)
               cdn21 = CMPLX(rho(i,ilh,iType,3),rho(i,ilh,iType,4))
               CALL rot_den_mat(noco%alph(iType),noco%beta(iType),cdn11,cdn22,cdn21)
               inRho(i,ilh,iType,1) = cdn11 + cdn22
               inRho(i,ilh,iType,2) = 2.0 * REAL(cdn21)
               ! Note: The minus sign in the following line is temporary to adjust for differences in the offdiagonal
               !       part of the density between this fleur version and ancient (v0.26) fleur.
               inRho(i,ilh,iType,3) = -2.0 * AIMAG(cdn21)
               inRho(i,ilh,iType,4) = cdn11 - cdn22
            END DO
         END DO
      END IF
   END DO

   ALLOCATE (rhoSphHarms(atoms%jmtd,(atoms%lmaxd+1)**2,atoms%ntype,4))
   ALLOCATE (rhoTempSphHarms(atoms%jmtd,(atoms%lmaxd+1)**2,atoms%ntype,4))
   ALLOCATE (rhoSphHarmsR(atoms%jmtd,(atoms%lmaxd+1)**2,atoms%ntype))

   rhoSphHarms = CMPLX(0.0,0.0)
   DO i = 1, 4
      DO iType = 1, atoms%ntype
         CALL lattHarmsRepToSphHarms(atoms,sphhar,iType,inRho(:,0:,iType,i),rhoSphHarms(:,:,iType,i))
      END DO
   END DO

   ! electric dipole moment (start)
   DO iType = 1, atoms%ntype
      DO i = 1, atoms%jri(iType)
         rhoSphHarmsR(i,:,iType) = rhoSphHarms(i,:,iType,1) * atoms%rmsh(i,iType)
      END DO
   END DO

   constA = SQRT(2.0*pi_const/3.0)
   rhoTempSphHarms = CMPLX(0.0,0.0)
   DO iType = 1, atoms%ntype
      DO lp = 0, MIN(2,atoms%lmax(iType))
         DO mp = -lp, lp
            DO l = 0, MIN(2,atoms%lmax(iType))
               DO m = -l, l
                  !note 1: For refinement maybe I could make use of selection rules.
                  !note 2: ls for r^\hat is 1, ms is -1..1.
                  gauntA = gaunt1(lp,1,l,mp,-1,m,atoms%lmaxd)
                  gauntB = gaunt1(lp,1,l,mp,0,m,atoms%lmaxd)
                  gauntC = gaunt1(lp,1,l,mp,1,m,atoms%lmaxd)
                  lmp = lp*(lp+1) + mp + 1
                  lm = l*(l+1) + m + 1
                  DO i = 1, atoms%jri(iType)
                     rhoTempSphHarms(i,lmp,iType,2) = rhoTempSphHarms(i,lmp,iType,2) +&
                        constA*(gauntA-gauntC)*rhoSphHarmsR(i,lm,iType)
                     rhoTempSphHarms(i,lmp,iType,3) = rhoTempSphHarms(i,lmp,iType,3) +&
                        constA*CMPLX(0.0,1.0)*(gauntA+gauntC)*rhoSphHarmsR(i,lm,iType)
                     rhoTempSphHarms(i,lmp,iType,4) = rhoTempSphHarms(i,lmp,iType,4) +&
                        constA*SQRT(2.0)*gauntB*rhoSphHarmsR(i,lm,iType)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO

   elecDipole = 0.0
   DO iType = 1, atoms%ntype
      CALL intgr3(REAL(rhoTempSphHarms(:,1,iType,2)),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),elecDipole(1))
      CALL intgr3(REAL(rhoTempSphHarms(:,1,iType,3)),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),elecDipole(2))
      CALL intgr3(REAL(rhoTempSphHarms(:,1,iType,4)),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),elecDipole(3))
      elecDipoles(:,iType) = elecDipole(:) * sfp_const
   END DO

   ! electric dipole moment (end)

!   WRITE(7534,*) '===================================================='
   DO iType = 1, atoms%ntype
      DO l = 0, 2
         DO m = -l, l
         magDipole(:) = 0.0
         myCharge = 0.0
         lm = l*(l+1) + m + 1
         CALL intgr3(REAL(rhoSphHarms(:,lm,iType,1)),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),myCharge)
         CALL intgr3(REAL(rhoSphHarms(:,lm,iType,2)),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),magDipole(1))
         CALL intgr3(REAL(rhoSphHarms(:,lm,iType,3)),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),magDipole(2))
         CALL intgr3(REAL(rhoSphHarms(:,lm,iType,4)),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),magDipole(3))
!         WRITE(7534,'(3i7,4e24.8)') iType, l, m, myCharge, magDipole(:)
         END DO
      END DO
      rhoSphHarms(:,:,iType,1) = CMPLX(0.0,0.0)
   END DO

   constA = SQRT(2.0*pi_const/3.0)
   DO iType = 1, atoms%ntype
      DO lp = 0, MIN(2,atoms%lmax(iType))
         DO mp = -lp, lp
            DO l = 0, MIN(2,atoms%lmax(iType))
               DO m = -l, l
                  !note 1: For refinement maybe I could make use of selection rules.
                  !note 2: ls for r^\hat is 1, ms is -1..1.
                  gauntA = gaunt1(lp,1,l,mp,-1,m,atoms%lmaxd)
                  gauntB = gaunt1(lp,1,l,mp,0,m,atoms%lmaxd)
                  gauntC = gaunt1(lp,1,l,mp,1,m,atoms%lmaxd)
                  lmp = lp*(lp+1) + mp + 1
                  lm = l*(l+1) + m + 1
                  DO i = 1, atoms%jri(iType)
                     rhoSphHarms(i,lmp,iType,1) = rhoSphHarms(i,lmp,iType,1) +&
                        constA*(gauntA-gauntC)*rhoSphHarms(i,lm,iType,2)
                     rhoSphHarms(i,lmp,iType,1) = rhoSphHarms(i,lmp,iType,1) +&
                        constA*CMPLX(0.0,1.0)*(gauntA+gauntC)*rhoSphHarms(i,lm,iType,3)
                     rhoSphHarms(i,lmp,iType,1) = rhoSphHarms(i,lmp,iType,1) +&
                        constA*SQRT(2.0)*gauntB*rhoSphHarms(i,lm,iType,4)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO

   rhoTempSphHarms = CMPLX(0.0,0.0)
   DO iType = 1, atoms%ntype
      DO lp = 0, MIN(2,atoms%lmax(iType))
         DO mp = -lp, lp
            DO l = 0, MIN(2,atoms%lmax(iType))
               DO m = -l, l
                  !note 1: For refinement maybe I could make use of selection rules.
                  !note 2: ls for r^\hat is 1, ms is -1..1.
                  gauntA = gaunt1(lp,1,l,mp,-1,m,atoms%lmaxd)
                  gauntB = gaunt1(lp,1,l,mp,0,m,atoms%lmaxd)
                  gauntC = gaunt1(lp,1,l,mp,1,m,atoms%lmaxd)
                  lmp = lp*(lp+1) + mp + 1
                  lm = l*(l+1) + m + 1
                  DO i = 1, atoms%jri(iType)
                     rhoTempSphHarms(i,lmp,iType,2) = rhoTempSphHarms(i,lmp,iType,2) +&
                        constA*(gauntA-gauntC)*rhoSphHarms(i,lm,iType,1)
                     rhoTempSphHarms(i,lmp,iType,3) = rhoTempSphHarms(i,lmp,iType,3) +&
                        constA*CMPLX(0.0,1.0)*(gauntA+gauntC)*rhoSphHarms(i,lm,iType,1)
                     rhoTempSphHarms(i,lmp,iType,4) = rhoTempSphHarms(i,lmp,iType,4) +&
                        constA*SQRT(2.0)*gauntB*rhoSphHarms(i,lm,iType,1)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO

   DO iType = 1, atoms%ntype
      DO l = 0, atoms%lmax(iType)
         DO m = -l, l
            lm = l*(l+1) + m + 1
            DO i = 1, atoms%jri(iType)
               rhoSphHarms(i,lm,iType,2) = rhoSphHarms(i,lm,iType,2) - 3.0 * rhoTempSphHarms(i,lm,iType,2)
               rhoSphHarms(i,lm,iType,3) = rhoSphHarms(i,lm,iType,3) - 3.0 * rhoTempSphHarms(i,lm,iType,3)
               rhoSphHarms(i,lm,iType,4) = rhoSphHarms(i,lm,iType,4) - 3.0 * rhoTempSphHarms(i,lm,iType,4)
            END DO
         END DO
      END DO
   END DO

   DO iType = 1, atoms%ntype
      DO i = 1, atoms%jri(iType)
         IF (ANY(AIMAG(rhoSphHarms(i,1,iType,:)).GT.1.0e-11)) THEN
            WRITE(6,*) 'imaginary part too large!'
         END IF
      END DO
   END DO

   magDipole = 0.0
   DO iType = 1, atoms%ntype
      CALL intgr3(REAL(rhoSphHarms(:,1,iType,2)),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),magDipole(1))
      CALL intgr3(REAL(rhoSphHarms(:,1,iType,3)),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),magDipole(2))
      CALL intgr3(REAL(rhoSphHarms(:,1,iType,4)),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),magDipole(3))
      magDipoles(:,iType) = magDipole(:) * sfp_const
   END DO

END SUBROUTINE magDiMom

END MODULE m_magDiMom
