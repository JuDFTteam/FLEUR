!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_hyperfine
   implicit none

   PRIVATE

   PUBLIC :: t_hyperfine

   TYPE t_hyperfine

      LOGICAL :: l_hyperfine

      CONTAINS

      PROCEDURE, PASS :: init => hyperfineInit
      PROCEDURE, PASS :: printValenceHyperfine => hyperfinePrintValence

      PROCEDURE, PASS :: calcPrintIsomerShifts => hyperfineCalcPrintIsomerShifts

   END TYPE t_hyperfine

   CONTAINS

   SUBROUTINE hyperfineInit(thisHF, input, atoms)

      USE m_types

      CLASS(t_hyperfine),  INTENT(INOUT) :: thisHF
      TYPE(t_input),       INTENT(IN)    :: input
      TYPE(t_atoms),       INTENT(IN)    :: atoms

      thisHF%l_hyperfine = (input%kcrel.EQ.1).AND.(input%jspins.EQ.2)

   END SUBROUTINE hyperfineInit

   SUBROUTINE hyperfinePrintValence(thisHF, input, atoms, fmpi, moments)
      USE m_constants
      USE m_types

      CLASS(t_hyperfine),  INTENT(INOUT) :: thisHF
      TYPE(t_input),       INTENT(IN)    :: input
      TYPE(t_atoms),       INTENT(IN)    :: atoms
      TYPE(t_mpi),         INTENT(IN)    :: fmpi
      TYPE(t_moments),     INTENT(INOUT) :: moments

      INTEGER :: iType

      REAL    :: a0, e0, cautog, bohrMagInCGS
      REAL    :: hyperfineResults(-1:3), hyperfineResultsTotal(-1:3)

      IF ((fmpi%irank.EQ.0).AND.thisHF%l_hyperfine) THEN
         ! Print out valence contributions to the hyperfine field
         a0 = bohr_to_angstrom_const * 1.0e-8
         e0 = 1.6021892e-19 * 2.997930e+09
         cautog = e0 / (a0*a0)
         bohrMagInCGS = 1.0/(2.0*c_light(1.0))
         WRITE(oUnit,*) ''
         WRITE(ounit,*) ' Hyperfine field valence contributions in kG '
         WRITE(ounit,*) ' ========================================================== '
         WRITE(ounit,*) ' atom type                          contribution'
         WRITE(ounit,*) '                total         s           p           d           f'
         DO iType = 1, atoms%ntype
            moments%hypFineContribs(:,iType,1,1) = moments%hypFineContribs(:,iType,1,1) - moments%hypFineContribs(:,iType,2,1)
            hyperfineResults(:) = moments%hypFineContribs(:,iType,1,1) * cautog * 0.001 * sfp_const * bohrMagInCGS * 8.0 * pi_const / 3.0
            WRITE(oUnit,'(i7,5x,5f12.5,5x,a)') iType, hyperfineResults(-1:3), 'contact term'
            hyperfineResultsTotal(:) = hyperfineResults(:)
            moments%hypFineContribs(:,iType,1,2) = moments%hypFineContribs(:,iType,1,2) - moments%hypFineContribs(:,iType,2,2)
            hyperfineResults(:) = moments%hypFineContribs(:,iType,1,2) * cautog * 0.001 * bohrMagInCGS
            WRITE(oUnit,'(i7,5x,5f12.5,5x,a)') iType, hyperfineResults(-1:3), 'dipolar term'
            hyperfineResultsTotal(:) = hyperfineResultsTotal(:) + hyperfineResults(:)
            moments%hypFineContribs(:,iType,1,3) = moments%hypFineContribs(:,iType,1,3) + moments%hypFineContribs(:,iType,2,3)
            hyperfineResults(:) = moments%hypFineContribs(:,iType,1,3) * cautog * 0.001 / c_light(1.0)
            WRITE(oUnit,'(i7,5x,5f12.5,5x,a)') iType, hyperfineResults(-1:3), 'orbital term'
            hyperfineResultsTotal(:) = hyperfineResultsTotal(:) + hyperfineResults(:)
            WRITE(oUnit,'(i7,5x,5f12.5,5x,a)') iType, hyperfineResultsTotal(-1:3), 'all terms'
         END DO
         WRITE(ounit,*) ' ========================================================== '
      END IF

   END SUBROUTINE hyperfinePrintValence

   SUBROUTINE hyperfineCalcPrintIsomerShifts(thisHF, input, atoms, fmpi, den)

      USE m_constants
      USE m_types
      USE m_intgr, ONLY : intgr2

      CLASS(t_hyperfine),  INTENT(INOUT) :: thisHF
      TYPE(t_input),       INTENT(IN)    :: input
      TYPE(t_atoms),       INTENT(IN)    :: atoms
      TYPE(t_mpi),         INTENT(IN)    :: fmpi
      TYPE(t_potden),      INTENT(IN)    :: den

      INTEGER :: iType, iRad

      REAL    :: nucRad, alpha, contactDenVal, smallRadDenVal, largeRadDenVal, averageDen
      REAL    :: indefInteg(atoms%jmtd), sphrDen(atoms%jmtd)

      LOGICAL :: isSmaller

      IF ((fmpi%irank.EQ.0) ) THEN !.AND.thisHF%l_hyperfine) THEN
         WRITE(oUnit,*) ''
         WRITE(ounit,*) ' Isomer shift output (contact charge density) '
         WRITE(ounit,*) ' ============================================================== '
         WRITE(ounit,*) ' atom type     radius of nucleus (Bohr)   smallest radial mesh point (Bohr)   charge density value (e/Bohr^3) '
         DO iType = 1, atoms%ntype
            indefInteg = 0.0
            nucRad = r0_const*(atomicMasses_const(atoms%nz(iType))**(1.0/3.0))
            IF (atoms%rmsh(1,iType).LT.nucRad) THEN
               iRad = 0
               isSmaller = .TRUE.
               DO WHILE (isSmaller)
                  iRad = iRad + 1
                  IF (atoms%rmsh(iRad,iType).GE.nucRad) isSmaller = .FALSE.
               END DO
               ! Poor man's solution: linear interpolation. :)
               alpha = (nucRad - atoms%rmsh(iRad-1,iType)) / (atoms%rmsh(iRad,iType) - atoms%rmsh(iRad-1,iType))
               IF (input%jspins.EQ.1) THEN
                  smallRadDenVal = den%mt(iRad-1,0,iType,1) / ((atoms%rmsh(iRad-1,iType)**2.0)*sfp_const)
                  largeRadDenVal = den%mt(iRad,0,iType,1) / ((atoms%rmsh(iRad,iType)**2.0)*sfp_const)
               ELSE
                  smallRadDenVal = (den%mt(iRad-1,0,iType,1) + den%mt(iRad-1,0,iType,2)) / ((atoms%rmsh(iRad-1,iType)**2.0)*sfp_const)
                  largeRadDenVal = (den%mt(iRad,0,iType,1) + den%mt(iRad,0,iType,2)) / ((atoms%rmsh(iRad,iType)**2.0)*sfp_const)
               END IF
               contactDenVal = (1.0-alpha)*smallRadDenVal + alpha*largeRadDenVal

               sphrDen = 0.0
               IF (input%jspins.EQ.1) THEN
                  sphrDen(:) = den%mt(:,0,iType,1)
               ELSE
                  sphrDen(:) = den%mt(:,0,iType,1) + den%mt(:,0,iType,2)
               END IF
               indefInteg = 0.0
               CALL intgr2(sphrDen(:),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),indefInteg)

               ! Again linear interpolation. :)
               smallRadDenVal = indefInteg(iRad-1)*sfp_const / ((4.0/3.0)*pi_const*(atoms%rmsh(iRad-1,iType)**3.0))
               largeRadDenVal = indefInteg(iRad)*sfp_const / ((4.0/3.0)*pi_const*(atoms%rmsh(iRad,iType)**3.0))
               averageDen = (1.0-alpha)*smallRadDenVal + alpha*largeRadDenVal
               averageDen = averageDen

               WRITE(oUnit,'(i7,10x,f15.8,15x,f15.8,20x,f17.8,5x,a)') iType, nucRad, atoms%rmsh(1,iType), contactDenVal, 'density at radius of nuncleus'
               WRITE(oUnit,'(i7,10x,15x,50x,f17.8,5x,a)') iType, averageDen, 'average density over nuncleus'
            ELSE
               WRITE(oUnit,'(i7,10x,f15.8,20x,a)') iType, nucRad, 'smallest radial mesh point outside radius of nucleus'
            END IF
         END DO
         WRITE(ounit,*) ' ============================================================== '
      END IF


   END SUBROUTINE hyperfineCalcPrintIsomerShifts

END MODULE m_types_hyperfine
