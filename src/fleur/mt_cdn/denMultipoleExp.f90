!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_denMultipoleExp
   implicit none


CONTAINS

SUBROUTINE denMultipoleExp(input, fmpi, atoms, sphhar, stars, sym, cell,   den)

   USE m_types
   USE m_constants
   USE m_mpmom

   TYPE(t_input),  INTENT(IN) :: input
   TYPE(t_mpi),    INTENT(IN) :: fmpi
   TYPE(t_atoms),  INTENT(IN) :: atoms
   TYPE(t_sphhar), INTENT(IN) :: sphhar
   TYPE(t_stars),  INTENT(IN) :: stars
   TYPE(t_sym),    INTENT(IN) :: sym
   TYPE(t_cell),   INTENT(IN) :: cell
    
   TYPE(t_potden), INTENT(IN) :: den

   type(t_potden)             :: workDen
   COMPLEX                    :: qlm(-atoms%lmaxd:atoms%lmaxd,0:atoms%lmaxd,atoms%ntype)

   IF(input%jspins == 2) THEN
      IF(fmpi%irank.EQ.0) THEN
         WRITE(oUnit,*) 'Multipole expansion for spin-up density:'
         WRITE(oUnit,*) '======================================='
      END IF
      qlm = CMPLX(0.0,0.0)
      workDen = den
      CALL mpmom(input,fmpi,atoms,sphhar,stars,sym,cell ,workDen,POTDEN_TYPE_DEN,qlm,1,l_coreCharge=.FALSE.)
      IF(fmpi%irank.EQ.0) THEN
         WRITE(oUnit,*) '======================================='
      END IF

      IF(fmpi%irank.EQ.0) THEN
         WRITE(oUnit,*) 'Multipole expansion for spin-down density:'
         WRITE(oUnit,*) '======================================='
      END IF
      qlm = CMPLX(0.0,0.0)
      CALL mpmom(input,fmpi,atoms,sphhar,stars,sym,cell ,workDen,POTDEN_TYPE_DEN,qlm,2,l_coreCharge=.FALSE.)
      IF(fmpi%irank.EQ.0) THEN
         WRITE(oUnit,*) '======================================='
      END IF
   END IF

   IF(fmpi%irank.EQ.0) THEN
      WRITE(oUnit,*) 'Multipole expansion for charge density:'
      WRITE(oUnit,*) '======================================='
   END IF
   qlm = CMPLX(0.0,0.0)
   workDen = den
   IF(input%jspins == 2) CALL workDen%SpinsToChargeAndMagnetisation()
   CALL mpmom(input,fmpi,atoms,sphhar,stars,sym,cell ,workDen,POTDEN_TYPE_DEN,qlm,1,l_coreCharge=.FALSE.)
   IF(fmpi%irank.EQ.0) THEN
      WRITE(oUnit,*) '======================================='
   END IF

   IF(input%jspins == 2) THEN
      IF(fmpi%irank.EQ.0) THEN
         WRITE(oUnit,*) 'Multipole expansion for magnetization density:'
         WRITE(oUnit,*) '======================================='
      END IF
      qlm = CMPLX(0.0,0.0)
      CALL mpmom(input,fmpi,atoms,sphhar,stars,sym,cell ,workDen,POTDEN_TYPE_DEN,qlm,2,l_coreCharge=.FALSE.)
      IF(fmpi%irank.EQ.0) THEN
         WRITE(oUnit,*) '======================================='
      END IF
   END IF

END SUBROUTINE denMultipoleExp

END MODULE m_denMultipoleExp
