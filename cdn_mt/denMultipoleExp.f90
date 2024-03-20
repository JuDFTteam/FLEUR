MODULE m_denMultipoleExp

IMPLICIT NONE

CONTAINS

SUBROUTINE denMultipoleExp(input, fmpi, atoms, sphhar, stars, sym, juphon, cell,   den)

   USE m_types
   USE m_constants
   USE m_mpmom

   TYPE(t_input),  INTENT(IN) :: input
   TYPE(t_mpi),    INTENT(IN) :: fmpi
   TYPE(t_atoms),  INTENT(IN) :: atoms
   TYPE(t_sphhar), INTENT(IN) :: sphhar
   TYPE(t_stars),  INTENT(IN) :: stars
   TYPE(t_sym),    INTENT(IN) :: sym
   TYPE(t_juphon), INTENT(IN) :: juphon
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
      CALL mpmom(input,fmpi,atoms,sphhar,stars,sym,juphon,cell ,workDen%pw(1:,1),workDen%mt(:,0:,1:,1),POTDEN_TYPE_DEN,qlm,l_coreCharge=.FALSE.)
      IF(fmpi%irank.EQ.0) THEN
         WRITE(oUnit,*) '======================================='
      END IF

      IF(fmpi%irank.EQ.0) THEN
         WRITE(oUnit,*) 'Multipole expansion for spin-down density:'
         WRITE(oUnit,*) '======================================='
      END IF
      qlm = CMPLX(0.0,0.0)
      CALL mpmom(input,fmpi,atoms,sphhar,stars,sym,juphon,cell ,workDen%pw(1:,2),workDen%mt(:,0:,1:,2),POTDEN_TYPE_DEN,qlm,l_coreCharge=.FALSE.)
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
   CALL mpmom(input,fmpi,atoms,sphhar,stars,sym,juphon,cell ,workDen%pw(1:,1),workDen%mt(:,0:,1:,1),POTDEN_TYPE_DEN,qlm,l_coreCharge=.FALSE.)
   IF(fmpi%irank.EQ.0) THEN
      WRITE(oUnit,*) '======================================='
   END IF

   IF(input%jspins == 2) THEN
      IF(fmpi%irank.EQ.0) THEN
         WRITE(oUnit,*) 'Multipole expansion for magnetization density:'
         WRITE(oUnit,*) '======================================='
      END IF
      qlm = CMPLX(0.0,0.0)
      CALL mpmom(input,fmpi,atoms,sphhar,stars,sym,juphon,cell ,workDen%pw(1:,2),workDen%mt(:,0:,1:,2),POTDEN_TYPE_DEN,qlm,l_coreCharge=.FALSE.)
      IF(fmpi%irank.EQ.0) THEN
         WRITE(oUnit,*) '======================================='
      END IF
   END IF

END SUBROUTINE denMultipoleExp

END MODULE m_denMultipoleExp
