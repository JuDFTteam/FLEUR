!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_regionCharges

IMPLICIT NONE

PRIVATE

   TYPE t_regionCharges

      REAL,    ALLOCATABLE :: sqal(:,:,:)
      REAL,    ALLOCATABLE :: ener(:,:,:)
      REAL,    ALLOCATABLE :: sqlo(:,:,:)
      REAL,    ALLOCATABLE :: enerlo(:,:,:)
      REAL,    ALLOCATABLE :: svac(:,:)
      REAL,    ALLOCATABLE :: pvac(:,:)

      CONTAINS
         PROCEDURE,PASS :: init => regionCharges_init
         PROCEDURE      :: sumBandsVac
   END TYPE t_regionCharges

PUBLIC t_regionCharges

CONTAINS

SUBROUTINE regionCharges_init(thisRegCharges,input,atoms)

   USE m_types_input
   USE m_types_atoms

   IMPLICIT NONE

   CLASS(t_regionCharges), INTENT(INOUT) :: thisRegCharges
   TYPE(t_input),          INTENT(IN)    :: input
   TYPE(t_atoms),          INTENT(IN)    :: atoms

   ALLOCATE(thisRegCharges%sqal(0:3,atoms%ntype,input%jspins))
   ALLOCATE(thisRegCharges%ener(0:3,atoms%ntype,input%jspins))

   ALLOCATE(thisRegCharges%sqlo(atoms%nlod,atoms%ntype,input%jspins))
   ALLOCATE(thisRegCharges%enerlo(atoms%nlod,atoms%ntype,input%jspins))
   ALLOCATE(thisRegCharges%svac(2,input%jspins))
   ALLOCATE(thisRegCharges%pvac(2,input%jspins))

   thisRegCharges%sqal = 0.0
   thisRegCharges%ener = 0.0
   thisRegCharges%sqlo = 0.0
   thisRegCharges%enerlo = 0.0
   thisRegCharges%svac = 0.0
   thisRegCharges%pvac = 0.0

END SUBROUTINE regionCharges_init

SUBROUTINE sumBandsVac(thisRegCharges,vacuum,dos,noccbd,ikpt,jsp_start,jsp_end,eig,we)

  USE m_types_vacuum
  USE m_types_dos
  
   USE m_types_dos

   IMPLICIT NONE

   CLASS(t_regionCharges), INTENT(INOUT) :: thisRegCharges
   TYPE(t_vacuum),         INTENT(IN)    :: vacuum
   TYPE(t_dos),            INTENT(IN)    :: dos
   INTEGER,                INTENT(IN)    :: noccbd
   INTEGER,                INTENT(IN)    :: ikpt
   INTEGER,                INTENT(IN)    :: jsp_start, jsp_end
   REAL,                   INTENT(IN)    :: eig(noccbd)
   REAL,                   INTENT(IN)    :: we(noccbd)

   INTEGER                               :: ispin, ivac

   ! perform Brillouin zone integration and summation over the bands in order to determine the vacuum energy parameters.
   DO ispin = jsp_start, jsp_end
      DO ivac = 1,vacuum%nvac
         thisRegCharges%pvac(ivac,ispin) = thisRegCharges%pvac(ivac,ispin)+&
            dot_product(eig(:noccbd)*dos%qvac(:noccbd,ivac,ikpt,ispin),we(:noccbd))
         thisRegCharges%svac(ivac,ispin)=thisRegCharges%svac(ivac,ispin)+&
            dot_product(dos%qvac(:noccbd,ivac,ikpt,ispin),we(:noccbd))
      END DO
   END DO

END SUBROUTINE sumBandsVac

END MODULE m_types_regionCharges
