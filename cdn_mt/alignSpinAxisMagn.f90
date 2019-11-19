!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and avhttps://gcc.gnu.org/onlinedocs/gfortran/SQRT.htmlailable as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!------------------------------------------------------------------------------
!  This routine 
!  
! 
!
! Robin Hilgers, Nov '19
MODULE m_alignSpinAxisMagn

USE m_magnMomFromDen
USE m_types
USE m_flipcdn

CONTAINS
SUBROUTINE rotateMagnetToSpinAxis(noco,vacuum,sphhar,stars&
,sym,oneD,cell,moments,input,atoms,den)
   TYPE(t_input), INTENT(INOUT)  :: input
   TYPE(t_atoms), INTENT(INOUT)  :: atoms
   TYPE(t_noco), INTENT(IN)      :: noco
   TYPE(t_stars),INTENT(IN)      :: stars
   TYPE(t_vacuum),INTENT(IN)     :: vacuum
   TYPE(t_sphhar),INTENT(IN)     :: sphhar
   TYPE(t_sym),INTENT(IN)        :: sym
   TYPE(t_oneD),INTENT(IN)       :: oneD
   TYPE(t_cell),INTENT(IN)       :: cell
   TYPE(t_potden), INTENT(INOUT) :: den 


   REAL                          :: moments(3)

   CALL magnMomFromDen(input,atoms,noco,den,moments)
   CALL flipcdn(atoms,input,vacuum,sphhar,stars,sym,noco,oneD,cell,atoms%flipSpinPhi,atoms%flipSpinTheta,den)

END SUBROUTINE rotateMagnetToSpinAxis


SUBROUTINE rotateMagnetFromSpinAxis(noco,vacuum,sphhar,stars&
,sym,oneD,cell,moments,input,atoms,den)
   TYPE(t_input), INTENT(INOUT)  :: input
   TYPE(t_atoms), INTENT(INOUT)  :: atoms
   TYPE(t_noco), INTENT(IN)	 :: noco
   TYPE(t_stars),INTENT(IN)	 :: stars
   TYPE(t_vacuum),INTENT(IN)     :: vacuum
   TYPE(t_sphhar),INTENT(IN)     :: sphhar
   TYPE(t_sym),INTENT(IN)        :: sym
   TYPE(t_oneD),INTENT(IN)	 :: oneD
   TYPE(t_cell),INTENT(IN)	 :: cell
   TYPE(t_potden), INTENT(INOUT) :: den 


   CALL flipcdn(atoms,input,vacuum,sphhar,stars,sym,noco,oneD,cell,-atoms%flipSpinPhi,-atoms%flipSpinTheta,den)
   atoms%flipSpinPhi=0
   atoms%flipSpinTheta=0


END SUBROUTINE rotateMagnetFromSpinAxis


END MODULE m_alignSpinAxisMagn

