!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_input

     TYPE t_input
      LOGICAL :: film=.false.
      INTEGER:: jspins=1
      LOGICAL:: total
      REAL :: rkmax
      REAL :: zelec
      LOGICAL :: strho =.FALSE.
      LOGICAL :: cdinf =.false.
      LOGICAL :: vchk =.false.
      LOGICAL :: l_f =.false.
      LOGICAL :: eonly =.false.
      LOGICAL :: ctail =.true.
      INTEGER :: coretail_lmax =0 
      INTEGER :: itmax =9 
      REAL    :: minDistance=1.0e-5
      INTEGER :: maxiter=99
      INTEGER :: imix=7
      INTEGER :: gw=0
      INTEGER :: gw_neigd=0
      INTEGER :: qfix=0
      REAL    :: forcealpha =1.0 !< mixing parameter for geometry optimzer
      REAL    :: epsdisp =0.00001!< minimal displacement. If all displacements are < epsdisp stop
      REAL    :: epsforce =0.00001!< minimal force. If all forces <epsforce stop
      REAL    :: force_converged=0.00001
      INTEGER :: forcemix=3
      REAL    :: delgau =0.001  !TODO = tkb?
      REAL    :: alpha=0.05
      REAL    :: preconditioning_param=0.0
      REAL    :: spinf=2.0
      REAL    :: tkb=0.001
      LOGICAL :: gauss=.false.
      LOGICAL :: l_bmt=.false.
      !INTEGER:: scale
      INTEGER:: kcrel =0
      LOGICAL:: frcor =.false. !frozen core
      LOGICAL:: lflip=.false.
      LOGICAL:: score=.false.
      LOGICAL:: swsp=.false.
      LOGICAL:: tria=.false.
      LOGICAL:: integ=.false.
      LOGICAL:: pallst=.false.
      LOGICAL:: l_coreSpec=.false.
      LOGICAL:: l_wann=.false.
      LOGICAL:: secvar=.false.
      LOGICAL:: evonly=.false.
      LOGICAL:: l_inpXML=.true.
      REAL :: scaleCell=1.0
      REAL :: scaleA1=1.0
      REAL :: scaleA2=1.0
      REAL :: scaleC=1.0
      REAL :: ellow=-1.8
      REAL :: elup=1.0
      REAL :: fixed_moment = 0.0
      CHARACTER(LEN=100) :: comment="FLEUR calculation without a title"
      REAL, POINTER :: sigma !this is the difference in charge due to the electric field it points to the value stored in t_efield
      LOGICAL :: l_core_confpot=.true. !Former CPP_CORE
      LOGICAL :: l_useapw=.false.
      LOGICAL :: ldauLinMix=.false.
      REAL    :: ldauMixParam=0.1
      REAL    :: ldauSpinf=2.0
      LOGICAL :: l_rdmft=.false.
      REAL    :: rdmftOccEps=0.0
      INTEGER :: rdmftStatesBelow=0
      INTEGER :: rdmftStatesAbove=0
      INTEGER :: rdmftFunctional=0
   END TYPE t_input

END MODULE m_types_input
  
