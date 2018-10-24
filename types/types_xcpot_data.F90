!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_xcpot_data
  !This module contains the xcpot-type used for the in-build xc-implementations
  IMPLICIT NONE

  TYPE t_xcpot_data
     !in the pbe case (exchpbe.F) lots of test are made
     !in addition some constants are set
     !to speed up this code precalculate things in init
     LOGICAL             :: is_rpbe !Rpbe
     LOGICAL             :: is_wc
     LOGICAL             :: is_hse !hse,lhse,vhse
     REAL                :: uk,um
     !many logicals to determine xcpot
     LOGICAL             :: is_pbes !is pbe-sol
     LOGICAL             :: is_pbe0 
     LOGICAL             :: is_bh 
     LOGICAL             :: is_mjw 
     REAL                :: exchange_weight
     INTEGER             :: krla !relativistic corrections
  END TYPE t_xcpot_data
END MODULE m_types_xcpot_data
