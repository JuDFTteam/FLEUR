!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_fleur_arguments
  IMPLICIT NONE
  PRIVATE
  TYPE t_param
     INTEGER             :: TYPE !can be 0,1,2 for a simple argument, an argument with a string or with a number
     CHARACTER(len=20)   :: name
     CHARACTER(len=200)  :: desc
     CHARACTER(len=200)  :: values
  END TYPE t_param

CONTAINS
 
END MODULE m_fleur_arguments
