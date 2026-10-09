!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
PROGRAM test
  USE xc_f03_lib_m
  TYPE(xc_f03_func_t) :: xc_func
  CALL xc_f03_func_init(xc_func, 2, XC_UNPOLARIZED)	  
END PROGRAM test
