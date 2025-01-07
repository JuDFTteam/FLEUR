!--------------------------------------------------------------------------------
! Copyright (c) 2024 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

PROGRAM test

  use iso_c_binding
  use chase_diag
  integer :: flag
  CALL pschase_finalize (flag)

END program
