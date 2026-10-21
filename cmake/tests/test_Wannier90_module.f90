!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
PROGRAM test_wannier90_module
  USE w90_library
  IMPLICIT NONE

  TYPE(lib_common_type) :: w90main

  ! Call a library routine so the test also needs the Wannier90 library to link
  CALL w90_set_option(w90main, 'num_wann', 1)
END PROGRAM test_wannier90_module
