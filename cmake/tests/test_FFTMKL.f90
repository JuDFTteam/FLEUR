!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
program test
  use MKL_DFTI
  type(DFTI_DESCRIPTOR), pointer :: dfti_handle
  integer :: dfti_status, length(3)
  length = 8
  dfti_status = DftiCreateDescriptor(dfti_handle, DFTI_DOUBLE, DFTI_COMPLEX, 3, length)
  dfti_status = DftiCommitDescriptor(dfti_handle)
  dfti_status = DftiFreeDescriptor(dfti_handle)
end program test
