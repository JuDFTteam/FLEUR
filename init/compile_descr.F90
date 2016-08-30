!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_compile_descr

  IMPLICIT NONE
!Use the Preprocessor variables to determine compile environment
  CONTAINS
  SUBROUTINE get_compile_desc(gitdesc,githash,compile_date,compile_user,compile_host)
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(OUT)::gitdesc,githash,compile_date,compile_user,compile_host

#include "compileinfo.h"
    
  END subroutine get_compile_desc
end MODULE m_compile_descr

