!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_compile_descr

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE get_compile_desc_string(info)
    USE m_constants
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(OUT):: info
    
    CHARACTER(len=50)::gitdesc,githash,compile_date,compile_user,compile_host,gitbranch

    call  get_compile_desc(gitdesc,githash,gitbranch,compile_date,compile_user,compile_host)
    info=new_line("a")// &
    "This is FLEUR version: "//trim(version_const)//new_line("a")// &
    "FLEUR was compiled:"//new_line("a")// &
    "   at: "//TRIM(compile_date)//new_line("a")// &
    "   by: "//TRIM(compile_user)//new_line("a")// &
    "   on: "//TRIM(compile_host)//new_line("a")// &
    "Its git version is:"//new_line("a")// &
    "   described by: "//TRIM(gitdesc)//new_line("a")// &
    "   from branch:  "//trim(gitbranch)//new_line("a")// &
    "   with hash:    "//TRIM(githash)
  end SUBROUTINE get_compile_desc_string


  SUBROUTINE get_compile_desc(gitdesc,githash,gitbranch,compile_date,compile_user,compile_host)
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(OUT)::gitdesc,githash,compile_date,compile_user,compile_host,gitbranch

!This file is created by cmake at time of configuration
#include "compileinfo.h"
    
  END subroutine get_compile_desc
end MODULE m_compile_descr

