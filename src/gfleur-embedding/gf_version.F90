!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_version
      use m_juDFT
      IMPLICIT NONE
!-----------------------------------------------
! DESC:Displays version information
!                 Daniel Wortmann, (05-08-12)
!-----------------------------------------------
      CONTAINS

      !<-- S: gf_hello()
      SUBROUTINE gf_hello()
!-----------------------------------------------
! Say hello!
!-----------------------------------------------
      USE m_gf_out
      USE m_constants, ONLY: oUnit, version_const
      USE m_compile_descr

      CHARACTER(LEN=:), ALLOCATABLE :: gitdesc, githash, gitbranch,      &
     &                                 compile_date, compile_user,       &
     &                                 compile_host, compile_flags,      &
     &                                 link_flags

      !OPEN files and say hello
      CALL gf_out_create()

      CALL get_compile_desc(gitdesc, githash, gitbranch, compile_date,   &
     &                      compile_user, compile_host, compile_flags,   &
     &                      link_flags)

      WRITE(oUnit,*) "              This is GFLEUR"
      WRITE(16,*)    "              This is GFLEUR"
      WRITE(oUnit,*)
      WRITE(16,*)
      WRITE(oUnit,*) "Settings:"
      WRITE(oUnit,"(a,a50)") "Based on   :",version_const
      WRITE(oUnit,"(a,a50)") "Branch     :",gitbranch
      WRITE(oUnit,"(a,a50)") "Git hash   :",githash
      WRITE(oUnit,"(a,a50)") "Compiled at:",compile_date
      WRITE(oUnit,"(a,a50)") "         by:",compile_user
      WRITE(oUnit,"(a,a50)") "         on:",compile_host
#ifdef CPP_MPI
      WRITE(oUnit,*) "          Parallel-version (CPP_MPI)"
#endif

      END SUBROUTINE
      !>

      END
