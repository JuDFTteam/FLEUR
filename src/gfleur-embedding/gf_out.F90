!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_out
      USE m_constants, ONLY: oUnit
      IMPLICIT NONE
!*****************************************************************
! DESC:Module provides subroutines to handle the output to the
!      out-file (oUnit, as in FLEUR) and the inf-file (unit 16)
!                          Daniel Wortmann, Tue May  6 17:33:58 2003
!*****************************************************************
      PRIVATE
      !opens the out&inf file
      PUBLIC gf_out_create
      !closes both files again
      PUBLIC gf_out_close
      !writes a new section header
      PUBLIC gf_out_newheader
      CONTAINS

      SUBROUTINE gf_out_create()
      LOGICAL :: l_opened
      !Open real files on Pe0, scratch files otherwise
      IF (isPe0()) THEN
         OPEN(16,FILE="inf",FORM="FORMATTED",ACTION="WRITE",             &
     &        STATUS="REPLACE")
         !the out-file might already have been opened by the juDFT setup
         INQUIRE(UNIT=oUnit,OPENED=l_opened)
         IF (.NOT.l_opened) OPEN(oUnit,FILE="out",FORM="FORMATTED",      &
     &        ACTION="WRITE",STATUS="REPLACE")
      ELSE
         CLOSE(oUnit,ERR=111)
         CLOSE(16,ERR=111)
  111    OPEN(16,FILE="inf.parallel",FORM="FORMATTED")
         OPEN(oUnit,FILE="out.parallel",FORM="FORMATTED")
      ENDIF

      END SUBROUTINE

      SUBROUTINE gf_out_close()
      CLOSE(oUnit)
      CLOSE(16)
      END SUBROUTINE


      SUBROUTINE gf_out_newheader(name)
      CHARACTER*(*),INTENT(IN)::name
      IF (.NOT.isPe0()) RETURN
      WRITE(oUnit,999) name
      WRITE(16,999) name
  999 FORMAT (/,/,'<***',/,10x,a,/,'***>')
      END SUBROUTINE

      FUNCTION isPe0()
!*****************************************************************
!  private function to check for pe=0
!*****************************************************************
#ifdef CPP_MPI
      USE mpi
#endif
      LOGICAL::isPe0
#ifdef CPP_MPI
      INTEGER ierr,irank
      CALL MPI_COMM_RANK (MPI_COMM_WORLD,irank,ierr)
      isPe0=irank==0
#else
      isPe0=.TRUE.
#endif
      END FUNCTION
      END
