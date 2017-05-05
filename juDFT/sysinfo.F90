MODULE m_judft_sysinfo
  IMPLICIT NONE
CONTAINS
  SUBROUTINE checkstack()
    CHARACTER(LEN=10):: l1,l2,l3,l4
    INTEGER          :: err
    LOGICAL          :: unlimited
    unlimited=.TRUE.  !set to true by default. 
                      !If /proc/self/limits does not exist
                      !or parsing fails no warning is issued
    OPEN(99,FILE="/proc/self/limits",ERR=999)
    DO
       READ(99,*,ERR=999,END=999) l1,l2,l3,l4
       IF (ALL((/INDEX(l1,"Max"),INDEX(l2,"stack"),INDEX(l3,"size")/)==1)) THEN
          unlimited=INDEX(l4,"unlim")==1
       ENDIF
    ENDDO
999 CLOSE(99,IOSTAT=err)
    IF (.NOT.unlimited) THEN
       WRITE(*,*) "*********** WARNING! ************"
       WRITE(*,*) "Your stacksize seems to be small"
       WRITE(*,*) "FLEUR might crash without further"
       WRITE(*,*) "notice. Try 'ulimit -s unlimited'"
       WRITE(*,*) "*********** WARNING! ************"
    ENDIF
  END SUBROUTINE checkstack
END MODULE m_judft_sysinfo

!program test
!  use m_judft_sysinfo
!  call checkstack()

!end program test
