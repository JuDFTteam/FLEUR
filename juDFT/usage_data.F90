!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_judft_usage
   IMPLICIT NONE
   PRIVATE
   CHARACTER(LEN=99),PARAMETER:: URL_STRING="www.flapw.de/collect.pl"
   INTEGER,PARAMETER :: MAX_NO_KEYS=20
   CHARACTER(LEN=200) :: keys(MAX_NO_KEYS)
   CHARACTER(LEN=200) :: values(MAX_NO_KEYS)
   INTEGER           :: no_keys=0

   INTERFACE add_usage_data
      MODULE PROCEDURE  add_usage_data_s,add_usage_data_i,add_usage_data_l,add_usage_data_r
   END INTERFACE add_usage_data

   PUBLIC :: add_usage_data,send_usage_data

CONTAINS
   SUBROUTINE add_usage_data_s(key,VALUE)
      IMPLICIT NONE
      CHARACTER(len=*),INTENT(IN)::key,VALUE
      no_keys=no_keys+1
      IF (no_keys>MAX_NO_KEYS) STOP "BUG, too many keys in usage_data"
      keys(no_keys)  =key
      values(no_keys)=VALUE
   END SUBROUTINE add_usage_data_s

   SUBROUTINE add_usage_data_i(key,VALUE)
      IMPLICIT NONE
      CHARACTER(len=*),INTENT(IN):: key
      INTEGER,intent(in)         :: value
      CHARACTER(len=20)::txt

      WRITE(txt,*) VALUE
      CALL add_usage_data_s(key,txt)
   END SUBROUTINE add_usage_data_i
   
   SUBROUTINE add_usage_data_r(key,VALUE)
      IMPLICIT NONE
      CHARACTER(len=*),INTENT(IN):: key
      REAL,intent(in)            :: value

      CHARACTER(len=20)::txt
      WRITE(txt,'(F18.2)') VALUE
      CALL add_usage_data_s(key,txt)
   END SUBROUTINE add_usage_data_r

   SUBROUTINE add_usage_data_l(key,VALUE)
      IMPLICIT NONE
      CHARACTER(len=*),INTENT(IN):: key
      LOGICAL,INTENT(in)         :: VALUE

      CHARACTER(len=20)::txt
      txt=MERGE("TRUE ","FALSE",value)
      CALL add_usage_data_s(key,txt)
   END SUBROUTINE add_usage_data_l

   SUBROUTINE send_usage_data()
      IMPLICIT NONE
      INTEGER :: i,ierr,pid,dt(8)
      INTEGER*8 :: r
#ifdef CPP_MPI
      INCLUDE 'mpif.h'

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,i,ierr)
      IF (i.NE.0) RETURN
#endif

!#ifdef CPP_ALLOW_USAGE_DATA
      r = random_id()

      !First write a json file
      OPEN(unit=961,file="usage.json",status='replace')
      WRITE(961,*) '{'
      WRITE(961,*) '  "url":"',TRIM(ADJUSTL(URL_STRING)),'",'
      WRITE(961,"(a,Z0.16,a)") '  "calculation-id":"',r,'",'
      WRITE(961,*) '  "data": {'
      DO i=1,no_keys
         WRITE(961,*) '       "',TRIM(ADJUSTL(keys(i))),'":"',TRIM(ADJUSTL(values(i))),'",'
      ENDDO
      WRITE(961,*) '     }'
      WRITE(961,*) '}'
      CLOSE(961)

#ifdef CPP_CURL
      IF (judft_was_argument("-no_send")) THEN
         PRINT *,"As requested by command line option usage data was not send, please send usage.json manually"
      ELSE
         !Send using curl
         CALL system('curl -H "Content-Type: application/json" --data @usage.json '\\URL_STRING)
         PRINT *,"Usage data send using curl: usage.json"
      ENDIF
#else
      PRINT *,'curl not found in compilation. Please send usage.json manually'
#endif
!#else
!    PRINT *,"No usage data collected"
!#endif
   END SUBROUTINE send_usage_data

   FUNCTION random_id() result(r)
      implicit none
      integer(8)     :: r, rand_i
      integer        :: dt(8), i
      real(8)       :: rand_r

      CALL DATE_AND_TIME(values = dt)
      r = (dt(1) - 1970) * 365 * 24 * 60 * 60 * 1000 + &
          dt(2) * 31 * 24 * 60 * 60 * 1000 + &
          dt(3) * 24 * 60 * 60  * 1000 + dt(5) * 60 * 60 * 1000 + &
          dt(6) * 60 * 1000 + dt(7) * 1000 + dt(8)

      ! 10 times xorshift64
      do i = 1,10
         r = ieor(r, ishft(r,13))
         r = ieor(r, ishft(r,-7))
         r = ieor(r, ishft(r,17))
      enddo
   END FUNCTION random_id
END MODULE m_judft_usage
