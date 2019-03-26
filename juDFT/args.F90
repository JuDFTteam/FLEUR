!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_juDFT_args
!This subroutine allows to query for command line arguments
  PRIVATE
  CHARACTER(len=5):: ENV_NAME="juDFT" !name of environment variable
  PUBLIC judft_was_argument,juDFT_string_for_argument
CONTAINS
  FUNCTION juDFT_was_argument(arg) RESULT(OK)
    USE m_check_arguments
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(IN)::arg
    LOGICAL ok

    INTEGER:: i
    CHARACTER(LEN=30)::str

    !check if argument is allowed
    IF (argument_type(arg)<0) THEN
       PRINT *,"Argument query invalid:",arg
       PRINT *,"Add specification for fleur_arguments"
       STOP "BUG in FLEUR, invalid argument query"
    END IF
    ok=.FALSE.
    DO i=1,COMMAND_ARGUMENT_COUNT()
       CALL GET_COMMAND_ARGUMENT(i,str)
       IF(ADJUSTL(str)==ADJUSTL(arg)) ok=.TRUE.
    ENDDO
    IF (ok) RETURN
    !Test for environment variable as well
    CALL GET_ENVIRONMENT_VARIABLE("juDFT",str,status=i)
    IF (i==0) ok=INDEX(str//' ',TRIM(ADJUSTL(arg))//' ')>0

  END FUNCTION juDFT_was_argument

  FUNCTION juDFT_string_for_argument(arg) RESULT(argstring)
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(IN)::arg
    CHARACTER(len=1000)::argstring

    INTEGER:: i
    CHARACTER(LEN=30)  ::str
    CHARACTER(LEN=1000)::env
    argstring=""
    IF (judft_was_argument(arg)) THEN
       DO i=1,COMMAND_ARGUMENT_COUNT()
          CALL GET_COMMAND_ARGUMENT(i,str)
          IF(ADJUSTL(str)==ADJUSTL(arg)) THEN
             IF (i<=COMMAND_ARGUMENT_COUNT()) CALL GET_COMMAND_ARGUMENT(i+1,argstring)
             RETURN !Argument found
          ENDIF
       ENDDO
       !Not found so must be in environment variable
       CALL GET_ENVIRONMENT_VARIABLE("juDFT",env,status=i)
       !cut of argument string
       env=ADJUSTL(env(INDEX(env//' ',TRIM(ADJUSTL(arg)))+1:))//' '
       i=INDEX(env,' ') !find first blank
       argstring=TRIM(env(:i))
    END IF
    
  END FUNCTION juDFT_string_for_argument
END MODULE m_juDFT_args
