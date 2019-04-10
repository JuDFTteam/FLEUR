!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_check_arguments
  IMPLICIT NONE
  PRIVATE
  TYPE t_param
     INTEGER             :: TYPE !can be 0,1,2 for a simple argument, an argument with a string or with a number
     CHARACTER(len=20)   :: name
     CHARACTER(len=200)  :: desc
     CHARACTER(len=200)  :: values
  END TYPE t_param

  TYPE(t_param),ALLOCATABLE:: params(:)

       
  PUBLIC argument_type,print_argument,check_arguments,new_argument
CONTAINS

  SUBROUTINE new_argument(argtype,arg,desc,values)
    INTEGER,INTENT(in)          :: argtype
    CHARACTER(len=*),INTENT(in) :: arg,desc,values
    
    TYPE(t_param),ALLOCATABLE::tmp(:)
    !extend the params array
    IF (ALLOCATED(params)) THEN
       CALL MOVE_ALLOC(params,tmp)
       ALLOCATE(params(SIZE(tmp)+1))
       params(:SIZE(tmp))=tmp
    ELSE
       ALLOCATE(params(1))
    ENDIF
    params(SIZE(params))=t_param(argtype,arg,desc,values)
  END SUBROUTINE new_argument
  
  FUNCTION argument_type(name)
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in):: name
    INTEGER                    :: n,argument_type
    argument_type=-1
    IF (.NOT.ALLOCATED(params)) THEN
       !nothing to check
       argument_type=0
       RETURN
    END IF
    DO n=1,SIZE(params)
       IF (TRIM(name)==params(n)%name) argument_type=params(n)%TYPE
    END DO
  END FUNCTION argument_type
    
  LOGICAL FUNCTION check_arguments()
    IMPLICIT NONE
    INTEGER :: i,n
    CHARACTER(len=200):: str

    check_arguments=.TRUE.
    IF (.NOT.ALLOCATED(params)) RETURN
    i=1
    DO WHILE(i<=COMMAND_ARGUMENT_COUNT())
       CALL GET_COMMAND_ARGUMENT(i,str)
       param_loop:DO n=1,SIZE(params)
          IF (TRIM(str)==params(n)%name) THEN
             SELECT CASE (params(n)%TYPE)
             CASE(1)
                i=i+1
                CALL GET_COMMAND_ARGUMENT(i,str)
                IF (TRIM(params(n)%values)/="") THEN
                   IF (INDEX(TRIM(params(n)%values),TRIM(str))==0) THEN
                      PRINT *,"Invalid value  :",TRIM(str)
                      PRINT *,"Possible values:",TRIM(params(n)%values)
                      check_arguments=.false.
                   END IF
                END IF
             CASE(2)
                i=i+1
             END SELECT
             EXIT param_loop
          END IF
       ENDDO param_loop
       IF (n>SIZE(params)) THEN
          PRINT *,"Unkown command line argument:"//str
          check_arguments=.FALSE.
       END IF
       i=i+1
    ENDDO
  END FUNCTION check_arguments
  
  SUBROUTINE print_argument(name)
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in):: name
  
    INTEGER :: n

    IF (.NOT.ALLOCATED(params)) RETURN
    DO n=1,size(params)
       IF (TRIM(name)==TRIM(params(n)%name)) THEN
          IF (params(n)%TYPE==0) THEN !parameter without option
             WRITE(*,1001) TRIM(params(n)%name),TRIM(params(n)%desc)
          ELSEIF (params(n)%TYPE==1) THEN
             IF (params(n)%values=="") THEN !parameter with string
                WRITE(*,1002) TRIM(params(n)%name),TRIM(params(n)%desc)
             ELSE !parameter with string and choice
                WRITE(*,1003) TRIM(params(n)%name),TRIM(params(n)%values),TRIM(params(n)%desc)
             END IF
          ELSE !parameter with number
             WRITE(*,1004) TRIM(params(n)%name),TRIM(params(n)%desc)
          ENDIF
          RETURN
       ENDIF
    END DO
1001 FORMAT(t5,a,t20,": ",a)
1002 FORMAT(t5,a," $$$",t20,": ",a)
1003 FORMAT(t5,a," [",a,"]",/,t20,": ",a)
1004 FORMAT(t5,a," #",t20,": ",a)
    
    PRINT *,"BUG, check handling of parameters in check_arguments.f90"
    PRINT *,name
  END SUBROUTINE print_argument

END MODULE m_check_arguments
