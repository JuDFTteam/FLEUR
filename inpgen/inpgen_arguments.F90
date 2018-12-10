!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_fleur_arguments
  IMPLICIT NONE
  PRIVATE
  TYPE t_fleur_param
     INTEGER             :: TYPE !can be 0,1,2 for a simple argument, an argument with a string or with a number
     CHARACTER(len=20)   :: name
     CHARACTER(len=200)  :: desc
     CHARACTER(len=200)  :: values
  END TYPE t_fleur_param
  
  INTEGER,PARAMETER:: no_params=7

  TYPE(t_fleur_param) :: fleur_param(no_params)=(/&
       t_fleur_param(0,"-old","Generate input file for old fleur versions",""),&
       t_fleur_param(0,"-explicit","write out k-point list, symmetry operations, and optional input to inp.xml",""),&
       t_fleur_param(0,"-genEnpara","generate an 'enpara' file",""),&
       t_fleur_param(0,"-electronConfig","explicitely write the electron configuration into inp.xml",""),&
       t_fleur_param(0,"-fast_defaults","generate more aggressive (and less stable) input parameters for faster calculations",""),&
       t_fleur_param(0,"-kpts_gw","add alternative k point set for GW",""),&
       t_fleur_param(0,"-h","print this help message","")&
       /)

       
  PUBLIC argument_type,print_argument,check_arguments
CONTAINS

  FUNCTION argument_type(name)
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in):: name
    INTEGER                    :: n,argument_type
    argument_type=-1
    DO n=1,SIZE(fleur_param)
       IF (TRIM(name)==fleur_param(n)%name) argument_type=fleur_param(n)%TYPE
    END DO
  END FUNCTION argument_type
    
  LOGICAL FUNCTION check_arguments()
    IMPLICIT NONE
    INTEGER :: i,n
    CHARACTER(len=200):: str

    check_arguments=.TRUE.
    i=1
    DO WHILE(i<=COMMAND_ARGUMENT_COUNT())
       CALL GET_COMMAND_ARGUMENT(i,str)
       param_loop:DO n=1,SIZE(fleur_param)
          IF (TRIM(str)==fleur_param(n)%name) THEN
             SELECT CASE (fleur_param(n)%TYPE)
             CASE(1)
                i=i+1
                CALL GET_COMMAND_ARGUMENT(i,str)
                IF (TRIM(fleur_param(n)%values)/="") THEN
                   IF (INDEX(TRIM(fleur_param(n)%values),TRIM(str))==0) THEN
                      PRINT *,"Invalid value  :",TRIM(str)
                      PRINT *,"Possible values:",TRIM(fleur_param(n)%values)
                      check_arguments=.false.
                   END IF
                END IF
             CASE(2)
                i=i+1
             END SELECT
             EXIT param_loop
          END IF
       ENDDO param_loop
       IF (n>SIZE(fleur_param)) THEN
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

    DO n=1,no_params
       IF (TRIM(name)==TRIM(fleur_param(n)%name)) THEN
          IF (fleur_param(n)%TYPE==0) THEN !parameter without option
             WRITE(*,1001) TRIM(fleur_param(n)%name),TRIM(fleur_param(n)%desc)
          ELSEIF (fleur_param(n)%TYPE==1) THEN
             IF (fleur_param(n)%values=="") THEN !parameter with string
                WRITE(*,1002) TRIM(fleur_param(n)%name),TRIM(fleur_param(n)%desc)
             ELSE !parameter with string and choice
                WRITE(*,1003) TRIM(fleur_param(n)%name),TRIM(fleur_param(n)%values),TRIM(fleur_param(n)%desc)
             END IF
          ELSE !parameter with number
             WRITE(*,1004) TRIM(fleur_param(n)%name),TRIM(fleur_param(n)%desc)
          ENDIF
          RETURN
       ENDIF
    END DO
1001 FORMAT(t5,a,t20,": ",a)
1002 FORMAT(t5,a," $$$",t20,": ",a)
1003 FORMAT(t5,a," [",a,"]",/,t20,": ",a)
1004 FORMAT(t5,a," #",t20,": ",a)
    
    PRINT *,"BUG IN FLEUR, check handling of parameters in fleur_help.f90"
    PRINT *,name
  END SUBROUTINE print_argument

END MODULE m_fleur_arguments
