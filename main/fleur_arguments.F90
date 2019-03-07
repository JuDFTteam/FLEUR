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
  
  INTEGER,PARAMETER:: no_params=25
  TYPE(t_fleur_param) :: fleur_param(no_params)=(/&
       !Input options
       t_fleur_param(0,"-toXML","Convert an old 'inp' file into the new XML format",""),&
       t_fleur_param(1,"-xmlXPath","modify the xml-xpath of the inp.xml file",""),&
       !Control the job
       t_fleur_param(0,"-check","run in check mode, i.e. stop after init",""),&
       t_fleur_param(0,"-info","Print out information on recommended parallelization and available charge densities",""),&
       t_fleur_param(2,"-wtime","run for # minutes (used to estimate if another iteration is started)",""),&
       t_fleur_param(1,"-j","Distribute MPI ranks to run subjobs (Format PE:DIR meaning run with PE in directory DIR)",""),&
       t_fleur_param(1,"-f","Obtain info on subjobs from file",""),&
       t_fleur_param(2,"-n_min_size","Try to use at least specified number of PE in eigenvalue parallelization",""),&
       t_fleur_param(1,"-diag","Choose method for diagonalization","lapack,debugout"&
#ifdef CPP_SCALAPACK
       //",scalapack"&
#endif
#ifdef CPP_ELPA_ONENODE
       //",elpa_1node"&
#endif
#ifdef CPP_ELPA
       //",elpa"&
#endif
#ifdef CPP_CHASE
       //",chase"&
#endif
#ifdef CPP_MAGMA
       //",magma"&
#endif
#ifdef CPP_GPU
       //",cusolver"&
#endif
       ),&
       t_fleur_param(1,"-eig","Method for storing the eigenvectors","mem,da"&
#ifdef CPP_MPI
       //",mpi"&
#endif
#ifdef CPP_HDF
       //",hdf"&
#endif
       ),&
       !Debugging 
       t_fleur_param(0,"-warn_only","Continue execution after a warning message",""),&
       t_fleur_param(0,"-trace","Try to generate a stacktrace in case of an error",""),&
       t_fleur_param(0,"-debugtime","Write the start/stop of all timers to the console",""),&
       !Output
       t_fleur_param(0,"-mix_io","Do not store mixing history in memory but do IO in each iteration",""),&
       t_fleur_param(0,"-no_out","Do not open the 'out' file but write to stdout",""),&
       t_fleur_param(0,"-genEnpara","Generate an 'enpara' file for the energy parameters",""),&
       t_fleur_param(0,"-kpts_gw","add alternative k point set for GW in all outputs for the XML input file",""),&
       t_fleur_param(0,"-noco","write out noco parameters in all outputs for inp.xml",""),&
       t_fleur_param(0,"-h","Print this message",""),&
       t_fleur_param(0,"-no_send","Do not send usage data","")&
       !HDF density
       ,t_fleur_param(0,"-no_cdn_hdf","Disable HDF charge density mode (activated by default if HDF5 is available)","")&
       ,t_fleur_param(0,"-last_extra","Generate an additional file cdn_last.hdf that contains only the last density","")&
       ,t_fleur_param(2,"-sd","use starting density N, where N is the index of the density according to -info","")&
       ,t_fleur_param(1,"-delden","delete densities (either an index N, a range N-M or the keyword 'allbutlast' should be given)","")&
       !GPU parameter
       ,t_fleur_param(0,"-gpu","Use GPU for computing","")&
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
