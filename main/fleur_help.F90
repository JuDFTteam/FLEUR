!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_fleur_help
  IMPLICIT NONE
  PRIVATE
  TYPE t_fleur_param
     INTEGER             :: TYPE
     CHARACTER(len=20)   :: name
     CHARACTER(len=200)  :: desc
     CHARACTER(len=200)  :: values
  END TYPE t_fleur_param
     
#ifdef CPP_HDF
  INTEGER,PARAMETER:: no_params=20
#else
  INTEGER,PARAMETER:: no_params=16
#endif
  TYPE(t_fleur_param) :: fleur_param(no_params)=(/&
       !Input options
       t_fleur_param(0,"-toXML","Convert an old 'inp' file into the new XML format",""),&
       t_fleur_param(1,"-xmlXpath","modify the xml-xpath of the inp.xml file",""),&
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
#ifdef CPP_ELPA
       //",elpa"&
#endif
#ifdef CPP_CHASE
       //",chase"&
#endif
#ifdef CPP_MAGMA
       //",magma"&
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
       t_fleur_param(0,"-no_out","Do not open the 'out' file but write to stdout",""),&
       t_fleur_param(0,"-gen_enpara","Generate an 'enpara' file for the energy parameters",""),&
       t_fleur_param(0,"-h","Print this message","")&
       !HDF density
#ifdef CPP_HDF       
       ,t_fleur_param(0,"-no_cdn_hdf","Disable HDF charge density mode (activated by default if HDF5 is available)","")&
       ,t_fleur_param(0,"-last_extra","Generate an additional file cdn_last.hdf that contains only the last density","")&
       ,t_fleur_param(2,"-sd","use starting density N, where N is the index of the density according to -info","")&
       ,t_fleur_param(1,"-delden","delete densities (either an index N, a range N-M or the keyword 'allbutlast' should be given)","")&
#endif       
       /)

       
  PUBLIC fleur_help
CONTAINS
  SUBROUTINE check_arguments()
    USE m_judft_stop
    IMPLICIT NONE
    INTEGER :: i,n
    CHARACTER(len=200):: str
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
                      CALL judft_warn("Invalid value to command line argument")
                   END IF
                END IF
             CASE(2)
                i=i+1
             END SELECT
             EXIT param_loop
          END IF
       ENDDO param_loop
       IF (n>SIZE(fleur_param)) CALL judft_warn("Unkown command line argument:"//str)
       i=i+1
    ENDDO

  END SUBROUTINE check_arguments
  
  SUBROUTINE print_param(name)
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
  END SUBROUTINE print_param
  
  SUBROUTINE fleur_help()
    USE m_compile_descr
    USE m_constants
    USE m_juDFT
    IMPLICIT NONE
    CHARACTER(LEN=500):: infostring

    PRINT *,"     Welcome to FLEUR        (www.flapw.de)   "
    PRINT *,"     MaX-Release 2.1       (www.max-centre.eu)"

    CALL check_arguments()
    
    IF (.NOT. (juDFT_was_argument("-h").OR.juDFT_was_argument("--help"))) RETURN

    !now print version info and help on command line arguments:
    CALL get_compile_desc_string(infostring)
    WRITE(*,'(a500)') infostring
    WRITE(*,'(a)')
    WRITE(*,'(a)')"------------------------------------------------------"
    WRITE(*,'(a)')"Usage info:"
    WRITE(*,'(a)')"The following command line options are known:"
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"Control the input:"
    CALL print_param("-toXML")
    CALL print_param("-xmlXpath")
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"Output options:"
    CALL print_param("-no_out")
    CALL print_param("-gen_enpara")
    CALL print_param("-h")
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"Control FLEUR job:"
    CALL print_param("-check")
    CALL print_param("-info")
    CALL print_param("-wtime")
    CALL print_param("-j")
    CALL print_param("-f")
    CALL print_param("-n_min_size")
    CALL print_param("-diag")
    CALL print_param("-eig")
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"Options useful for debugging:"
    CALL print_param("-warn_only")
    CALL print_param("-trace")
    CALL print_param("-debugtime")
#ifdef CPP_HDF
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"HDF density file relevant options:"

    CALL print_param("-no_cdn_hdf")
    CALL print_param("-last_extra")
    CALL print_param("-sd")
    CALL print_param("-delden")
#endif  
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"Please check the documentation on www.flapw.de for more details."

    CALL juDFT_end("",l_endXML=.FALSE.) !No message so do a not print more on exit
  END SUBROUTINE fleur_help
END MODULE m_fleur_help
