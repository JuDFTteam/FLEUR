!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_fleur_info
  IMPLICIT NONE
CONTAINS
  SUBROUTINE fleur_info()
    USE m_compile_descr
    USE m_constants
    USE m_juDFT
    IMPLICIT NONE
    CHARACTER(LEN=50):: gitdesc,githash,compile_date,compile_user,compile_host
    
    PRINT *,"     Welcome to FLEUR        (www.flapw.de)   "
    PRINT *,"     MaX-Release 1         (www.max-centre.eu)"

    IF (.NOT. (juDFT_was_argument("-h").OR.juDFT_was_argument("--help"))) RETURN

    !now print version info and help on command line arguments:
    CALL get_compile_desc(gitdesc,githash,compile_date,compile_user,compile_host)
    WRITE(*,*) "This is version: ",version_const
    WRITE(*,*) "FLEUR was compiled at ",TRIM(compile_date)," by ",TRIM(compile_user)," on ",TRIM(compile_host)
    WRITE(*,*) "Its git version is ",TRIM(gitdesc)," with hash: ",TRIM(githash)
    WRITE(*,*)
    WRITE(*,*)"------------------------------------------------------"
    WRITE(*,*)"Usage info:"
    WRITE(*,*)"The following command line options are known:"
    WRITE(*,*)"-da,-mem,-mpi,-hdf: choose a storage for the eigenvalues"
    WRITE(*,*)"                    and eigenvectors. The default will depend"
    WRITE(*,*)"                    be -mem for serial and -mpi for parallel builds" 
    WRITE(*,*)""
    WRITE(*,*)"-lapack,-lapack2,"
    write(*,*)"-elpa,-scalapack,"
    WRITE(*,*)"-elemental,-magma : choose diagonalization. Not all might be available"
    WRITE(*,*)""
    WRITE(*,*)"-debugtime        : write out the start/stop of all timers to STDOUT"
    WRITE(*,*)""
    WRITE(*,*)"-genEnpara        : write enpara file"
    WRITE(*,*)""
    WRITE(*,*)"-xmlInput         : use inp.xml instead of inp"
    WRITE(*,*)""
    WRITE(*,*)"-j #:DIR          : run subjob in directory DIR using # PEs"
    WRITE(*,*)"-f FILENAME       : obtain info on subjobs from file FILENAME"
    WRITE(*,*)""
    WRITE(*,*)"-h, --help        : print this text :-)"
    WRITE(*,*)""
    WRITE(*,*)"Please check the documentation on www.flapw.de for more details"

    CALL juDFT_error("help was written")
  END SUBROUTINE fleur_info
END MODULE m_fleur_info
