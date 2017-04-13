!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_fleur_help
  IMPLICIT NONE
CONTAINS
  SUBROUTINE fleur_help()
    USE m_compile_descr
    USE m_constants
    USE m_juDFT
    IMPLICIT NONE
    CHARACTER(LEN=500):: infostring

    PRINT *,"     Welcome to FLEUR        (www.flapw.de)   "
    PRINT *,"     MaX-Release 1         (www.max-centre.eu)"

    IF (.NOT. (juDFT_was_argument("-h").OR.juDFT_was_argument("--help"))) RETURN

    !now print version info and help on command line arguments:
    CALL get_compile_desc_string(infostring)
    WRITE(*,'(a500)') infostring
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
    WRITE(*,*)"-xmlInput or -xml : use inp.xml instead of inp"
    WRITE(*,*)""
    WRITE(*,*)"-check            : run in check mode, i.e. stop after init"
    WRITE(*,*)""
    WRITE(*,*)"-n_min_size  XXX  : try to use at least XXX PE in Eigenvalue parallelization" 
    WRITE(*,*)""
    WRITE(*,*)"-wtime XXXXX      : run for XXXX minutes (used to estimate if another iteration is started"
    WRITE(*,*)""
    WRITE(*,*)"-j #:DIR          : run subjob in directory DIR using # PEs"
    WRITE(*,*)"-f FILENAME       : obtain info on subjobs from file FILENAME"
    WRITE(*,*)""
    WRITE(*,*)"-info             : Print out information on recommended parallelization and available charge densities"
    WRITE(*,*)""
    WRITE(*,*)"-h, --help        : print this text :-)"
    WRITE(*,*)""
    WRITE(*,*)"HDF density file relevant options:"
    WRITE(*,*)""
    WRITE(*,*)"-hdf_cdn          : activate HDF charge density mode"
    WRITE(*,*)"-sd N             : use starting density N, where N is the index of the density according to -info"
    WRITE(*,*)"-delden N-M       : delete densities N to M"
    WRITE(*,*)"-delden N         : delete density N"
    WRITE(*,*)""
    WRITE(*,*)"Please check the documentation on www.flapw.de for more details"

    CALL juDFT_error("help was written")
  END SUBROUTINE fleur_help
END MODULE m_fleur_help
