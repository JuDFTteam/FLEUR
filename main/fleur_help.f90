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
    PRINT *,"     MaX-Release 3.0       (www.max-centre.eu)"

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
    WRITE(*,'(a)')"-xmlInput or -xml : use inp.xml instead of inp"
    WRITE(*,'(a)')"-toXML            : convert inp file to XML input file"
    WRITE(*,'(a)')"-xmlXPath XXX=YYY : modify the xml-xpath of the inp.xml file"
    WRITE(*,'(a)')"-n_min_size  XXX  : try to use at least XXX PE in Eigenvalue parallelization" 
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"Output options:"
    WRITE(*,'(a)')"-no_out           : do not open out file but write to stdout"
    WRITE(*,'(a)')"-genEnpara        : write enpara file"   
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"Control FLEUR job:"
    WRITE(*,'(a)')"-check            : run in check mode, i.e. stop after init"
    WRITE(*,'(a)')"-wtime XXXXX      : run for XXXX minutes (used to estimate if another iteration is started)"
    WRITE(*,'(a)')"-j #:DIR          : run subjob in directory DIR using # PEs"
    WRITE(*,'(a)')"-f FILENAME       : obtain info on subjobs from file FILENAME"
    WRITE(*,'(a)')"-info             : Print out information on recommended parallelization and available charge densities"
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"Options useful for debugging:"
    WRITE(*,'(a)')"-warn_only        : Do not stop in case of warnings"
    WRITE(*,'(a)')"-trace            : Try to generate stacktrace in case of an error"
    WRITE(*,'(a)')"-debugtime        : write out the start/stop of all timers to STDOUT"
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"Storage for eigenvalues:"
    WRITE(*,'(a)')"-da,-mem,-mpi,-hdf: choose a storage for the eigenvalues"
    WRITE(*,'(a)')"                    and eigenvectors. The default will"
    WRITE(*,'(a)')"                    be -mem for serial and -mpi for parallel builds" 
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"Method used for diagonalization:"
    WRITE(*,'(a)')"-lapack,-lapack2,"
    WRITE(*,'(a)')"-elpa,-scalapack,"
    WRITE(*,'(a)')"-elemental,-magma : choose diagonalization. Not all might be available"
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"HDF density file relevant options:"
    WRITE(*,'(a)')"-no_cdn_hdf        : disable HDF charge density mode (activated by default if HDF5 is available)"
    WRITE(*,'(a)')"-sd N              : use starting density N, where N is the index of the density according to -info"
    WRITE(*,'(a)')"-delden N-M        : delete densities N to M"
    WRITE(*,'(a)')"-delden N          : delete density N"
    WRITE(*,'(a)')"-delden allbutlast : delete all but the last density"
    WRITE(*,'(a)')""   
    WRITE(*,'(a)')"-h, --help        : print this text :-)"
    WRITE(*,'(a)')"Please check the documentation on www.flapw.de for more details"

    CALL juDFT_end("help was written",l_endXML=.FALSE.)
  END SUBROUTINE fleur_help
END MODULE m_fleur_help
