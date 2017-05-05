!--------------------------------------------------------------------------------
! Copyright (c) 2017 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_inpgen_help
  IMPLICIT NONE
CONTAINS
  SUBROUTINE inpgen_help()
    USE m_compile_descr
    USE m_constants
    USE m_juDFT
    IMPLICIT NONE
    CHARACTER(LEN=500):: infostring

    PRINT *,"     Welcome to FLEUR - inpgen   (www.flapw.de)   "
    PRINT *,"     MaX-Release 1.2           (www.max-centre.eu)"

    IF (.NOT. (juDFT_was_argument("-h").OR.juDFT_was_argument("--help"))) RETURN

    !now print version info and help on command line arguments:
    CALL get_compile_desc_string(infostring)
    WRITE(*,'(a500)') infostring
    WRITE(*,'(a)')
    WRITE(*,'(a)')"------------------------------------------------------"
    WRITE(*,'(a)')"inpgen usage info:"
    WRITE(*,'(a)')"The following command line options are known:"
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"-old              : generate input files for old fleur versions"
    WRITE(*,'(a)')"-genEnpara        : write enpara file"
    WRITE(*,'(a)')"-explicit         : write out k-point list, symmetry operations,"
    WRITE(*,'(a)')"                    and optional input to inp.xml"
    WRITE(*,'(a)')"-fast_defaults    : generate more aggressive (and less stable)"
    WRITE(*,'(a)')"                    input parameters for faster calculations"
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"-h, --help        : print this text :-)"
    WRITE(*,'(a)')""
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"Please check the documentation on www.flapw.de for more details"

    CALL juDFT_error("help was written")
  END SUBROUTINE inpgen_help
END MODULE m_inpgen_help
