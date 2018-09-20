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
    USE m_fleur_arguments
    IMPLICIT NONE
    CHARACTER(LEN=500):: infostring

    PRINT *,"     Welcome to FLEUR - inpgen   (www.flapw.de)   "
    PRINT *,"     MaX-Release 2.1          (www.max-centre.eu)"

    IF (.NOT. juDFT_was_argument("-h")) RETURN

    !now print version info and help on command line arguments:
    CALL get_compile_desc_string(infostring)
    WRITE(*,'(a500)') infostring
    WRITE(*,'(a)')
    WRITE(*,'(a)')"------------------------------------------------------"
    WRITE(*,'(a)')"inpgen usage info:"
    WRITE(*,'(a)')"The following command line options are known:"
    WRITE(*,'(a)')""
    CALL print_argument("-old")
    CALL print_argument("-genEnpara")
    CALL print_argument("-explicit")
    CALL print_argument("-electronConfig")
    CALL print_argument("-fast_defaults")
    CALL print_argument("-h")
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"Please check the documentation on www.flapw.de for more details"

    CALL juDFT_error("help was written")
  END SUBROUTINE inpgen_help
END MODULE m_inpgen_help
