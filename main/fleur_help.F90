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
    USE m_fleur_arguments
    IMPLICIT NONE
 
    CHARACTER(:), ALLOCATABLE:: infostring

    PRINT *,"     Welcome to FLEUR        (www.flapw.de)   "
    PRINT *,"     MaX-Release 2.1       (www.max-centre.eu)"

    IF (.NOT.check_arguments()) CALL judft_warn("Invalid command line arguments",hint="Use -h option to see valid choices")
    
    
    IF (.NOT. judft_was_argument("-h")) RETURN

    !now print version info and help on command line arguments:
    CALL get_compile_desc_string(infostring)
    WRITE(*,'(a)') infostring
    WRITE(*,'(a)')
    WRITE(*,'(a)')"------------------------------------------------------"
    WRITE(*,'(a)')"Usage info:"
    WRITE(*,'(a)')"The following command line options are known:"
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"Control the input:"
    CALL print_argument("-toXML")
    CALL print_argument("-xmlXpath")
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"Output options:"
    CALL print_argument("-no_out")
    CALL print_argument("-gen_enpara")
    CALL print_argument("-h")
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"Control FLEUR job:"
#ifdef CPP_GPU
    CALL print_argument("-gpu")
#endif
    CALL print_argument("-check")
    CALL print_argument("-info")
    CALL print_argument("-wtime")
    CALL print_argument("-j")
    CALL print_argument("-f")
    CALL print_argument("-n_min_size")
    CALL print_argument("-diag")
    CALL print_argument("-eig")
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"Options useful for debugging:"
    CALL print_argument("-warn_only")
    CALL print_argument("-trace")
    CALL print_argument("-debugtime")
#ifdef CPP_HDF
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"HDF density file relevant options:"

    CALL print_argument("-no_cdn_hdf")
    CALL print_argument("-last_extra")
    CALL print_argument("-sd")
    CALL print_argument("-delden")
#endif
    WRITE (*,'(a)') "Options for privacy sensitive users"
    CALL print_argument("-no_send")
    
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"Please check the documentation on www.flapw.de for more details."

    CALL juDFT_end("",l_endXML=.FALSE.) !No message so do a not print more on exit
  END SUBROUTINE fleur_help
END MODULE m_fleur_help
