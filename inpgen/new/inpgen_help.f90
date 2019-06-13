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
    USE m_check_arguments
    IMPLICIT NONE
    CHARACTER(:), ALLOCATABLE:: infostring

    PRINT *,"     Welcome to FLEUR - inpgen   (www.flapw.de)   "
    PRINT *,"     MaX-Release 3.0          (www.max-centre.eu)"
    
    CALL new_argument(0,"-genEnpara","Generate an 'enpara' file for the energy parameters","") 
    CALL new_argument(0,"-explicit","Write out k-point list, symmetry operations, and optional input to inp.xml","") 
    CALL new_argument(0,"-kpts_gw","add alternative k point set for GW in all outputs for the XML input file","")
    CALL new_argument(0,"-noco","write out noco parameters into inp.xml","")
    CALL new_argument(0,"-electronConfig","explicitely write the electron configuration into inp.xml","")
    CALL new_argument(0,"-fast_defaults","generate more aggressive (and less stable) input parameters for faster calculations","")
    CALL new_argument(0,"-inp.xml","modify existing inp.xml file","")
    CALL new_argument(0,"-inp","convert old inp file to inp.xml file","")
    CALL new_argument(1,"-f","filename to process","")
    CALL new_argument(0,"-warn_only","do not stop for warnings","")
    
    CALL new_argument(0,"-h","Print this help message","")
    
    IF (.NOT.check_arguments()) CALL judft_warn("Invalid command line arguments",hint="Use -h option to see valid choices")
    IF (.NOT. juDFT_was_argument("-h")) RETURN

    !now print version info and help on command line arguments:
    CALL get_compile_desc_string(infostring)
    WRITE(*,'(a)') infostring
    WRITE(*,'(a)')
    WRITE(*,'(a)')"------------------------------------------------------"
    WRITE(*,'(a)')"inpgen usage info:"
    WRITE(*,'(a)')"The following command line options are known:"
    WRITE(*,'(a)')""
    CALL print_argument("-genEnpara")
    CALL print_argument("-explicit")
    CALL print_argument("-noco")
    CALL print_argument("-electronConfig")
    CALL print_argument("-fast_defaults")
    CALL print_argument("-kpts_gw")
    CALL print_argument("-h")
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"Please check the documentation on www.flapw.de for more details"

    CALL juDFT_error("help was written")
  END SUBROUTINE inpgen_help
END MODULE m_inpgen_help
