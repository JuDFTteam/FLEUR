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
    PRINT *,version_const_MaX

    CALL new_argument(0,"-genEnpara","Generate an 'enpara' file for the energy parameters","")
    CALL new_argument(0,"-explicit","Write out k-point list, symmetry operations, and optional input to inp.xml","")
    CALL new_argument(0,"-kpts_gw","add alternative k point set for GW in all outputs for the XML input file","")
    CALL new_argument(0,"-noco","write out noco parameters into inp.xml","")
    CALL new_argument(0,"-greensf","write out green's function parameters into inp.xml","")
    CALL new_argument(0,"-electronConfig","explicitely write the electron configuration into inp.xml","")
    CALL new_argument(0,"-fast_defaults","generate more aggressive (and less stable) input parameters for faster calculations","")
    CALL new_argument(0,"-inp.xml","modify existing inp.xml file","")
    CALL new_argument(0,"-inp","convert old inp file to inp.xml file","")
    CALL new_argument(1,"-f","filename to process","")
    CALL new_argument(1,"-o","output filename (default inp.xml) (Use only if you convert an old inp.xml file)","")
    CALL new_argument(0,"-warn_only","do not stop for warnings","")
    CALL new_argument(1,"-inc","which data to include in inp.xml, e.g. +all,-species,+operations,-kpts","")
    CALL new_argument(1,"-kptsPath","define a special-k point e.g. -kptsPath 'X=0.5,0.5,0.5;g=0.0,0.0,0.0'","")

    CALL new_argument(1,"-xmlXPath","specify an xml path and value to overwrite inp.xml","")
    CALL new_argument(1,"-kpt","String to define k-point set e.g. -kpt name#nk=2","")
    call new_argument(0,"-no_send","Do not send usage data","")
    CALL new_argument(0,"-overwrite","Overwrite inp.xml if present","")
    CALL new_argument(0,"-h","Print this help message","")
    CALL new_argument(0,"-version","Show important version information about the inpgen executable","")
    CALL new_argument(0,"-dropXMLSchema","Write out the default XML schema files","")
    CALL new_argument(0,"-trace","Try to generate a traceback in case of an error","")
    CALL new_argument(0,"-debugtime","Write the start/stop of all timers to the console","")
    CALL new_argument(0,"-all_times","Write json files of timing for all PE, not only for PE=0","")

    CALL new_argument(0,"-nosym","In general: Only use identity as symmetry operation","")
    CALL new_argument(0,"-noKsym","For the generation of the k-point set: Only use identity as symmetry operation","")
    CALL new_argument(0,"-noInpgenComment","Disable printing inpgen input as a comment to the inp.xml file","")

    CALL new_argument(0,"-precise","Short for '-profile precise'","")
    CALL new_argument(1,"-profile","Generate default parameters according to provided profile name. This can be one of 'precise'.","")

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
    !CALL print_argument("-genEnpara")
    !CALL print_argument("-explicit")
    !CALL print_argument("-noco")
    !CALL print_argument("-electronConfig")
    !CALL print_argument("-fast_defaults")
    CALL print_argument("-f")
    CALL print_argument("-inp.xml")
    CALL print_argument("-o")
    CALL print_argument("-inp")
    CALL print_argument("-kpt")
    CALL print_argument("-inc")
    CALL print_argument("-overwrite")
    CALL print_argument("-noco")
    CALL print_argument("-greensf")
    CALL print_argument("-warn_only")
    CALL print_argument("-kptsPath")
    CALL print_argument("-nosym")
    CALL print_argument("-noKsym")
    CALL print_argument("-trace")
    CALL print_argument("-version")
    CALL print_argument("-dropXMLSchema")
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"Please check the documentation on www.flapw.de for more details"

    CALL juDFT_error("help was written")
  END SUBROUTINE inpgen_help
END MODULE m_inpgen_help
