!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_fleur_help
  IMPLICIT NONE
  PRIVATE
  PUBLIC fleur_help
CONTAINS
  SUBROUTINE fleur_help()
    USE m_compile_descr
    USE m_constants
    USE m_juDFT
    USE m_check_arguments
    IMPLICIT NONE

    CHARACTER(:), ALLOCATABLE:: infostring

    PRINT *,"     Welcome to FLEUR        (www.flapw.de)   "
    PRINT *, version_const_MaX
    CALL add_fleur_arguments()
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
    CALL print_argument("-xmlXPath")
    CALL print_argument("-dropXMLSchema")
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"Output options:"
    CALL print_argument("-no_out")
    CALL print_argument("-minimalOutput")
    CALL print_argument("-genEnpara")
    CALL print_argument("-h")
    call print_argument("-version")
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
    CALL print_argument("-add_name")
    CALL print_argument("-fft")
    CALL print_argument("-diag")
    CALL print_argument("-eig")
    CALL print_argument("-disable_progress_thread")
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"Options useful for debugging:"
    CALL print_argument("-warn_only")
    CALL print_argument("-trace")
    CALL print_argument("-debugtime")
    CALL print_argument("-all_times")
#ifdef CPP_HDF
    WRITE(*,'(a)')""
    WRITE(*,'(a)')"HDF density file relevant options:"

    CALL print_argument("-no_cdn_hdf")
    CALL print_argument("-last_extra")
    CALL print_argument("-sd")
    CALL print_argument("-delden")
    CALL print_argument("-storeSF")
#endif
    WRITE (*,'(a)') "Options for privacy sensitive users"
    CALL print_argument("-no_send")

    WRITE(*,'(a)')""
    WRITE(*,'(a)')"Please check the documentation on www.flapw.de for more details."

    CALL juDFT_end("",l_endXML=.FALSE.) !No message so do a not print more on exit
  END SUBROUTINE fleur_help

  SUBROUTINE add_fleur_arguments()
    USE m_check_arguments

    CALL new_argument(1,"-xmlXPath","modify the xml-xpath of the inp.xml file","")
    CALL new_argument(0,"-dropXMLSchema","Write out the default XML schema files","")
    !Control the job
    CALL new_argument(0,"-check","run in check mode, i.e. stop after init","")
    CALL new_argument(0,"-info","Print out information on recommended parallelization and available charge densities","")
    CALL new_argument(2,"-wtime","run for # minutes (used to estimate if another iteration is started)","")
    CALL new_argument(1,"-j","Distribute MPI ranks to run subjobs (Format PE:DIR meaning run with PE in directory DIR)","")
    CALL new_argument(1,"-f","Obtain info on subjobs from file","")
    CALL new_argument(2,"-n_min_size","Try to use at least specified number of PE in eigenvalue parallelization","")
    CALL new_argument(1,"-add_name","file descriptor that will be put before the output file names","")
    CALL new_argument(1,"-fft","library used for Fast Fourier Transformations","inbuilt"&
#ifdef CPP_FFTW
    //",fftw"&
#endif
#ifdef CPP_FFT_MKL
         //",mkl"&
#endif
#ifdef CPP_SPFFT
         //",spfft"&
#endif
         )
    CALL new_argument(1,"-diag","Choose method for diagonalization","lapack,debugout,lapack_singlePrec,"&
#ifdef CPP_SCALAPACK
         //",scalapack"&
#endif
#ifdef CPP_ELPA_ONENODE
         //",elpa_1node"&
#endif
#ifdef CPP_ELSI
      //",elsielpa"&
#endif
#ifdef CPP_ELSI
      //",elsichase"&
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
#ifdef CPP_CUSOLVER
       //",cusolver"&
#endif
       )
    CALL new_argument(1,"-eig","Method for storing the eigenvectors","mem,da"&
#ifdef CPP_MPI
         //",mpi"&
#endif
#ifdef CPP_HDF
         //",hdf"&
#endif
         )
    CALL new_argument(0,"-serial_diag","Use serial diagonalization in parallel code","")
    CALL new_argument(0,"-disable_progress_thread","Do not use progress_thread","")
    !Debugging
    CALL new_argument(0,"-warn_only","Continue execution after a warning message","")
    CALL new_argument(0,"-trace","Try to generate a stacktrace in case of an error","")
    CALL new_argument(0,"-debugtime","Write the start/stop of all timers to the console","")
    CALL new_argument(0,"-all_times","Write json files of timing for all PE, not only for PE=0","")
    !Output
    CALL new_argument(0,"-mix_io","Do not store mixing history in memory but do IO in each iteration","")
    CALL new_argument(0,"-no_out","Do not open the 'out' file but write to stdout","")
    CALL new_argument(0,"-minimalOutput","Reduce the amount of output in the out.xml file","")
    CALL new_argument(0,"-genEnpara","Generate an 'enpara' file for the energy parameters","")
    CALL new_argument(0,"-noco","write out noco parameters in all outputs for inp.xml","")
    CALL new_argument(0,"-greensf","write out green's function parameters in all outputs for inp.xml","")
    CALL new_argument(0,"-h","Print this message","")
    CALL new_argument(0,"-version","Show important version information about the fleur executable","")
    CALL new_argument(0,"-no_send","Do not send usage data","")
    !HDF density
    CALL new_argument(0,"-no_cdn_hdf","Disable HDF charge density mode (activated by default if HDF5 is available)","")
    CALL new_argument(0,"-last_extra","Generate an additional file cdn_last.hdf that contains only the last density","")
    CALL new_argument(2,"-sd","use starting density N, where N is the index of the density according to -info","")
    CALL new_argument(1,"-delden","delete densities (either an index N, a range N-M or the keyword 'allbutlast' should be given)","")
    CALL new_argument(0,"-storeSF","Store the stepfunction in the charge density file to avoid its recalculation when restarting Fleur.","")
    !GPU parameter
    CALL new_argument(0,"-gpu","Use GPU for computing","")
  END SUBROUTINE add_fleur_arguments
END MODULE m_fleur_help
