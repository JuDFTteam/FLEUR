module m_fleur_version
    implicit none
    private
    public fleur_version
  contains
    subroutine fleur_version()
      use m_compile_descr
      use m_constants
      use m_juDFT
      use m_check_arguments
      use m_types_xml
!$    use omp_lib
  
      character(:), allocatable:: infostring, additional_info, omp_string
      character(len=10) :: outputVersionString
      integer :: omp
      type(t_xml)::xml

      omp = 1
      if (.NOT. judft_was_argument("-version")) return

      !now print version info and help on command line arguments:
      call get_compile_desc_string(infostring)
      write(*,'(a)') infostring
      write(outputVersionString,'(a,i0)') '0.', xml%currentversionNumber

      !$ omp=omp_get_max_threads()
      if (omp==1) then
        omp_string = "   OpenMP  : False"
      else
        omp_string = "   OpenMP  : True"
      endif
      additional_info = "Libraries:"//new_LINE("a")// &
#ifdef CPP_MPI
                      "   MPI     : True"//new_LINE("a")// &
#else
                      "   MPI     : False"//new_LINE("a")// &
#endif

#ifdef _OPENACC
                      "   OpenMP  : False"//new_LINE("a")// &
                      "   OpenACC : True"//new_LINE("a")// &
#else
                      trim(omp_string)//new_LINE("a")// &
                      "   OpenACC : False"//new_LINE("a")// &
#endif
#ifdef CPP_HDF
                      "   HDF5    : True"//new_LINE("a")// &
#else
                      "   HDF5    : False"//new_LINE("a")// &
#endif
#ifdef CPP_LIBXC
                      "   LibXC   : True"//new_LINE("a")// &
#else
                      "   LibXC   : False"//new_LINE("a")// &
#endif
                      "   Diagonalization: lapack" &
#ifdef CPP_SCALAPACK
                      //",scalapack"&
#endif
#ifdef CPP_ELPA_ONENODE
                      //",elpa_1node"&
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
#ifdef CPP_GPU
                      //",cusolver"&
#endif
                      //new_LINE("a")//&
                      "Default XML file version :  "//TRIM(outputVersionString)

      write(*,'(a)') additional_info

      CALL juDFT_end("Version output completed")
    end subroutine

end module