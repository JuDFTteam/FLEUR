module m_inpgen_version
    implicit none
    private
    public inpgen_version
  contains
    subroutine inpgen_version()
      use m_compile_descr
      use m_constants
      use m_juDFT
      use m_check_arguments
      use m_types_xml
  
      character(:), allocatable:: infostring, additional_info
      character(len=10) :: outputVersionString
      integer :: omp
      type(t_xml)::xml

      if (.NOT. judft_was_argument("-version")) return

      !now print version info and help on command line arguments:
      call get_compile_desc_string(infostring)
      write(*,'(a)') infostring
      write(outputVersionString,'(a,i0)') '0.', xml%currentversionNumber
      additional_info = "Default XML file version :  "//TRIM(outputVersionString)
      write(*,'(a)') additional_info

      CALL juDFT_end("Version output completed")
    end subroutine

end module