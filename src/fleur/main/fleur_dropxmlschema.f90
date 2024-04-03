module m_fleur_dropxmlschema
    implicit none
    private
    public fleur_dropxmlschema
  contains
    subroutine fleur_dropxmlschema()
      use m_juDFT
      use m_types_xml
      use iso_c_binding

      CHARACTER(LEN=200, KIND=c_char) :: versionString
      type(t_xml)::xml

      if (.NOT. judft_was_argument("-dropXMLSchema")) return
      write(versionString,'(a,i0)') '0.', xml%currentversionNumber
      call drop_schema_files(versionString, versionString)

      CALL juDFT_end("XML Schema files written")
    end subroutine

end module