!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_fleur_dropxmlschema
    implicit none
    private
    public fleur_dropxmlschema, fleur_write_xmlschema_files
  contains
    subroutine fleur_write_xmlschema_files() BIND(C, name="drop_xmlschema")
      use m_types_xml
      use iso_c_binding

      CHARACTER(LEN=200, KIND=c_char) :: versionString
      type(t_xml)::xml

      write(versionString,'(a,i0)') '0.', xml%currentversionNumber
      call drop_schema_files(versionString, versionString)
    end subroutine

    subroutine fleur_dropxmlschema()
      use m_juDFT

      if (.NOT. judft_was_argument("-dropXMLSchema")) return
      call fleur_write_xmlschema_files()

      CALL juDFT_end("XML Schema files written")
    end subroutine

end module