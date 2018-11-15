!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------


MODULE m_inp_xml
  USE m_judft
CONTAINS
  CHARACTER (len=200) FUNCTION inp_xml_speciesxpath_for_group(n)
    IMPLICIT NONE
    INTEGER,INTENT(in)::n

    CHARACTER(len=200)::xpath,species
    !First determine name of species from group
    WRITE(xPath,*) '/fleurInput/atomGroups/atomGroup[',n,']/@species'
    species=TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPath)))))

    DO i=1,xmlGetNumberOfNodes('/fleurInput/atomSpecies/species')
       WRITE(xPath,*) '/fleurInput/atomSpecies/species[',i,']'
       IF (TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@name')))==TRIM(species)) THEN
          inp_xml_speciesxpath_for_group=xpath
          RETURN
       END IF
    END DO
    WRITE(xpath,*) n
    CALL judft_error("No species found for name "//TRIM(species)//" used in atom group "//TRIM(xpath))
  END FUNCTION inp_xml_speciesxpath_for_group

  CHARACTER (len=200) FUNCTION inp_xml_xpath_for_group(n)
    IMPLICIT NONE
    INTEGER,INTENT(in)::n
    WRITE(inp_xml_xpath_for_group,*) '/fleurInput/atomGroups/atomGroup[',n,']'
  END FUNCTION inp_xml_xpath_for_group
END MODULE m_inp_xml
