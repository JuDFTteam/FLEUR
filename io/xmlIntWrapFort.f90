!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Wrapper routines for XML IO - Fortran side
!
!                                   GM'16
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE m_xmlIntWrapFort

USE m_juDFT

CONTAINS
SUBROUTINE xmlInitInterface()

   USE iso_c_binding
   USE m_types

   IMPLICIT NONE

   INTEGER :: errorStatus

   interface
      function initializeXMLInterface() bind(C, name="initializeXMLInterface")
         use iso_c_binding
         INTEGER(c_int) initializeXMLInterface
      end function initializeXMLInterface
   end interface

   errorStatus = 0
   errorStatus = initializeXMLInterface()
   IF(errorStatus.NE.0) THEN
      CALL juDFT_error("Could not initialize XML interface.",calledby="xmlInitInterface")
   END IF

END SUBROUTINE xmlInitInterface

SUBROUTINE xmlParseSchema(schemaFilename)

   USE iso_c_binding
   USE m_types

   IMPLICIT NONE

   CHARACTER(LEN=200,KIND=c_char), INTENT(IN) :: schemaFilename

   INTEGER :: errorStatus

   interface
      function parseXMLSchema(schemaFilename) bind(C, name="parseXMLSchema")
         use iso_c_binding
         INTEGER(c_int) parseXMLSchema
         character(kind=c_char) :: schemaFilename(*)
      end function parseXMLSchema
   end interface

   errorStatus = 0
   errorStatus = parseXMLSchema(schemaFilename)
   IF(errorStatus.NE.0) THEN
      CALL juDFT_error("XML Schema file not parsable: "//TRIM(ADJUSTL(schemaFilename)),calledby="xmlParseSchema")
   END IF

END SUBROUTINE xmlParseSchema

SUBROUTINE xmlParseDoc(docFilename)

   USE iso_c_binding
   USE m_types

   IMPLICIT NONE

   CHARACTER(LEN=200,KIND=c_char), INTENT(IN) :: docFilename

   INTEGER :: errorStatus

   interface
      function parseXMLDocument(docFilename) bind(C, name="parseXMLDocument")
         use iso_c_binding
         INTEGER(c_int) parseXMLDocument
         character(kind=c_char) :: docFilename(*)
      end function parseXMLDocument
   end interface

   errorStatus = 0
   errorStatus = parseXMLDocument(docFilename)
   IF(errorStatus.NE.0) THEN
      CALL juDFT_error("XML document file not parsable: "//TRIM(ADJUSTL(docFilename)),calledby="xmlParseDoc")
   END IF

END SUBROUTINE xmlParseDoc

SUBROUTINE xmlValidateDoc()

   USE iso_c_binding
   USE m_types

   IMPLICIT NONE

   INTEGER :: errorStatus

   interface
      function validateXMLDocument() bind(C, name="validateXMLDocument")
         use iso_c_binding
         INTEGER(c_int) validateXMLDocument
      end function validateXMLDocument
   end interface

   errorStatus = 0
   errorStatus = validateXMLDocument()
   IF(errorStatus.NE.0) THEN
      CALL juDFT_error("XML document cannot be validated against Schema.",calledby="xmlValidateDoc")
   END IF

END SUBROUTINE xmlValidateDoc

SUBROUTINE xmlInitXPath()

   USE iso_c_binding
   USE m_types

   IMPLICIT NONE

   INTEGER :: errorStatus

   interface
      function initializeXPath() bind(C, name="initializeXPath")
         use iso_c_binding
         INTEGER(c_int) initializeXPath
      end function initializeXPath
   end interface

   errorStatus = 0
   errorStatus = initializeXPath()
   IF(errorStatus.NE.0) THEN
      CALL juDFT_error("Could not initialize XPath.",calledby="xmlInitXPath")
   END IF

END SUBROUTINE xmlInitXPath

FUNCTION xmlGetNumberOfNodes(xPath)

   USE iso_c_binding
   USE m_types

   IMPLICIT NONE

   INTEGER :: xmlGetNumberOfNodes
   CHARACTER(LEN=*,KIND=c_char), INTENT(IN) :: xPath

   interface
      function getNumberOfXMLNodes(xPathExpression) bind(C, name="getNumberOfXMLNodes")
         use iso_c_binding
         INTEGER(c_int) getNumberOfXMLNodes
         character(kind=c_char) :: xPathExpression(*)
      end function getNumberOfXMLNodes
   end interface

   xmlGetNumberOfNodes = getNumberOfXMLNodes(TRIM(ADJUSTL(xPath))//C_NULL_CHAR)

END FUNCTION xmlGetNumberOfNodes

FUNCTION xmlGetAttributeValue(xPath)

   USE iso_c_binding
   USE m_types

   IMPLICIT NONE

   CHARACTER(LEN=255) :: xmlGetAttributeValue

   CHARACTER(LEN=*, KIND=c_char), INTENT(IN) :: xPath

   CHARACTER (LEN=1, KIND=c_char), POINTER, DIMENSION (:) :: valueFromC => null()
   CHARACTER*255 :: value
   INTEGER :: length, errorStatus, i
   TYPE(c_ptr) :: c_string

   interface
      function getXMLAttributeValue(xPathExpression) bind(C, name="getXMLAttributeValue")
         use iso_c_binding
         CHARACTER(KIND=c_char) :: xPathExpression(*)
         TYPE(c_ptr) :: getXMLAttributeValue
      end function getXMLAttributeValue
   end interface

   c_string = getXMLAttributeValue(TRIM(ADJUSTL(xPath))//C_NULL_CHAR)

   CALL C_F_POINTER(c_string, valueFromC, [ 255 ])
   IF (.NOT.c_associated(c_string)) THEN
      WRITE(*,*) 'Error in trying to obtain attribute value from XPath:'
      WRITE(*,*) TRIM(ADJUSTL(xPath))
      CALL juDFT_error("Attribute value could not be obtained.",calledby="xmlGetAttributeValue")
   END IF
   DO i=1, 255
     value(i:i) = valueFromC(i)
   END DO
   length = LEN_TRIM(value(1:INDEX(value, CHAR(0))))

   xmlGetAttributeValue = value(1:length-1)

END FUNCTION xmlGetAttributeValue

SUBROUTINE xmlFreeResources()

   USE iso_c_binding
   USE m_types

   IMPLICIT NONE

   INTEGER :: errorStatus

   interface
      function freeXMLResources() bind(C, name="freeXMLResources")
         use iso_c_binding
         INTEGER freeXMLResources
      end function freeXMLResources
   end interface

   errorStatus = 0
   errorStatus = freeXMLResources()
   IF(errorStatus.NE.0) THEN
      CALL juDFT_error("Could not free XML resources.",calledby="xmlFreeResources")
      STOP 'Error!'
   END IF

END SUBROUTINE xmlFreeResources

END MODULE m_xmlIntWrapFort
