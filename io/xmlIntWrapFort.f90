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


PRIVATE :: init_from_command_line
CONTAINS

SUBROUTINE init_from_command_line()
  IMPLICIT NONE
  CHARACTER(len=1000)::xpath
  INTEGER:: i,ii
  
  IF (judft_was_argument("-xmlXPath")) THEN
     xpath=judft_string_for_argument("-xmlXPath")
     DO WHILE(INDEX(xpath,"=")>0)
        i=INDEX(xpath,"=")
        ii=INDEX(xpath,":")
        IF (ii==0) ii=LEN(TRIM(xpath))+1
        IF (i>100.OR.ii-i>100) CALL judft_error("Too long xmlXPath argument",calledby="xmlIntWarpFort.f90")
        CALL xmlSetAttributeValue(xpath(:i-1),xpath(i+1:ii-1))
        WRITE(*,*) "Set from command line:",TRIM(xpath(:i-1)),"=",TRIM(xpath(i+1:ii-1))
        IF (ii+1<len(xpath))THEN
           xpath=xpath(ii+1:)
        ELSE
           xpath=""
        ENDIF
     END DO
  END IF
END SUBROUTINE init_from_command_line
     

SUBROUTINE xmlInitInterface()

   USE iso_c_binding

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
   CALL init_from_command_line()
END SUBROUTINE xmlInitXPath

FUNCTION xmlGetNumberOfNodes(xPath)

   USE iso_c_binding

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

   value = ''
   i = 1
   DO WHILE ((valueFromC(i).NE.C_NULL_CHAR).AND.(i.LE.255))
      value(i:i) = valueFromC(i)
      i = i + 1
   END DO
   length = i-1

   xmlGetAttributeValue = value(1:length)

END FUNCTION xmlGetAttributeValue


SUBROUTINE xmlSetAttributeValue(xPath,VALUE)

   USE iso_c_binding
 
   IMPLICIT NONE
  
   CHARACTER(LEN=*, KIND=c_char), INTENT(IN) :: xPath
   CHARACTER(len=*, KIND=c_char), INTENT(IN) :: value

   INTEGER :: errorStatus

   INTERFACE
      FUNCTION setXMLAttributeValue(xPathExpression,valueExpression) BIND(C, name="setXMLAttributeValue")
         use iso_c_binding
         CHARACTER(KIND=c_char) :: xPathExpression(*)
         CHARACTER(KIND=c_char) :: valueExpression(*)
         INTEGER(c_int) :: setXMLAttributeValue
       END FUNCTION setXMLAttributeValue
    END INTERFACE

    errorStatus = setXMLAttributeValue(TRIM(ADJUSTL(xPath))//C_NULL_CHAR,TRIM(ADJUSTL(VALUE))//C_NULL_CHAR)
    IF (errorStatus.ne.0) THEN
      WRITE(*,*) 'Error in trying to setting attribute value from XPath:'
      WRITE(*,*) TRIM(ADJUSTL(xPath))
      WRITE(*,*) TRIM(ADJUSTL(VALUE))
      CALL juDFT_error("Attribute value could not be set.",calledby="xmlSetAttributeValue")
   END IF

 END SUBROUTINE xmlSetAttributeValue

SUBROUTINE xmlFreeResources()

   USE iso_c_binding
 
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
