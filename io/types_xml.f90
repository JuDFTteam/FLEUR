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
MODULE m_types_xml

  USE m_juDFT

  PRIVATE

  TYPE t_xml
     INTEGER:: id
   CONTAINS
     PROCEDURE,NOPASS :: InitInterface
     PROCEDURE,NOPASS :: ParseSchema
     PROCEDURE,NOPASS :: ParseDoc 
     PROCEDURE,NOPASS :: ValidateDoc 
     PROCEDURE,NOPASS :: InitXPath
     PROCEDURE,NOPASS :: GetNumberOfNodes  
     PROCEDURE,NOPASS :: SetAttributeValue  
     PROCEDURE,NOPASS :: GetAttributeValue  
     PROCEDURE,NOPASS :: FreeResources
  END TYPE t_xml
  PUBLIC t_xml

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
          CALL SetAttributeValue(xpath(:i-1),xpath(i+1:ii-1))
          WRITE(*,*) "Set from command line:",TRIM(xpath(:i-1)),"=",TRIM(xpath(i+1:ii-1))
          IF (ii+1<LEN(xpath))THEN
             xpath=xpath(ii+1:)
          ELSE
             xpath=""
          ENDIF
       END DO
    END IF
  END SUBROUTINE init_from_command_line


  SUBROUTINE InitInterface()

    USE iso_c_binding
    USE m_types

    IMPLICIT NONE

    INTEGER :: errorStatus

    INTERFACE
       FUNCTION initializeXMLInterface() BIND(C, name="initializeXMLInterface")
         USE iso_c_binding
         INTEGER(c_int) initializeXMLInterface
       END FUNCTION initializeXMLInterface
    END INTERFACE

    errorStatus = 0
    errorStatus = initializeXMLInterface()
    IF(errorStatus.NE.0) THEN
       CALL juDFT_error("Could not initialize XML interface.",calledby="xmlInitInterface")
    END IF

  END SUBROUTINE InitInterface

  SUBROUTINE ParseSchema(schemaFilename)

    USE iso_c_binding
    USE m_types

    IMPLICIT NONE

    CHARACTER(LEN=200,KIND=c_char), INTENT(IN) :: schemaFilename

    INTEGER :: errorStatus

    INTERFACE
       FUNCTION parseXMLSchema(schemaFilename) BIND(C, name="parseXMLSchema")
         USE iso_c_binding
         INTEGER(c_int) parseXMLSchema
         CHARACTER(kind=c_char) :: schemaFilename(*)
       END FUNCTION parseXMLSchema
    END INTERFACE

    errorStatus = 0
    errorStatus = parseXMLSchema(schemaFilename)
    IF(errorStatus.NE.0) THEN
       CALL juDFT_error("XML Schema file not parsable: "//TRIM(ADJUSTL(schemaFilename)),calledby="xmlParseSchema")
    END IF

  END SUBROUTINE ParseSchema

  SUBROUTINE ParseDoc(docFilename)

    USE iso_c_binding
    USE m_types

    IMPLICIT NONE

    CHARACTER(LEN=200,KIND=c_char), INTENT(IN) :: docFilename

    INTEGER :: errorStatus

    INTERFACE
       FUNCTION parseXMLDocument(docFilename) BIND(C, name="parseXMLDocument")
         USE iso_c_binding
         INTEGER(c_int) parseXMLDocument
         CHARACTER(kind=c_char) :: docFilename(*)
       END FUNCTION parseXMLDocument
    END INTERFACE

    errorStatus = 0
    errorStatus = parseXMLDocument(docFilename)
    IF(errorStatus.NE.0) THEN
       CALL juDFT_error("XML document file not parsable: "//TRIM(ADJUSTL(docFilename)),calledby="xmlParseDoc")
    END IF

  END SUBROUTINE ParseDoc

  SUBROUTINE ValidateDoc()

    USE iso_c_binding
    USE m_types

    IMPLICIT NONE

    INTEGER :: errorStatus

    INTERFACE
       FUNCTION validateXMLDocument() BIND(C, name="validateXMLDocument")
         USE iso_c_binding
         INTEGER(c_int) validateXMLDocument
       END FUNCTION validateXMLDocument
    END INTERFACE

    errorStatus = 0
    errorStatus = validateXMLDocument()
    IF(errorStatus.NE.0) THEN
       CALL juDFT_error("XML document cannot be validated against Schema.",calledby="xmlValidateDoc")
    END IF

  END SUBROUTINE ValidateDoc

  SUBROUTINE InitXPath()

    USE iso_c_binding
    USE m_types

    IMPLICIT NONE

    INTEGER :: errorStatus

    INTERFACE
       FUNCTION initializeXPath() BIND(C, name="initializeXPath")
         USE iso_c_binding
         INTEGER(c_int) initializeXPath
       END FUNCTION initializeXPath
    END INTERFACE

    errorStatus = 0
    errorStatus = initializeXPath()
    IF(errorStatus.NE.0) THEN
       CALL juDFT_error("Could not initialize XPath.",calledby="InitXPath")
    END IF
    CALL init_from_command_line()
  END SUBROUTINE InitXPath

  FUNCTION GetNumberOfNodes(xPath)

    USE iso_c_binding
    USE m_types

    IMPLICIT NONE

    INTEGER :: GetNumberOfNodes
    CHARACTER(LEN=*,KIND=c_char), INTENT(IN) :: xPath

    INTERFACE
       FUNCTION getNumberOfXMLNodes(xPathExpression) BIND(C, name="getNumberOfXMLNodes")
         USE iso_c_binding
         INTEGER(c_int) getNumberOfXMLNodes
         CHARACTER(kind=c_char) :: xPathExpression(*)
       END FUNCTION getNumberOfXMLNodes
    END INTERFACE

    GetNumberOfNodes = getNumberOfXMLNodes(TRIM(ADJUSTL(xPath))//C_NULL_CHAR)

  END FUNCTION GetNumberOfNodes

  FUNCTION GetAttributeValue(xPath)

    USE iso_c_binding
    USE m_types

    IMPLICIT NONE

    CHARACTER(LEN=:),ALLOCATABLE :: GetAttributeValue

    CHARACTER(LEN=*, KIND=c_char), INTENT(IN) :: xPath

    CHARACTER (LEN=1, KIND=c_char), POINTER, DIMENSION (:) :: valueFromC => NULL()
    CHARACTER*255 :: VALUE
    INTEGER :: length, errorStatus, i
    TYPE(c_ptr) :: c_string

    INTERFACE
       FUNCTION getXMLAttributeValue(xPathExpression) BIND(C, name="getXMLAttributeValue")
         USE iso_c_binding
         CHARACTER(KIND=c_char) :: xPathExpression(*)
         TYPE(c_ptr) :: getXMLAttributeValue
       END FUNCTION getXMLAttributeValue
    END INTERFACE



    c_string = getXMLAttributeValue(TRIM(ADJUSTL(xPath))//C_NULL_CHAR)

    CALL C_F_POINTER(c_string, valueFromC, [ 255 ])
    IF (.NOT.C_ASSOCIATED(c_string)) THEN
       WRITE(*,*) 'Error in trying to obtain attribute value from XPath:'
       WRITE(*,*) TRIM(ADJUSTL(xPath))
       CALL juDFT_error("Attribute value could not be obtained.",calledby="xmlGetAttributeValue")
    END IF

    VALUE = ''
    i = 1
    DO WHILE ((valueFromC(i).NE.C_NULL_CHAR).AND.(i.LE.255))
       VALUE(i:i) = valueFromC(i)
       i = i + 1
    END DO
    length = i-1

    GetAttributeValue = TRIM(ADJUSTL(VALUE(1:length)))

  END FUNCTION GetAttributeValue


  SUBROUTINE SetAttributeValue(xPath,VALUE)

    USE iso_c_binding
    USE m_types

    IMPLICIT NONE

    CHARACTER(LEN=*, KIND=c_char), INTENT(IN) :: xPath
    CHARACTER(len=*, KIND=c_char), INTENT(IN) :: VALUE

    INTEGER :: errorStatus

    INTERFACE
       FUNCTION setXMLAttributeValue(xPathExpression,valueExpression) BIND(C, name="setXMLAttributeValue")
         USE iso_c_binding
         CHARACTER(KIND=c_char) :: xPathExpression(*)
         CHARACTER(KIND=c_char) :: valueExpression(*)
         INTEGER(c_int) :: setXMLAttributeValue
       END FUNCTION setXMLAttributeValue
    END INTERFACE

    errorStatus = setXMLAttributeValue(TRIM(ADJUSTL(xPath))//C_NULL_CHAR,TRIM(ADJUSTL(VALUE))//C_NULL_CHAR)
    IF (errorStatus.NE.0) THEN
       WRITE(*,*) 'Error in trying to setting attribute value from XPath:'
       WRITE(*,*) TRIM(ADJUSTL(xPath))
       WRITE(*,*) TRIM(ADJUSTL(VALUE))
       CALL juDFT_error("Attribute value could not be set.",calledby="xmlSetAttributeValue")
    END IF

  END SUBROUTINE SetAttributeValue

  SUBROUTINE FreeResources()

    USE iso_c_binding
    USE m_types

    IMPLICIT NONE

    INTEGER :: errorStatus

    INTERFACE
       FUNCTION freeXMLResources() BIND(C, name="freeXMLResources")
         USE iso_c_binding
         INTEGER freeXMLResources
       END FUNCTION freeXMLResources
    END INTERFACE

    errorStatus = 0
    errorStatus = freeXMLResources()
    IF(errorStatus.NE.0) THEN
       CALL juDFT_error("Could not free XML resources.",calledby="xmlFreeResources")
       STOP 'Error!'
    END IF

  END SUBROUTINE FreeResources

END MODULE m_types_xml
