MODULE m_xmlOutput

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!   XML output service routines
!!!
!!!   This module provides several subroutines that simplify the
!!!   generation of the out.xml file.
!!!                                         GM'16
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   IMPLICIT NONE

   PRIVATE
   INTEGER, SAVE :: currentElementIndex
   INTEGER, SAVE :: maxNumElements
   INTEGER, SAVE :: xmlOutputUnit
   CHARACTER(LEN= 40), ALLOCATABLE :: elementList(:)

   PUBLIC startXMLOutput, endXMLOutput
   PUBLIC writeXMLElementFormPoly, writeXMLElementPoly
   PUBLIC writeXMLElementForm, writeXMLElement
   PUBLIC openXMLElementPoly, openXMLElementNoAttributes, openXMLElement, closeXMLElement
   PUBLIC getXMLOutputUnitNumber

   CONTAINS

   FUNCTION getXMLOutputUnitNumber()
      IMPLICIT NONE
      INTEGER getXMLOutputUnitNumber
      getXMLOutputUnitNumber = xmlOutputUnit
   END FUNCTION getXMLOutputUnitNumber

   SUBROUTINE startXMLOutput()
      IMPLICIT NONE
      maxNumElements = 10
      ALLOCATE(elementList(maxNumElements))
      elementList = ''
      currentElementIndex = 0
      xmlOutputUnit = 53
      OPEN (xmlOutputUnit,file='out.xml',form='formatted',status='unknown')
      WRITE (xmlOutputUnit,'(a)') '<?xml version="1.0" encoding="UTF-8" standalone="no"?>'
      WRITE (xmlOutputUnit,'(a)') '<fleurOutput fleurOutputVersion="0.27">'
      ! TODO: Write out standard stuff
   END SUBROUTINE startXMLOutput

   SUBROUTINE endXMLOutput()
      IMPLICIT NONE
      DO WHILE (currentElementIndex.NE.0)
         CALL closeXMLElement(elementList(currentElementIndex))
      END DO
      DEALLOCATE(elementList)
      WRITE (xmlOutputUnit,'(a)') '</fleurOutput>'
      CLOSE(xmlOutputUnit)
   END SUBROUTINE endXMLOutput

   SUBROUTINE writeXMLElementFormPoly(elementName,attributeNames,attributeValues,lengths,contentList)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)           :: elementName
      CHARACTER(LEN=*), INTENT(IN)           :: attributeNames(:)
      CLASS(*),         INTENT(IN)           :: attributeValues(:)
      INTEGER,          INTENT(IN)           :: lengths(:,:)
      CLASS(*),         INTENT(IN), OPTIONAL :: contentList(:)

      CHARACTER(LEN= 30), ALLOCATABLE :: charAttributeValues(:)
      CHARACTER(LEN= 30), ALLOCATABLE :: charContentList(:)
      INTEGER                         :: i

      ALLOCATE(charAttributeValues(SIZE(attributeValues)))
      IF (PRESENT(contentList)) THEN
         ALLOCATE(charContentList(SIZE(contentList)))
      END IF

      charAttributeValues = ''
      charContentList = ''

      DO i = 1, SIZE(attributeValues)
         SELECT TYPE (attributeValues)
            TYPE IS(INTEGER)
               WRITE(charAttributeValues(i),'(i0)'), attributeValues(i)
            TYPE IS(REAL)
               WRITE(charAttributeValues(i),'(f20.10)'), attributeValues(i)
            TYPE IS(LOGICAL)
               WRITE(charAttributeValues(i),'(l1)'), attributeValues(i)
            CLASS DEFAULT
               STOP 'Type of attributeValues not allowed'
         END SELECT
      END DO

      IF (PRESENT(contentList)) THEN
         DO i = 1, SIZE(contentList)
            SELECT TYPE(contentList)
               TYPE IS(INTEGER)
                  WRITE(charContentList(i),'(i0)'), contentList(i)
               TYPE IS(REAL)
                  WRITE(charContentList(i),'(f20.10)'), contentList(i)
               TYPE IS(LOGICAL)
                  WRITE(charContentList(i),'(l1)'), contentList(i)
               CLASS DEFAULT
                  STOP 'Type of contentList not allowed'
            END SELECT
         END DO
         CALL writeXMLElementForm(elementName,attributeNames,charAttributeValues,lengths,charContentList)
         DEALLOCATE(charContentList,charAttributeValues)
      ELSE
         CALL writeXMLElementForm(elementName,attributeNames,charAttributeValues,lengths)
      END IF

   END SUBROUTINE writeXMLElementFormPoly

   SUBROUTINE writeXMLElementPoly(elementName,attributeNames,attributeValues,contentList)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)           :: elementName
      CHARACTER(LEN=*), INTENT(IN)           :: attributeNames(:)
      CLASS(*),         INTENT(IN)           :: attributeValues(:)
      CLASS(*),         INTENT(IN), OPTIONAL :: contentList(:)

      INTEGER, ALLOCATABLE :: lengths(:,:)
      INTEGER              :: contentListSize

      contentListSize = 0
      IF (PRESENT(contentList)) THEN
         contentListSize = SIZE(contentList)
      END IF

      ALLOCATE(lengths(MAX(SIZE(attributeNames),contentListSize),3))
      lengths = 0
      CALL writeXMLElementFormPoly(elementName,attributeNames,attributeValues,lengths,contentList)
      DEALLOCATE(lengths)
   END SUBROUTINE writeXMLElementPoly

   SUBROUTINE writeXMLElementForm(elementName,attributeNames,attributeValues,lengths,contentList)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)           :: elementName
      CHARACTER(LEN=*), INTENT(IN)           :: attributeNames(:)
      CHARACTER(LEN=*), INTENT(IN)           :: attributeValues(:)
      INTEGER,          INTENT(IN)           :: lengths(:,:)
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: contentList(:)

      CHARACTER(LEN=200), ALLOCATABLE :: contentLineList(:)
      INTEGER, ALLOCATABLE            :: bigLengths(:,:)
      INTEGER :: i, j, contentLineLength, contentLineListSize
      CHARACTER(LEN=70)               :: format
      CHARACTER(LEN=200)              :: outputString
      INTEGER                         :: contentListSize, overallListSize
      INTEGER                         :: lengthsShape(2)

      IF(SIZE(attributeNames).NE.SIZE(attributeValues)) THEN
         WRITE(*,*) 'attributeNames', attributeNames
         WRITE(*,*) 'attributeValues', attributeValues
         STOP 'ERROR: SIZE(attributeNames).NE.SIZE(attributeValues)'
      END IF

      lengthsShape = SHAPE(lengths)
      contentListSize = 0
      IF (PRESENT(contentList)) THEN
         contentListSize = SIZE(contentList)
      END IF
      overallListSize = MAX(SIZE(attributeNames),contentListSize)
      ALLOCATE(bigLengths(overallListSize,2))
      bigLengths = 0
      DO j = 1, 2
         DO i = 1, MIN(overallListSize,lengthsShape(1))
            bigLengths(i,j) = lengths(i,j)
         END DO
      END DO

      outputString = '<'//TRIM(ADJUSTL(elementName))
      DO i = 1, SIZE(attributeNames)
         WRITE(format,'(a)') "(a,a,a"
         IF (bigLengths(i,1).GT.0) WRITE(format,'(a,i0)') TRIM(ADJUSTL(format)),bigLengths(i,1)
         WRITE(format,'(a,a)') TRIM(ADJUSTL(format)),",a,a"
         IF (bigLengths(i,2).GT.0) WRITE(format,'(a,i0)') TRIM(ADJUSTL(format)),bigLengths(i,2)
         WRITE(format,'(a,a)') TRIM(ADJUSTL(format)),",a)"
         WRITE(*,*) TRIM(ADJUSTL(format))
         WRITE(outputString,format) TRIM(ADJUSTL(outputString)), ' ', TRIM(ADJUSTL(attributeNames(i))),&
                                    '="', TRIM(ADJUSTL(attributeValues(i))), '"'
!         outputString = TRIM(ADJUSTL(outputString))//' '//TRIM(ADJUSTL(attributeNames(i)))//'="'//TRIM(ADJUSTL(attributeValues(i)))//'"'
      END DO
      WRITE(format,'(a,i0,a)') "(a",3*(currentElementIndex+1),",a)"
      IF (PRESENT(contentList)) THEN
         contentLineLength = 5 ! At most 5 data elements per line
         contentLineListSize = SIZE(contentList) / contentLineLength
         IF(contentLineListSize*contentLineLength.NE.SIZE(contentList)) THEN
            contentLineListSize = contentLineListSize + 1
         END IF
         ALLOCATE(contentLineList(contentLineListSize))
         CALL fillContentLineList(contentList,contentLineList,contentLineLength)
         IF(SIZE(contentLineList).LE.1) THEN
            outputString = TRIM(ADJUSTL(outputString))//'>'//contentLineList(1)//'</'//TRIM(ADJUSTL(elementName))//'>'
         ELSE
            outputString = TRIM(ADJUSTL(outputString))//'>'
         END IF
         WRITE(xmlOutputUnit,format) ' ', TRIM(ADJUSTL(outputString))
         IF(SIZE(contentLineList).GT.1) THEN
            WRITE(format,'(a,i0,a)') "(a",3*(currentElementIndex+2),",a)"
            DO i = 1, SIZE(contentLineList)
               WRITE(xmlOutputUnit,format) ' ', TRIM(ADJUSTL(contentLineList(i)))
            END DO
            WRITE(format,'(a,i0,a)') "(a",3*(currentElementIndex+1),",a)"
            outputString = '</'//TRIM(ADJUSTL(elementName))//'>'
            WRITE(xmlOutputUnit,format) ' ', TRIM(ADJUSTL(outputString))
         END IF
      ELSE
         outputString = TRIM(ADJUSTL(outputString))//'/>'
         WRITE(xmlOutputUnit,format) ' ', TRIM(ADJUSTL(outputString))
      END IF
      DEALLOCATE (bigLengths)
   END SUBROUTINE writeXMLElementForm

   SUBROUTINE writeXMLElement(elementName,attributeNames,attributeValues,contentList)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)           :: elementName
      CHARACTER(LEN=*), INTENT(IN)           :: attributeNames(:)
      CHARACTER(LEN=*), INTENT(IN)           :: attributeValues(:)
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: contentList(:)

      INTEGER, ALLOCATABLE :: lengths(:,:)
      INTEGER              :: contentListSize

      contentListSize = 0
      IF (PRESENT(contentList)) THEN
         contentListSize = SIZE(contentList)
      END IF

      ALLOCATE(lengths(MAX(SIZE(attributeNames),contentListSize),2))
      lengths = 0
      CALL writeXMLElementForm(elementName,attributeNames,attributeValues,lengths,contentList)
      DEALLOCATE(lengths)
   END SUBROUTINE writeXMLElement
   

   SUBROUTINE fillContentLineList(contentList,contentLineList,contentLineLength)

      IMPLICIT NONE

      CHARACTER(LEN= 30), INTENT(IN)    :: contentList(:)
      CHARACTER(LEN=200), INTENT(INOUT) :: contentLineList(:)
      INTEGER,            INTENT(IN)    :: contentLineLength
      INTEGER :: i, j

      contentLineList = ''
      DO i = 1, SIZE(contentLineList)
         DO j = 1, contentLineLength
            IF((i-1)*contentLineLength+j.GT.SIZE(contentList)) THEN
               RETURN
            END IF
            contentLineList(i) = TRIM(ADJUSTL(contentLineList(i)))//' '//TRIM(ADJUSTL(contentList((i-1)*contentLineLength+j)))
         END DO
      END DO
   END SUBROUTINE fillContentLineList

   SUBROUTINE openXMLElementPoly(elementName,attributeNames,attributeValues)

      IMPLICIT NONE

      CHARACTER(LEN=*)             :: elementName
      CHARACTER(LEN=*), INTENT(IN) :: attributeNames(:)
      CLASS(*),         INTENT(IN) :: attributeValues(:)

      CHARACTER(LEN= 30), ALLOCATABLE :: charAttributeValues(:)
      INTEGER                         :: i

      ALLOCATE(charAttributeValues(SIZE(attributeValues)))

      charAttributeValues = ''

      DO i = 1, SIZE(attributeValues)
         SELECT TYPE (attributeValues)
            TYPE IS(INTEGER)
               WRITE(charAttributeValues(i),'(i0)'), attributeValues(i)
            TYPE IS(REAL)
               WRITE(charAttributeValues(i),'(f20.10)'), attributeValues(i)
            TYPE IS(LOGICAL)
               WRITE(charAttributeValues(i),'(l1)'), attributeValues(i)
            CLASS DEFAULT
               STOP 'Type of attributeValues not allowed'
         END SELECT
      END DO

      CALL openXMLElement(elementName,attributeNames,charAttributeValues)

      DEALLOCATE(charAttributeValues)

   END SUBROUTINE openXMLElementPoly

   SUBROUTINE openXMLElementNoAttributes(elementName)

      IMPLICIT NONE

      CHARACTER(LEN=*)             :: elementName

      INTEGER :: i
      CHARACTER(LEN=70)  :: format
      CHARACTER(LEN=200) :: openingString

      IF(currentElementIndex.EQ.maxNumElements) THEN
         WRITE(*,*) 'elementName', TRIM(ADJUSTL(elementName))
         STOP 'ERROR: xml hierarchy too deep!'
      END IF
      currentElementIndex = currentElementIndex + 1
      elementList(currentElementIndex) = TRIM(ADJUSTL(elementName))
      openingString = '<'//TRIM(ADJUSTL(elementName))//'>'
      WRITE(format,'(a,i0,a)') "(a",3*currentElementIndex,",a)"
      WRITE(xmlOutputUnit,format) ' ', TRIM(ADJUSTL(openingString))
   END SUBROUTINE openXMLElementNoAttributes

   SUBROUTINE openXMLElement(elementName,attributeNames,attributeValues)

      IMPLICIT NONE

      CHARACTER(LEN=*)             :: elementName
      CHARACTER(LEN=*), INTENT(IN) :: attributeNames(:)
      CHARACTER(LEN=*), INTENT(IN) :: attributeValues(:)

      INTEGER :: i
      CHARACTER(LEN=70)  :: format
      CHARACTER(LEN=200) :: openingString

      IF(SIZE(attributeNames).NE.SIZE(attributeValues)) THEN
         WRITE(*,*) 'elementName', TRIM(ADJUSTL(elementName))
         WRITE(*,*) 'attributeNames', attributeNames
         WRITE(*,*) 'attributeValues', attributeValues
         STOP 'ERROR: SIZE(attributeNames).NE.SIZE(attributeValues)'
      END IF
      IF(currentElementIndex.EQ.maxNumElements) THEN
         WRITE(*,*) 'elementName', TRIM(ADJUSTL(elementName))
         WRITE(*,*) 'attributeNames', attributeNames
         WRITE(*,*) 'attributeValues', attributeValues
         STOP 'ERROR: xml hierarchy too deep!'
      END IF
      currentElementIndex = currentElementIndex + 1
      elementList(currentElementIndex) = TRIM(ADJUSTL(elementName))
      openingString = '<'//TRIM(ADJUSTL(elementName))
      DO i = 1, SIZE(attributeNames)
         openingString = TRIM(ADJUSTL(openingString))//' '//TRIM(ADJUSTL(attributeNames(i)))//'="'//TRIM(ADJUSTL(attributeValues(i)))//'"'
      END DO
      openingString = TRIM(ADJUSTL(openingString))//'>'
      WRITE(format,'(a,i0,a)') "(a",3*currentElementIndex,",a)"
      WRITE(xmlOutputUnit,format) ' ', TRIM(ADJUSTL(openingString))
   END SUBROUTINE openXMLElement

   SUBROUTINE closeXMLElement(elementName)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: elementName
      INTEGER :: i
      CHARACTER(LEN=70) :: format
      CHARACTER(LEN=70) :: closingString

      IF(TRIM(ADJUSTL(elementList(currentElementIndex))).NE.TRIM(ADJUSTL(elementName))) THEN
         WRITE(*,*) 'elementList(currentElementIndex):', TRIM(ADJUSTL(elementList(currentElementIndex)))
         WRITE(*,*) 'elementName', TRIM(ADJUSTL(elementName))
         STOP 'ERROR: Closing xml element inconsistency!'
      END IF
      closingString = '</'//TRIM(ADJUSTL(elementName))//'>'
      WRITE(format,'(a,i0,a)') "(a",3*currentElementIndex,",a)"
      WRITE(xmlOutputUnit,format) ' ', TRIM(ADJUSTL(closingString))
      currentElementIndex = currentElementIndex - 1
   END SUBROUTINE closeXMLElement

END MODULE m_xmlOutput
