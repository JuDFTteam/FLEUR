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

   PUBLIC startXMLOutput, endXMLOutput, writeXMLElementPoly, writeXMLElement
   PUBLIC openXMLElementPoly, openXMLElementNoAttributes, openXMLElement, closeXMLElement

   CONTAINS

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

   SUBROUTINE writeXMLElementPoly(elementName,attributeNames,attributeValues,contentList)

      IMPLICIT NONE

      CHARACTER(LEN=*)             :: elementName
      CHARACTER(LEN= 40), INTENT(IN) :: attributeNames(:)
      CLASS(*),           INTENT(IN) :: attributeValues(:)
      CLASS(*),           INTENT(IN) :: contentList(:)

      CHARACTER(LEN= 30), ALLOCATABLE :: charAttributeValues(:)
      CHARACTER(LEN= 30), ALLOCATABLE :: charContentList(:)
      INTEGER                         :: i

      ALLOCATE(charAttributeValues(SIZE(attributeValues)))
      ALLOCATE(charContentList(SIZE(contentList)))

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

      CALL writeXMLElement(elementName,attributeNames,charAttributeValues,charContentList)

      DEALLOCATE(charContentList,charAttributeValues)
   END SUBROUTINE writeXMLElementPoly

   SUBROUTINE writeXMLElement(elementName,attributeNames,attributeValues,contentList)

      IMPLICIT NONE

      !! Overloading for different types of attributes, data elements?
      CHARACTER(LEN=*)             :: elementName
      CHARACTER(LEN= 40), INTENT(IN) :: attributeNames(:)
      CHARACTER(LEN= 30), INTENT(IN) :: attributeValues(:)
      CHARACTER(LEN= 30), INTENT(IN) :: contentList(:)

      CHARACTER(LEN=200), ALLOCATABLE :: contentLineList(:)
      INTEGER :: i, contentLineLength, contentLineListSize
      CHARACTER(LEN=70)               :: format
      CHARACTER(LEN=200)              :: outputString

      IF(SIZE(attributeNames).NE.SIZE(attributeValues)) THEN
         WRITE(*,*) 'attributeNames', attributeNames
         WRITE(*,*) 'attributeValues', attributeValues
         STOP 'ERROR: SIZE(attributeNames).NE.SIZE(attributeValues)'
      END IF
      contentLineLength = 6 ! At most 6 data elements per line
      contentLineListSize = SIZE(contentList) / contentLineLength
      IF(contentLineListSize*contentLineLength.NE.SIZE(contentList)) THEN
         contentLineListSize = contentLineListSize + 1
      END IF
      ALLOCATE(contentLineList(contentLineListSize))
      CALL fillContentLineList(contentList,contentLineList,contentLineLength)
      outputString = '<'//TRIM(ADJUSTL(elementName))
      DO i = 1, SIZE(attributeNames)
         outputString = TRIM(ADJUSTL(outputString))//' '//TRIM(ADJUSTL(attributeNames(i)))//'="'//TRIM(ADJUSTL(attributeValues(i)))//'"'
      END DO
      IF(SIZE(contentLineList).EQ.0) THEN
         outputString = TRIM(ADJUSTL(outputString))//'/>'
      ELSE IF(SIZE(contentLineList).EQ.1) THEN
         outputString = TRIM(ADJUSTL(outputString))//'>'//contentLineList(1)//'</'//TRIM(ADJUSTL(elementName))//'>'
      ELSE
         outputString = TRIM(ADJUSTL(outputString))//'>'
      END IF
      WRITE(format,'(a,i0,a)') "(a",3*(currentElementIndex+1),",a)"
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
