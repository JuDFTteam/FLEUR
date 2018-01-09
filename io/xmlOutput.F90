!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

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
   PUBLIC writeXMLElementMatrixPoly, writeXMLElementMatrixFormPoly
   PUBLIC writeXMLElementMatrixForm
   PUBLIC openXMLElementFormPoly, openXMLElementPoly
   PUBLIC openXMLElementForm, openXMLElement
   PUBLIC openXMLElementNoAttributes, closeXMLElement
   PUBLIC getXMLOutputUnitNumber, isCurrentXMLElement

   CONTAINS

   FUNCTION getXMLOutputUnitNumber()
      IMPLICIT NONE
      INTEGER getXMLOutputUnitNumber
      getXMLOutputUnitNumber = xmlOutputUnit
   END FUNCTION getXMLOutputUnitNumber

   SUBROUTINE startXMLOutput()

      USE m_juDFT_args
      USE m_constants
      USE m_utility
      USE m_compile_descr
      
      IMPLICIT NONE

#ifdef CPP_MPI
      include "mpif.h"
      INTEGER           :: err, isize
#endif
      INTEGER           :: numFlags
      CHARACTER(LEN=8)  :: date
      CHARACTER(LEN=10) :: time
      CHARACTER(LEN=10) :: zone
      CHARACTER(LEN=10) :: dateString
      CHARACTER(LEN=10) :: timeString
      CHARACTER(LEN=6)  :: precisionString
      CHARACTER(LEN=9)  :: flags(11)
      CHARACTER(LEN=20) :: structureSpecifiers(11)
      CHARACTER(LEN=50) :: gitdesc,githash,gitbranch,compile_date,compile_user,compile_host
      CHARACTER(LEN=50) :: compile_flags,link_flags
      
      maxNumElements = 10
      ALLOCATE(elementList(maxNumElements))
      elementList = ''
      currentElementIndex = 0
      xmlOutputUnit = 53
      CALL DATE_AND_TIME(date,time,zone)
      WRITE(dateString,'(a4,a1,a2,a1,a2)') date(1:4),'/',date(5:6),'/',date(7:8)
      WRITE(timeString,'(a2,a1,a2,a1,a2)') time(1:2),':',time(3:4),':',time(5:6)
      IF (juDFT_was_argument("-info")) THEN
         OPEN (xmlOutputUnit,status='scratch')
      ELSE
         OPEN (xmlOutputUnit,file='out.xml',form='formatted',status='unknown')
      ENDIF
      WRITE (xmlOutputUnit,'(a)') '<?xml version="1.0" encoding="UTF-8" standalone="no"?>'
      WRITE (xmlOutputUnit,'(a)') '<fleurOutput fleurOutputVersion="0.27">'
      CALL openXMLElement('programVersion',(/'version'/),(/version_const/))
      CALL get_compile_desc(gitdesc,githash,gitbranch,compile_date,compile_user,compile_host,compile_flags,link_flags)
      CALL writeXMLElement('compilationInfo',(/'date','user','host','flag','link'/),(/compile_date,compile_user,compile_host,compile_flags,link_flags/))
      CALL writeXMLElement('gitInfo',(/'version       ','branch        ','lastCommitHash'/),(/gitdesc,gitbranch,githash/))
      CALL getComputerArchitectures(flags, numFlags)
      IF (numFlags.EQ.0) THEN
         numFlags = 1
         flags(numFlags) = 'GEN'
      END IF
      CALL writeXMLElementNoAttributes('targetComputerArchitectures',flags(1:numFlags))
      IF (numFlags.GT.1) THEN 
         STOP "ERROR: Define only one system architecture! (called by xmlOutput)"
      END IF
      CALL getPrecision(precisionString)
      CALL writeXMLElement('precision',(/'type'/),(/precisionString/))
      CALL getTargetStructureProperties(structureSpecifiers, numFlags)
      CALL writeXMLElementNoAttributes('targetStructureClass',structureSpecifiers(1:numFlags))
      CALL getAdditionalCompilationFlags(flags, numFlags)
      IF (numFlags.GE.1) THEN
         CALL writeXMLElementNoAttributes('additionalCompilerFlags',flags(1:numFlags))
      END IF
      CALL closeXMLElement('programVersion')
#ifdef CPP_MPI
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,isize,err)
      CALL writeXMLElementPoly('parallelizationParameters',(/'mpiPEs'/),(/isize/))
#endif
      CALL writeXMLElement('startDateAndTime',(/'date','time','zone'/),(/dateString,timeString,zone/))
   END SUBROUTINE startXMLOutput

   SUBROUTINE endXMLOutput()

      IMPLICIT NONE

      CHARACTER(LEN=8)  :: date
      CHARACTER(LEN=10) :: time
      CHARACTER(LEN=10)  :: zone
      CHARACTER(LEN=10) :: dateString
      CHARACTER(LEN=10)  :: timeString

      DO WHILE (currentElementIndex.NE.0)
         CALL closeXMLElement(elementList(currentElementIndex))
      END DO
      DEALLOCATE(elementList)
      CALL DATE_AND_TIME(date,time,zone)
      WRITE(dateString,'(a4,a1,a2,a1,a2)') date(1:4),'/',date(5:6),'/',date(7:8)
      WRITE(timeString,'(a2,a1,a2,a1,a2)') time(1:2),':',time(3:4),':',time(5:6)
      CALL writeXMLElement('endDateAndTime',(/'date','time','zone'/),(/dateString,timeString,zone/))
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
      ELSE
         ALLOCATE(charContentList(1))
      END IF

      charAttributeValues = ''
      charContentList = ''

      DO i = 1, SIZE(attributeValues)
         SELECT TYPE (attributeValues)
           TYPE IS(INTEGER)
               WRITE(charAttributeValues(i),'(i0)') attributeValues(i)
            TYPE IS(REAL)
               WRITE(charAttributeValues(i),'(f19.10)') attributeValues(i)
            TYPE IS(LOGICAL)
               WRITE(charAttributeValues(i),'(l1)') attributeValues(i)
#ifndef __PGI
            TYPE IS(CHARACTER(LEN=*))
               WRITE(charAttributeValues(i),'(a)') TRIM(ADJUSTL(attributeValues(i)))
            CLASS DEFAULT
               STOP 'Type of attributeValues not allowed'
#endif
         END SELECT
      END DO

      IF (PRESENT(contentList)) THEN
         DO i = 1, SIZE(contentList)
            SELECT TYPE(contentList)
               TYPE IS(INTEGER)
                  WRITE(charContentList(i),'(i0)') contentList(i)
               TYPE IS(REAL)
                  WRITE(charContentList(i),'(f19.10)') contentList(i)
               TYPE IS(LOGICAL)
                  WRITE(charContentList(i),'(l1)') contentList(i)
#ifndef __PGI
               TYPE IS(CHARACTER(LEN=*))
                  WRITE(charContentList(i),'(a)') TRIM(ADJUSTL(contentList(i)))
               CLASS DEFAULT
                  STOP 'Type of contentList not allowed'
#endif
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
      INTEGER                         :: contentListSize, overallListSize, numContentLineChars
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
         WRITE(outputString,format) TRIM(ADJUSTL(outputString)), ' ', TRIM(ADJUSTL(attributeNames(i))),&
                                    '="', TRIM(ADJUSTL(attributeValues(i))), '"'
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
            outputString = TRIM(ADJUSTL(outputString))//'>'//TRIM(ADJUSTL(contentLineList(1)))//'</'//&
                 TRIM(ADJUSTL(elementName))//'>'
         ELSE
            outputString = TRIM(ADJUSTL(outputString))//'>'
         END IF
         WRITE(xmlOutputUnit,format) ' ', TRIM(ADJUSTL(outputString))
         IF(SIZE(contentLineList).GT.1) THEN
            DO i = 1, SIZE(contentLineList)
               IF (i.EQ.SIZE(contentLineList)) THEN
                  numContentLineChars = 20*MOD(SIZE(contentList),contentLineLength)
                  IF(numContentLineChars.EQ.0) numContentLineChars = 20 * contentLineLength
                  WRITE(format,'(a,i0,a,i0,a)') "(a",3*(currentElementIndex+2),",a",numContentLineChars,")"
               ELSE
                  WRITE(format,'(a,i0,a)') "(a",3*(currentElementIndex+2),",a100)"
               END IF
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

   SUBROUTINE writeXMLElementMatrixPoly(elementName,attributeNames,attributeValues,matrix)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)           :: elementName
      CHARACTER(LEN=*), INTENT(IN)           :: attributeNames(:)
      CLASS(*),         INTENT(IN)           :: attributeValues(:)
      CLASS(*),         INTENT(IN)           :: matrix(:,:)

      INTEGER, ALLOCATABLE :: lengths(:,:)
      INTEGER              :: numDifferentMatrixElementLengths

      numDifferentMatrixElementLengths = 1 ! At the moment nothing else is implemented.
                                           ! In principle one could think about an own
                                           ! length for each matrix element.

      ALLOCATE(lengths(MAX(SIZE(attributeNames),numDifferentMatrixElementLengths),3))
      lengths = 0
      CALL writeXMLElementMatrixFormPoly(elementName,attributeNames,attributeValues,lengths,matrix)
      DEALLOCATE(lengths)

   END SUBROUTINE writeXMLElementMatrixPoly

   SUBROUTINE writeXMLElementMatrixFormPoly(elementName,attributeNames,attributeValues,lengths,matrix)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)           :: elementName
      CHARACTER(LEN=*), INTENT(IN)           :: attributeNames(:)
      CLASS(*),         INTENT(IN)           :: attributeValues(:)
      INTEGER,          INTENT(IN)           :: lengths(:,:)
      CLASS(*),         INTENT(IN)           :: matrix(:,:)

      CHARACTER(LEN= 30), ALLOCATABLE :: charAttributeValues(:)
      CHARACTER(LEN= 50), ALLOCATABLE :: charMatrix(:,:)
      INTEGER                         :: i, j

      ALLOCATE(charAttributeValues(SIZE(attributeValues)))
      ALLOCATE(charMatrix(SIZE(matrix,1),SIZE(matrix,2)))

      charAttributeValues = ''
      charMatrix = ''

      DO i = 1, SIZE(attributeValues)
         SELECT TYPE (attributeValues)
           TYPE IS(INTEGER)
               WRITE(charAttributeValues(i),'(i0)') attributeValues(i)
            TYPE IS(REAL)
               WRITE(charAttributeValues(i),'(f19.10)') attributeValues(i)
            TYPE IS(LOGICAL)
               WRITE(charAttributeValues(i),'(l1)') attributeValues(i)
#ifndef __PGI
            TYPE IS(CHARACTER(LEN=*))
               WRITE(charAttributeValues(i),'(a)') TRIM(ADJUSTL(attributeValues(i)))
            CLASS DEFAULT
               STOP 'Type of attributeValues not allowed'
#endif
         END SELECT
      END DO

      DO j = 1, SIZE(matrix,2)
         DO i = 1, SIZE(matrix,1)
            SELECT TYPE(matrix)
               TYPE IS(INTEGER)
                  WRITE(charMatrix(i,j),'(i0)') matrix(i,j)
               TYPE IS(REAL)
                  WRITE(charMatrix(i,j),'(f19.10)') matrix(i,j)
               TYPE IS(COMPLEX)
                  WRITE(charMatrix(i,j),'(a,SP,f8.5,a,f8.5,a)') '(',REAL(matrix(i,j)),',', AIMAG(matrix(i,j)),')'
               TYPE IS(LOGICAL)
                  WRITE(charMatrix(i,j),'(l1)') matrix(i,j)
#ifndef __PGI
               TYPE IS(CHARACTER(LEN=*))
                  WRITE(charMatrix(i,j),'(a)') TRIM(ADJUSTL(matrix(i,j)))
               CLASS DEFAULT
                  STOP 'Type of matrix not allowed'
#endif
            END SELECT
         END DO
      END DO

      CALL writeXMLElementMatrixForm(elementName,attributeNames,charAttributeValues,lengths,charMatrix)
      DEALLOCATE(charMatrix,charAttributeValues)
   END SUBROUTINE writeXMLElementMatrixFormPoly

   SUBROUTINE writeXMLElementMatrixForm(elementName,attributeNames,attributeValues,lengths,matrix)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: elementName
      CHARACTER(LEN=*), INTENT(IN) :: attributeNames(:)
      CHARACTER(LEN=*), INTENT(IN) :: attributeValues(:)
      INTEGER,          INTENT(IN) :: lengths(:,:)
      CHARACTER(LEN=*), INTENT(IN) :: matrix(:,:)

      CHARACTER(LEN=400), ALLOCATABLE :: matrixRows(:)
      INTEGER, ALLOCATABLE            :: bigLengths(:,:)
      INTEGER                         :: i, j
      CHARACTER(LEN=70)               :: format
      CHARACTER(LEN=400)              :: outputString
      INTEGER                         :: matrixRowLength, numMatrixRows
      INTEGER                         :: overallListSize
      INTEGER                         :: lengthsShape(2)

      IF(SIZE(attributeNames).NE.SIZE(attributeValues)) THEN
         WRITE(*,*) 'attributeNames', attributeNames
         WRITE(*,*) 'attributeValues', attributeValues
         STOP 'ERROR in writeXMLElementMatrixForm: SIZE(attributeNames).NE.SIZE(attributeValues)'
      END IF

      lengthsShape = SHAPE(lengths)

      overallListSize = MAX(SIZE(attributeNames),1)
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
         WRITE(outputString,format) TRIM(ADJUSTL(outputString)), ' ', TRIM(ADJUSTL(attributeNames(i))),&
                                    '="', TRIM(ADJUSTL(attributeValues(i))), '"'
      END DO
      WRITE(format,'(a,i0,a)') "(a",3*(currentElementIndex+1),",a)"

      matrixRowLength = SIZE(matrix,1)
      numMatrixRows = SIZE(matrix,2)

      ALLOCATE(matrixRows(numMatrixRows))

      matrixRows = ''
      DO i = 1, numMatrixRows
         DO j = 1, matrixRowLength
            WRITE(matrixRows(i),'(a,a1,a)') TRIM(ADJUSTL(matrixRows(i))),' ',TRIM(ADJUSTL(matrix(j,i)))
         END DO
      END DO

      outputString = TRIM(ADJUSTL(outputString))//'>'
      WRITE(xmlOutputUnit,format) ' ', TRIM(ADJUSTL(outputString))

      DO i = 1, numMatrixRows
         WRITE(format,'(a,i0,a)') "(a",3*(currentElementIndex+2),",a)"
         WRITE(xmlOutputUnit,format) ' ', TRIM(ADJUSTL(matrixRows(i)))
      END DO

      WRITE(format,'(a,i0,a)') "(a",3*(currentElementIndex+1),",a)"
      outputString = '</'//TRIM(ADJUSTL(elementName))//'>'
      WRITE(xmlOutputUnit,format) ' ', TRIM(ADJUSTL(outputString))

      DEALLOCATE (bigLengths)

   END SUBROUTINE writeXMLElementMatrixForm

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

   SUBROUTINE writeXMLElementNoAttributes(elementName,contentList)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)           :: elementName
      CHARACTER(LEN=*), INTENT(IN)           :: contentList(:)

      CHARACTER(LEN=200), ALLOCATABLE :: contentLineList(:)
      INTEGER :: i, j, contentLineLength, contentLineListSize
      CHARACTER(LEN=70)               :: format
      CHARACTER(LEN=200)              :: outputString
      INTEGER                         :: contentListSize, overallListSize, numContentLineChars

      outputString = '<'//TRIM(ADJUSTL(elementName))

      WRITE(format,'(a,i0,a)') "(a",3*(currentElementIndex+1),",a)"

      contentLineLength = 5 ! At most 5 data elements per line
      contentLineListSize = SIZE(contentList) / contentLineLength
      IF(contentLineListSize*contentLineLength.NE.SIZE(contentList)) THEN
         contentLineListSize = contentLineListSize + 1
      END IF
      ALLOCATE(contentLineList(contentLineListSize))
      CALL fillContentLineList(contentList,contentLineList,contentLineLength)
      IF(contentLineListSize.EQ.0) THEN
         outputString = TRIM(ADJUSTL(outputString))//'> </'//TRIM(ADJUSTL(elementName))//'>'
      ELSE IF(contentLineListSize.EQ.1) THEN
         outputString = TRIM(ADJUSTL(outputString))//'>'//TRIM(ADJUSTL(contentLineList(1)))//'</'//&
              TRIM(ADJUSTL(elementName))//'>'
      ELSE
         outputString = TRIM(ADJUSTL(outputString))//'>'
      END IF
      WRITE(xmlOutputUnit,format) ' ', TRIM(ADJUSTL(outputString))
      IF(contentLineListSize.GT.1) THEN
         DO i = 1, SIZE(contentLineList)
            IF (i.EQ.SIZE(contentLineList)) THEN
               numContentLineChars = 20*MOD(SIZE(contentList),contentLineLength)
               IF(numContentLineChars.EQ.0) numContentLineChars = 20 * contentLineLength
               WRITE(format,'(a,i0,a,i0,a)') "(a",3*(currentElementIndex+2),",a",numContentLineChars,")"
            ELSE
               WRITE(format,'(a,i0,a)') "(a",3*(currentElementIndex+2),",a100)"
            END IF
            WRITE(xmlOutputUnit,format) ' ', TRIM(ADJUSTL(contentLineList(i)))
         END DO
         WRITE(format,'(a,i0,a)') "(a",3*(currentElementIndex+1),",a)"
         outputString = '</'//TRIM(ADJUSTL(elementName))//'>'
         WRITE(xmlOutputUnit,format) ' ', TRIM(ADJUSTL(outputString))
      END IF

   END SUBROUTINE writeXMLElementNoAttributes

   SUBROUTINE fillContentLineList(contentList,contentLineList,contentLineLength)

      IMPLICIT NONE

      CHARACTER(LEN= *), INTENT(IN)    :: contentList(:)
      CHARACTER(LEN=200), INTENT(INOUT) :: contentLineList(:)
      INTEGER,            INTENT(IN)    :: contentLineLength
      INTEGER :: i, j

      contentLineList = ''
      DO i = 1, SIZE(contentLineList)
         DO j = 1, contentLineLength
            IF((i-1)*contentLineLength+j.GT.SIZE(contentList)) THEN
               RETURN
            END IF
            WRITE(contentLineList(i),'(a,a20)') TRIM(ADJUSTL(contentLineList(i))),&
                        TRIM(ADJUSTL(contentList((i-1)*contentLineLength+j)))
         END DO
      END DO
   END SUBROUTINE fillContentLineList

   SUBROUTINE openXMLElementPoly(elementName,attributeNames,attributeValues)

      IMPLICIT NONE

      CHARACTER(LEN=*)             :: elementName
      CHARACTER(LEN=*), INTENT(IN) :: attributeNames(:)
      CLASS(*),         INTENT(IN) :: attributeValues(:)

      INTEGER, ALLOCATABLE :: lengths(:,:)

      ALLOCATE(lengths(SIZE(attributeNames),2))
      lengths = 0
      CALL openXMLElementFormPoly(elementName,attributeNames,attributeValues,lengths)
      DEALLOCATE(lengths)

   END SUBROUTINE openXMLElementPoly

   SUBROUTINE openXMLElementFormPoly(elementName,attributeNames,attributeValues,lengths)

      IMPLICIT NONE

      CHARACTER(LEN=*)             :: elementName
      CHARACTER(LEN=*), INTENT(IN) :: attributeNames(:)
      CLASS(*),         INTENT(IN) :: attributeValues(:)
      INTEGER,          INTENT(IN) :: lengths(:,:)

      CHARACTER(LEN= 30), ALLOCATABLE :: charAttributeValues(:)
      INTEGER                         :: i

      ALLOCATE(charAttributeValues(SIZE(attributeValues)))

      charAttributeValues = ''

      DO i = 1, SIZE(attributeValues)
         SELECT TYPE (attributeValues)
            TYPE IS(INTEGER)
               WRITE(charAttributeValues(i),'(i0)') attributeValues(i)
            TYPE IS(REAL)
               WRITE(charAttributeValues(i),'(f19.10)') attributeValues(i)
            TYPE IS(LOGICAL)
               WRITE(charAttributeValues(i),'(l1)') attributeValues(i)
#ifndef __PGI
            TYPE IS(CHARACTER(LEN=*))
               WRITE(charAttributeValues(i),'(a)') TRIM(ADJUSTL(attributeValues(i)))
            CLASS DEFAULT
               STOP 'Type of attributeValues not allowed'
#endif
         END SELECT
      END DO

      CALL openXMLElementForm(elementName,attributeNames,charAttributeValues,lengths)
      DEALLOCATE(charAttributeValues)

   END SUBROUTINE openXMLElementFormPoly

   SUBROUTINE openXMLElementNoAttributes(elementName)

      IMPLICIT NONE

      CHARACTER(LEN=*)             :: elementName

      INTEGER :: i
      CHARACTER(LEN=70)  :: format
      CHARACTER(LEN=200) :: openingString

      IF(currentElementIndex.EQ.maxNumElements) THEN
         WRITE(*,*) 'elementName ', TRIM(ADJUSTL(elementName))
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

      CHARACTER(LEN=*), INTENT(IN)           :: elementName
      CHARACTER(LEN=*), INTENT(IN)           :: attributeNames(:)
      CHARACTER(LEN=*), INTENT(IN)           :: attributeValues(:)

      INTEGER, ALLOCATABLE :: lengths(:,:)

      ALLOCATE(lengths(SIZE(attributeNames),2))
      lengths = 0
      CALL openXMLElementForm(elementName,attributeNames,attributeValues,lengths)
      DEALLOCATE(lengths)

   END SUBROUTINE

   SUBROUTINE openXMLElementForm(elementName,attributeNames,attributeValues,lengths)

      IMPLICIT NONE

      CHARACTER(LEN=*)             :: elementName
      CHARACTER(LEN=*), INTENT(IN) :: attributeNames(:)
      CHARACTER(LEN=*), INTENT(IN) :: attributeValues(:)
      INTEGER,          INTENT(IN) :: lengths(:,:)

      CHARACTER(LEN=70)  :: format
      CHARACTER(LEN=200) :: openingString

      INTEGER, ALLOCATABLE            :: bigLengths(:,:)
      INTEGER                         :: i, j
      INTEGER                         :: overallListSize
      INTEGER                         :: lengthsShape(2)

      IF(SIZE(attributeNames).NE.SIZE(attributeValues)) THEN
         WRITE(*,*) 'elementName ', TRIM(ADJUSTL(elementName))
         WRITE(*,*) 'attributeNames ', attributeNames
         WRITE(*,*) 'attributeValues ', attributeValues
         STOP 'ERROR: SIZE(attributeNames).NE.SIZE(attributeValues)'
      END IF
      IF(currentElementIndex.EQ.maxNumElements) THEN
         WRITE(*,*) 'elementName ', TRIM(ADJUSTL(elementName))
         WRITE(*,*) 'attributeNames ', attributeNames
         WRITE(*,*) 'attributeValues ', attributeValues
         WRITE(*,*) 'elementList ', elementList
         STOP 'ERROR: xml hierarchy too deep!'
      END IF

      lengthsShape = SHAPE(lengths)
      overallListSize = SIZE(attributeNames)
      ALLOCATE(bigLengths(overallListSize,2))
      bigLengths = 0
      DO j = 1, 2
         DO i = 1, MIN(overallListSize,lengthsShape(1))
            bigLengths(i,j) = lengths(i,j)
         END DO
      END DO

      openingString = '<'//TRIM(ADJUSTL(elementName))
      DO i = 1, SIZE(attributeNames)
         WRITE(format,'(a)') "(a,a,a"
         IF (bigLengths(i,1).GT.0) WRITE(format,'(a,i0)') TRIM(ADJUSTL(format)),bigLengths(i,1)
         WRITE(format,'(a,a)') TRIM(ADJUSTL(format)),",a,a"
         IF (bigLengths(i,2).GT.0) WRITE(format,'(a,i0)') TRIM(ADJUSTL(format)),bigLengths(i,2)
         WRITE(format,'(a,a)') TRIM(ADJUSTL(format)),",a)"
         WRITE(openingString,format) TRIM(ADJUSTL(openingString)), ' ', TRIM(ADJUSTL(attributeNames(i))),&
                                    '="', TRIM(ADJUSTL(attributeValues(i))), '"'
      END DO
      openingString = TRIM(ADJUSTL(openingString))//'>'

      currentElementIndex = currentElementIndex + 1
      elementList(currentElementIndex) = TRIM(ADJUSTL(elementName))
      WRITE(format,'(a,i0,a)') "(a",3*currentElementIndex,",a)"
      WRITE(xmlOutputUnit,format) ' ', TRIM(ADJUSTL(openingString))
      DEALLOCATE (bigLengths)

   END SUBROUTINE openXMLElementForm

   SUBROUTINE closeXMLElement(elementName)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: elementName
      INTEGER :: i
      CHARACTER(LEN=70) :: format
      CHARACTER(LEN=70) :: closingString

      IF(TRIM(ADJUSTL(elementList(currentElementIndex))).NE.TRIM(ADJUSTL(elementName))) THEN
         WRITE(*,*) 'elementList(currentElementIndex): ', TRIM(ADJUSTL(elementList(currentElementIndex)))
         WRITE(*,*) 'elementName ', TRIM(ADJUSTL(elementName))
         STOP 'ERROR: Closing xml element inconsistency!'
      END IF
      closingString = '</'//TRIM(ADJUSTL(elementName))//'>'
      WRITE(format,'(a,i0,a)') "(a",3*currentElementIndex,",a)"
      WRITE(xmlOutputUnit,format) ' ', TRIM(ADJUSTL(closingString))
      currentElementIndex = currentElementIndex - 1

   END SUBROUTINE closeXMLElement

   LOGICAL FUNCTION isCurrentXMLElement(elementName)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: elementName
      if (currentElementIndex==0.or.currentElementIndex>SIZE(elementList)) THEN
                  isCurrentXMLElement=.false.
                  return
      endif
      isCurrentXMLElement = (TRIM(ADJUSTL(elementList(currentElementIndex))).EQ.TRIM(ADJUSTL(elementName)))

   END FUNCTION isCurrentXMLElement

END MODULE m_xmlOutput
