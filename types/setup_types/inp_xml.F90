!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------


MODULE m_inp_xml
  USE m_xmlIntWrapFort
  USE m_judft
CONTAINS
  CHARACTER (len=200) FUNCTION inp_xml_speciesxpath_for_group(n)
    IMPLICIT NONE
    INTEGER,INTENT(in):: n
    
    INTEGER           :: i
    CHARACTER(len=200)::xpath,species
    !First determine name of species from group
    WRITE(xPath,*) '/fleurInput/atomGroups/atomGroup[',n,']/@species'
    species=TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPath)))))

    DO i=1,xmlGetNumberOfNodes('/fleurInput/atomSpecies/species')
       WRITE(xPath,*) '/fleurInput/atomSpecies/species[',i,']'
       IF (TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@name')))==TRIM(species)) THEN
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

   SUBROUTINE getIntegerSequenceFromString(string, sequence, count)

      IMPLICIT NONE

      CHARACTER(*),         INTENT(IN)  :: string
      INTEGER, ALLOCATABLE, INTENT(OUT) :: sequence(:)
      INTEGER,              INTENT(OUT) :: count

      INTEGER :: i, length, start, lastNumber, currentNumber, index
      LOGICAL singleNumber, comma, dash

      ! 3 cases: 1. a single number
      !          2. number - number
      !          3. comma separated numbers

      length = LEN(string)
      count = 0
      start = 1
      singleNumber = .TRUE.
      comma = .FALSE.
      dash = .FALSE.
      lastNumber = 0
      count = 0

      ! 1. Determine number count

      DO i = 1, length
         SELECT CASE (string(i:i))
         CASE ('0':'9')
         CASE (',')
            IF ((start.EQ.i).OR.(dash)) THEN
               STOP 'String has wrong syntax (in getIntegerSequenceFromString)'
            END IF
            singleNumber = .FALSE.
            comma = .TRUE.
            READ(string(start:i-1),*) lastNumber
            count = count + 1
            start = i+1
         CASE ('-')
            IF ((start.EQ.i).OR.(dash).OR.(comma)) THEN
               STOP 'String has wrong syntax (in getIntegerSequenceFromString)'
            END IF
            singleNumber = .FALSE.
            dash = .TRUE.
            READ(string(start:i-1),*) lastNumber
            start = i+1
         CASE DEFAULT
            STOP 'String has wrong syntax (in getIntegerSequenceFromString)'
         END SELECT
      END DO
      IF(start.GT.length) THEN
         STOP 'String has wrong syntax (in getIntegerSequenceFromString)'
      END IF
      READ(string(start:length),*) currentNumber
      IF (dash) THEN
         count = currentNumber - lastNumber + 1
      ELSE
         count = count + 1
      END IF

      IF (ALLOCATED(sequence)) THEN
         DEALLOCATE(sequence)
      END IF
      ALLOCATE(sequence(count))

      ! 2. Read in numbers iff comma separation ...and store numbers in any case

      IF (singleNumber) THEN
         sequence(1) = currentNumber
      ELSE IF (dash) THEN
         DO i = 1, count
            sequence(i) = lastNumber + i - 1
         END DO
      ELSE
         index = 1
         start = 1
         DO i = 1, length
            SELECT CASE (string(i:i))
            CASE (',')
               comma = .TRUE.
               READ(string(start:i-1),*) lastNumber
               start = i+1
               sequence(index) = lastNumber
               index = index + 1
            END SELECT
         END DO
         sequence(index) = currentNumber
      END IF

   END SUBROUTINE getIntegerSequenceFromString

END MODULE m_inp_xml
