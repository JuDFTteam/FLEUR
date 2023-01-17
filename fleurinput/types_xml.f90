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
   USE m_calculator
   IMPLICIT NONE
   PRIVATE

   LOGICAL :: INITIALIZED = .false.
   TYPE t_xml
      INTEGER:: id
      character(len=200):: basepath = ""
      integer           :: versionNumber = 0
      INTEGER           :: currentversionNumber = 36 !parameters are not allowed here
   CONTAINS
      PROCEDURE        :: init
      PROCEDURE        :: GetNumberOfNodes
      PROCEDURE, NOPASS :: SetAttributeValue
      PROCEDURE        :: GetAttributeValue
      PROCEDURE        :: GetAttributeValue_List
      PROCEDURE, NOPASS :: getIntegerSequenceFromString
      PROCEDURE        :: read_q_list
      PROCEDURE, NOPASS :: popFirstStringToken
      PROCEDURE, NOPASS :: countStringTokens
      PROCEDURE        :: speciesPath
      PROCEDURE, NOPASS :: groupPath
      PROCEDURE :: get_nat
      PROCEDURE :: get_lmaxd
      PROCEDURE :: get_nlo
      PROCEDURE :: get_ntype
      PROCEDURE :: posPath
      PROCEDURE :: set_basepath
      procedure, NOPASS :: writexml
      PROCEDURE, NOPASS :: FreeResources
   END TYPE t_xml
   PUBLIC t_xml, evaluateFirstOnly, EvaluateFirst, evaluateFirstBoolOnly, evaluateFirstIntOnly, &
      evaluateList, drop_schema_files

CONTAINS
   subroutine set_basepath(xml, xpath)
      CLASS(t_xml), INTENT(INOUT)::xml
      character(len=*)::xpath

      xml%basepath = xpath
   end subroutine

   subroutine drop_schema_files(version, output_version)
      USE iso_c_binding
      character(len=4, kind=c_char), intent(in):: version
      character(len=4, kind=c_char), intent(in):: output_version
      integer :: errorStatus
      INTERFACE
         FUNCTION dropInputSchema(version) BIND(C, name="dropInputSchema")
            USE iso_c_binding
            INTEGER(c_int) ::dropInputSchema
            CHARACTER(kind=c_char) ::version
         END FUNCTION dropInputSchema

      END INTERFACE
      INTERFACE

         FUNCTION dropOutputSchema(version) BIND(C, name="dropOutputSchema")
            USE iso_c_binding
            INTEGER(c_int) ::dropOutputSchema
            CHARACTER(kind=c_char) ::version
         END FUNCTION dropOutputSchema
      END INTERFACE
      
      !Now validate with schema
      errorStatus = 0
      errorStatus = dropInputSchema(version//C_NULL_CHAR)

      IF (errorStatus .NE. 0) THEN
         CALL juDFT_error('Error: Cannot print out FleurInputSchema.xsd for version '//version)
      END IF

      errorStatus = 0
      errorStatus = dropOutputSchema(output_version//C_NULL_CHAR)

      IF (errorStatus .NE. 0) THEN
         WRITE (*, *) 'Cannot print out FleurOutputSchema.xsd for version '//version
      END IF

   end subroutine

   subroutine validate_with_schema(version, output_version)
      use iso_c_binding
      character(len=4, kind=c_char), intent(in):: version
      character(len=4, kind=c_char), intent(in):: output_version
      integer :: errorStatus
      character(len=200, KIND=c_char):: schemaFilename
     
      call drop_schema_files(version, output_version)
      schemaFilename = "FleurInputSchema.xsd"//C_NULL_CHAR
      CALL ParseSchema(schemaFilename)
      CALL ValidateDoc()
      CALL InitXPath()
   end

   SUBROUTINE init(xml, filename_add, old_version)
      USE iso_c_binding

      CLASS(t_xml), INTENT(INOUT) :: xml
      CHARACTER(len=100), INTENT(IN) :: filename_add
      LOGICAL, OPTIONAL, INTENT(inout):: old_version

      LOGICAL                        :: l_allow_old
      INTEGER                        :: errorStatus
      CHARACTER(LEN=200, KIND=c_char) :: docFilename, versionString, outputVersionString
      INTEGER                        :: i, numberNodes
      CHARACTER(LEN=255)             :: xPathA, xPathB, valueString
      REAL                           :: tempReal

!   I comment this initialization test out because it causes trouble when generating additional k-point sets.
!   Don't know for what it is needed, anyway. I leave rest of this mechanism in the code. ...for now. G.M.
!    if (INITIALIZED) RETURN
      INITIALIZED = .true.

      l_allow_old = .false.
      if (present(old_version)) then
         l_allow_old = old_version
         old_version = .false.
      endif

      !Open inp.xml
      docFilename = TRIM(filename_add)//"inp.xml"//C_NULL_CHAR
      CALL InitInterface()
      CALL ParseDoc(docFilename)
      CALL InitXPath()

      ! Check version of inp.xml
      versionString = adjustl(xml%GetAttributeValue('/fleurInput/@fleurInputVersion'))
      read (versionString, *) tempReal
      write(outputVersionString,'(a,i0)') '0.', xml%currentversionNumber
      xml%versionNumber = nint(tempReal*100)
      IF (xml%versionNumber .NE. xml%currentversionNumber) THEN
         if (.not. l_allow_old .and. xml%versionNumber<33) CALL juDFT_error('Version number of '//TRIM(filename_add)//'inp.xml file is not compatible with this fleur version')
         if (present(old_version)) old_version = .true.
      END IF

      call validate_with_Schema(versionString, outputVersionString)

      ! Read in constants
      xPathA = '/fleurInput/constants/constant'
      numberNodes = xml%GetNumberOfNodes(xPathA)
      DO i = 1, numberNodes
         WRITE (xPathB, '(a,a,i0,a)') TRIM(ADJUSTL(xPathA)), '[', i, ']'
         tempReal = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@value'))
         valueString = xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@name')
         CALL ASSIGN_var(valueString, tempReal)
      END DO
   END SUBROUTINE init

   INTEGER FUNCTION get_lmaxd(xml)
      CLASS(t_xml), INTENT(IN)::xml
      INTEGER :: n
      get_lmaxd = 0
      DO n = 1, xml%get_ntype()
         get_lmaxd = MAX(get_lmaxd, evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(xml%speciesPath(n))//'/atomicCutoffs/@lmax')))
      ENDDO
   END FUNCTION get_lmaxd

   FUNCTION get_nlo(xml)
      CLASS(t_xml), INTENT(IN)::xml
      INTEGER, ALLOCATABLE::get_nlo(:)

      INTEGER n, numNodes, lNumCount, nNumCount, iLO, nLO
      CHARACTER(LEN=200) :: xPathA, xPathB, xPathC
      CHARACTER(LEN=200) :: lString, nString

      INTEGER, ALLOCATABLE :: lNumbers(:), nNumbers(:)

      ALLOCATE (get_nlo(xml%get_ntype()))
      DO n = 1, xml%get_ntype()
         xPathA = TRIM(ADJUSTL(xml%speciesPath(n))//'/lo')
         numNodes = xml%GetNumberOfNodes(TRIM(xPathA))
         nLO = 0
         DO iLO = 1, numNodes
            WRITE(xPathB,*) TRIM(xPathA),'[',iLO,']/@l'
            WRITE(xPathC,*) TRIM(xPathA),'[',iLO,']/@n'
            lString = xml%getAttributeValue(TRIM(ADJUSTL(xPathB)))
            nString = xml%getAttributeValue(TRIM(ADJUSTL(xPathC)))
            CALL getIntegerSequenceFromString(TRIM(ADJUSTL(lString)), lNumbers, lNumCount)
            CALL getIntegerSequenceFromString(TRIM(ADJUSTL(nString)), nNumbers, nNumCount)
            IF(lNumCount.NE.nNumCount) THEN
               CALL judft_error('Error in LO input: l quantum number count does not equal n quantum number count')
            END IF
            nLO = nLO + lNumCount
            DEALLOCATE (lNumbers, nNumbers)
         END DO

         get_nlo(n) = nLO
      ENDDO
   END FUNCTION get_nlo

   SUBROUTINE getIntegerSequenceFromString(string, sequence, count)
      use m_juDFT_stop
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
               CALL judft_error('String has wrong syntax (in getIntegerSequenceFromString)')
            END IF
            singleNumber = .FALSE.
            comma = .TRUE.
            READ(string(start:i-1),*) lastNumber
            count = count + 1
            start = i+1
         CASE ('-')
            IF ((start.EQ.i).OR.(dash).OR.(comma)) THEN
               CALL judft_error('String has wrong syntax (in getIntegerSequenceFromString)')
            END IF
            singleNumber = .FALSE.
            dash = .TRUE.
            READ(string(start:i-1),*) lastNumber
            start = i+1
         CASE DEFAULT
            CALL judft_error('String has wrong syntax (in getIntegerSequenceFromString)')
         END SELECT
      END DO
      IF(start.GT.length) THEN
         CALL judft_error('String has wrong syntax (in getIntegerSequenceFromString)')
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
               SEQUENCE(index) = lastNumber
               index = index + 1
            END SELECT
         END DO
         sequence(index) = currentNumber
      END IF
   END SUBROUTINE getIntegerSequenceFromString

   FUNCTION speciesPath(xml, itype)
      CLASS(t_xml), INTENT(IN)::xml
      INTEGER::itype
      CHARACTER(len=:), ALLOCATABLE::speciesPath

      INTEGER           :: i
      CHARACTER(len=200)::xpath, species
      !First determine name of species from group
      WRITE (xPath, '(a,i0,a)') '/fleurInput/atomGroups/atomGroup[', itype, ']/@species'
      species = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPath)))))

      DO i = 1, xml%GetNumberOfNodes('/fleurInput/atomSpecies/species')
         WRITE (xPath, '(a,i0,a)') '/fleurInput/atomSpecies/species[', i, ']'
         IF (TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPath))//'/@name'))) == TRIM(species)) THEN
            speciesPath = TRIM(xpath)
            RETURN
         END IF
      END DO
      WRITE (xpath, *) itype
      CALL judft_error("No species found for name "//TRIM(species)//" used in atom group "//TRIM(xpath))
   END FUNCTION speciesPath

   FUNCTION groupPath(itype)
      IMPLICIT NONE
      INTEGER, INTENT(in)::itype
      CHARACTER(len=:), allocatable::groupPath
      CHARACTER(len=100)::str
      WRITE (str, "(a,i0,a)") '/fleurInput/atomGroups/atomGroup[', itype, ']'
      groupPath = TRIM(str)
   END FUNCTION groupPath

   INTEGER FUNCTION get_nat(xml)
      CLASS(t_xml), INTENT(IN)::xml

      INTEGER ntype, n
      CHARACTER(len=100)::xpath
      get_nat = 0
      DO n = 1, xml%get_ntype()
         xpath = xml%groupPath(n)
         get_nat = get_nat + xml%getNumberOfNodes(TRIM(ADJUSTL(xPath))//'/relPos')
         get_nat = get_nat + xml%getNumberOfNodes(TRIM(ADJUSTL(xPath))//'/absPos')
         get_nat = get_nat + xml%getNumberOfNodes(TRIM(ADJUSTL(xPath))//'/filmPos')
      END DO
   END FUNCTION get_nat

   INTEGER FUNCTION get_ntype(xml)
      CLASS(t_xml), INTENT(IN)::xml

      get_ntype = xml%getNumberOfNodes('/fleurInput/atomGroups/atomGroup')
   END FUNCTION get_ntype

   FUNCTION posPath(xml, nat)
      CLASS(t_xml), INTENT(IN):: xml
      INTEGER, intent(in)     :: nat
      CHARACTER(len=:), ALLOCATABLE::posPath

      INTEGER na, n
      CHARACTER(len=100)::xpath, xpath2
      na = nat
      DO n = 1, xml%get_ntype()
         xpath = xml%groupPath(n)
         IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPath))//'/relPos') > 0) xpath = TRIM(ADJUSTL(xPath))//'/relPos'
         IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPath))//'/absPos') > 0) xpath = TRIM(ADJUSTL(xPath))//'/absPos'
         IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPath))//'/filmPos') > 0) xpath = TRIM(ADJUSTL(xPath))//'/filmPos'
         IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPath))) < na) THEN
            na = na - xml%getNumberOfNodes(TRIM(ADJUSTL(xPath)))
         ELSE
            WRITE (xpath2, "(a,a,i0,a)") TRIM(ADJUSTL(xpath)), '[', na, ']'
            posPath = TRIM(xpath2)
            RETURN
         END IF
      END DO
      CALL judft_error("Not so many positions found in inp.xml")
   END FUNCTION posPath
   FUNCTION popFirstStringToken(line) RESULT(firstToken)

      IMPLICIT NONE

      CHARACTER(*), INTENT(INOUT) :: line
      CHARACTER(LEN=LEN(line)) :: firstToken

      INTEGER separatorIndex

      separatorIndex = 0
      line = TRIM(ADJUSTL(line))

      separatorIndex = INDEX(line, ' ')
      IF (separatorIndex .LE. 1) THEN
         firstToken = ' '
      ELSE
         firstToken = line(1:separatorIndex - 1)
         line = line(separatorIndex + 1:)
      END IF

   END FUNCTION popFirstStringToken

   FUNCTION countStringTokens(line) RESULT(tokenCount)

      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: line
      INTEGER                  :: tokenCount

      CHARACTER(LEN=LEN(line)) :: tempLine
      INTEGER separatorIndex

      tokenCount = 0

      tempLine = TRIM(ADJUSTL(line))

      DO WHILE (tempLine .NE. '')
         separatorIndex = 0
         separatorIndex = INDEX(tempLine, ' ')
         IF (separatorIndex .EQ. 0) THEN
            tempLine = ''
         ELSE
            tempLine = TRIM(ADJUSTL(tempLine(separatorIndex + 1:)))
         END IF
         tokenCount = tokenCount + 1
      END DO

   END FUNCTION countStringTokens

   FUNCTION read_q_list(xml, path) RESULT(q)
      IMPLICIT NONE
      CLASS(t_xml), INTENT(IN)::xml
      CHARACTER(len=*), INTENT(in):: path
      REAL, ALLOCATABLE::q(:, :)

      INTEGER:: n, i
      CHARACTER(len=256):: xpatha, valueString

      n = xml%GetNumberOfNodes(TRIM(ADJUSTL(path))//'/q')
      ALLOCATE (q(3, n))
      DO i = 1, n
         !PRINT *, path,'/q[',i,']'
         WRITE (xPathA, "(a,a,i0,a)") TRIM(ADJUSTL(path)), '/q[', i, ']'
         !PRINT *,xpatha
         valueString = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA)))))
         !PRINT *,"Q:",valueString
         q(1,i)=evaluateFirst(valueString)
         q(2,i)=evaluateFirst(valueString)
         q(3,i)=evaluateFirst(valueString)
         !READ (valueString, *) q(1, i), q(2, i), q(3, i)
      END DO
   END FUNCTION read_q_list

   SUBROUTINE init_from_command_line()
      IMPLICIT NONE
      CHARACTER(len=1000)::xpath
      INTEGER:: i, ii

      IF (judft_was_argument("-xmlXPath")) THEN
         xpath = judft_string_for_argument("-xmlXPath")
         DO WHILE (INDEX(xpath, "=") > 0)
            i = INDEX(xpath, "=")
            ii = INDEX(xpath, ":")
            IF (ii == 0) ii = LEN(TRIM(xpath)) + 1
            IF (i > 100 .OR. ii - i > 100) CALL judft_error("Too long xmlXPath argument", calledby="xmlIntWarpFort.f90")
            CALL SetAttributeValue(xpath(:i - 1), xpath(i + 1:ii - 1))
            WRITE (*, *) "Set from command line:", TRIM(xpath(:i - 1)), "=", TRIM(xpath(i + 1:ii - 1))
            IF (ii + 1 < LEN(xpath)) THEN
               xpath = xpath(ii + 1:)
            ELSE
               xpath = ""
            ENDIF
         END DO
      END IF
   END SUBROUTINE init_from_command_line

   SUBROUTINE InitInterface()

      USE iso_c_binding

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
      IF (errorStatus .NE. 0) THEN
         CALL juDFT_error("Could not initialize XML interface.", calledby="xmlInitInterface")
      END IF

   END SUBROUTINE InitInterface

   SUBROUTINE ParseSchema(schemaFilename)

      USE iso_c_binding

      IMPLICIT NONE

      CHARACTER(LEN=200, KIND=c_char), INTENT(IN) :: schemaFilename

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
      IF (errorStatus .NE. 0) THEN
         CALL juDFT_error("XML Schema file not parsable: "//TRIM(ADJUSTL(schemaFilename)), calledby="xmlParseSchema")
      END IF

   END SUBROUTINE ParseSchema

   SUBROUTINE ParseDoc(docFilename)

      USE iso_c_binding

      IMPLICIT NONE

      CHARACTER(LEN=200, KIND=c_char), INTENT(IN) :: docFilename

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
      IF (errorStatus .NE. 0) THEN
         CALL juDFT_error("XML document file not parsable: "//TRIM(ADJUSTL(docFilename)), calledby="xmlParseDoc", &
                          hint="If there are no details listed above this message use 'xmllint --xinclude "//TRIM(ADJUSTL(docFilename))//"'")
      END IF

   END SUBROUTINE ParseDoc

   SUBROUTINE ValidateDoc()

      USE iso_c_binding

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
      IF (errorStatus .NE. 0) THEN
         CALL juDFT_error("XML document cannot be validated against Schema.", calledby="xml%ValidateDoc", &
                          hint="If there are no details listed above this message use 'xmllint --xinclude --schema FleurInputSchema.xsd inp.xml'.")
      END IF

   END SUBROUTINE ValidateDoc

   SUBROUTINE InitXPath()

      USE iso_c_binding

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
      IF (errorStatus .NE. 0) THEN
         CALL juDFT_error("Could not initialize XPath.", calledby="InitXPath")
      END IF
      CALL init_from_command_line()
   END SUBROUTINE InitXPath

   FUNCTION GetNumberOfNodes(xml, xPath)

      USE iso_c_binding

      IMPLICIT NONE

      INTEGER :: GetNumberOfNodes
      CLASS(t_xml), INTENT(IN):: xml
      CHARACTER(LEN=*, KIND=c_char), INTENT(IN) :: xPath

      INTERFACE
         FUNCTION getNumberOfXMLNodes(xPathExpression) BIND(C, name="getNumberOfXMLNodes")
            USE iso_c_binding
            INTEGER(c_int) getNumberOfXMLNodes
            CHARACTER(kind=c_char) :: xPathExpression(*)
         END FUNCTION getNumberOfXMLNodes
      END INTERFACE

      GetNumberOfNodes = getNumberOfXMLNodes(trim(adjustl(xml%basepath))//TRIM(ADJUSTL(xPath))//C_NULL_CHAR)

   END FUNCTION GetNumberOfNodes

   FUNCTION GetAttributeValue(xml, xPath, l_nocheck)

      USE iso_c_binding

      IMPLICIT NONE

      CHARACTER(LEN=:), ALLOCATABLE :: GetAttributeValue
      CLASS(t_xml), INTENT(IN):: xml
      CHARACTER(LEN=*, KIND=c_char), INTENT(IN) :: xPath
      LOGICAL, INTENT(IN), OPTIONAL :: l_nocheck

      CHARACTER(LEN=1, KIND=c_char), POINTER, DIMENSION(:) :: valueFromC => NULL()
      CHARACTER*255 :: VALUE
      INTEGER :: length, errorStatus, i
      TYPE(c_ptr) :: c_string
      logical :: l_docheck

      INTERFACE
         FUNCTION getXMLAttributeValue(xPathExpression) BIND(C, name="getXMLAttributeValue")
            USE iso_c_binding
            CHARACTER(KIND=c_char) :: xPathExpression(*)
            TYPE(c_ptr) :: getXMLAttributeValue
         END FUNCTION getXMLAttributeValue
      END INTERFACE

      l_docheck = .not. present(l_nocheck)
      if (.not. l_docheck) l_docheck = .not. l_nocheck

      if (l_docheck) then
         IF (xml%GetNumberOfNodes(xPath) < 1) THEN
            call judft_warn("Invalid xPath:"//xPath)
            GetAttributeValue = ""
            RETURN
         ENDIF
      endif
      c_string = getXMLAttributeValue(trim(adjustl(xml%basepath))//TRIM(ADJUSTL(xPath))//C_NULL_CHAR)

      CALL C_F_POINTER(c_string, valueFromC, [255])
      IF (.NOT. C_ASSOCIATED(c_string)) THEN
         WRITE (*, *) 'Error in trying to obtain attribute value from XPath:'
         WRITE (*, *) TRIM(ADJUSTL(xPath))
         CALL juDFT_error("Attribute value could not be obtained.", calledby="xml%getAttributeValue")
      END IF

      VALUE = ''
      i = 1
      DO WHILE ((valueFromC(i) .NE. C_NULL_CHAR) .AND. (i .LE. 255))
         VALUE(i:i) = valueFromC(i)
         i = i + 1
      END DO
      length = i - 1

      GetAttributeValue = TRIM(ADJUSTL(VALUE(1:length)))

   END FUNCTION GetAttributeValue

   subroutine GetAttributeValue_List(xml, xPath, list)

      USE iso_c_binding

      IMPLICIT NONE

      CHARACTER(LEN=255), intent(out) :: List(:)
      CLASS(t_xml), INTENT(IN):: xml
      CHARACTER(LEN=*, KIND=c_char), INTENT(IN) :: xPath

      CHARACTER(LEN=1, KIND=c_char), POINTER, DIMENSION(:) :: valueFromC => NULL()
      CHARACTER*255 :: VALUE
      INTEGER :: length, errorStatus, i, n
      TYPE(c_ptr) :: c_string, node
      logical :: l_docheck

      INTERFACE
         FUNCTION getXMLAttributeValueNode(node) BIND(C, name="getXMLAttributeValueNode")
            USE iso_c_binding
            type(c_ptr), VALUE :: node
            TYPE(c_ptr) :: getXMLAttributeValueNode
         END FUNCTION getXMLAttributeValueNode
      END INTERFACE

      INTERFACE
         FUNCTION getXMLNextNode(node) BIND(C, name="getXMLNextNode")
            USE iso_c_binding
            type(c_ptr), VALUE :: node
            TYPE(c_ptr) :: getXMLNextNode
         END FUNCTION getXMLNextNode
      END INTERFACE
      INTERFACE
         FUNCTION getXMLNode(xPathExpression) BIND(C, name="getXMLNode")
            USE iso_c_binding
            CHARACTER(KIND=c_char) :: xPathExpression(*)
            TYPE(c_ptr) :: getXMLNode
         END FUNCTION getXMLNode
      END INTERFACE
      INTERFACE
         FUNCTION getXMLAttributeValue(xPathExpression) BIND(C, name="getXMLAttributeValue")
            USE iso_c_binding
            CHARACTER(KIND=c_char) :: xPathExpression(*)
            TYPE(c_ptr) :: getXMLAttributeValue
         END FUNCTION getXMLAttributeValue
      END INTERFACE

      DO n = 1, size(List)
         if (n == 1) THEN
            node = getXMLNode(xPath//'[1]'//C_NULL_CHAR)
         else
            node = getXMLNextNode(node)
         endif
         c_string = getXMLAttributeValueNode(node)
         !c_string = getXMLAttributeValue(xPath//'[1]'//C_NULL_CHAR)
         CALL C_F_POINTER(c_string, valueFromC, [255])
         IF (.NOT. C_ASSOCIATED(c_string)) THEN
            WRITE (*, *) 'Error in trying to obtain attribute value from XPath:'
            WRITE (*, *) TRIM(ADJUSTL(xPath)), ":", n
            CALL juDFT_error("Attribute value could not be obtained.", calledby="xml%getAttributeValue_List")
         END IF
         VALUE = ''
         i = 1
         DO WHILE ((valueFromC(i) .NE. C_NULL_CHAR) .AND. (i .LE. 255))
            VALUE(i:i) = valueFromC(i)
            i = i + 1
         END DO
         length = i - 1

         List(n) = TRIM(ADJUSTL(VALUE(1:length)))
      enddo

   END subroutine GetAttributeValue_List

   SUBROUTINE SetAttributeValue(xPath, VALUE)

      USE iso_c_binding

      IMPLICIT NONE

      CHARACTER(LEN=*, KIND=c_char), INTENT(IN) :: xPath
      CHARACTER(len=*, KIND=c_char), INTENT(IN) :: VALUE

      INTEGER :: errorStatus

      INTERFACE
         FUNCTION setXMLAttributeValue(xPathExpression, valueExpression) BIND(C, name="setXMLAttributeValue")
            USE iso_c_binding
            CHARACTER(KIND=c_char) :: xPathExpression(*)
            CHARACTER(KIND=c_char) :: valueExpression(*)
            INTEGER(c_int) :: setXMLAttributeValue
         END FUNCTION setXMLAttributeValue
      END INTERFACE

      errorStatus = setXMLAttributeValue(TRIM(ADJUSTL(xPath))//C_NULL_CHAR, TRIM(ADJUSTL(VALUE))//C_NULL_CHAR)
      IF (errorStatus .NE. 0) THEN
         WRITE (*, *) 'Error in trying to setting attribute value from XPath:'
         WRITE (*, *) TRIM(ADJUSTL(xPath))
         WRITE (*, *) TRIM(ADJUSTL(VALUE))
         CALL juDFT_error("Attribute value could not be set.", calledby="xmlSetAttributeValue")
      END IF

   END SUBROUTINE SetAttributeValue

   SUBROUTINE FreeResources()

      USE iso_c_binding

      IMPLICIT NONE

      INTEGER :: errorStatus

      INTERFACE
         FUNCTION freeXMLResources() BIND(C, name="freeXMLResources")
            USE iso_c_binding
            INTEGER freeXMLResources
         END FUNCTION freeXMLResources
      END INTERFACE
      return
      errorStatus = 0
      errorStatus = freeXMLResources()
      IF (errorStatus .NE. 0) THEN
         CALL juDFT_error("Could not free XML resources.", calledby="xmlFreeResources")
         STOP 'Error!'
      END IF

   END SUBROUTINE FreeResources

   subroutine writexml(fileNum)
      USE iso_c_binding
      interface
         subroutine write_xml_file() bind(C, name="write_xml_file")
         end subroutine
      end interface
      integer, intent(in):: filenum

      integer :: err
      logical :: firstline = .true.
      character(len=500):: line

      if (fileNum == 0) return !no dump if called without proper fileNum, i.e. from inpgen

      !This will dump the inp.xml into inp_dump.xml
      call write_xml_file()
      print *, "Now copying inp_dump.xml"
      open (99, file="inp_dump.xml")
      firstline = .true. !Do not copy first line of file
      write (fileNum, *) "  <!-- Now follows a dump of the inp.xml file after evaluating the Schema -->"
      DO
         READ (99, '(a)', iostat=err) line
         if (firstline) then
            firstline = .false.
            CYCLE
         endif
         if (err .ne. 0) exit
         write (fileNum, "(a,a)") "  ", trim(line)
      end do
      write (fileNum, *) "  <!-- END of dump of the inp.xml file -->"
      close (99, status='delete')
   end subroutine

END MODULE m_types_xml
