MODULE m_generic_txtio

   !------------------------------------------------------------------------------
   !
   ! MODULE:  m_generic_txtio
   !
   !> @author
   !> Henning Janßen
   !
   ! DESCRIPTION: 
   !>  Contains relatively general subroutines for writing txt files with keywords
   !>  written for hubbard 1 solver 
   !
   !------------------------------------------------------------------------------
   USE m_juDFT

   IMPLICIT NONE
   
   PUBLIC :: writeValue,comment,header,startSection,endSection

   PRIVATE

   INTERFACE writeValue
      PROCEDURE writeInt, writeReal, writeRealArray, writeKeyword, writeCharacter
   END INTERFACE

   INTERFACE comment 
      PROCEDURE add_comment
   END INTERFACE

   INTERFACE header 
      PROCEDURE add_header
   END INTERFACE
   
   
   !Format specifiers:
   INTEGER, PARAMETER            :: indent_before_key = 3
   INTEGER, PARAMETER            :: pos_numbers       = 18
   INTEGER, PARAMETER            :: float_width       = 15
   INTEGER, PARAMETER            :: decimal_places    = 8
   INTEGER, PARAMETER            :: int_width         = 6

CONTAINS
   SUBROUTINE add_comment(iounit,comment_str,indent)

      INTEGER,           INTENT(IN)  :: iounit 
      CHARACTER(len=*),  INTENT(IN)  :: comment_str
      INTEGER, OPTIONAL, INTENT(IN)  :: indent

      INTEGER ind
      CHARACTER(len=300) fmt
      ind = 0
      IF(PRESENT(indent)) ind = indent
      WRITE(fmt,'("(A1,TR",I1,",(A))")') ind

      WRITE(iounit,fmt) "#",TRIM(ADJUSTL(comment_str))

   END SUBROUTINE add_comment

   SUBROUTINE add_header(iounit,header_str,indent)

      INTEGER,           INTENT(IN)  :: iounit 
      CHARACTER(len=*),  INTENT(IN)  :: header_str
      INTEGER, OPTIONAL, INTENT(IN)  :: indent

      INTEGER ind
      ind = 0
      IF(PRESENT(indent)) ind = indent

      WRITE(iounit,"(A)") "#**********************************************************"
      CALL add_comment(iounit,header_str,ind)
      WRITE(iounit,"(A)") "#**********************************************************"

   END SUBROUTINE add_header

   SUBROUTINE writeReal(iounit,key,value)

      IMPLICIT NONE

      INTEGER,             INTENT(IN)  :: iounit
      CHARACTER(len=*),    INTENT(IN)  :: key 
      REAL,                INTENT(IN)  :: value

      INTEGER indent_after_key
      CHARACTER(len=300) fmt

      indent_after_key = pos_numbers-LEN(key)-indent_before_key

      IF(indent_after_key.LT.1) CALL juDFT_error("indent_after_key<0",calledby="writeReal")
      !Define Format specifier
      WRITE(fmt,'("(TR",I2.2,",A",I2.2,",TR",I2.2,",f",I2.2,".",I2.2,")")') &
            indent_before_key,LEN(key),indent_after_key,float_width,decimal_places

      WRITE(iounit,fmt) TRIM(ADJUSTL(key)), value

   END SUBROUTINE writeReal

   SUBROUTINE writeRealArray(iounit,key,value)

      IMPLICIT NONE

      INTEGER,             INTENT(IN)  :: iounit
      CHARACTER(len=*),    INTENT(IN)  :: key 
      REAL,                INTENT(IN)  :: value(:)

      INTEGER indent_after_key
      CHARACTER(len=300) fmt
      indent_after_key = pos_numbers-LEN(key)-indent_before_key

      IF(indent_after_key.LT.1) CALL juDFT_error("indent_after_key<0",calledby="writeRealArray")
      !Define Format specifier
      WRITE(fmt,'("(TR",I2.2,",A",I2.2,",TR",I2.2,",",I2.2,"f",I2.2,".",I2.2,")")') &
            indent_before_key,LEN(key),indent_after_key,SIZE(value,1),float_width,decimal_places

      WRITE(iounit,fmt) TRIM(ADJUSTL(key)), value

   END SUBROUTINE writeRealArray

   SUBROUTINE writeInt(iounit,key,value)

      IMPLICIT NONE

      INTEGER,             INTENT(IN)  :: iounit
      CHARACTER(len=*),    INTENT(IN)  :: key 
      INTEGEr,             INTENT(IN)  :: value

      INTEGER indent_after_key
      CHARACTER(len=300) fmt

      indent_after_key = pos_numbers-LEN(key)-indent_before_key

      IF(indent_after_key.LT.1) CALL juDFT_error("indent_after_key<0",calledby="writeInt")
      !Define Format specifier
      WRITE(fmt,'("(TR",I2.2,",A",I2.2,",TR",I2.2,",I",I2.2,")")') &
            indent_before_key,LEN(key),indent_after_key,int_width

      WRITE(iounit,fmt) TRIM(ADJUSTL(key)), value

   END SUBROUTINE writeInt

   SUBROUTINE writeKeyword(iounit,key)

      INTEGER,          INTENT(IN)  :: iounit
      CHARACTER(len=*), INTENT(IN)  :: key

      CHARACTER(len=300) fmt

      WRITE(fmt,'("(TR",I2.2,",A",I2.2,")")') indent_before_key,LEN(key)

      WRITE(iounit,fmt) key

   END SUBROUTINE writeKeyword

   SUBROUTINE writeCharacter(iounit,key,value)

      INTEGER,          INTENT(IN)  :: iounit
      CHARACTER(len=*), INTENT(IN)  :: key
      CHARACTER(len=*), INTENT(IN)  :: value

      CHARACTER(len=300) fmt
      INTEGER indent_after_key

      indent_after_key = pos_numbers-LEN(key)-indent_before_key + 5

      WRITE(fmt,'("(TR",I2.2,",A",I2.2,",TR",I2.2,",A",I2.2,")")') &
            indent_before_key,LEN(value),indent_after_key,LEN(value)
         
      WRITE(iounit,fmt) TRIM(ADJUSTL(key)),TRIM(ADJUSTL(value))
   END SUBROUTINE writeCharacter
      
   SUBROUTINE startSection(iounit,sectionname)

      IMPLICIT NONE

      INTEGER,             INTENT(IN)  :: iounit 
      CHARACTER(len=*),    INTENT(IN)  :: sectionname

      CHARACTER(len=300) fmt

      WRITE(fmt,'("(A",I2.2,",TR1,A1)")') LEN(sectionname)
      WRITE(iounit,fmt) sectionname , "{" 

   END SUBROUTINE startSection

   SUBROUTINE endSection(iounit)

      IMPLICIT NONE

      INTEGER,             INTENT(IN)  :: iounit 

      WRITE(iounit,"(A1)") "}" 

   END SUBROUTINE endSection


END MODULE m_generic_txtio