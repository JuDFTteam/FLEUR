MODULE m_generic_txtio

   USE m_juDFT

   IMPLICIT NONE
   
   PUBLIC :: writeValue,comment,header

   PRIVATE

   INTERFACE writeValue
      PROCEDURE writeInt, writeReal, writeRealArray
   END INTERFACE

   INTERFACE comment 
      PROCEDURE add_comment
   END INTERFACE

   INTERFACE header 
      PROCEDURE add_header
   END INTERFACE
   
   
   !Format specifiers:
   INTEGER, PARAMETER            :: indent_before_key = 3
   INTEGER, PARAMETER            :: pos_numbers       = 20
   INTEGER, PARAMETER            :: float_width       = 15
   INTEGER, PARAMETER            :: decimal_places    = 8
   INTEGER, PARAMETER            :: int_width         = 6

CONTAINS
   SUBROUTINE add_comment(iounit,comment,indent)

      INTEGER,           INTENT(IN)  :: iounit 
      CHARACTER(len=*),  INTENT(IN)  :: comment
      INTEGER, OPTIONAL, INTENT(IN)  :: indent

      INTEGER ind
      CHARACTER(len=300) fmt
      ind = 0
      IF(PRESENT(indent)) ind = indent
      WRITE(fmt,'("(A1,TR",I1,",(A))")') ind

      WRITE(iounit,fmt) "#",TRIM(ADJUSTL(comment))

   END SUBROUTINE add_comment

   SUBROUTINE add_header(iounit,header,indent)

      INTEGER,           INTENT(IN)  :: iounit 
      CHARACTER(len=*),  INTENT(IN)  :: header
      INTEGER, OPTIONAL, INTENT(IN)  :: indent

      INTEGER ind
      ind = 0
      IF(PRESENT(indent)) ind = indent

      WRITE(iounit,"(A)") "#**********************************************************"
      CALL add_comment(iounit,header,ind)
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

END MODULE m_generic_txtio