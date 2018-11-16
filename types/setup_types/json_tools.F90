!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_json_tools
  IMPLICIT NONE
  PRIVATE
  interface json_print
     module procedure json_rprint,json_iprint,json_lprint
     module procedure json_print_rarray,json_print_iarray,json_print_larray
  end interface json_print

  interface json_read
     module procedure json_rread,json_iread,json_lread
     module procedure json_read_rarray,json_read_iarray,json_read_larray
  end interface json_read
  PUBLIC:: json_read, json_write,json_open_class,json_close_class

CONTAINS
  SUBROUTINE json_open_class(name,unit,iostat)
    CHARACTER(len=*)    ::name
    INTEGER,INTENT(in)  :: unit
    INTEGER,INTENT(OUT) ::iostat

    CHARACTER(len=200)::line

    READ(unit,*,iostat=iostat) line

    IF (iostat==0.AND.(INDEX(line,name)==0)) iostat=1
    RETURN
  END SUBROUTINE json_open_class

  SUBROUTINE json_rprint(unit,name,val,last)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    REAL,INTENT(in)     :: val
    LOGICAL,OPTIONAL,INTENT(in)::last

    character :: endchar
    endchar=","
    IF (PRESENT(last)) THEN
       IF (last) endchar=" "
    END IF
    WRITE(unit,'(a,f0.15,a)') '"'//TRIM(name)//'":',val,endchar
  END SUBROUTINE json_rprint

 

  SUBROUTINE json_lprint(unit,name,val,last)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    LOGICAL,INTENT(in)     :: val
    LOGICAL,OPTIONAL,INTENT(in)::last

    character :: endchar
    endchar=","
    IF (PRESENT(last)) THEN
       IF (last) endchar=" "
    END IF
    WRITE(unit,'(a,l1,a)') '"'//TRIM(name)//'":',val,endchar
  END SUBROUTINE json_lprint
  SUBROUTINE json_iprint(unit,name,val,last)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    INTEGER,INTENT(in)     :: val
    LOGICAL,OPTIONAL,INTENT(in)::last

    character :: endchar
    endchar=","
    IF (PRESENT(last)) THEN
       IF (last) endchar=" "
    END IF
    WRITE(unit,'(a,i0,a)') '"'//TRIM(name)//'":',val,endchar
  END SUBROUTINE json_iprint

  SUBROUTINE json_print_rarray(unit,name,val,last)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    REAL,INTENT(in)     :: val(:)
    LOGICAL,OPTIONAL,INTENT(in)::last

    character :: endchar
    endchar=","
    IF (PRESENT(last)) THEN
       IF (last) endchar=" "
    END IF
    WRITE(unit,'(a)') '"'//TRIM(name)//'": ['
    WRITE(unit,"(5(f0.15,','))") val
    WRITE(unit,'(a)') ' ]'//endchar
    
  END SUBROUTINE json_print_rarray

  SUBROUTINE json_print_larray(unit,name,val,last)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    logical,INTENT(in)     :: val(:)
    LOGICAL,OPTIONAL,INTENT(in)::last

    character :: endchar
    endchar=","
    IF (PRESENT(last)) THEN
       IF (last) endchar=" "
    END IF
    WRITE(unit,'(a)') '"'//TRIM(name)//'": ['
    WRITE(unit,"(5(l1,','))") val
    WRITE(unit,'(a)') ' ]'//endchar
    
  END SUBROUTINE json_print_larray

  SUBROUTINE json_print_iarray(unit,name,val,last)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    INTEGER,INTENT(in)     :: val(:)
    LOGICAL,OPTIONAL,INTENT(in)::last

    character :: endchar
    endchar=","
    IF (PRESENT(last)) THEN
       IF (last) endchar=" "
    END IF
    WRITE(unit,'(a)') '"'//TRIM(name)//'": ['
    WRITE(unit,"(5(i0,','))") val
    WRITE(unit,'(a)') ' ]'//endchar    
  END SUBROUTINE json_print_iarray

  SUBROUTINE json_rread(unit,name,val)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    REAL,INTENT(out)    :: val

    character(len=200):: line
    read(unit,*)line
    if (index(line,name)==0) call judft_warn("Key not found in json file:"//name)
    if (index(line,":")==0) call judft_warn("Malformed json file:"//name)
    line=line(index(line,":")+1:)

    read(line,*) val
  END SUBROUTINE json_rread

  SUBROUTINE json_iread(unit,name,val)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    INTEGER,INTENT(out)    :: val

    character(len=200):: line
    read(unit,*)line
    if (index(line,name)==0) call judft_warn("Key not found in json file:"//name)
    if (index(line,":")==0) call judft_warn("Malformed json file:"//name)
    line=line(index(line,":")+1:)

    read(line,*) val
  END SUBROUTINE json_iread

  SUBROUTINE json_lread(unit,name,val)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    logical,INTENT(out)    :: val

    character(len=200):: line
    read(unit,*)line
    if (index(line,name)==0) call judft_warn("Key not found in json file:"//name)
    if (index(line,":")==0) call judft_warn("Malformed json file:"//name)
    line=line(index(line,":")+1:)

    read(line,*) val
  END SUBROUTINE json_lread

  SUBROUTINE json_read_rarray(unit,name,val)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    real,INTENT(out)    :: val(:)

    character(len=200):: line
    read(unit,*)line
    if (index(line,name)==0) call judft_warn("Key not found in json file:"//name)
    read(unit,*) val
    read(unit,*)line
    if (index(line,"]")==0) call judft_warn("Malformed json file:"//name)
  END SUBROUTINE json_read_rarray

  SUBROUTINE json_read_iarray(unit,name,val)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    integer,INTENT(out)    :: val(:)

    character(len=200):: line
    read(unit,*)line
    if (index(line,name)==0) call judft_warn("Key not found in json file:"//name)
    read(unit,*) val
    read(unit,*)line
    if (index(line,"]")==0) call judft_warn("Malformed json file:"//name)
  END SUBROUTINE json_read_iarray

  SUBROUTINE json_read_larray(unit,name,val)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    logical,INTENT(out)    :: val(:)

    character(len=200):: line
    read(unit,*)line
    if (index(line,name)==0) call judft_warn("Key not found in json file:"//name)
    read(unit,*) val
    read(unit,*)line
    if (index(line,"]")==0) call judft_warn("Malformed json file:"//name)
  END SUBROUTINE json_read_larray

end MODULE m_json_tools
