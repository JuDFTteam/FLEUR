!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_json_tools
CONTAINS
  SUBROUTINE json_open_class(name,unit,iostat)
    CHARACTER(len=*)    ::name
    INTEGER,INTENT(in)  :: unit
    INTEGER,INTENT(OUT) ::iostat

    CHARACTER(len=200)::line

    READ(unit,*,iostat) line

    IF (iostat==0.AND.(INDEX(line,name)==0)) iostat=1
    RETURN
  END SUBROUTINE json_open_class

  SUBROUTINE json_rprint(unit,name,val,last)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    REAL,INTENT(in)     :: val
    LOGICAL,OPTIONAL,INTENT(in)::last

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

    endchar=","
    IF (PRESENT(last)) THEN
       IF (last) endchar=" "
    END IF
    WRITE(unit,'(a)') '"'//TRIM(name)//'": ['
    WRITE(unit,"5(f0.15,','))") val
    WRITE(unit,'(a)') ' ]'//endchar
    
  END SUBROUTINE json_print_rarray

  SUBROUTINE json_print_larray(unit,name,val,last)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    logical,INTENT(in)     :: val(:)
    LOGICAL,OPTIONAL,INTENT(in)::last

    endchar=","
    IF (PRESENT(last)) THEN
       IF (last) endchar=" "
    END IF
    WRITE(unit,'(a)') '"'//TRIM(name)//'": ['
    WRITE(unit,"5(l1,','))") val
    WRITE(unit,'(a)') ' ]'//endchar
    
  END SUBROUTINE json_print_larray

  SUBROUTINE json_print_iarray(unit,name,val,last)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    INTEGER,INTENT(in)     :: val(:)
    LOGICAL,OPTIONAL,INTENT(in)::last

    endchar=","
    IF (PRESENT(last)) THEN
       IF (last) endchar=" "
    END IF
    WRITE(unit,'(a)') '"'//TRIM(name)//'": ['
    WRITE(unit,"5(i0,','))") val
    WRITE(unit,'(a)') ' ]'//endchar
    
  END SUBROUTINE json_print_iarray
