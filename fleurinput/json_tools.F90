!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_json_tools
  USE m_judft
  IMPLICIT NONE
  PRIVATE
  interface json_print
     module procedure json_rprint,json_iprint,json_lprint
     module procedure json_print_rarray,json_print_iarray,json_print_larray
     module procedure json_print_rrarray,json_print_iiarray,json_print_llarray
     module procedure json_print_rrrarray,json_print_iiiarray,json_print_lllarray
  end interface json_print

  interface json_read
     module procedure json_rread,json_iread,json_lread
     module procedure json_read_rarray,json_read_iarray,json_read_larray
     module procedure json_read_rrarray,json_read_iiarray,json_read_llarray
     module procedure json_read_rrrarray,json_read_iiiarray,json_read_lllarray
  end interface json_read
  PUBLIC:: json_read, json_print,json_open_class,json_close_class

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
  SUBROUTINE json_close_class(unit,iostat)
    INTEGER,INTENT(in)  :: unit
    INTEGER,INTENT(OUT) ::iostat
    
    CHARACTER(len=200)::line
    
    READ(unit,*,iostat=iostat) line

    IF (iostat==0.AND.(INDEX(line,"}")==0)) iostat=1
    RETURN
  END SUBROUTINE json_close_class

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

  !1-d array
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
    WRITE(unit,'(a,i0,a)') '"'//TRIM(name)//'_shape":',shape(val),','
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
    WRITE(unit,'(a,i0,a)') '"'//TRIM(name)//'_shape":',shape(val),','
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
    WRITE(unit,'(a,i0,a)') '"'//TRIM(name)//'_shape":',shape(val),','
    WRITE(unit,'(a)') '"'//TRIM(name)//'": ['
    WRITE(unit,"(5(i0,','))") val
    WRITE(unit,'(a)') ' ]'//endchar    
  END SUBROUTINE json_print_iarray

    !2-d array
  SUBROUTINE json_print_rrarray(unit,name,val,last)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    REAL,INTENT(in)     :: val(:,:)
    LOGICAL,OPTIONAL,INTENT(in)::last

    character :: endchar
    endchar=","
    IF (PRESENT(last)) THEN
       IF (last) endchar=" "
    END IF
    WRITE(unit,'(a,2i0,a)') '"'//TRIM(name)//'_shape":',shape(val),','
    WRITE(unit,'(a)') '"'//TRIM(name)//'": ['
    WRITE(unit,"(5(f0.15,','))") val
    WRITE(unit,'(a)') ' ]'//endchar
    
  END SUBROUTINE json_print_rrarray

  SUBROUTINE json_print_llarray(unit,name,val,last)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    LOGICAL,INTENT(in)     :: val(:,:)
    LOGICAL,OPTIONAL,INTENT(in)::last

    character :: endchar
    endchar=","
    IF (PRESENT(last)) THEN
       IF (last) endchar=" "
    END IF
    WRITE(unit,'(a,2i0,a)') '"'//TRIM(name)//'_shape":',shape(val),','
    WRITE(unit,'(a)') '"'//TRIM(name)//'": ['
    WRITE(unit,"(5(l1,','))") val
    WRITE(unit,'(a)') ' ]'//endchar
    
  END SUBROUTINE json_print_llarray

  SUBROUTINE json_print_iiarray(unit,name,val,last)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    INTEGER,INTENT(in)     :: val(:,:)
    LOGICAL,OPTIONAL,INTENT(in)::last

    character :: endchar
    endchar=","
    IF (PRESENT(last)) THEN
       IF (last) endchar=" "
    END IF
    WRITE(unit,'(a,2i0,a)') '"'//TRIM(name)//'_shape":',shape(val),','
    WRITE(unit,'(a)') '"'//TRIM(name)//'": ['
    WRITE(unit,"(5(i0,','))") val
    WRITE(unit,'(a)') ' ]'//endchar    
  END SUBROUTINE json_print_iiarray
      !3-d array
  SUBROUTINE json_print_rrrarray(unit,name,val,last)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    REAL,INTENT(in)     :: val(:,:,:)
    LOGICAL,OPTIONAL,INTENT(in)::last

    character :: endchar
    endchar=","
    IF (PRESENT(last)) THEN
       IF (last) endchar=" "
    END IF
    WRITE(unit,'(a,3i0,a)') '"'//TRIM(name)//'_shape":',shape(val),','
    WRITE(unit,'(a)') '"'//TRIM(name)//'": ['
    WRITE(unit,"(5(f0.15,','))") val
    WRITE(unit,'(a)') ' ]'//endchar
    
  END SUBROUTINE json_print_rrrarray

  SUBROUTINE json_print_lllarray(unit,name,val,last)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    LOGICAL,INTENT(in)     :: val(:,:,:)
    LOGICAL,OPTIONAL,INTENT(in)::last

    character :: endchar
    endchar=","
    IF (PRESENT(last)) THEN
       IF (last) endchar=" "
    END IF
    WRITE(unit,'(a,3i0,a)') '"'//TRIM(name)//'_shape":',SHAPE(val),','
    WRITE(unit,'(a)') '"'//TRIM(name)//'": ['
    WRITE(unit,"(5(l1,','))") val
    WRITE(unit,'(a)') ' ]'//endchar
    
  END SUBROUTINE json_print_lllarray

  SUBROUTINE json_print_iiiarray(unit,name,val,last)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    INTEGER,INTENT(in)     :: val(:,:,:)
    LOGICAL,OPTIONAL,INTENT(in)::last

    character :: endchar
    endchar=","
    IF (PRESENT(last)) THEN
       IF (last) endchar=" "
    END IF
    WRITE(unit,'(a,3i0,a)') '"'//TRIM(name)//'_shape":',shape(val),','
    WRITE(unit,'(a)') '"'//TRIM(name)//'": ['
    WRITE(unit,"(5(i0,','))") val
    WRITE(unit,'(a)') ' ]'//endchar    
  END SUBROUTINE json_print_iiiarray

  !Now array reading routines

  
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
  !1-D
  SUBROUTINE json_read_rarray(unit,name,val)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    REAL,ALLOCATABLE,INTENT(out)    :: val(:)
    INTEGER           :: n(1)
    CHARACTER(len=200):: line
    IF (ALLOCATED(val)) DEALLOCATE(val)
    read(unit,*)line
    IF (INDEX(line,name//'_shape')==0) CALL judft_warn("Key not found in json file:"//name)
    line=line(INDEX(line,":")+1:)
    READ(line,*) n
    ALLOCATE(val(n(1)))
    READ(unit,*)line
    IF (INDEX(line,name//'_shape')==0) CALL judft_warn("Key not found in json file:"//name)
    read(unit,*) val
    read(unit,*)line
    if (index(line,"]")==0) call judft_warn("Malformed json file:"//name)
  END SUBROUTINE json_read_rarray

  SUBROUTINE json_read_iarray(unit,name,val)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    integer,ALLOCATABLE,INTENT(out)    :: val(:)

    CHARACTER(len=200):: line
    integer           :: n(1)
    IF (ALLOCATED(val)) DEALLOCATE(val)
    read(unit,*)line
    IF (INDEX(line,name//'_shape')==0) CALL judft_warn("Key not found in json file:"//name)
    line=line(INDEX(line,":")+1:)
    READ(line,*) n
    ALLOCATE(val(n(1)))
    READ(unit,*)line
    if (index(line,name)==0) call judft_warn("Key not found in json file:"//name)
    read(unit,*) val
    read(unit,*)line
    if (index(line,"]")==0) call judft_warn("Malformed json file:"//name)
  END SUBROUTINE json_read_iarray

  SUBROUTINE json_read_larray(unit,name,val)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    logical,ALLOCATABLE,INTENT(out)    :: val(:)

    character(len=200):: line
    INTEGER           :: n(1)
    IF (ALLOCATED(val)) DEALLOCATE(val)
    read(unit,*)line
    IF (INDEX(line,name//'_shape')==0) CALL judft_warn("Key not found in json file:"//name)
    line=line(INDEX(line,":")+1:)
    READ(line,*) n
    ALLOCATE(val(n(1)))
    READ(unit,*)line
    if (index(line,name)==0) call judft_warn("Key not found in json file:"//name)
    read(unit,*) val
    read(unit,*)line
    if (index(line,"]")==0) call judft_warn("Malformed json file:"//name)
  END SUBROUTINE json_read_larray
 !2-D
  SUBROUTINE json_read_rrarray(unit,name,val)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    REAL,ALLOCATABLE,INTENT(out)    :: val(:,:)
    INTEGER           :: n(2)
    CHARACTER(len=200):: line
    IF (ALLOCATED(val)) DEALLOCATE(val)
    read(unit,*)line
    IF (INDEX(line,name//'_shape')==0) CALL judft_warn("Key not found in json file:"//name)
    line=line(INDEX(line,":")+1:)
    READ(line,*) n
    ALLOCATE(val(n(1),n(2)))
    READ(unit,*)line
    IF (INDEX(line,name//'_shape')==0) CALL judft_warn("Key not found in json file:"//name)
    read(unit,*) val
    read(unit,*)line
    if (index(line,"]")==0) call judft_warn("Malformed json file:"//name)
  END SUBROUTINE json_read_rrarray

  SUBROUTINE json_read_iiarray(unit,name,val)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    INTEGER,ALLOCATABLE,INTENT(out)    :: val(:,:)

    CHARACTER(len=200):: line
    integer           :: n(2)
    IF (ALLOCATED(val)) DEALLOCATE(val)
    read(unit,*)line
    IF (INDEX(line,name//'_shape')==0) CALL judft_warn("Key not found in json file:"//name)
    line=line(INDEX(line,":")+1:)
    READ(line,*) n
    ALLOCATE(val(n(1),n(2)))
    READ(unit,*)line
    if (index(line,name)==0) call judft_warn("Key not found in json file:"//name)
    read(unit,*) val
    read(unit,*)line
    if (index(line,"]")==0) call judft_warn("Malformed json file:"//name)
  END SUBROUTINE json_read_iiarray

  SUBROUTINE json_read_llarray(unit,name,val)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    LOGICAL,ALLOCATABLE,INTENT(out)    :: val(:,:)

    character(len=200):: line
    INTEGER           :: n(2)
    IF (ALLOCATED(val)) DEALLOCATE(val)
    read(unit,*)line
    IF (INDEX(line,name//'_shape')==0) CALL judft_warn("Key not found in json file:"//name)
    line=line(INDEX(line,":")+1:)
    READ(line,*) n
    ALLOCATE(val(n(1),n(2)))
    READ(unit,*)line
    if (index(line,name)==0) call judft_warn("Key not found in json file:"//name)
    read(unit,*) val
    read(unit,*)line
    if (index(line,"]")==0) call judft_warn("Malformed json file:"//name)
  END SUBROUTINE json_read_llarray
   !3-D
  SUBROUTINE json_read_rrrarray(unit,name,val)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    REAL,ALLOCATABLE,INTENT(out)    :: val(:,:,:)
    INTEGER           :: n(3)
    CHARACTER(len=200):: line
    IF (ALLOCATED(val)) DEALLOCATE(val)
    read(unit,*)line
    IF (INDEX(line,name//'_shape')==0) CALL judft_warn("Key not found in json file:"//name)
    line=line(INDEX(line,":")+1:)
    READ(line,*) n
    ALLOCATE(val(n(1),n(2),n(3)))
    READ(unit,*)line
    IF (INDEX(line,name//'_shape')==0) CALL judft_warn("Key not found in json file:"//name)
    read(unit,*) val
    read(unit,*)line
    if (index(line,"]")==0) call judft_warn("Malformed json file:"//name)
  END SUBROUTINE json_read_rrrarray

  SUBROUTINE json_read_iiiarray(unit,name,val)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    INTEGER,ALLOCATABLE,INTENT(out)    :: val(:,:,:)

    CHARACTER(len=200):: line
    integer           :: n(3)
    IF (ALLOCATED(val)) DEALLOCATE(val)
    read(unit,*)line
    IF (INDEX(line,name//'_shape')==0) CALL judft_warn("Key not found in json file:"//name)
    line=line(INDEX(line,":")+1:)
    READ(line,*) n
    ALLOCATE(val(n(1),n(2),n(3)))
    READ(unit,*)line
    if (index(line,name)==0) call judft_warn("Key not found in json file:"//name)
    read(unit,*) val
    read(unit,*)line
    if (index(line,"]")==0) call judft_warn("Malformed json file:"//name)
  END SUBROUTINE json_read_iiiarray

  SUBROUTINE json_read_lllarray(unit,name,val)
    CHARACTER(len=*)    :: name
    INTEGER,INTENT(in)  :: unit
    LOGICAL,ALLOCATABLE,INTENT(out)    :: val(:,:,:)

    character(len=200):: line
    INTEGER           :: n(3)
    IF (ALLOCATED(val)) DEALLOCATE(val)
    read(unit,*)line
    IF (INDEX(line,name//'_shape')==0) CALL judft_warn("Key not found in json file:"//name)
    line=line(INDEX(line,":")+1:)
    READ(line,*) n
    ALLOCATE(val(n(1),n(2),n(3)))
    READ(unit,*)line
    if (index(line,name)==0) call judft_warn("Key not found in json file:"//name)
    read(unit,*) val
    read(unit,*)line
    if (index(line,"]")==0) call judft_warn("Malformed json file:"//name)
  END SUBROUTINE json_read_lllarray
end MODULE m_json_tools
