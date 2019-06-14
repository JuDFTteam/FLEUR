!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_make_kpoints
  USE m_juDFT
  IMPLICIT NONE
CONTAINS
  SUBROUTINE make_kpoints(kpts,cell,sym,film,l_socorss,str)
    USE m_types_kpts
    USE m_types_cell
    USE m_types_sym
    TYPE(t_kpts),INTENT(out)::kpts
    TYPE(t_cell),INTENT(in) ::cell
    TYPE(t_sym),INTENT(in)  ::sym
    LOGICAL,INTENT(in)::l_socorss,film
    CHARACTER(len=*),INTENT(inout)::str

    LOGICAL:: tria,l_gamma,l_soc_or_ss
    REAL   :: den
    INTEGER:: nk,grid(3)
    character(len=20)::name=""
    l_soc_or_ss=l_socorss

    IF (judft_was_argument("-k")) THEN
       IF (LEN_TRIM(str)>1) CALL judft_error("Do not specify k-points in file and on command line")
       str=judft_string_for_argument("-k")
    END IF

    !set name
    IF (INDEX(str,"#")>0) THEN
       name=trim(adjustl(str(INDEX(str,"#")+1:)))
       str=str(:INDEX(str,"#")-1)
    END IF

    str=ADJUSTL(str)
    DO WHILE(INDEX(str,'@')>0)
       IF (INDEX(str,'tria@')==1) THEN
          tria=.TRUE.
          str=str(6:)
       ENDIF
       IF (INDEX(str,'gamma@')==1) THEN
          l_gamma=.TRUE.
          str=str(7:)
       ENDIF
       IF (INDEX(str,'soc@')==1) THEN
          l_soc_or_ss=.TRUE.
          str=str(5:)
       ENDIF
    END DO
    IF (INDEX(str,'den=')==1) THEN
       str=str(5:)
       READ(str,*) den
       CALL kpts%init_by_density(den,cell,sym,film,tria,l_soc_or_ss,l_gamma)
    ELSEIF(INDEX(str,'nk=')==1) THEN
       str=str(4:)
       READ(str,*) nk
       CALL kpts%init_by_number(nk,cell,sym,film,tria,l_soc_or_ss,l_gamma)
    ELSEIF(INDEX(str,'grid=')==1) THEN
       str=str(6:)
       READ(str,*) grid
       CALL kpts%init_by_grid(grid,cell,sym,film,tria,l_soc_or_ss,l_gamma)
    ELSEIF(INDEX(str,'file')==1) THEN
       CALL kpts%init_by_kptsfile(film)
    ELSEIF(LEN_TRIM(str)<1.OR.INDEX(ADJUSTL(str),'#')==1) THEN
       CALL kpts%init_defaults(cell,sym,film,tria,l_soc_or_ss,l_gamma)
    ELSE
       CALL judft_error(("Could not process -k argument:"//str))
    ENDIF
    
    if (len_trim(name)>0) kpts%name=name
  END SUBROUTINE make_kpoints
END MODULE m_make_kpoints
