!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_make_kpoints
  use m_juDFT
  implicit none
contains
  subroutine make_kpoints(kpts,cell,sym,film,l_socorss,str)
    use m_types_kpts
    use m_types_cell
    use m_types_sym
    type(t_kpts),intent(out)::kpts
    type(t_cell),intent(in) ::cell
    type(t_sym),intent(in)  ::sym
    logical,intent(in)::l_socorss
    character(len=*),intent(inout)::str
  
    logical::tria,l_gamma,l_soc_or_ss
    l_soc_or_ss=l_socorss

    if (judft_was_argument("-k")) then
       if (len_trim(str)>1) call judft_error("Do not specify k-points in file and on command line")
       str=judft_string_for_argument("-k")
    end if

    str=adjustl(str)
    do while(index(str,'@')>0)
       if (index(str,'tria@')==1) then
          tria=.true.
          str=str(6:)
       endif
       if (index(str,'gamma@')==1) then
          l_gamma=.true.
          str=str(7:)
       endif
       if (index(str,'soc@')==1) then
          l_soc_or_ss=.true.
          str=str(5:)
       endif
    end do
    if (index(str,'den=')==1) THEN
       str=str(5:)
       read(str,*) den
       call kpts%init_by_density(den,cell,sym,film,tria,l_soc_or_ss,l_gamma)
    elseif(index(str,'nk=')==1) then
       str=str(4:)
       read(str,*) nk
       call kpts%init_by_number(nk,cell,sym,film,tria,l_soc_or_ss,l_gamma)
    elseif(index(str,'grid=')==1) then
       str=str(6:)
       read(str,*) grid
       call kpts%init_by_grid(grid,cell,sym,film,tria,l_soc_or_ss,l_gamma)
    elseif(index(str,'file')==1) then
       call kpts%init_by_kptsfile(cell,sym,film,tria,l_soc_or_ss,l_gamma)
    elseif(index(str,'special')==1) then
       call kpts%init_by_kptsfile(cell,sym,film,tria,l_soc_or_ss,l_gamma)
    elseif(len_trim(str)<1.or.index(adjustl(str),'#')==1) then
       call kpts%init_defaults(cell,sym,film,tria,l_soc_or_ss,l_gamma)
    else
       call judft_error(("Could not process -k argument:"//str))
    endif
    !set name
    if (index(str,"#")>0) then
       str=str(index(str,"#")+1:)
    else
       str='default'
    end if
  end subroutine make_kpoints
