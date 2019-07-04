!--------------------------------------------------------------------------------
! Copyright (c) 2017 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_make_forcetheo
  implicit none
 
contains
  subroutine make_forcetheo(forcetheo_data,cell,sym,atoms,forcetheo)
    use m_types
    USE m_types_forcetheo_extended
    TYPE(t_sym),      INTENT(IN)     :: sym
    TYPE(t_atoms),    INTENT(IN)     :: atoms
    TYPE(t_cell),     INTENT(IN)     :: cell
    TYPE(t_forcetheo_data),INTENT(IN):: forcetheo_data
    CLASS(t_forcetheo),ALLOCATABLE,INTENT(OUT):: forcetheo
    
    !Finish setup of forcetheorem
    SELECT CASE (forcetheo_data%mode)
    CASE(1)
       ALLOCATE(t_forcetheo_mae::forcetheo)
    CASE(2)
       ALLOCATE(t_forcetheo_dmi::forcetheo)
    CASE(3)
       ALLOCATE(t_forcetheo_jij::forcetheo)
    CASE(4)
       ALLOCATE(t_forcetheo_ssdisp::forcetheo)
    CASE default
       ALLOCATE(t_forcetheo::forcetheo)
    END SELECT
    
    SELECT TYPE(forcetheo)
    TYPE IS(t_forcetheo_mae)
       CALL forcetheo%init(forcetheo_data%theta,forcetheo_data%phi,cell,sym)
    TYPE IS(t_forcetheo_dmi)
       CALL forcetheo%init(forcetheo_data%qvec,forcetheo_data%theta,forcetheo_data%phi)
    TYPE IS(t_forcetheo_jij)
       CALL forcetheo%init(forcetheo_data%qvec,forcetheo_data%theta(1),atoms)
    TYPE IS(t_forcetheo_ssdisp)
       CALL forcetheo%init(forcetheo_data%qvec)
    END SELECT
  end subroutine make_forcetheo
end MODULE m_make_forcetheo
