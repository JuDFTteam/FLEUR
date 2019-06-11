!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_noco
  
 TYPE t_noco
    LOGICAL:: l_ss= .FALSE.
    LOGICAL:: l_soc= .FALSE.
    LOGICAL:: l_noco = .FALSE.
    LOGICAL:: l_mperp = .FALSE.
    LOGICAL:: l_constr = .FALSE.
    LOGICAL:: l_mtNocoPot = .FALSE.
    REAL   :: qss(3)=[0.,0.,0.]
    REAL   :: mix_b=0.0
    LOGICAL:: l_spav= .FALSE.
    REAL   :: theta=0.0
    REAL   :: phi=0.0
    
    LOGICAL, ALLOCATABLE :: l_relax(:)
    REAL, ALLOCATABLE    :: alphInit(:)
    REAL, ALLOCATABLE    :: alph(:)
    REAL, ALLOCATABLE    :: beta(:)
    REAL, ALLOCATABLE    :: b_con(:, :)
    REAL, ALLOCATABLE    :: socscale(:)
   END TYPE t_noco

 END MODULE m_types_noco
