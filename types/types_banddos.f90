!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_banddos

   TYPE t_banddos
      LOGICAL :: dos =.false. 
      LOGICAL :: band =.false.
      LOGICAL :: l_mcd =.false.
      LOGICAL :: l_orb =.false.
      LOGICAL :: vacdos =.false.
      INTEGER :: ndir =0
      INTEGER :: orbCompAtom=0
      REAL    :: e1_dos=0.5
      REAL    :: e2_dos=-0.5
      REAL    :: sig_dos=0.015
      REAL    :: e_mcd_lo =-10.0
      REAL    :: e_mcd_up= 0.0
      LOGICAL :: unfoldband =.false.
      INTEGER :: s_cell_x
      INTEGER :: s_cell_y
      INTEGER :: s_cell_z
      REAL    :: alpha,beta,gamma !For orbital decomp. (was orbcomprot)
   END TYPE t_banddos
 END MODULE m_types_banddos
  
