!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_vecforlo
  USE m_juDFT
  !----------------------------------------------------------------------------
  ! For a given atom (na) set up 2*llo(lo)+1 k-vectors for each lo on this atom.
  ! if it is the first of two inversion-related atoms, set up twice as much.
  !
  ! nkvec(nlod,1) = number of k-vectors that were set up
  ! kvec(2*(2*llod+1),nlod) = index of these k-vectors. Stored as kveclo on the
  !                           eig-file, for later use in charge-density part.
  !----------------------------------------------------------------------------
CONTAINS

  
      !MOVED to types_lapw


END MODULE m_vecforlo
