!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This module is used to write out data that can be used to generate augmented
! DOS or bandstructure plots. For the augmentation additional data, e.g., weights
! for the orbital character of states, is stored.
!
!                                          GM' 2018
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE m_banddos_io

#ifdef CPP_HDF

   PUBLIC openBandDOSFile, closeBandDOSFile, writeBandDOSData

   CONTAINS

   SUBROUTINE openBandDOSFile()
   END SUBROUTINE

   SUBROUTINE closeBandDOSFile()
   END SUBROUTINE

   SUBROUTINE writeBandDOSData()
   END SUBROUTINE

#endif

END MODULE m_banddos_io
