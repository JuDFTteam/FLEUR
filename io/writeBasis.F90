!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_writeBasis

CONTAINS

SUBROUTINE writeBasis()

   USE m_types
   USE m_juDFT

   IMPLICIT NONE

#ifdef CPP_HDF

#else
   CALL juDFT_error("writeBasis called without HDF5! ",calledby="writeBasis")
#endif

END SUBROUTINE writeBasis

END MODULE m_writeBasis
