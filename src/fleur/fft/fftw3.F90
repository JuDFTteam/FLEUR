!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module FFTW3
#ifdef CPP_FFTW
   use, intrinsic :: iso_c_binding
   implicit none
   include 'fftw3.f03'
#endif
   implicit none
end module FFTW3