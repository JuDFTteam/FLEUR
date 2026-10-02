!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module FFTW3
#ifdef CPP_FFTW
   use, intrinsic :: iso_c_binding
#endif
   implicit none
   private
#ifdef CPP_FFTW
   public :: fftw_plan_dft_3d, fftw_execute_dft, fftw_destroy_plan, fftw_alloc_complex, fftw_free, &
             fftw_forward, fftw_backward, fftw_measure
   include 'fftw3.f03'
#endif
end module FFTW3
