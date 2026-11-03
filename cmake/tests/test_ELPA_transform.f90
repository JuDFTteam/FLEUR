!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      program test  !Basically test if the elpa module is available and if can be used to do a tranform to standard problems
      use elpa
      class(elpa_t),pointer :: elpa_obj
      integer:: err
      integer, parameter:: dp = selected_real_kind(15)
      real(kind=dp):: h(2,2),s(2,2),tmp(2,2)

      call elpa_obj%elpa_transform_generalized_double(h,s,tmp,.false.,err)

      end program
