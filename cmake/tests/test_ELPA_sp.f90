!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      program test  !Basically test if the elpa module is available and if it can diagonalize a single precision matrix
      use elpa
      class(elpa_t),pointer :: elpa_obj
      integer:: err
      integer, parameter:: sp = selected_real_kind(6)
      real(kind=sp):: h(2,2),e(2),z(2,2)

      call elpa_obj%eigenvectors(h,e,z,err)
      
      end program
