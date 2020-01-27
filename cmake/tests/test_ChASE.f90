!--------------------------------------------------------------------------------
! Copyright (c) 2017 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
! Author: Miriam Hinzen

PROGRAM test

  use iso_c_binding

  interface 
    subroutine chase_r( h, n, v, ritzv, nev, nex, deg, tol, mode, opt ) bind( c, name = 'dchase_' )
      use iso_c_binding
      real(c_double)                :: h(n,*), v(n,*)
      integer(c_int)                :: n, deg, nev, nex
      real(c_double)                :: ritzv(*), tol
      character(len=1,kind=c_char)  :: mode, opt
    end subroutine chase_r
  end interface

  real(c_double), ALLOCATABLE   :: h(:,:), v(:,:), ritzv(:)
  integer(c_int)                :: n, deg, nev, nex
  real(c_double)                :: tol
  character(len=1,kind=c_char)  :: mode, opt

  CALL chase_r( h, n, v, ritzv, nev, nex, deg, tol, mode, opt )

END
