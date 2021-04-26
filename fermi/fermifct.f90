!--------------------------------------------------------------------------------
! Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_fermifct
CONTAINS
  real elemental function fermifct(e,e_f,tkb)
    implicit NONE
    REAL, INTENT(IN) :: e,e_f,tkb

    REAL:: expo
    expo= EXP(-ABS(e-e_f)/tkb)
    IF (e<e_f) THEN
       fermifct = 1./ (expo+1.)
    ELSE
       fermifct= expo/ (expo+1.)
    ENDIF
  end function
end module
