!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
program test_prec
    integer,PARAMETER:: d=selected_real_kind(14)

    if(kind(1.0)/=kind(1.0_d)) stop 1
    if(kind(1.0d0)/=kind(1.0_d)) stop 2
end
