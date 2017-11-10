!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_tlmplm
    TYPE t_tlmplm
     COMPLEX,ALLOCATABLE :: tdd(:,:,:)
     COMPLEX,ALLOCATABLE :: tdu(:,:,:)
     !(0:lmplmd,ntypd,tspin)
     COMPLEX,ALLOCATABLE :: tud(:,:,:)
     COMPLEX,ALLOCATABLE :: tuu(:,:,:)
     !(0:lmplmd,ntypd,tspin)
     INTEGER,ALLOCATABLE :: ind(:,:,:,:)
     !(0:lmd,0:lmd,ntypd,tspin)
     COMPLEX,ALLOCATABLE :: tdulo(:,:,:,:)
     !(0:lmd,-llod:llod,mlotot,tspin)
     COMPLEX,ALLOCATABLE :: tuulo(:,:,:,:)
     !(0:lmd,-llod:llod,mlotot,tspin)
     COMPLEX,ALLOCATABLE :: tuloulo(:,:,:,:)
     !(-llod:llod,-llod:llod,mlolotot,tspin)
     COMPLEX,ALLOCATABLE :: h_loc(:,:,:,:)
     REAL                :: e_shift
  END TYPE t_tlmplm
END MODULE m_types_tlmplm
