!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

! NOTE: this contains only the interface, the actual code is included by the
! proprocessor from the file zsymsecloc_cpp.F90
!

MODULE m_zsymsecloc
  use m_juDFT
!*******************************************************
!  Solve the generalized secular equation. 
!  For film-systems exhibiting
!  z-reflexion symmetry, the basis is transformed to
!  even and odd functions and the even-even and odd-odd 
!  blocks are diagonalized separately.
!  If local orbitals are present in a film with z-reflection,
!  locrectify is used to construct linear combinations of
!  the local orbitals that are eigenfunctions of the z-
!  reflexion operation.
!  Frank Freimuth, January 2006
!*******************************************************
  INTERFACE zsymsecloc
     MODULE procedure zsymsecloc_r,zsymsecloc_c
  END INTERFACE zsymsecloc
CONTAINS
  SUBROUTINE zsymsecloc_r(jsp,input,lapw,bkpt,atoms, kveclo, sym,cell, dimension,matsize, nsize, jij,matind,nred,eig,ne, a,b, z)

#define CPP_REALDATA
#include "zsymsecloc_cpp.F90"
  END SUBROUTINE zsymsecloc_r

  SUBROUTINE zsymsecloc_c(jsp,input,lapw,bkpt,atoms, kveclo, sym,cell, dimension,matsize, nsize, jij,matind,nred,eig,ne, a,b, z)

#undef CPP_REALDATA
#include "zsymsecloc_cpp.F90"
  END SUBROUTINE zsymsecloc_c

END MODULE m_zsymsecloc
