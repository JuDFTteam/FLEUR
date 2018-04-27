!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_denCoeffsOffdiag

IMPLICIT NONE

PRIVATE

   TYPE t_denCoeffsOffdiag
      LOGICAL              :: l_fmpl

      ! spherical
      COMPLEX, ALLOCATABLE :: uu21(:,:)
      COMPLEX, ALLOCATABLE :: dd21(:,:)
      COMPLEX, ALLOCATABLE :: du21(:,:)
      COMPLEX, ALLOCATABLE :: ud21(:,:)

      ! nonspherical
      COMPLEX, ALLOCATABLE :: uunmt21(:,:,:)
      COMPLEX, ALLOCATABLE :: ddnmt21(:,:,:)
      COMPLEX, ALLOCATABLE :: dunmt21(:,:,:)
      COMPLEX, ALLOCATABLE :: udnmt21(:,:,:)

      ! spherical - LOs
      COMPLEX, ALLOCATABLE :: uulo21(:,:)
      COMPLEX, ALLOCATABLE :: dulo21(:,:)
      COMPLEX, ALLOCATABLE :: ulou21(:,:)
      COMPLEX, ALLOCATABLE :: ulod21(:,:)

      COMPLEX, ALLOCATABLE :: uloulop21(:,:,:)

      ! norms
      REAL, ALLOCATABLE     :: uu21n(:,:)
      REAL, ALLOCATABLE     :: ud21n(:,:)
      REAL, ALLOCATABLE     :: du21n(:,:)
      REAL, ALLOCATABLE     :: dd21n(:,:)

      REAL, ALLOCATABLE     :: uulo21n(:,:)
      REAL, ALLOCATABLE     :: dulo21n(:,:)
      REAL, ALLOCATABLE     :: ulou21n(:,:)
      REAL, ALLOCATABLE     :: ulod21n(:,:)

      REAL, ALLOCATABLE     :: uloulop21n(:,:,:)

      CONTAINS

      PROCEDURE,PASS :: init => denCoeffsOffdiag_init
      PROCEDURE      :: addRadFunScalarProducts
      PROCEDURE      :: calcCoefficients

   END TYPE t_denCoeffsOffdiag

PUBLIC t_denCoeffsOffdiag

CONTAINS

SUBROUTINE denCoeffsOffdiag_init(thisDenCoeffsOffdiag, atoms, noco, sphhar, l_fmpl)

   USE m_types_setup

   IMPLICIT NONE

   CLASS(t_denCoeffsOffdiag), INTENT(INOUT) :: thisDenCoeffsOffdiag
   TYPE(t_atoms),      INTENT(IN)    :: atoms
   TYPE(t_noco),       INTENT(IN)    :: noco
   TYPE(t_sphhar),     INTENT(IN)    :: sphhar
   LOGICAL,            INTENT(IN)    :: l_fmpl

   thisDenCoeffsOffdiag%l_fmpl = l_fmpl

   IF (noco%l_mperp) THEN
      ALLOCATE (thisDenCoeffsOffdiag%uu21(0:atoms%lmaxd,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%ud21(0:atoms%lmaxd,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%du21(0:atoms%lmaxd,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%dd21(0:atoms%lmaxd,atoms%ntype))

      ALLOCATE (thisDenCoeffsOffdiag%uulo21(atoms%nlod,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%dulo21(atoms%nlod,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%ulou21(atoms%nlod,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%ulod21(atoms%nlod,atoms%ntype))

      ALLOCATE (thisDenCoeffsOffdiag%uloulop21(atoms%nlod,atoms%nlod,atoms%ntype))

      ALLOCATE (thisDenCoeffsOffdiag%uu21n(0:atoms%lmaxd,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%ud21n(0:atoms%lmaxd,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%du21n(0:atoms%lmaxd,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%dd21n(0:atoms%lmaxd,atoms%ntype))

      ALLOCATE (thisDenCoeffsOffdiag%uulo21n(atoms%nlod,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%dulo21n(atoms%nlod,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%ulou21n(atoms%nlod,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%ulod21n(atoms%nlod,atoms%ntype))

      ALLOCATE (thisDenCoeffsOffdiag%uloulop21n(atoms%nlod,atoms%nlod,atoms%ntype))
   ELSE
      ALLOCATE (thisDenCoeffsOffdiag%uu21(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%ud21(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%du21(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%dd21(1,1))

      ALLOCATE (thisDenCoeffsOffdiag%uulo21(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%dulo21(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%ulou21(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%ulod21(1,1))

      ALLOCATE (thisDenCoeffsOffdiag%uloulop21(1,1,1))

      ALLOCATE (thisDenCoeffsOffdiag%uu21n(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%ud21n(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%du21n(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%dd21n(1,1))

      ALLOCATE (thisDenCoeffsOffdiag%uulo21n(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%dulo21n(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%ulou21n(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%ulod21n(1,1))

      ALLOCATE (thisDenCoeffsOffdiag%uloulop21n(1,1,1))
   END IF

   IF (noco%l_mperp.AND.l_fmpl) THEN
      ALLOCATE (thisDenCoeffsOffdiag%uunmt21((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%udnmt21((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%dunmt21((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%ddnmt21((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntype))
   ELSE
      ALLOCATE (thisDenCoeffsOffdiag%uunmt21(1,1,1))
      ALLOCATE (thisDenCoeffsOffdiag%udnmt21(1,1,1))
      ALLOCATE (thisDenCoeffsOffdiag%dunmt21(1,1,1))
      ALLOCATE (thisDenCoeffsOffdiag%ddnmt21(1,1,1))
   END IF

   thisDenCoeffsOffdiag%uu21 = CMPLX(0.0,0.0)
   thisDenCoeffsOffdiag%ud21 = CMPLX(0.0,0.0)
   thisDenCoeffsOffdiag%du21 = CMPLX(0.0,0.0)
   thisDenCoeffsOffdiag%dd21 = CMPLX(0.0,0.0)

   thisDenCoeffsOffdiag%uulo21 = CMPLX(0.0,0.0)
   thisDenCoeffsOffdiag%dulo21 = CMPLX(0.0,0.0)
   thisDenCoeffsOffdiag%ulou21 = CMPLX(0.0,0.0)
   thisDenCoeffsOffdiag%ulod21 = CMPLX(0.0,0.0)

   thisDenCoeffsOffdiag%uloulop21 = CMPLX(0.0,0.0)

   thisDenCoeffsOffdiag%uu21n = 0.0
   thisDenCoeffsOffdiag%ud21n = 0.0
   thisDenCoeffsOffdiag%du21n = 0.0
   thisDenCoeffsOffdiag%dd21n = 0.0

   thisDenCoeffsOffdiag%uulo21n = 0.0
   thisDenCoeffsOffdiag%dulo21n = 0.0
   thisDenCoeffsOffdiag%ulou21n = 0.0
   thisDenCoeffsOffdiag%ulod21n = 0.0

   thisDenCoeffsOffdiag%uloulop21n = 0.0

   thisDenCoeffsOffdiag%uunmt21 = CMPLX(0.0,0.0)
   thisDenCoeffsOffdiag%udnmt21 = CMPLX(0.0,0.0)
   thisDenCoeffsOffdiag%dunmt21 = CMPLX(0.0,0.0)
   thisDenCoeffsOffdiag%ddnmt21 = CMPLX(0.0,0.0)

END SUBROUTINE denCoeffsOffdiag_init

SUBROUTINE addRadFunScalarProducts(thisDenCoeffsOffdiag, atoms, f, g, flo, iType)

   USE m_types_setup
   USE m_int21       ! integrate (spin) off-diagonal radial functions
   USE m_int21lo     ! -"- for u_lo

   IMPLICIT NONE

   CLASS(t_denCoeffsOffdiag), INTENT(INOUT) :: thisDenCoeffsOffdiag
   TYPE(t_atoms),             INTENT(IN)    :: atoms
   REAL,                      INTENT(IN)    :: f(:,:,0:,:)!(atoms%jmtd,2,0:atoms%lmaxd,dimension%jspd)
   REAL,                      INTENT(IN)    :: g(:,:,0:,:)!(atoms%jmtd,2,0:atoms%lmaxd,dimension%jspd)
   REAL,                      INTENT(IN)    :: flo(:,:,:,:)!(atoms%jmtd,2,atoms%nlod,dimension%jspd)
   INTEGER,                   INTENT(IN)    :: iType

   INTEGER :: l, ilo

   DO l = 0,atoms%lmax(iType)
      CALL int_21(f,g,atoms,iType,l,thisDenCoeffsOffdiag%uu21n,thisDenCoeffsOffdiag%ud21n,&
                                    thisDenCoeffsOffdiag%du21n,thisDenCoeffsOffdiag%dd21n)
   END DO
   DO ilo = 1, atoms%nlo(iType)
      CALL int_21lo(f,g,atoms,iType,flo,ilo,thisDenCoeffsOffdiag%uulo21n,thisDenCoeffsOffdiag%ulou21n,&
                                            thisDenCoeffsOffdiag%dulo21n,thisDenCoeffsOffdiag%ulod21n,&
                                            thisDenCoeffsOffdiag%uloulop21n)
   END DO

END SUBROUTINE addRadFunScalarProducts

SUBROUTINE calcCoefficients(thisDenCoeffsOffdiag,atoms,sphhar,sym,eigVecCoeffs,we,noccbd)

   USE m_types_setup
   USE m_types_cdnval
   USE m_rhomt21     ! calculate (spin) off-diagonal MT-density coeff's
   USE m_rhonmt21    ! -"-                       non-MT-density coeff's

   IMPLICIT NONE

   CLASS(t_denCoeffsOffdiag), INTENT(INOUT) :: thisDenCoeffsOffdiag
   TYPE(t_atoms),             INTENT(IN)    :: atoms
   TYPE(t_sphhar),            INTENT(IN)    :: sphhar
   TYPE(t_sym),               INTENT(IN)    :: sym
   TYPE(t_eigVecCoeffs),      INTENT(IN)    :: eigVecCoeffs
   REAL,                      INTENT(IN)    :: we(noccbd)
   INTEGER,                   INTENT(IN)    :: noccbd

   CALL rhomt21(atoms,we,noccbd,eigVecCoeffs,thisDenCoeffsOffdiag%uu21,thisDenCoeffsOffdiag%ud21,&
                thisDenCoeffsOffdiag%du21,thisDenCoeffsOffdiag%dd21,thisDenCoeffsOffdiag%uulo21,&
                thisDenCoeffsOffdiag%dulo21,thisDenCoeffsOffdiag%ulou21,thisDenCoeffsOffdiag%ulod21,&
                thisDenCoeffsOffdiag%uloulop21)
   IF (thisDenCoeffsOffdiag%l_fmpl) THEN
      CALL rhonmt21(atoms,sphhar,we,noccbd,sym,eigVecCoeffs,thisDenCoeffsOffdiag%uunmt21,thisDenCoeffsOffdiag%udnmt21,&
                                                            thisDenCoeffsOffdiag%dunmt21,thisDenCoeffsOffdiag%ddnmt21)
   END IF

END SUBROUTINE calcCoefficients

END MODULE m_types_denCoeffsOffdiag
