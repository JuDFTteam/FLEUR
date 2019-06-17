!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_oneD
    ! types for 1D calculations
   TYPE od_dim
      LOGICAL :: d1=.false.
      INTEGER :: mb=0, M=0, k3=0, m_cyl=0
      INTEGER :: chi=0, rot=0
      LOGICAL :: invs=.false., zrfs=.false.
      INTEGER :: n2d=0, nq2=0, nn2d=0
      INTEGER :: kimax2
      INTEGER :: nop=0, nat=0
   END TYPE od_dim

   TYPE od_inp
      LOGICAL :: d1
      INTEGER :: mb, M, k3, m_cyl
      INTEGER :: chi, rot
      LOGICAL :: invs, zrfs
      INTEGER :: n2d, nq2, nn2d
      INTEGER :: kimax2
      INTEGER, POINTER :: ig(:, :)  !(-k3:k3,-M:M)
      INTEGER, POINTER :: kv(:, :)        !(2,n2d)
      INTEGER, POINTER :: nst2(:)        !(n2d)
   END TYPE od_inp

   TYPE od_sym
      INTEGER :: nop, nat
      INTEGER, POINTER :: ngopr(:)     !(nat)
      REAL, POINTER :: mrot(:, :, :)  !(3,3,nop)
      REAL, POINTER :: tau(:, :)     !(3,nop)
      INTEGER, POINTER :: invtab(:)    !(nop)
      INTEGER, POINTER :: multab(:, :)  !(nop,nop)
   END TYPE od_sym

   TYPE od_lda
      INTEGER :: nn2d
      INTEGER, POINTER :: igf(:, :)  !(0:nn2d-1,2)
      REAL, POINTER :: pgf(:)    !(0:nn2d-1)
   END TYPE od_lda

   TYPE od_gga
      INTEGER          :: nn2d
      REAL, POINTER    :: pgfx(:)  ! (0:nn2d-1)
      REAL, POINTER    :: pgfy(:)
      REAL, POINTER    :: pgfxx(:)
      REAL, POINTER    :: pgfyy(:)
      REAL, POINTER    :: pgfxy(:)
   END TYPE od_gga
 TYPE t_oneD
      TYPE(od_dim) :: odd
      TYPE(od_inp) :: odi
      TYPE(od_sym) :: ods
      TYPE(od_lda) :: odl
      TYPE(od_gga) :: odg
      INTEGER, POINTER :: ig1(:, :)
      INTEGER, POINTER :: kv1(:, :)
      INTEGER, POINTER :: nstr1(:)
      INTEGER, POINTER :: ngopr1(:)
      REAL, POINTER :: mrot1(:, :, :)
      REAL, POINTER :: tau1(:, :)
      INTEGER, POINTER :: invtab1(:)
      INTEGER, POINTER :: multab1(:, :)
      INTEGER, POINTER :: igfft1(:, :)
      REAL, POINTER :: pgfft1(:)
      REAL, POINTER :: pgft1x(:)
      REAL, POINTER :: pgft1y(:)
      REAL, POINTER :: pgft1xx(:)
      REAL, POINTER :: pgft1yy(:)
      REAL, POINTER :: pgft1xy(:)
    contains
      procedure :: read_xml
   END TYPE t_oneD
 contains
   subroutine read_xml(oneD,xml)
     use m_types_xml
     class(t_oned),intent(out)::oneD
     type(t_xml),intent(in)   ::xml

     
      ! Read in optional 1D parameters if present
     character(len=100):: xpathA
     integer :: numberNodes
     
      xPathA = '/fleurInput/calculationSetup/oneDParams'
      numberNodes = xml%GetNumberOfNodes(xPathA)

      oneD%odd%d1 = .FALSE.

      IF (numberNodes.EQ.1) THEN
         oneD%odd%d1 = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@d1'))
         oneD%odd%M = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@MM'))
         oneD%odd%mb = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@vM'))
         oneD%odd%m_cyl = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@m_cyl'))
         oneD%odd%chi = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@chi'))
         oneD%odd%rot = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@rot'))
         oneD%odd%invs = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@invs1'))
         oneD%odd%zrfs = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@zrfs1'))
      END IF

    end subroutine read_xml
   
 END MODULE m_types_oneD
