!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_hybrid
  use m_judft
  IMPLICIT NONE
  PRIVATE

     TYPE t_hybrid
      LOGICAL               ::  l_hybrid = .false.
      LOGICAL               ::  l_subvxc = .false.
      LOGICAL               ::  l_calhf = .false.
      LOGICAL               ::  l_addhf = .false.
      INTEGER               ::  ewaldlambda =3
      INTEGER               ::  lexp =16
      INTEGER               ::  bands1 !Only read in
      INTEGER               ::  nbasp
      INTEGER               ::  maxlcutm1
      INTEGER               ::  maxindxm1
      INTEGER               ::  maxbasm1
      INTEGER               ::  maxindxp1
      INTEGER               ::  maxgptm
      INTEGER               ::  maxgptm1
      INTEGER               ::  maxindx
      INTEGER               ::  maxlmindx
      INTEGER               ::  gptmd
      INTEGER, ALLOCATABLE   ::  nindx(:, :)
      INTEGER, ALLOCATABLE   ::  select1(:, :)
      INTEGER, ALLOCATABLE   ::  lcutm1(:)
      INTEGER, ALLOCATABLE   ::  nindxm1(:, :)
      INTEGER, ALLOCATABLE   ::  gptm(:, :)
      INTEGER, ALLOCATABLE   ::  ngptm1(:)
      INTEGER, ALLOCATABLE   ::  pgptm1(:, :)
      INTEGER, ALLOCATABLE   ::  ngptm(:)
      INTEGER, ALLOCATABLE   ::  pgptm(:, :)
      INTEGER, ALLOCATABLE   ::  lcutwf(:)
      INTEGER, ALLOCATABLE   ::  map(:, :)
      INTEGER, ALLOCATABLE   ::  tvec(:, :, :)
      INTEGER, ALLOCATABLE ::  nbasm(:)
      REAL                  ::  gcutm1
      REAL                  ::  tolerance1 = 1e-4 !only read in
      REAL, ALLOCATABLE   ::  basm1(:, :, :, :)
      COMPLEX, ALLOCATABLE   ::  d_wgn2(:, :, :, :)
      INTEGER, ALLOCATABLE   ::  ne_eig(:), nbands(:), nobd(:)                   !alloc in eigen_HF_init
      REAL, ALLOCATABLE   ::  div_vv(:, :, :)
    CONTAINS
      PROCEDURE :: read_xml =>read_xml_hybrid

   END TYPE t_hybrid
   PUBLIC t_hybrid

 CONTAINS
   SUBROUTINE read_xml_hybrid(hybrid,xml)
     USE m_types_xml
     CLASS(t_hybrid),INTENT(out):: hybrid
     TYPE(t_xml),INTENT(in)     :: xml
     
     
     INTEGER::numberNodes,ntype,itype
     CHARACTER(len=100)::xPathA

     ntype=xml%GetNumberOfNodes('/fleurInput/atomGroups/atomGroup')
     ALLOCATE(hybrid%lcutm1(ntype),hybrid%lcutwf(ntype),hybrid%select1(4,ntype))
     numberNodes = xml%GetNumberOfNodes('/fleurInput/calculationSetup/prodBasis')
     IF (numberNodes==1) THEN
         hybrid%gcutm1=evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/prodBasis/@gcutm'))
         hybrid%tolerance1=evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/prodBasis/@tolerance'))
         hybrid%ewaldlambda=evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/prodBasis/@ewaldlambda'))
         hybrid%lexp=evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/prodBasis/@lexp'))
         hybrid%bands1=evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/prodBasis/@bands'))
      ENDIF

      DO itype=1,ntype
         xpatha=xml%SpeciesPath(itype)//'/prodBasis'
         numberNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA)))
         IF (numberNodes==1) THEN
            hybrid%lcutm1(iType) =evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@lcutm'))
            hybrid%lcutwf(iType) =evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@lcutwf'))
            xPathA=xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@select')
            hybrid%select1(1,iType) = NINT(evaluateFirst(xPathA))
            hybrid%select1(2,iType) = NINT(evaluateFirst(xPathA))
            hybrid%select1(3,iType) = NINT(evaluateFirst(xPathA))
            hybrid%select1(4,iType) = NINT(evaluateFirst(xPathA))
         ELSE
            hybrid%lcutm1(iType) =-1
         ENDIF
      END DO
    END SUBROUTINE read_xml_hybrid
 END MODULE m_types_hybrid
