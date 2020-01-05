!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_hybinp
   USE m_judft
   USE m_types_fleurinput_base
   IMPLICIT NONE
   PRIVATE

   TYPE, EXTENDS(t_fleurinput_base):: t_hybinp
      LOGICAL                ::  l_hybrid = .false.
      LOGICAL                ::  l_subvxc = .false.
      LOGICAL                ::  l_calhf = .false.
      LOGICAL                ::  l_addhf = .false.
      INTEGER                ::  ewaldlambda = -1
      INTEGER                ::  lexp = -1
      INTEGER                ::  bands1 = -1 !Only read in
      INTEGER, ALLOCATABLE   ::  select1(:, :)
      INTEGER, ALLOCATABLE   ::  lcutm1(:)
      INTEGER, ALLOCATABLE   ::  lcutwf(:)
      INTEGER, ALLOCATABLE   ::  map(:, :)
      INTEGER, ALLOCATABLE   ::  tvec(:, :, :)
      !REAL, ALLOCATABLE      ::  radbasfn_mt(:,:,:,:)
      COMPLEX, ALLOCATABLE   ::  d_wgn2(:, :, :, :)

   CONTAINS
      PROCEDURE :: read_xml => read_xml_hybinp
      PROCEDURE :: mpi_bc => mpi_bc_hybinp
   END TYPE t_hybinp
   PUBLIC t_hybinp

CONTAINS

   SUBROUTINE mpi_bc_hybinp(this, mpi_comm, irank)
      USE m_mpi_bc_tool
      CLASS(t_hybinp), INTENT(INOUT)::this
      INTEGER, INTENT(IN):: mpi_comm
      INTEGER, INTENT(IN), OPTIONAL::irank
      INTEGER ::rank
      IF (PRESENT(irank)) THEN
         rank = irank
      ELSE
         rank = 0
      END IF
      CALL mpi_bc(this%l_hybrid, rank, mpi_comm)
      CALL mpi_bc(this%l_subvxc, rank, mpi_comm)
      CALL mpi_bc(this%l_calhf, rank, mpi_comm)
      CALL mpi_bc(this%l_addhf, rank, mpi_comm)
      CALL mpi_bc(this%ewaldlambda, rank, mpi_comm)
      CALL mpi_bc(this%lexp, rank, mpi_comm)
      CALL mpi_bc(this%bands1, rank, mpi_comm)
      CALL mpi_bc(this%select1, rank, mpi_comm)
      CALL mpi_bc(this%lcutm1, rank, mpi_comm)
      CALL mpi_bc(this%lcutwf, rank, mpi_comm)
      CALL mpi_bc(this%map, rank, mpi_comm)
      CALL mpi_bc(this%tvec, rank, mpi_comm)
      CALL mpi_bc(this%d_wgn2, rank, mpi_comm)
   END SUBROUTINE mpi_bc_hybinp
   SUBROUTINE read_xml_hybinp(this, xml)
      USE m_types_xml
      CLASS(t_hybinp), INTENT(INout):: this
      TYPE(t_xml), INTENT(in)     :: xml

      INTEGER::numberNodes, ntype, itype
      CHARACTER(len=100)  :: xPathA
      CHARACTER(len=4), allocatable  :: xc_name

      ntype = xml%GetNumberOfNodes('/fleurInput/atomGroups/atomGroup')
      ALLOCATE (this%lcutm1(ntype), this%lcutwf(ntype), this%select1(4, ntype), source=0)
      numberNodes = xml%GetNumberOfNodes('/fleurInput/calculationSetup/prodBasis')
      IF (numberNodes == 1) THEN
         ! this%g_cutoff=evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/prodBasis/@gcutm'))
         ! this%linear_dep_tol=evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/prodBasis/@tolerance'))
         this%ewaldlambda = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/prodBasis/@ewaldlambda'))
         this%lexp = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/prodBasis/@lexp'))
         this%bands1 = evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/prodBasis/@bands'))
      ENDIF

      DO itype = 1, ntype
         xpatha = xml%SpeciesPath(itype)//'/prodBasis'
         numberNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA)))
         IF (numberNodes == 1) THEN
            this%lcutm1(iType) = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@lcutm'))
            this%lcutwf(iType) = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@lcutwf'))
            xPathA = xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@select')
            this%select1(1, iType) = NINT(evaluateFirst(xPathA))
            this%select1(2, iType) = NINT(evaluateFirst(xPathA))
            this%select1(3, iType) = NINT(evaluateFirst(xPathA))
            this%select1(4, iType) = NINT(evaluateFirst(xPathA))
         ELSE
            this%lcutm1(iType) = -1
         ENDIF
      END DO

      xc_name = trim(xml%GetAttributeValue('/fleurInput/xcFunctional/@name'))
      if(trim(xc_name) == "pbe0") then
         this%l_hybrid = .True.
      else
         this%l_hybrid = .False.
      endif
   END SUBROUTINE read_xml_hybinp
END MODULE m_types_hybinp
