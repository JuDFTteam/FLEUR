!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_hybrid
  USE m_judft
  USE m_types_fleurinput_base
  IMPLICIT NONE
  PRIVATE

  TYPE,EXTENDS(t_fleurinput_base):: t_hybrid
      LOGICAL                ::  l_hybrid = .false.
      LOGICAL                ::  l_subvxc = .false.
      LOGICAL                ::  l_calhf = .false.
      LOGICAL                ::  l_addhf = .false.
      INTEGER                ::  ewaldlambda
      INTEGER                ::  lexp = 0
      INTEGER                ::  bands1 !Only read in
      INTEGER                ::  nbasp
      INTEGER                ::  maxbasm1
      INTEGER                ::  max_indx_p_1 !new
      INTEGER                ::  maxlmindx
      INTEGER, ALLOCATABLE   ::  select1(:,:)
      INTEGER, ALLOCATABLE   ::  lcutm1(:)
      INTEGER, ALLOCATABLE   ::  lcutwf(:)
      INTEGER, ALLOCATABLE   ::  map(:,:)
      INTEGER, ALLOCATABLE   ::  tvec(:,:,:)
      INTEGER, ALLOCATABLE   ::  nbasm(:)
      !REAL, ALLOCATABLE      ::  radbasfn_mt(:,:,:,:)
      COMPLEX, ALLOCATABLE   ::  d_wgn2(:,:,:,:)
      INTEGER, ALLOCATABLE   ::  ne_eig(:)
      INTEGER, ALLOCATABLE   ::  nbands(:)
      INTEGER, ALLOCATABLE   ::  nobd(:,:)
      REAL, ALLOCATABLE      ::  div_vv(:,:,:)
    CONTAINS
      PROCEDURE :: read_xml =>read_xml_hybrid
      PROCEDURE :: mpi_bc =>mpi_bc_hybrid
   END TYPE t_hybrid
   PUBLIC t_hybrid

 CONTAINS

   SUBROUTINE mpi_bc_hybrid(this,mpi_comm,irank)
     USE m_mpi_bc_tool
     CLASS(t_hybrid),INTENT(INOUT)::this
     INTEGER,INTENT(IN):: mpi_comm
     INTEGER,INTENT(IN),OPTIONAL::irank
     INTEGER ::rank
     IF (PRESENT(irank)) THEN
        rank=irank
     ELSE
        rank=0
     END IF
     CALL mpi_bc(this%l_hybrid,rank,mpi_comm)
     CALL mpi_bc(this%l_subvxc ,rank,mpi_comm)
     CALL mpi_bc(this%l_calhf ,rank,mpi_comm)
     CALL mpi_bc(this%l_addhf ,rank,mpi_comm)
     CALL mpi_bc(this%ewaldlambda ,rank,mpi_comm)
     CALL mpi_bc(this%lexp ,rank,mpi_comm)
     CALL mpi_bc(this%bands1,rank,mpi_comm)
     CALL mpi_bc(this%nbasp,rank,mpi_comm)
     !CALL mpi_bc(this%maxlcutm1,rank,mpi_comm)
     !CALL mpi_bc(this%maxindxm1,rank,mpi_comm)
     CALL mpi_bc(this%maxbasm1,rank,mpi_comm)
     CALL mpi_bc(this%max_indx_p_1,rank,mpi_comm)
     !CALL mpi_bc(this%maxindxp1,rank,mpi_comm)
     !CALL mpi_bc(this%maxgptm,rank,mpi_comm)
     !CALL mpi_bc(this%maxgptm1,rank,mpi_comm)
     !CALL mpi_bc(this%maxindx,rank,mpi_comm)
     CALL mpi_bc(this%maxlmindx,rank,mpi_comm)
     !CALL mpi_bc(this%gptmd,rank,mpi_comm)
     !CALL mpi_bc(this%nindx,rank,mpi_comm)
     CALL mpi_bc(this%select1,rank,mpi_comm)
     CALL mpi_bc(this%lcutm1,rank,mpi_comm)
     !CALL mpi_bc(this%nindxm1,rank,mpi_comm)
     !CALL mpi_bc(this%gptm,rank,mpi_comm)
     !CALL mpi_bc(this%ngptm1,rank,mpi_comm)
     !CALL mpi_bc(this%pgptm1,rank,mpi_comm)
     !CALL mpi_bc(this%ngptm,rank,mpi_comm)
     !CALL mpi_bc(this%pgptm,rank,mpi_comm)
     CALL mpi_bc(this%lcutwf,rank,mpi_comm)
     CALL mpi_bc(this%map,rank,mpi_comm)
     CALL mpi_bc(this%tvec,rank,mpi_comm)
     CALL mpi_bc(this%nbasm,rank,mpi_comm)
     !CALL mpi_bc(this%gcutm1,rank,mpi_comm)
     !CALL mpi_bc(this%tolerance1 ,rank,mpi_comm)
     !CALL mpi_bc(this%basm1,rank,mpi_comm)
     CALL mpi_bc(this%d_wgn2,rank,mpi_comm)
     CALL mpi_bc(this%ne_eig,rank,mpi_comm)
     CALL mpi_bc(this%nbands,rank,mpi_comm)
     CALL mpi_bc(this%nobd,rank,mpi_comm)
     CALL mpi_bc(this%div_vv,rank,mpi_comm)


   END SUBROUTINE mpi_bc_hybrid
   SUBROUTINE read_xml_hybrid(this,xml)
     USE m_types_xml
     CLASS(t_hybrid),INTENT(INout):: this
     TYPE(t_xml),INTENT(in)     :: xml


     INTEGER::numberNodes,ntype,itype
     CHARACTER(len=100)::xPathA

     ntype=xml%GetNumberOfNodes('/fleurInput/atomGroups/atomGroup')
     ALLOCATE(this%lcutm1(ntype),this%lcutwf(ntype),this%select1(4,ntype))
     numberNodes = xml%GetNumberOfNodes('/fleurInput/calculationSetup/prodBasis')
     IF (numberNodes==1) THEN
         !this%gcutm1=evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/prodBasis/@gcutm'))
         !this%tolerance1=evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/prodBasis/@tolerance'))
         this%ewaldlambda=evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/prodBasis/@ewaldlambda'))
         this%lexp=evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/prodBasis/@lexp'))
         this%bands1=evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/prodBasis/@bands'))
      ENDIF

      DO itype=1,ntype
         xpatha=xml%SpeciesPath(itype)//'/prodBasis'
         numberNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA)))
         IF (numberNodes==1) THEN
            this%lcutm1(iType) =evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@lcutm'))
            this%lcutwf(iType) =evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@lcutwf'))
            xPathA=xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@select')
            this%select1(1,iType) = NINT(evaluateFirst(xPathA))
            this%select1(2,iType) = NINT(evaluateFirst(xPathA))
            this%select1(3,iType) = NINT(evaluateFirst(xPathA))
            this%select1(4,iType) = NINT(evaluateFirst(xPathA))
         ELSE
            this%lcutm1(iType) =-1
         ENDIF
      END DO
    END SUBROUTINE read_xml_hybrid
 END MODULE m_types_hybrid
