!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_noco
  USE m_judft
  USE m_types_fleurinput_base
  IMPLICIT NONE
  PRIVATE
  TYPE,EXTENDS(t_fleurinput_base):: t_noco
    LOGICAL:: l_ss= .FALSE.
    LOGICAL:: l_soc= .FALSE.
    LOGICAL:: l_noco = .FALSE.
    LOGICAL:: l_mperp = .FALSE.
    LOGICAL:: l_constr = .FALSE.
    LOGICAL:: l_mtNocoPot = .FALSE.
    LOGICAL:: l_alignMT = .FALSE.
    LOGICAL:: l_sourceFree = .FALSE.
    REAL   :: mix_b=0.0
    LOGICAL:: l_spav= .FALSE.
    REAL   :: theta_inp=0.0
    REAL   :: phi_inp=0.0
    REAL   :: qss_inp(3)=[0.,0.,0.]

    LOGICAL, ALLOCATABLE :: l_relax(:)
    REAL, ALLOCATABLE    :: alph_inp(:)
    REAL, ALLOCATABLE    :: beta_inp(:)
    REAL, ALLOCATABLE    :: socscale(:)

  CONTAINS
    PROCEDURE :: read_xml=>read_xml_noco
    PROCEDURE :: mpi_bc =>mpi_bc_noco
   END TYPE t_noco

   PUBLIC t_noco

 CONTAINS
   SUBROUTINE mpi_bc_noco(this,mpi_comm,irank)
     USE m_mpi_bc_tool
     CLASS(t_noco),INTENT(INOUT)::this
     INTEGER,INTENT(IN):: mpi_comm
     INTEGER,INTENT(IN),OPTIONAL::irank
     INTEGER ::rank
     IF (PRESENT(irank)) THEN
        rank=irank
     ELSE
        rank=0
     END IF

     CALL mpi_bc(this%l_ss,rank,mpi_comm)
     CALL mpi_bc(this%l_soc,rank,mpi_comm)
     CALL mpi_bc(this%l_noco ,rank,mpi_comm)
     CALL mpi_bc(this%l_mperp ,rank,mpi_comm)
     CALL mpi_bc(this%l_constr ,rank,mpi_comm)
     CALL mpi_bc(this%l_mtNocoPot ,rank,mpi_comm)
     CALL mpi_bc(this%l_sourceFree ,rank,mpi_comm)
     CALL mpi_bc(rank,mpi_comm,this%qss_inp)
     CALL mpi_bc(this%mix_b,rank,mpi_comm)
     CALL mpi_bc(this%l_spav,rank,mpi_comm)
     CALL mpi_bc(this%theta_inp,rank,mpi_comm)
     CALL mpi_bc(this%phi_inp,rank,mpi_comm)

     CALL mpi_bc(this%l_relax,rank,mpi_comm)
     CALL mpi_bc(this%alph_inp,rank,mpi_comm)
     CALL mpi_bc(this%beta_inp,rank,mpi_comm)
     !CALL mpi_bc(this%b_con,rank,mpi_comm)
     CALL mpi_bc(this%socscale,rank,mpi_comm)


   END SUBROUTINE mpi_bc_noco

   SUBROUTINE read_xml_noco(this,xml)
     USE m_types_xml
     CLASS(t_noco),INTENT(inout):: this
     TYPE(t_xml),INTENT(INOUT) ::xml

     INTEGER:: numberNodes,ntype,itype
     CHARACTER(len=100)::xpathA,xpathB,valueString

      this%l_noco = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/magnetism/@l_noco'))

      ! Read in optional SOC parameters if present
      xPathA = '/fleurInput/calculationSetup/soc'
      numberNodes = xml%GetNumberOfNodes(xPathA)
      IF (numberNodes.EQ.1) THEN
         this%theta_inp=evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@theta'))
         this%phi_inp=evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@phi'))
         this%l_soc = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_soc'))
         this%l_spav = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@spav'))
      END IF

      ! Read in optional noco parameters if present
      xPathA = '/fleurInput/calculationSetup/nocoParams'
      numberNodes = xml%GetNumberOfNodes(xPathA)
      IF ((this%l_noco).AND.(numberNodes.EQ.0)) THEN
         CALL juDFT_error('Error: l_noco is true but no noco parameters set in xml input file!')
      END IF

      IF (numberNodes.EQ.1) THEN
         this%l_ss = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_ss'))
         this%l_mperp = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_mperp'))
         this%l_constr = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_constr'))
         this%l_mtNocoPot = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_mtNocoPot'))
         this%l_sourceFree = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_sourceFree'))
         this%mix_b = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@mix_b'))
         valueString = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/qss')))
         READ(valueString,*) this%qss_inp
      END IF

      ntype=xml%GetNumberOfNodes('/fleurInput/atomGroups/atomGroup')
      ALLOCATE(this%l_relax(ntype))
      this%l_relax=.false.
      ALLOCATE(this%alph_inp(ntype),this%beta_inp(ntype))
      this%alph_inp=0.0;this%beta_inp=0.0
      ALLOCATE(this%socscale(ntype))
      this%socscale=0.0

      DO itype=1,ntype
        this%socscale(Itype)=1.0
         IF (xml%GetNumberOfNodes(TRIM(ADJUSTL(xml%speciesPath(itype)))//'/special/@socscale')>0) &
              this%socscale(Itype)=evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xml%speciesPath(itype)))//'/special/@socscale'))
         !Read in atom group specific noco parameters
         xPathB = TRIM(ADJUSTL(xml%groupPath(itype)))//'/nocoParams'
         numberNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathB)))
         IF (numberNodes.GE.1) THEN
            this%l_relax(iType) = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@l_relax'))
            this%alph_inp(iType) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@alpha'))
            this%beta_inp(iType) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@beta'))
            !this%b_con(1,iType) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@b_cons_x'))
            !this%b_con(2,iType) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@b_cons_y'))
         END IF
      ENDDO
    END SUBROUTINE read_xml_noco





 END MODULE m_types_noco
