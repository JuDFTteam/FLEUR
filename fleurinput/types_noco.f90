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
    LOGICAL:: l_soc= .FALSE.
    LOGICAL:: l_spav= .FALSE.

    LOGICAL:: l_noco = .FALSE.
    LOGICAL:: l_ss= .FALSE.
    LOGICAL:: l_mperp = .FALSE.
    INTEGER:: mag_mixing_scheme=0

    LOGICAL:: l_sourceFree = .FALSE.
    LOGICAL:: l_scaleMag = .FALSE.
    REAL   :: mag_scale=1.0
    REAL   :: mix_b=1.0


    REAL   :: theta_inp=0.0
    REAL   :: phi_inp=0.0
    REAL   :: qss_inp(3)=[0.0,0.0,0.0]

    REAL, ALLOCATABLE   :: mix_RelaxWeightOffD(:)
    LOGICAL, ALLOCATABLE :: l_constrained(:)
    LOGICAL, ALLOCATABLE :: l_unrestrictMT(:)
    LOGICAL, ALLOCATABLE :: l_spinoffd_ldau(:)
    LOGICAL,allocatable  :: l_alignMT(:)
    REAL, ALLOCATABLE    :: alph_inp(:)
    REAL, ALLOCATABLE    :: beta_inp(:)
    REAL, ALLOCATABLE    :: socscale(:)

  CONTAINS
    PROCEDURE :: read_xml=>read_xml_noco
    PROCEDURE :: mpi_bc =>mpi_bc_noco
    PROCEDURE :: init => init_noco
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
     CALL mpi_bc(this%l_alignMT ,rank,mpi_comm)
     CALL mpi_bc(this%l_sourceFree ,rank,mpi_comm)
     CALL mpi_bc(this%l_scaleMag ,rank,mpi_comm)
     CALL mpi_bc(this%mag_scale ,rank,mpi_comm)
     CALL mpi_bc(rank,mpi_comm,this%qss_inp)
     call mpi_bc(this%mag_mixing_scheme,rank,mpi_comm)
     CALL mpi_bc(this%mix_RelaxWeightOffD,rank,mpi_comm)
     CALL mpi_bc(this%l_spav,rank,mpi_comm)
     CALL mpi_bc(this%theta_inp,rank,mpi_comm)
     CALL mpi_bc(this%phi_inp,rank,mpi_comm)
     CALL mpi_bc(this%alph_inp,rank,mpi_comm)
     CALL mpi_bc(this%beta_inp,rank,mpi_comm)
     CALL mpi_bc(this%l_constrained,rank,mpi_comm)
     CALL mpi_bc(this%l_unrestrictMT,rank,mpi_comm)
     CALL mpi_bc(this%mix_b,rank,mpi_comm)
     CALL mpi_bc(this%socscale,rank,mpi_comm)
     CALL mpi_bc(this%l_spinoffd_ldau,rank,mpi_comm)

   END SUBROUTINE mpi_bc_noco

   SUBROUTINE read_xml_noco(this,xml)
     USE m_types_xml
     CLASS(t_noco),INTENT(inout):: this
     TYPE(t_xml),INTENT(INOUT) ::xml

     INTEGER:: numberNodes,ntype,itype
     CHARACTER(len=100)::xpathA,xpathB,valueString

     if (xml%versionNumber<33) THEN
       call read_xml_noco_old(this,xml)
       RETURN
     endif

      this%l_noco = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/magnetism/@l_noco'))
      this%l_ss = .FALSE.
      this%l_ss = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/magnetism/@l_ss'))

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
      !IF (xml%versionNumber > 31) THEN
      !   xPathA = '/fleurInput/calculationSetup/magnetism/nocoParams'
      !ELSE
      !   xPathA = '/fleurInput/calculationSetup/nocoParams'
      !END IF
      !numberNodes = xml%GetNumberOfNodes(xPathA)
      !IF ((this%l_noco).AND.(numberNodes.EQ.0)) THEN
      !   CALL juDFT_error('Error: l_noco is true but no noco parameters set in xml input file!')
      !END IF


      this%qss_inp=0.0
      ntype=xml%GetNumberOfNodes('/fleurInput/atomGroups/atomGroup')
      ALLOCATE(this%l_constrained(ntype),this%mix_RelaxWeightOffD(ntype))
      ALLOCATE(this%l_unrestrictMT(ntype),this%l_alignMT(ntype))
      ALLOCATE(this%l_spinoffd_ldau(ntype))
      this%l_unrestrictMT=.false.
      this%l_constrained=.false.
      this%mix_RelaxWeightOffD=1.0
      this%l_alignMT=.false.
      this%l_spinoffd_ldau=.false.
      ALLOCATE(this%alph_inp(ntype),this%beta_inp(ntype))
      this%alph_inp=0.0;this%beta_inp=0.0
      ALLOCATE(this%socscale(ntype))
      this%socscale=0.0


      xPathA = '/fleurInput/calculationSetup/magnetism/qss'
      numberNodes = xml%GetNumberOfNodes(xPathA)
      IF (numberNodes.EQ.1) THEN
        valueString = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA)))))
        this%qss_inp(1) = evaluatefirst(valueString)
        this%qss_inp(2) = evaluatefirst(valueString)
        this%qss_inp(3) = evaluatefirst(valueString)
      END IF

      ! Read in noco MT parameters
      xPathA = '/fleurInput/calculationSetup/magnetism/mtNocoParams'
      numberNodes = xml%GetNumberOfNodes(xPathA)
      IF (numberNodes.EQ.1) THEN
         this%l_alignMT   = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_relaxSQA'))
         this%l_constrained = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_constrained'))
         this%l_mperp= evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_mperp'))
         this%mix_b = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@mix_constr'))
         this%mag_mixing_scheme= evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@mag_mixing_scheme'))
         this%mix_RelaxWeightOffD = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@mix_RelaxWeightOffD'))
         this%l_unrestrictMT=evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_mtNocoPot'))
      END IF

      ! Read in optional source free magnetism parameters if present
      xPathA = '/fleurInput/calculationSetup/magnetism/sourceFreeMag'
      numberNodes = xml%GetNumberOfNodes(xPathA)
      IF(numberNodes.EQ.1) THEN
         this%l_sourceFree = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_sourceFree'))
         this%l_scaleMag = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_scaleMag'))
         this%mag_scale = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@mag_scale'))
      END IF



      DO itype=1,ntype
        this%socscale(Itype)=1.0
         IF (xml%GetNumberOfNodes(TRIM(ADJUSTL(xml%speciesPath(itype)))//'/special/@socscale')>0) &
              this%socscale(Itype)=evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xml%speciesPath(itype)))//'/special/@socscale'))
         !Read in atom group specific noco parameters
         xPathB = TRIM(ADJUSTL(xml%groupPath(itype)))//'/nocoParams'
         numberNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathB)))
         IF (numberNodes.GE.1) THEN
           if (xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathB))//'/@l_constrained')>0) this%l_constrained(iType) = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@l_constrained'))
           if (xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathB))//'/@l_mtNocoPot')>0) this%l_unrestrictMT(iType) = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@l_mtNocoPot'))
           if (xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathB))//'/@l_relaxSQA')>0) this%l_alignMT(iType) = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@l_relaxSQA'))
           this%alph_inp(iType) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@alpha'))
           this%beta_inp(iType) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@beta'))
         END IF
      ENDDO
      if (any(this%l_unrestrictMT)) this%l_mperp=.true.

    END SUBROUTINE read_xml_noco

   SUBROUTINE read_xml_noco_old(this,xml)
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

      ntype=xml%GetNumberOfNodes('/fleurInput/atomGroups/atomGroup')
      ALLOCATE(this%l_alignMT(ntype),this%l_constrained(ntype),this%mix_RelaxWeightOffD(ntype))
      ALLOCATE(this%l_unrestrictMT(ntype),this%l_spinoffd_ldau(ntype))
      this%l_alignMT=.false.
      this%l_unrestrictMT=.false.
      this%mix_RelaxWeightOffD=1.0
      this%l_constrained=.false.
      this%l_spinoffd_ldau=.false.
      ALLOCATE(this%alph_inp(ntype),this%beta_inp(ntype))
      this%alph_inp=0.0;this%beta_inp=0.0
      ALLOCATE(this%socscale(ntype))
      this%socscale=1.0

      IF (numberNodes.EQ.1) THEN
         this%l_ss = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_ss'))
         this%l_mperp = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_mperp'))
         this%l_constrained = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_constr'))
         this%l_unrestrictMT = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_mtNocoPot'))
         this%l_sourceFree = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_sourceFree'))
         this%l_scaleMag = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_scaleMag'))
         this%mag_scale = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@mag_scale'))
         this%mix_b = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@mix_b'))
         this%l_alignMT = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@l_RelaxMT'))
         this%mix_RelaxWeightOffD = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@mix_RelaxWeightOffD'))
         valueString = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/qss')))
         READ(valueString,*) this%qss_inp
      END IF



      DO itype=1,ntype
        this%socscale(Itype)=1.0
         IF (xml%GetNumberOfNodes(TRIM(ADJUSTL(xml%speciesPath(itype)))//'/special/@socscale')>0) &
              this%socscale(Itype)=evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xml%speciesPath(itype)))//'/special/@socscale'))
         !Read in atom group specific noco parameters
         xPathB = TRIM(ADJUSTL(xml%groupPath(itype)))//'/nocoParams'
         numberNodes = xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathB)))
         IF (numberNodes.GE.1) THEN
            this%l_alignMt(iType) = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@l_relax'))
            this%alph_inp(iType) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@alpha'))
            this%beta_inp(iType) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@beta'))
            !this%b_con(1,iType) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@b_cons_x'))
            !this%b_con(2,iType) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@b_cons_y'))
         END IF
      ENDDO


    END SUBROUTINE read_xml_noco_old

  SUBROUTINE init_noco(noco,atoms,l_spinoffd_ldau)
     USE m_types_atoms
     CLASS(t_noco),INTENT(inout):: noco
     TYPE(t_atoms),INTENT(in)   :: atoms
     LOGICAL, INTENT(in)        :: l_spinoffd_ldau
     INTEGER :: i_u, atomType

     IF (noco%l_noco .AND. (ABS(noco%theta_inp)>1E-5.OR.ABS(noco%phi_inp)>1E-5)) &
           CALL judft_warn("You specified a theta/phi angle for SOC in a noco-calculation. These angles are used only in non-Noco SOC calculations.")

     IF (.NOT.noco%l_noco) THEN
        IF (ANY(noco%alph_inp(1:atoms%ntype).NE.0.0).OR.ANY(noco%beta_inp(1:atoms%ntype).NE.0.0)) THEN
           CALL judft_warn("You specified noco-calculation alpha/beta angles unequal to 0.0 but don't perform a noco-calculation.")
        END IF
     END IF

 
     IF(l_spinoffd_ldau) THEN
        DO i_u = 1, atoms%n_u+atoms%n_hia
           atomType = atoms%lda_u(i_u)%atomType
           noco%l_spinoffd_ldau(atomType) = .TRUE.
        ENDDO
     ENDIF

  END SUBROUTINE init_noco

 END MODULE m_types_noco
