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
    LOGICAL:: l_sourceFree = .FALSE.
    REAL   :: qss(3)=[0.,0.,0.]
    REAL   :: mix_b=0.0
    LOGICAL:: l_spav= .FALSE.
    REAL   :: theta=0.0
    REAL   :: phi=0.0

    LOGICAL, ALLOCATABLE :: l_relax(:)
    REAL, ALLOCATABLE    :: alphInit(:)
    REAL, ALLOCATABLE    :: alph(:)
    REAL, ALLOCATABLE    :: beta(:)
    REAL, ALLOCATABLE    :: b_con(:, :)
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
     CALL mpi_bc(this%l_constr ,rank,mpi_comm)
     CALL mpi_bc(this%l_mtNocoPot ,rank,mpi_comm)
     CALL mpi_bc(this%l_sourceFree ,rank,mpi_comm)
     CALL mpi_bc(rank,mpi_comm,this%qss)
     CALL mpi_bc(this%mix_b,rank,mpi_comm)
     CALL mpi_bc(this%l_spav,rank,mpi_comm)
     CALL mpi_bc(this%theta,rank,mpi_comm)
     CALL mpi_bc(this%phi,rank,mpi_comm)

     CALL mpi_bc(this%l_relax,rank,mpi_comm)
     CALL mpi_bc(this%alphInit,rank,mpi_comm)
     CALL mpi_bc(this%alph,rank,mpi_comm)
     CALL mpi_bc(this%beta,rank,mpi_comm)
     CALL mpi_bc(this%b_con,rank,mpi_comm)
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
         this%theta=evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@theta'))
         this%phi=evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@phi'))
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
         READ(valueString,*) this%qss(1), this%qss(2), this%qss(3)
      END IF

      ntype=xml%GetNumberOfNodes('/fleurInput/atomGroups/atomGroup')
      ALLOCATE(this%l_relax(ntype),this%b_con(2,ntype))
      this%l_relax=.false.; this%b_con=0.0
      ALLOCATE(this%alphInit(ntype),this%alph(ntype),this%beta(ntype))
      this%alphInit=0.0;this%alph=0.0;this%beta=0.0
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
            this%alphInit(iType) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@alpha'))
            this%alph(iType) = this%alphInit(iType)
            this%beta(iType) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@beta'))
            this%b_con(1,iType) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@b_cons_x'))
            this%b_con(2,iType) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@b_cons_y'))
         END IF
      ENDDO
    END SUBROUTINE read_xml_noco

    subroutine init_noco(noco,atoms)
      use m_types_atoms
      use m_constants
      CLASS(t_noco),INTENT(inout):: noco
      TYPE(t_atoms),INTENT(IN)::atoms


      integer :: na,itype

      ! Check noco stuff and calculate missing noco parameters
      IF (noco%l_noco) THEN
         IF (noco%l_ss) THEN
            !--->    the angle beta is relative to the spiral in a spin-spiral
            !--->    calculation, i.e. if beta = 0 for all atoms in the unit cell
            !--->    that means that the moments are "in line" with the spin-spiral
            !--->    (beta = qss * taual). note: this means that only atoms within
            !--->    a plane perpendicular to qss can be equivalent!
            na = 1
            DO iType = 1,atoms%ntype
               noco%alph(iType) = noco%alphInit(iType) + tpi_const*dot_product(noco%qss,atoms%taual(:,na))
               na = na + atoms%neq(iType)
            END DO
         END IF
      ELSE
         IF (noco%l_ss) THEN
            CALL judft_warn("l_noco=F and l_ss=T is meaningless. Setting l_ss to F.")
            noco%l_ss = .FALSE.
         END IF
      END IF



    end subroutine init_noco



 END MODULE m_types_noco
