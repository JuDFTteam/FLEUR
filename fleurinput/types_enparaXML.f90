!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_enparaXML
  USE m_judft
  USE m_types_fleurinput_base
  IMPLICIT NONE
  PRIVATE
  TYPE,EXTENDS(t_fleurinput_base):: t_enparaXML
     INTEGER,ALLOCATABLE  :: qn_el(:,:,:)    !if these are .ne.0 they are understood as
     INTEGER,ALLOCATABLE  :: qn_ello(:,:,:)  !quantum numbers
     REAL                 :: evac0(2,2)
   CONTAINS
     PROCEDURE :: init
     procedure :: set_quantum_numbers
     PROCEDURE :: read_xml=>read_xml_enpara
     procedure :: mpi_bc => mpi_bc_enpara
  END TYPE t_enparaXML



  PUBLIC:: t_enparaXML

CONTAINS

  SUBROUTINE mpi_bc_enpara(this,mpi_comm,irank)
    USE m_mpi_bc_tool
    CLASS(t_enparaXML),INTENT(INOUT)::this
    INTEGER,INTENT(IN):: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL::irank
    INTEGER ::rank
    IF (PRESENT(irank)) THEN
       rank=irank
    ELSE
       rank=0
    END IF

    CALL mpi_bc(this%qn_el,rank,mpi_comm)
    CALL mpi_bc(this%qn_ello,rank,mpi_comm)
    CALL mpi_bc(rank,mpi_comm,this%evac0)

  END SUBROUTINE mpi_bc_enpara

  SUBROUTINE read_xml_enpara(this,xml)
    use m_types_xml
    CLASS(t_enparaXML),INTENT(INOUT):: this
    TYPE(t_xml),INTENT(INOUT)   ::xml

    LOGICAL :: l_enpara,film
    INTEGER :: jspins,ntype,n,iLO,i, nLOLocal, lNumCount, nNumCount
    CHARACTER(len=200) :: xPath, xPathB
    CHARACTER(LEN=200) :: lString, nString
    INTEGER, ALLOCATABLE :: nlo(:)
    INTEGER, ALLOCATABLE :: lNumbers(:), nNumbers(:)

    ntype=xml%get_ntype()
    nlo=xml%get_nlo()
    jspins=evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/magnetism/@jspins'))
    film=xml%GetNumberOfNodes('/fleurInput/cell/filmLattice')==1

    CALL this%init(ntype,MAXVAL(nlo),jspins)

    DO n=1,ntype
       xPath = TRIM(xml%speciesPath(n))//'/energyParameters'
       this%qn_el(0,n,:)=evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(xPath)//'/@s'))
       this%qn_el(1,n,:)=evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(xPath)//'/@p'))
       this%qn_el(2,n,:)=evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(xPath)//'/@d'))
       this%qn_el(3,n,:)=evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(xPath)//'/@f'))
       nLOLocal = 0
       DO iLO = 1, xml%getNumberOfNodes(TRIM(ADJUSTL(xml%speciesPath(n))//'/lo'))
          WRITE(xPathB,"(a,a,i0,a)") TRIM(xml%speciesPath(n)),'/lo[',iLO,']'
          lString = xml%getAttributeValue(TRIM(ADJUSTL(xPathB))//'/@l')
          nString = xml%getAttributeValue(TRIM(ADJUSTL(xPathB))//'/@n')
          CALL xml%getIntegerSequenceFromString(TRIM(ADJUSTL(lString)), lNumbers, lNumCount)
          CALL xml%getIntegerSequenceFromString(TRIM(ADJUSTL(nString)), nNumbers, nNumCount)
          IF(lNumCount.NE.nNumCount) THEN
             CALL judft_error('Error in LO input: l quantum number count does not equal n quantum number count')
          END IF
          DO i = 1, lNumCount
             this%qn_ello(nLOLocal+i,n,:) = nNumbers(i)
             IF (TRIM(ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(xPathB))//'/@type')))=='HELO') THEN
                this%qn_ello(nLOLocal+i,n,:)=-1*this%qn_ello(nLOLocal+i,n,:)
             END IF
          END DO
          nLOLocal = nLOLocal +lnumcount
          DEALLOCATE (lNumbers, nNumbers)
       END DO
    END DO
    !Read vacuum
    IF (xml%GetNumberOfNodes("/fleurInput/cell/filmLattice")==1) THEN
       DO i=1,xml%GetNumberOfNodes("/fleurInput/cell/filmLattice/vacuumEnergyParameters")
          WRITE(xPath,"(a,i0,a)") '/fleurInput/cell/filmLattice/vacuumEnergyParameters[',i,']'
          n = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPath))//'/@vacuum'))
          this%evac0(n,:) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPath))//'/@spinUp'))
          IF (xml%GetNumberOfNodes(TRIM(ADJUSTL(xPath))//'/@spinDown').GE.1) THEN
             this%evac0(n,2) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPath))//'/@spinDown'))
          END IF
          IF (i==1.AND.n==1) this%evac0(2,:)=this%evac0(1,:)
       END DO
    END IF


  END SUBROUTINE read_xml_enpara

  SUBROUTINE set_quantum_numbers(enpara,ntype,atoms,str,lo,addLOs)
    USE m_types_atoms
    !sets the energy parameters according to simple electronic config string and lo string
    CLASS(t_enparaXML),INTENT(inout):: enpara
    TYPE(t_atoms),INTENT(IN)     :: atoms
    INTEGER,INTENT(in)           :: ntype
    CHARACTER(len=*),INTENT(in)  :: str,lo
    INTEGER, INTENT(IN)          :: addLOs(atoms%ntype)

    CHARACTER(len=100):: val_str
    character         :: ch
    INTEGER           :: qn,n,i,l

    !Process lo's
    DO i=1,LEN_TRIM(lo)/2
       READ(lo(2*i-1:2*i),"(i1,a1)") qn,ch
       enpara%qn_ello(i,ntype,:)=qn
    ENDDO

    !Valence string
    val_str=ADJUSTL(str(INDEX(str,"|")+1:))
    DO WHILE(LEN_TRIM(val_str)>1)
       READ(val_str,"(i1,a1)") qn,ch
       l=INDEX("spdf",ch)-1
       !check if we have an lo for this l-channel
       DO i=1,atoms%nlo(ntype) - addLOs(ntype)
          IF (l==atoms%llo(i,ntype).AND.qn==enpara%qn_ello(i,ntype,1)) EXIT
       ENDDO
       IF (atoms%nlo(ntype)==0.OR.i>atoms%nlo(ntype)) THEN
          !set the lapw parameter
          IF (enpara%qn_el(l,ntype,1)<0) CALL judft_error("Electronic configuration needs more LOs")
          enpara%qn_el(l,ntype,:)=-qn
       ENDIF
       IF (LEN_TRIM(val_str)>5) THEN
          val_str=ADJUSTL(val_str(5:))
       ELSE
          val_str=""
       ENDIF
    ENDDO
    enpara%qn_el(:,ntype,:)=ABS(enpara%qn_el(:,ntype,:))
  END SUBROUTINE set_quantum_numbers

  SUBROUTINE Init(This,Ntype,Nlod,Jspins,L_defaults,Nz)
    USE m_constants
    CLASS(t_enparaXML),INTENT(inout):: this
    INTEGER,INTENT(IN)           :: jspins,nlod,ntype
    LOGICAL,INTENT(IN),OPTIONAL  :: l_defaults
    INTEGER,INTENT(IN),OPTIONAL  :: nz(:)

    INTEGER :: n,i,jsp,l

    this%evac0=-1E99
    if (allocated(this%qn_el)) deallocate(this%qn_el)
    if (allocated(this%qn_ello)) deallocate(this%qn_ello)

    ALLOCATE(this%qn_el(0:3,ntype,jspins))
    ALLOCATE(this%qn_ello(nlod,ntype,jspins))

    IF (PRESENT(l_defaults)) THEN
       IF (.NOT.l_defaults) RETURN
    ELSE
       RETURN
    ENDIF
    !Set most simple defaults
    DO jsp=1,jspins
       DO n = 1,ntype
          IF ( nz(n) < 3 ) THEN
             this%qn_el(0:3,n,jsp) =  (/1,2,3,4/)
          ELSEIF ( nz(n) < 11 ) THEN
             this%qn_el(0:3,n,jsp) =  (/2,2,3,4/)
          ELSEIF ( nz(n) < 19 ) THEN
             this%qn_el(0:3,n,jsp) =  (/3,3,3,4/)
          ELSEIF ( nz(n) < 31 ) THEN
             this%qn_el(0:3,n,jsp) =  (/4,4,3,4/)
          ELSEIF ( nz(n) < 35 ) THEN
             this%qn_el(0:3,n,jsp) =  (/4,4,3,4/)
          ELSEIF ( nz(n) < 37 ) THEN
             this%qn_el(0:3,n,jsp) =  (/4,4,4,4/)
          ELSEIF ( nz(n) < 49 ) THEN
             this%qn_el(0:3,n,jsp) =  (/5,5,4,4/)
          ELSEIF ( nz(n) < 55 ) THEN
             this%qn_el(0:3,n,jsp) =  (/5,5,5,4/)
          ELSEIF ( nz(n) < 72 ) THEN
             this%qn_el(0:3,n,jsp) =  (/6,6,5,4/)
          ELSEIF ( nz(n) < 81 ) THEN
             this%qn_el(0:3,n,jsp) =  (/6,6,5,5/)
          ELSEIF ( nz(n) < 87 ) THEN
             this%qn_el(0:3,n,jsp) =  (/6,6,6,5/)
          ELSE
             this%qn_el(0:3,n,jsp) =  (/7,7,6,5/)
          ENDIF

          this%qn_ello(:,n,jsp) = 0
       ENDDO
    ENDDO

    this%evac0=eVac0Default_const

  END SUBROUTINE init

END MODULE m_types_enparaXML
