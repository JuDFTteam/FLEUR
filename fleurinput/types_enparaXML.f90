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
     INTEGER,ALLOCATABLE  :: qn_el(:,:,:)    !>if these are .ne.0 they are understood as
     INTEGER,ALLOCATABLE  :: qn_ello(:,:,:)  !>quantum numbers
     REAL                 :: evac0(2,2)
   CONTAINS
     PROCEDURE :: init
     procedure :: set_quantum_numbers
     PROCEDURE :: read_xml=>read_xml_enpara
  END TYPE t_enparaXML



  PUBLIC:: t_enparaXML

CONTAINS

  SUBROUTINE read_xml_enpara(this,xml)
    use m_types_xml
    CLASS(t_enparaXML),INTENT(INOUT):: this
    TYPE(t_xml),INTENT(IN)     :: xml

    LOGICAL :: l_enpara,film
    INTEGER :: jspins,ntype,lmaxd,n,lo,i
    CHARACTER(len=100)::xpath,xpath2
    INTEGER, ALLOCATABLE :: nlo(:),neq(:)

    ntype=xml%get_ntype()
    nlo=xml%get_nlo()
    lmaxd=xml%get_lmaxd()
    jspins=evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/magnetism/@jspins'))
    film=xml%GetNumberOfNodes('/fleurInput/cell/filmLattice')==1
    
    CALL this%init(ntype,MAXVAL(nlo),jspins)

    DO n=1,ntype
       xPath=TRIM(xml%speciesPath(n))//'/energyParameters'
       this%qn_el(0,n,:)=evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(xpath)//'/@s'))
       this%qn_el(1,n,:)=evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(xpath)//'/@p'))
       this%qn_el(2,n,:)=evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(xpath)//'/@d'))
       this%qn_el(3,n,:)=evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(xpath)//'/@f'))
       DO lo=1,nlo(n)
          WRITE(xpath2,"(a,a,i0,a)") TRIM(xml%speciesPath(n)),'/lo[',lo,']'
          this%qn_ello(lo,n,:)=evaluateFirstINTOnly(xml%GetAttributeValue(TRIM(xpath)//'/@n'))
          IF (TRIM(ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'@type')))=='HELO') &
               this%qn_ello(lo,n,:)=-1*this%qn_ello(lo,n,:)
       END DO
    END DO
    !Read vacuum
    IF (xml%GetNumberOfNodes("/fleurInput/cell/filmLattice")==1) THEN
       DO i=1,xml%GetNumberOfNodes("/fleurInput/cell/filmLattice/vacuumEnergyParameters")
          WRITE(xpath,"(a,i0,a)") '/fleurInput/cell/filmLattice/vacuumEnergyParameters[',i,']'
          n = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPath))//'/@vacuum'))
          this%evac0(n,:) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPath))//'/@spinUp'))
          IF (xml%GetNumberOfNodes(TRIM(ADJUSTL(xPath))//'/@spinDown').GE.1) THEN
             this%evac0(n,2) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPath))//'/@spinDown'))
          END IF
          IF (i==1.AND.n==1) this%evac0(2,:)=this%evac0(1,:)
       END DO
    END IF
    

  END SUBROUTINE read_xml_enpara

  SUBROUTINE set_quantum_numbers(enpara,ntype,atoms,str,lo)
    use m_types_atoms
    !sets the energy parameters according to simple electronic config string and lo string
    CLASS(t_enparaXML),INTENT(inout):: enpara
    TYPE(t_atoms),INTENT(IN)     :: atoms
    INTEGER,INTENT(in)           :: ntype
    CHARACTER(len=*),INTENT(in)  :: str,lo
    
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
       DO i=1,atoms%nlo(ntype)
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
  
  SUBROUTINE init(this,ntype,nlod,jspins,l_defaults,nz)
    USE m_constants
    CLASS(t_enparaXML),INTENT(inout):: this
    INTEGER,INTENT(IN)           :: jspins,nlod,ntype
    LOGICAL,INTENT(IN),OPTIONAL  :: l_defaults
    INTEGER,INTENT(IN),OPTIONAL  :: nz(:)

    INTEGER :: n,i,jsp,l

    this%evac0=-1E99
    ALLOCATE(this%qn_el(0:3,ntype,jspins))
    ALLOCATE(this%qn_ello(nlod,ntype,jspins))

    IF (PRESENT(l_defaults)) THEN
       IF (.NOT.l_defaults) RETURN
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
