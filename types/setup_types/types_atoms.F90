!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_atoms
  USE m_judft
  USE m_types_fleur_setup
  use m_json_tools
  IMPLICIT NONE
  
  TYPE t_utype !for LDA+U
     SEQUENCE
     REAL u,j         ! the actual U and J parameters
     INTEGER l        ! the l quantum number to which this U parameter belongs
     INTEGER atomType ! The atom type to which this U parameter belongs
     LOGICAL :: l_amf ! logical switch to choose the "around mean field" LDA+U limit
  END TYPE t_utype
  
  TYPE,EXTENDS(t_fleursetup):: t_atoms
        !<no of types
     INTEGER :: ntype
     !<total-no of atoms
     INTEGER :: nat
     !<dimensions of LO's
     INTEGER ::nlod
     INTEGER ::llod
     INTEGER ::nlotot
     !lmaxd=maxval(lmax)
     INTEGER:: lmaxd
     ! no of lda+us
     INTEGER ::n_u
     ! dimensions
     INTEGER :: jmtd
     !No of element
     INTEGER,ALLOCATABLE ::nz(:)
     !atoms per type
     INTEGER,ALLOCATABLE::neq(:)
     !radial grid points
     INTEGER,ALLOCATABLE::jri(:)
     !core states
     INTEGER,ALLOCATABLE::ncst(:)
     !How many states are explicitely provided?
     INTEGER,ALLOCATABLE::numStatesProvided(:)
     !core state occupations
     REAL,ALLOCATABLE::coreStateOccs(:,:,:)
     !core state nprnc
     INTEGER,ALLOCATABLE::coreStateNprnc(:,:)
     !core state kappa
     INTEGER,ALLOCATABLE::coreStateKappa(:,:)
     !lmax
     INTEGER,ALLOCATABLE::lmax(:)
     !lmax non-spherical
     INTEGER,ALLOCATABLE::lnonsph(:)
     !expansion of pseudo-charge
     INTEGER,ALLOCATABLE::ncv(:)
     !no of LO
     INTEGER,ALLOCATABLE::nlo(:)
     !l of LO (nlo,ntype)
     INTEGER,ALLOCATABLE::llo(:,:)
     !lmax for lapw (ntype)
     INTEGER,ALLOCATABLE::lapw_l(:)
     !first LO with a given l (max(nlo
     INTEGER,ALLOCATABLE::lo1l(:,:)
     !??
     INTEGER,ALLOCATABLE::ulo_der(:,:)
     !no of LOs per l (max(nlo1),ntype
     INTEGER,ALLOCATABLE::nlol(:,:)
     !true if LO is formed by \dot u (
     LOGICAL,ALLOCATABLE::l_dulo(:,:)
     !no of op that maps atom into
     INTEGER,ALLOCATABLE::ngopr(:)
     !symetry of atom (nat)
     INTEGER,ALLOCATABLE::ntypsy(:)
     !no of sphhar for atom type(ntype
     INTEGER,ALLOCATABLE ::nlhtyp(:)
     !atom mapped to by inversion (nat
     INTEGER,ALLOCATABLE ::invsat(:)
     !Calaculate forces for this atom?
     LOGICAL,ALLOCATABLE :: l_geo(:)
     !MT-Radius (ntype)
     REAL,ALLOCATABLE CPP_MANAGED::rmt(:)
     !log increment(ntype)
     REAL,ALLOCATABLE::dx(:)
     !vol of MT(ntype)
     REAL,ALLOCATABLE::volmts(:)
     !radial grid points(max(jri),ntyp
     REAL,ALLOCATABLE::rmsh(:,:)
     !charge of nucleus(ntype)
     REAL,ALLOCATABLE::zatom(:)
     !initial mag moment(ntype)
     REAL,ALLOCATABLE::bmu(:)
     !pos of atom (absol) (3,nat)
     REAL,ALLOCATABLE::pos(:,:)
     !pos of atom (relat)(3,nat)
     REAL,ALLOCATABLE CPP_MANAGED::taual(:,:)  
     !labels
     CHARACTER(LEN=20), ALLOCATABLE :: label(:)          !Not saved or broadcasted
     CHARACTER(len=20), ALLOCATABLE :: speciesName(:)    !Not saved or broadcasted
     !lda_u information(ntype)
     TYPE(t_utype),ALLOCATABLE::lda_u(:)
     INTEGER,ALLOCATABLE :: relax(:,:) !<(3,ntype)
     INTEGER, ALLOCATABLE :: nflip(:) !<flip magnetisation of this atom
  
   CONTAINS
     PROCEDURE,PASS :: broadcast=>broadcast_atoms
     PROCEDURE,PASS :: write=>WRITE_atoms
     PROCEDURE,PASS :: read=>READ_atoms
     PROCEDURE,PASS :: read_xml=>read_xml_atoms
  END TYPE t_atoms

CONTAINS
  SUBROUTINE broadcast_atoms(tt,mpi_comm,origin)
#ifdef CPP_MPI
    USE mpi_bc_tool
#endif    
    IMPLICIT NONE
    CLASS(t_atoms),INTENT(INOUT):: tt
    INTEGER,INTENT(IN)               :: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL      :: origin

#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER :: pe,ierr,irank,n
    CALL MPI_COMM_RANK(mpi_comm,irank,ierr)
    IF (PRESENT(origin)) THEN
       pe=origin
    ELSE
       pe=0
    ENDIF

    CALL MPI_BCAST(tt%ntype,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%nat,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%nlod,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%llod,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%nlotot,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%lmaxd,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%n_u,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%jmtd,1,MPI_INTEGER,pe,mpi_comm,ierr)
 
    CALL MPI_BC(tt%nz,pe,mpi_comm)
    CALL mpi_bc(tt%neq,pe,mpi_comm)
    CALL mpi_bc(tt%jri,pe,mpi_comm)
    CALL mpi_bc(tt%ncst,pe,mpi_comm)
    CALL mpi_bc(tt%numStatesProvided,pe,mpi_comm)
    CALL mpi_bc(tt%coreStateOccs,pe,mpi_comm)
    CALL mpi_bc(tt%coreStateNprnc,pe,mpi_comm)
    CALL mpi_bc(tt%coreStateKappa,pe,mpi_comm):,:)
    CALL mpi_bc(tt%lmax,pe,mpi_comm)
    
    CALL mpi_bc(tt%lnonsph,pe,mpi_comm)
    CALL mpi_bc(tt%ncv,pe,mpi_comm)
    CALL mpi_bc(tt%nlo,pe,mpi_comm)
    CALL mpi_bc(tt%llo,pe,mpi_comm)
    CALL mpi_bc(tt%lapw_l,pe,mpi_comm)
    CALL mpi_bc(tt%lo1l,pe,mpi_comm)
    CALL mpi_bc(tt%ulo_der,pe,mpi_comm)
    CALL mpi_bc(tt%nlol,pe,mpi_comm)
    CALL mpi_bc(tt%l_dulo,pe,mpi_comm)
    CALL mpi_bc(tt%ngopr,pe,mpi_comm)
    CALL mpi_bc(tt%ntypsy,pe,mpi_comm)
    CALL mpi_bc(tt%nlhtyp,pe,mpi_comm)
    CALL mpi_bc(tt%invsat,pe,mpi_comm)
    CALL mpi_bc(tt%l_geo,pe,mpi_comm)
    CALL mpi_bc(tt%rmt,pe,mpi_comm)
    CALL mpi_bc(tt%dx,pe,mpi_comm)
    CALL mpi_bc(tt%volmts,pe,mpi_comm)
    CALL mpi_bc(tt%rmsh,pe,mpi_comm)
    CALL mpi_bc(tt%zatom,pe,mpi_comm)
    CALL mpi_bc(tt%bmu,pe,mpi_comm)
    CALL mpi_bc(tt%pos,pe,mpi_comm)
    CALL mpi_bc(tt%taual,pe,mpi_comm)  
    CALL mpi_bc(tt%relax,pe,mpi_comm)
     CALL mpi_bc(tt%nflip,pe,mpi_comm)
     !lda_u information(ntype)
     IF (irank.NE.pe) THEN
        IF (ALLOCATED(atoms%lda_u)) DEALLOCATE(atoms%lda_u)
        ALLOCATE(atoms%lda_u(4*atoms%ntype))
     END IF
     DO n=1,atoms%ntype
        CALL MPI_BCAST(tt%lda_u(n)%u,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
        CALL MPI_BCAST(tt%lda_u(n)%j,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
        CALL MPI_BCAST(tt%lda_u(n)%l,MPI_INTEGER,pe,mpi_comm,ierr)
        CALL MPI_BCAST(tt%lda_u(n)%atomtype,MPI_INTEGER,pe,mpi_comm,ierr)
        CALL MPI_BCAST(tt%lda_u(n)%l_amf,MPI_LOGICAL,pe,mpi_comm,ierr)
     END DO
#endif
      
  END SUBROUTINE broadcast_atoms

  SUBROUTINE write_atoms(tt, unit, iotype, v_list, iostat, iomsg)
    IMPLICIT NONE
    CLASS(t_atoms),INTENT(IN):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    CHARACTER(len=20)::s
    INTEGER :: n

    WRITE(unit,*,IOSTAT=iostat) '"atoms":{'

CALL json_print(unit,"ntype",tt%ntype)
CALL json_print(unit,"nat",tt%nat)
CALL json_print(unit,"nlod",tt%nlod)
CALL json_print(unit,"llod",tt%llod)
CALL json_print(unit,"nlotot",tt%nlotot)
CALL json_print(unit,"lmaxd",tt%lmaxd)
CALL json_print(unit,"n_u",tt%n_u)
CALL json_print(unit,"jmtd",tt%jmtd)
 
CALL json_print(unit,"nz",tt%nz)
CALL json_print(unit,"neq",tt%neq)
CALL json_print(unit,"jri",tt%jri)
CALL json_print(unit,"ncst",tt%ncst)
CALL json_print(unit,"numStatesProvided",tt%numStatesProvided)
CALL json_print(unit,"coreStateOccs",tt%coreStateOccs)
CALL json_print(unit,"coreStateNprnc",tt%coreStateNprnc)
CALL json_print(unit,"coreStateKappa",tt%coreStateKappa)
CALL json_print(unit,"lmax",tt%lmax)
    
CALL json_print(unit,"lnonsph",tt%lnonsph)
CALL json_print(unit,"ncv",tt%ncv)
CALL json_print(unit,"nlo",tt%nlo)
CALL json_print(unit,"llo",tt%llo)
CALL json_print(unit,"lapw_l",tt%lapw_l)
CALL json_print(unit,"lo1l",tt%lo1l)
CALL json_print(unit,"ulo_der",tt%ulo_der)
CALL json_print(unit,"nlol",tt%nlol)
CALL json_print(unit,"l_dulo",tt%l_dulo)
CALL json_print(unit,"ngopr",tt%ngopr)
CALL json_print(unit,"ntypsy",tt%ntypsy)
CALL json_print(unit,"nlhtyp",tt%nlhtyp)
CALL json_print(unit,"invsat",tt%invsat)
CALL json_print(unit,"l_geo",tt%l_geo)
CALL json_print(unit,"rmt",tt%rmt)
CALL json_print(unit,"dx",tt%dx)
CALL json_print(unit,"volmts",tt%volmts)
CALL json_print(unit,"rmsh",tt%rmsh)
CALL json_print(unit,"zatom",tt%zatom)
CALL json_print(unit,"bmu",tt%bmu)
CALL json_print(unit,"pos",tt%pos)
CALL json_print(unit,"taual",tt%taual)
CALL json_print(unit,"relax",tt%relax)
CALL json_print(unit,"nflip",tt%nflip)

DO n=1,tt%ntype
   WRITE(s,"(a,i0,a)") "ldau(",n,")%"
   CALL json_print(unit,trim(s)//"u",tt%lda_u(n)%u)
   CALL json_print(unit,trim(s)//"j",tt%lda_u(n)%j)
   CALL json_print(unit,trim(s)//"l",tt%lda_u(n)%l)
   CALL json_print(unit,trim(s)//"atomtype",tt%lda_u(n)%atomtype)
   CALL json_print(unit,trim(s)//"l_amf",tt%lda_u(n)%l_amf)
END DO
 
    WRITE(unit,*,IOSTAT=iostat) '}'
    
  END SUBROUTINE write_atoms
  SUBROUTINE read_atoms(tt, unit, iotype, v_list, iostat, iomsg)
    IMPLICIT NONE
    CLASS(t_atoms),INTENT(INOUT):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    CHARACTER(len=40):: s
    INTEGER          :: n
    CALL json_open_class("atoms",unit,iostat)
    IF (iostat.NE.0)   RETURN

   CALL json_read(unit,"ntype",tt%ntype)
CALL json_read(unit,"nat",tt%nat)
CALL json_read(unit,"nlod",tt%nlod)
CALL json_read(unit,"llod",tt%llod)
CALL json_read(unit,"nlotot",tt%nlotot)
CALL json_read(unit,"lmaxd",tt%lmaxd)
CALL json_read(unit,"n_u",tt%n_u)
CALL json_read(unit,"jmtd",tt%jmtd)
 
CALL json_read(unit,"nz",tt%nz)
CALL json_read(unit,"neq",tt%neq)
CALL json_read(unit,"jri",tt%jri)
CALL json_read(unit,"ncst",tt%ncst)
CALL json_read(unit,"numStatesProvided",tt%numStatesProvided)
CALL json_read(unit,"coreStateOccs",tt%coreStateOccs)
CALL json_read(unit,"coreStateNprnc",tt%coreStateNprnc)
CALL json_read(unit,"coreStateKappa",tt%coreStateKappa)
CALL json_read(unit,"lmax",tt%lmax)
    
CALL json_read(unit,"lnonsph",tt%lnonsph)
CALL json_read(unit,"ncv",tt%ncv)
CALL json_read(unit,"nlo",tt%nlo)
CALL json_read(unit,"llo",tt%llo)
CALL json_read(unit,"lapw_l",tt%lapw_l)
CALL json_read(unit,"lo1l",tt%lo1l)
CALL json_read(unit,"ulo_der",tt%ulo_der)
CALL json_read(unit,"nlol",tt%nlol)
CALL json_read(unit,"l_dulo",tt%l_dulo)
CALL json_read(unit,"ngopr",tt%ngopr)
CALL json_read(unit,"ntypsy",tt%ntypsy)
CALL json_read(unit,"nlhtyp",tt%nlhtyp)
CALL json_read(unit,"invsat",tt%invsat)
CALL json_read(unit,"l_geo",tt%l_geo)
CALL json_read(unit,"rmt",tt%rmt)
CALL json_read(unit,"dx",tt%dx)
CALL json_read(unit,"volmts",tt%volmts)
CALL json_read(unit,"rmsh",tt%rmsh)
CALL json_read(unit,"zatom",tt%zatom)
CALL json_read(unit,"bmu",tt%bmu)
CALL json_read(unit,"pos",tt%pos)
CALL json_read(unit,"taual",tt%taual)
CALL json_read(unit,"relax",tt%relax)
CALL json_read(unit,"nflip",tt%nflip)

IF (ALLOCATED(tt%lda_u)) DEALLOCATE(tt%lda_u)
ALLOCATE(tt%lda_u(4*tt%ntype))
DO n=1,tt%ntype
   WRITE(s,"(a,i0,a)") "ldau(",n,")%"
   CALL json_read(unit,trim(s)//"u",tt%lda_u(n)%u)
   CALL json_read(unit,trim(s)//"j",tt%lda_u(n)%j)
   CALL json_read(unit,trim(s)//"l",tt%lda_u(n)%l)
   CALL json_read(unit,trim(s)//"atomtype",tt%lda_u(n)%atomtype)
   CALL json_read(unit,s//"l_amf",tt%lda_u(n)%l_amf)
END DO

    CALL json_close_class(unit,iostat)
    
  END SUBROUTINE read_atoms

 


  SUBROUTINE read_xml_atoms(tt)
    USE m_xmlIntWrapFort
    USE m_constants
    USE m_calculator
    USE m_econfig
    use m_inp_xml
    IMPLICIT NONE
    CLASS(t_atoms),INTENT(OUT):: tt


    CHARACTER(len=200):: xpaths,xpathg,xpath,valueString,lstring,nstring
    INTEGER           :: i,n,j,numberNodes,ilo,lNumCount,nNumCount
    INTEGER,ALLOCATABLE::lNumbers(:),nNumbers(:)
    LOGICAL           :: relaxx,relaxy,relaxz
    INTEGER,ALLOCATABLE :: itmp(:,:)
    
    tt%ntype= xmlGetNumberOfNodes('/fleurInput/atomGroups/atomGroup')
    tt%nat =  xmlGetNumberOfNodes('/fleurInput/atomGroups/atomGroup/relPos')&
         + xmlGetNumberOfNodes('/fleurInput/atomGroups/atomGroup/absPos')&
         + xmlGetNumberOfNodes('/fleurInput/atomGroups/atomGroup/filmPos')
    ALLOCATE(tt%nz(tt%ntype))     !nz and zatom have the same content!
    ALLOCATE(tt%zatom(tt%ntype))  !nz and zatom have the same content!
    ALLOCATE(tt%jri(tt%ntype))
    ALLOCATE(tt%dx(tt%ntype))
    ALLOCATE(tt%lmax(tt%ntype))
    ALLOCATE(tt%nlo(tt%ntype))
    ALLOCATE(tt%ncst(tt%ntype))
    ALLOCATE(tt%lnonsph(tt%ntype))
    ALLOCATE(tt%nflip(tt%ntype))
    ALLOCATE(tt%l_geo(tt%ntype))
    ALLOCATE(tt%lda_u(4*tt%ntype))
    ALLOCATE(tt%bmu(tt%ntype))
    ALLOCATE(tt%relax(3,tt%ntype))
    ALLOCATE(tt%neq(tt%ntype))
    ALLOCATE(tt%taual(3,tt%nat))
    ALLOCATE(tt%label(tt%nat))
    ALLOCATE(tt%pos(3,tt%nat))
    ALLOCATE(tt%rmt(tt%ntype))
    ALLOCATE(tt%numStatesProvided(tt%ntype))
    ALLOCATE(tt%ncv(tt%ntype)) ! For what is this?
    ALLOCATE(tt%ngopr(tt%nat)) ! For what is this?
    ALLOCATE(tt%lapw_l(tt%ntype)) ! Where do I put this?
    ALLOCATE(tt%invsat(tt%nat)) ! Where do I put this?
    tt%numStatesProvided = 0
    tt%lapw_l(:) = -1
    tt%n_u = 0
    DO n = 1, tt%ntype
       !in Species:
       !@name,element,atomicNumber,coreStates
       !@optional: magMom,flipSpin,magField,vcaAddCharge
       !mtSphere,atomicCutoffs
       !optional: energyParameters,prodBasis,special,force,electronConfig,nocoParams,ldaU(up to 4),lo(as many as needed)
       xpathg=inp_xml_xpath_for_group(n)
       xpaths=inp_xml_speciesxpath_for_group(n)
       tt%nz(n)=evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPaths))//'/@atomicNumber'))
       IF (tt%nz(n).EQ.0) THEN
          WRITE(*,*) 'Note: Replacing atomic number 0 by 1.0e-10 on atom type ', n
          tt%zatom(n) = 1.0e-10
       END IF
       tt%zatom(n) = tt%nz(n)
       tt%ncst(n) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPaths))//'/@coreStates'))
      
       
       IF (evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPaths))//'/@flipSpin'))) THEN
          tt%nflip(n) = 1
       ELSE
          tt%nflip(n) = 0
       ENDIF
       tt%bmu(n) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPaths))//'/@magMom'))
       !Now the xml elements
       !mtSphere
       xpath=xpaths
       IF (xmlGetNumberOfNodes(TRIM(xpathg)//'/mtSphere')==1) xpath=xpathg
       tt%rmt(n) =  evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/mtSphere/@radius'))
       tt%jri(n) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/mtSphere/@gridPoints'))
       tt%dx(n) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/mtSphere/@logIncrement'))
       !atomicCuttoffs
       xpath=xpaths
       IF (xmlGetNumberOfNodes(TRIM(xpathg)//'/atomicCutoffs')==1) xpath=xpathg
       tt%lmax(n) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/atomicCutoffs/@lmax'))
       tt%lnonsph(n) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/atomicCutoffs/@lnonsphr'))
       tt%lapw_l(n) = -1
       IF (xmlGetNumberOfNodes(TRIM(ADJUSTL(xPath))//'/atomicCutoffs/@lmaxAPW').EQ.1) &
            tt%lapw_l(n) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/atomicCutoffs/@lmaxAPW'))
       !force type
       xpath=''
       IF(xmlGetNumberOfNodes(TRIM(ADJUSTL(xPaths))//'/force')==1) xpath=xpaths
       IF (xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathg))//'/force')==1) xpath=xpathg
       if (xpath.ne.'') THEN
          tt%l_geo(n) = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathg))//'/force/@calculate'))
          valueString = xmlGetAttributeValue(TRIM(ADJUSTL(xPathg))//'force/@relaxXYZ')
          READ(valueString,'(3l1)') relaxX, relaxY, relaxZ
          IF (relaxX) tt%relax(1,n) = 1
          IF (relaxY) tt%relax(2,n) = 1
          IF (relaxZ) tt%relax(3,n) = 1
       ELSE
          tt%l_geo(n) = .FALSE.
           tt%relax(:,n) = 0
        END IF
        !LO stuff  
        !(As we have not determined tt%nlod yet, we have to reallocate tt%llo if it becomes too small)
        tt%nlo(n) = 0
        DO ilo = 1,xmlGetNumberOfNodes(TRIM(ADJUSTL(xpaths))//'/lo')+xmlGetNumberOfNodes(TRIM(ADJUSTL(xpathg))//'/lo')
           IF (ilo>xmlGetNumberOfNodes(TRIM(ADJUSTL(xpaths))//'/lo')) THEN
              WRITE(xpath,*) TRIM(ADJUSTL(xpathg))//'/lo[',ilo-xmlGetNumberOfNodes(TRIM(ADJUSTL(xpaths))//'/lo'),']'
           ELSE
              WRITE(xpath,*) TRIM(ADJUSTL(xpaths))//'/lo[',ilo,']'
           END IF
          lString = xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@l')
          nString = xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@n')
          CALL getIntegerSequenceFromString(TRIM(ADJUSTL(lString)), lNumbers, lNumCount)
          CALL getIntegerSequenceFromString(TRIM(ADJUSTL(nString)), nNumbers, nNumCount)
          IF(lNumCount.NE.nNumCount) THEN
             CALL judft_error('Error in LO input: l quantum number count does not equal n quantum number count')
          END IF
          IF (.NOT. ALLOCATED(tt%llo)) ALLOCATE(tt%llo(1,tt%ntype),tt%ulo_der(1,tt%ntype))
          IF (SIZE(tt%llo,1)<lNumCount+tt%nlo(n)) THEN
             CALL move_ALLOC(tt%llo,itmp)
             ALLOCATE(tt%llo(lNumCount+tt%nlo(n),tt%ntype))
             tt%llo(:SIZE(itmp,1),:)=itmp
             CALL move_ALLOC(tt%ulo_der,itmp)
             ALLOCATE(tt%ulo_der(lNumCount+tt%nlo(n),tt%ntype))
             tt%ulo_der(:SIZE(itmp,1),:)=itmp
          ENDIF
          DO i=1,lNumCount
             tt%llo(tt%nlo(n)+i,n) = lNumbers(i)
             tt%ulo_der(tt%nlo(n)+i,n) = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@eDeriv'))
          ENDDO
          tt%nlo(n) = tt%nlo(n) +lnumcount 
          DEALLOCATE (lNumbers, nNumbers)
       END DO
       !LDA+U
       DO i = 1,xmlGetNumberOfNodes(TRIM(ADJUSTL(xPaths))//'/ldaU')+xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathg))//'/ldaU')
          IF (i>xmlGetNumberOfNodes(TRIM(ADJUSTL(xPaths))//'/ldaU')) THEN
             WRITE(xpath,*) TRIM(ADJUSTL(xPathg))//'/ldaU[',i-xmlGetNumberOfNodes(TRIM(ADJUSTL(xPaths))//'/ldaU'),']'
          ELSE
             WRITE(xpath,*) TRIM(ADJUSTL(xPaths))//'/ldaU[',i,']'
          END IF
          IF (i.GT.4) CALL juDFT_error("Too many U parameters provided for a certain species (maximum is 4).",calledby ="types_atoms")
          tt%n_u = tt%n_u + 1
          tt%lda_u(tt%n_u)%l = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@l'))
          
          tt%lda_u(tt%n_u)%u =  evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@U'))
          tt%lda_u(tt%n_u)%j = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@J'))
          tt%lda_u(tt%n_u)%l_amf =  evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/@l_amf'))
          tt%lda_u(tt%n_u)%atomType = n
       END DO
       !electron config
       tt%numStatesProvided(n)=0
       IF (xmlGetNumberOfNodes(TRIM(ADJUSTL(xmlGetAttributeValue(xPaths)))//'/electronConfig')==1) THEN
          CALL parse_econfig(&
               xmlGetAttributeValue(TRIM(ADJUSTL(xmlGetAttributeValue(xPaths)))//'/electronConfig/coreConfig'),&
               .TRUE.,tt%coreStateOccs(:,:,n),tt%coreStateNprnc(:,n),tt%coreStateKappa(:,n)&
               ,tt%numStatesProvided(n))
          IF (xmlGetNumberOfNodes(TRIM(ADJUSTL(xmlGetAttributeValue(xPaths)))//'/electronConfig/valenceConfig')==1) &
               CALL parse_econfig(&
               xmlGetAttributeValue(TRIM(ADJUSTL(xmlGetAttributeValue(xPaths)))//'/electronConfig/valenceConfig'),&
               .FALSE.,tt%coreStateOccs(:,:,n),tt%coreStateNprnc(:,n),tt%coreStateKappa(:,n),&
               tt%numStatesProvided(n))
          
          numberNodes = xmlGetNumberOfNodes(TRIM(ADJUSTL(xmlGetAttributeValue(xPaths)))//'/electronConfig/stateOccupation')
          IF (numberNodes.GE.1) THEN
             DO i = 1, numberNodes
                WRITE(xpath,*) TRIM(ADJUSTL(xPath)),'/electronConfig/stateOccupation[',i,']'
                CALL parse_occupation(xpath,tt%coreStateOccs(:,:,n),&
                     tt%coreStateNprnc(:,n),tt%coreStateKappa(:,n),&
                     tt%numStatesProvided(n))
             END DO
          END IF
       END IF
     
    END DO
      
  END SUBROUTINE read_xml_atoms


  SUBROUTINE init_atoms(atoms)
    IMPLICIT NONE
    CLASS(t_atoms),INTENT(OUT):: atoms
    
  END SUBROUTINE init_atoms

    
END MODULE m_types_atoms
