!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
#ifndef CPP_MANAGED
#define CPP_MANAGED 
#endif
MODULE m_types_atoms
  use m_juDFT
  USE m_types_econfig
  IMPLICIT NONE
  PRIVATE
  TYPE t_utype
      SEQUENCE
      REAL :: u, j         ! the actual U and J parameters
      REAL :: theta,phi   !the rotation angles by which the density metrics is rotated
      INTEGER :: l        ! the l quantum number to which this U parameter belongs
      INTEGER :: atomType ! The atom type to which this U parameter belongs
      LOGICAL :: l_amf ! logical switch to choose the "around mean field" LDA+U limit
   END TYPE t_utype


   TYPE t_atoms
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
      INTEGER, ALLOCATABLE ::nz(:)
      !atoms per type
      INTEGER, ALLOCATABLE::neq(:)
      !radial grid points
      INTEGER, ALLOCATABLE::jri(:)
      !core states
      TYPE(t_econfig),ALLOCATABLE::econf(:)
      !lmax
      INTEGER, ALLOCATABLE::lmax(:)
      !lmax non-spherical
      INTEGER, ALLOCATABLE::lnonsph(:)
      !expansion of pseudo-charge
      INTEGER, ALLOCATABLE::ncv(:)
      !no of LO
      INTEGER, ALLOCATABLE::nlo(:)
      !l of LO (nlo,ntype)
      INTEGER, ALLOCATABLE::llo(:, :)
      !lmax for lapw (ntype)
      INTEGER, ALLOCATABLE::lapw_l(:)
      !first LO with a given l (max(nlo
      INTEGER, ALLOCATABLE::lo1l(:, :)
      !??
      INTEGER, ALLOCATABLE::ulo_der(:, :)
      !no of LOs per l (max(nlo1),ntype
      INTEGER, ALLOCATABLE::nlol(:, :)
      !true if LO is formed by \dot u (
      LOGICAL, ALLOCATABLE::l_dulo(:, :)
      !no of op that maps atom into
      INTEGER, ALLOCATABLE::ngopr(:)
      !symetry of atoms(nat)
      INTEGER, ALLOCATABLE::ntypsy(:)
      !no of sphhar for atom type(ntype
      INTEGER, ALLOCATABLE ::nlhtyp(:)
      !atom mapped to by inversion (nat
      INTEGER, ALLOCATABLE ::invsat(:)
      !Calaculate forces for this atom?
      LOGICAL, ALLOCATABLE :: l_geo(:)
      !MT-Radius (ntype)
      REAL, ALLOCATABLE CPP_MANAGED::rmt(:)
      !log increment(ntype)
      REAL, ALLOCATABLE::dx(:)
      !vol of MT(ntype)
      REAL, ALLOCATABLE::volmts(:)
      !radial grid points(max(jri),ntyp
      REAL, ALLOCATABLE::rmsh(:, :)
      !charge of nucleus(ntype)
      REAL, ALLOCATABLE::zatom(:)
      !initial mag moment(ntype)
      REAL, ALLOCATABLE::bmu(:)
      !pos of atom (absol) (3,nat)
      REAL, ALLOCATABLE::pos(:, :)
      !pos of atom (relat)(3,nat)
      REAL, ALLOCATABLE CPP_MANAGED::taual(:, :)
      !labels
      CHARACTER(LEN=20), ALLOCATABLE :: label(:)
      CHARACTER(len=20), ALLOCATABLE :: speciesName(:)
      !name and other data of explicitely provided xc functional
      CHARACTER(len=4), ALLOCATABLE :: namex(:)
      INTEGER, ALLOCATABLE :: icorr(:)
      INTEGER, ALLOCATABLE :: igrd(:)
      INTEGER, ALLOCATABLE :: krla(:)
      LOGICAL, ALLOCATABLE :: relcor(:)
      !lda_u information(ntype)
      TYPE(t_utype), ALLOCATABLE::lda_u(:)
      INTEGER, ALLOCATABLE :: relax(:, :) !<(3,ntype)
      INTEGER, ALLOCATABLE :: nflip(:) !<flip magnetisation of this atom
   CONTAINS
      procedure :: nsp => calc_nsp_atom
      PROCEDURE :: same_species
      PROCEDURE :: read_xml => read_xml_atoms
   END TYPE t_atoms

   public :: t_atoms

 contains
   LOGICAL function same_species(atoms,n,nn)
     use m_judft
     implicit none
     class(t_atoms),INTENT(IN)::atoms
     integer,intent(in)::n,nn

     if (n>atoms%ntype.or.nn>atoms%ntype) call judft_error("Same species checked for non-existing atom")

     same_species=atoms%nz(n)==atoms%nz(nn)
     same_species=same_species.and.atoms%jri(n)==atoms%jri(nn)
     same_species=same_species.and.atoms%dx(n)==atoms%dx(nn)
     same_species=same_species.and.atoms%rmt(n)==atoms%rmt(nn)
     same_species=same_species.and.atoms%lmax(n)==atoms%lmax(nn)
     same_species=same_species.and.atoms%lnonsph(n)==atoms%lnonsph(nn)
     same_species=same_species.and.atoms%nlo(n)==atoms%nlo(nn)
     if (atoms%nlo(n)==atoms%nlo(nn)) same_species=same_species.and.all(atoms%llo(:,n)==atoms%llo(:,nn))
     same_species=same_species.and.atoms%lapw_l(n)==atoms%lapw_l(nn)
     same_species=same_species.and.atoms%l_geo(n)==atoms%l_geo(nn)
     same_species=same_species.and.trim(atoms%econf(n)%coreconfig)==trim(atoms%econf(nn)%coreconfig)
     same_species=same_species.and.trim(atoms%econf(n)%valenceconfig)==trim(atoms%econf(nn)%valenceconfig)
     same_species=same_species.and.trim(atoms%econf(n)%valenceconfig)==trim(atoms%econf(nn)%valenceconfig)
  end function
  pure function calc_nsp_atom(self) result(nsp)
    implicit none
    CLASS(t_atoms),INTENT(IN)      :: self
    INTEGER                        :: nsp
      
    nsp = (self%lmaxd+1+MOD(self%lmaxd+1,2))*(2*self%lmaxd+1)
   end function

   SUBROUTINE read_xml_atoms(atoms,xml)
    USE m_types_xml
    IMPLICIT NONE
    CLASS(t_atoms),INTENT(OUT):: atoms
    TYPE(t_xml),INTENT(IN)    :: xml

    CHARACTER(len=200):: xpaths,xpathg,xpath,valueString,lstring,nstring,core,valence
    INTEGER           :: i,j,numberNodes,ilo,lNumCount,nNumCount,l,n
    INTEGER,ALLOCATABLE::lNumbers(:),nNumbers(:)
    LOGICAL           :: relaxx,relaxy,relaxz
    INTEGER,ALLOCATABLE :: itmp(:,:)
    REAL                :: down,up
    CHARACTER(len=20)   :: state
    atoms%ntype= xml%get_ntype()
    atoms%nat =  xml%get_nat()
    ALLOCATE(atoms%nz(atoms%ntype))     !nz and zatom have the same content!
    ALLOCATE(atoms%zatom(atoms%ntype))  !nz and zatom have the same content!
    ALLOCATE(atoms%jri(atoms%ntype))
    ALLOCATE(atoms%dx(atoms%ntype))
    ALLOCATE(atoms%lmax(atoms%ntype))
    ALLOCATE(atoms%nlo(atoms%ntype))
    ALLOCATE(atoms%lnonsph(atoms%ntype))
    ALLOCATE(atoms%nflip(atoms%ntype))
    ALLOCATE(atoms%l_geo(atoms%ntype))
    ALLOCATE(atoms%lda_u(4*atoms%ntype))
    ALLOCATE(atoms%bmu(atoms%ntype))
    ALLOCATE(atoms%relax(3,atoms%ntype))
    ALLOCATE(atoms%neq(atoms%ntype))
    ALLOCATE(atoms%taual(3,atoms%nat))
    ALLOCATE(atoms%label(atoms%nat))
    ALLOCATE(atoms%pos(3,atoms%nat))
    ALLOCATE(atoms%rmt(atoms%ntype))
    ALLOCATE(atoms%econf(atoms%ntype))
    ALLOCATE(atoms%ncv(atoms%ntype)) ! For what is this?
    ALLOCATE(atoms%ngopr(atoms%nat)) ! For what is this?
    ALLOCATE(atoms%lapw_l(atoms%ntype)) ! Where do I put this?
    ALLOCATE(atoms%invsat(atoms%nat)) ! Where do I put this?
    ALLOCATE(atoms%llo(MAXVAL(xml%get_nlo()),atoms%ntype))
    atoms%lapw_l(:) = -1
    atoms%n_u = 0
    DO n = 1, atoms%ntype
       !in Species:
       !@name,element,atomicNumber,coreStates
       !@optional: magMom,flipSpin,magField,vcaAddCharge
       !mtSphere,atomicCutoffs
       !optional: energyParameters,prodBasis,special,force,electronConfig,nocoParams,ldaU(up to 4),lo(as many as needed)
       xpathg=xml%groupPath(n)
       xpaths=xml%speciesPath(n)
       atoms%nz(n)=evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPaths))//'/@atomicNumber'))
       IF (atoms%nz(n).EQ.0) THEN
          WRITE(*,*) 'Note: Replacing atomic number 0 by 1.0e-10 on atom type ', n
          atoms%zatom(n) = 1.0e-10
       END IF
       atoms%zatom(n) = atoms%nz(n)
       
       IF (evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPaths))//'/@flipSpin'))) THEN
          atoms%nflip(n) = 1
       ELSE
          atoms%nflip(n) = 0
       ENDIF
       atoms%bmu(n) = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPaths))//'/@magMom'))
       !Now the xml elements
       !mtSphere
       xpath=xpaths
       atoms%rmt(n) =  evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/mtSphere/@radius'))
       atoms%jri(n) = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/mtSphere/@gridPoints'))
       atoms%dx(n) = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/mtSphere/@logIncrement'))
       !atomicCuttoffs
       xpath=xpaths
       atoms%lmax(n) = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/atomicCutoffs/@lmax'))
       atoms%lnonsph(n) = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/atomicCutoffs/@lnonsphr'))
       atoms%lapw_l(n) = -1
       IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPath))//'/atomicCutoffs/@lmaxAPW').EQ.1) &
            atoms%lapw_l(n) = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/atomicCutoffs/@lmaxAPW'))
       !force type
       xpath=''
       IF(xml%getNumberOfNodes(TRIM(ADJUSTL(xPaths))//'/force')==1) xpath=xpaths
       IF (xpath.NE.'') THEN
          atoms%l_geo(n) = evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathg))//'/force/@calculate'))
          valueString = xml%getAttributeValue(TRIM(ADJUSTL(xPathg))//'force/@relaxXYZ')
          READ(valueString,'(3l1)') relaxX, relaxY, relaxZ
          IF (relaxX) atoms%relax(1,n) = 1
          IF (relaxY) atoms%relax(2,n) = 1
          IF (relaxZ) atoms%relax(3,n) = 1
       ELSE
          atoms%l_geo(n) = .FALSE.
          atoms%relax(:,n) = 0
       END IF
       !LO stuff  
       atoms%nlo(n) = 0
        DO ilo = 1,xml%getNumberOfNodes(TRIM(ADJUSTL(xpaths))//'/lo')+xml%getNumberOfNodes(TRIM(ADJUSTL(xpathg))//'/lo')
           IF (ilo>xml%getNumberOfNodes(TRIM(ADJUSTL(xpaths))//'/lo')) THEN
              WRITE(xpath,*) TRIM(ADJUSTL(xpathg))//'/lo[',ilo-xml%getNumberOfNodes(TRIM(ADJUSTL(xpaths))//'/lo'),']'
           ELSE
              WRITE(xpath,*) TRIM(ADJUSTL(xpaths))//'/lo[',ilo,']'
           END IF
          lString = xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@l')
          nString = xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@n')
          CALL getIntegerSequenceFromString(TRIM(ADJUSTL(lString)), lNumbers, lNumCount)
          CALL getIntegerSequenceFromString(TRIM(ADJUSTL(nString)), nNumbers, nNumCount)
          IF(lNumCount.NE.nNumCount) THEN
             CALL judft_error('Error in LO input: l quantum number count does not equal n quantum number count')
          END IF
          IF (.NOT. ALLOCATED(atoms%llo)) ALLOCATE(atoms%llo(1,atoms%ntype),atoms%ulo_der(1,atoms%ntype))
          DO i=1,lNumCount
             atoms%llo(atoms%nlo(n)+i,n) = lNumbers(i)
             atoms%ulo_der(atoms%nlo(n)+i,n) = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@eDeriv'))
          ENDDO
          atoms%nlo(n) = atoms%nlo(n) +lnumcount 
          DEALLOCATE (lNumbers, nNumbers)
       END DO
       !LDA+U
       DO i = 1,xml%getNumberOfNodes(TRIM(ADJUSTL(xPaths))//'/ldaU')
          WRITE(xpath,*) TRIM(ADJUSTL(xPaths))//'/ldaU[',i,']'
          IF (i.GT.4) CALL juDFT_error("Too many U parameters provided for a certain species (maximum is 4).",calledby ="types_atoms")
          atoms%n_u = atoms%n_u + 1
          atoms%lda_u(atoms%n_u)%l = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@l'))
          
          atoms%lda_u(atoms%n_u)%u =  evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@U'))
          atoms%lda_u(atoms%n_u)%j = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@J'))
          atoms%lda_u(atoms%n_u)%l_amf =  evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@l_amf'))
          atoms%lda_u(atoms%n_u)%atomType = n
       END DO
       !electron config
       IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xml%getAttributeValue(xPaths)))//'/electronConfig')==1) THEN
          core=xml%getAttributeValue(TRIM(ADJUSTL(xml%getAttributeValue(xPaths)))//'/electronConfig/coreConfig')
          valence=xml%getAttributeValue(TRIM(ADJUSTL(xml%getAttributeValue(xPaths)))//'/electronConfig/valenceConfig')
          CALL atoms%econf(n)%init(core,valence)
          numberNodes = xml%getNumberOfNodes(TRIM(ADJUSTL(xml%getAttributeValue(xPaths)))//'/electronConfig/stateOccupation')
          IF (numberNodes.GE.1) THEN
             DO i = 1, numberNodes
                WRITE(xpath,"(a,a,i0,a)") TRIM(ADJUSTL(xPath)),'/electronConfig/stateOccupation[',i,']'
                state=xml%getAttributeValue(TRIM(xpath)//'/@state')
                up=evaluateFirstOnly(xml%getAttributeValue(TRIM(xpath)//'/@spinUp'))
                down=evaluateFirstOnly(xml%getAttributeValue(TRIM(xpath)//'/@spinDown'))
                CALL atoms%econf(n)%set_occupation(state,up,down)
             END DO
          END IF
       END IF
     
    END DO

    atoms%nlotot = 0
    DO n = 1, atoms%ntype
       DO l = 1,atoms%nlo(n)
          atoms%nlotot = atoms%nlotot + atoms%neq(n) * ( 2*atoms%llo(l,n) + 1 )
       ENDDO
    ENDDO
    
    ! Check the LO stuff and call setlomap (from inped):
    
    ALLOCATE(atoms%lo1l(0:atoms%llod,atoms%ntype))
    ALLOCATE(atoms%nlol(0:atoms%llod,atoms%ntype))
    
    
    DO n = 1, atoms%ntype
       IF (atoms%nlo(n).GE.1) THEN
          IF (atoms%nlo(n).GT.atoms%nlod) THEN
             WRITE (6,*) 'nlo(n) =',atoms%nlo(n),' > nlod =',atoms%nlod
             CALL juDFT_error("nlo(n)>nlod",calledby ="postprocessInput")
          END IF
          DO j=1,atoms%nlo(n)
             IF ( (atoms%llo(j,n).GT.atoms%llod).OR.(MOD(-atoms%llod,10)-1).GT.atoms%llod ) THEN
                WRITE (6,*) 'llo(j,n) =',atoms%llo(j,n),' > llod =',atoms%llod
                CALL juDFT_error("llo(j,n)>llod",calledby ="postprocessInput")
             END IF
          END DO
          
          ! Replace call to setlomap with the following 3 loops (preliminary).
          ! atoms%nlol and atoms%lo1l arrays are strange. This should be solved differently.
          DO l = 0,atoms%llod
             atoms%nlol(l,n) = 0
             atoms%lo1l(l,n) = 0
          END DO
          
          DO ilo = 1,atoms%nlod
             atoms%l_dulo(ilo,n) = .FALSE.
          END DO
          
          DO ilo = 1,atoms%nlo(n)
             WRITE(6,'(A,I2,A,I2)') 'I use',atoms%ulo_der(ilo,n),'. derivative of l =',atoms%llo(ilo,n)
             IF (atoms%llo(ilo,n)>atoms%llod) CALL juDFT_error(" l > llod!!!",calledby="postprocessInput")
             l = atoms%llo(ilo,n)
             IF (ilo.EQ.1) THEN
                atoms%lo1l(l,n) = ilo
             ELSE
                IF (l.NE.atoms%llo(ilo-1,n)) THEN
                   atoms%lo1l(l,n) = ilo
                END IF
             END IF
             atoms%nlol(l,n) = atoms%nlol(l,n) + 1
          END DO
          WRITE (6,*) 'atoms%lapw_l(n) = ',atoms%lapw_l(n)
       END IF
       
    END DO
    
    ! Check lda+u stuff (from inped)
    
    DO i = 1, atoms%n_u
       n = atoms%lda_u(i)%atomType
       IF (atoms%nlo(n).GE.1) THEN
          DO j = 1, atoms%nlo(n)
             IF ((ABS(atoms%llo(j,n)).EQ.atoms%lda_u(i)%l) .AND. (.NOT.atoms%l_dulo(j,n)) ) &
                  WRITE (*,*) 'LO and LDA+U for same l not implemented'
          END DO
       END IF
    END DO
    
    
  END SUBROUTINE read_xml_atoms
  

 END MODULE m_types_atoms
