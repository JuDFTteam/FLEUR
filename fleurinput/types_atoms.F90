!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
#ifndef CPP_MANAGED
#define CPP_MANAGED
#endif
MODULE m_types_atoms
  USE m_juDFT
  USE m_types_econfig
  USE m_types_fleurinput_base
  IMPLICIT NONE
  PRIVATE

  TYPE t_utype
      SEQUENCE
      REAL :: u=0.0, j=0.0         ! the actual U and J parameters
      REAL :: theta=0.0,phi=0.0   !the rotation angles by which the density metrics is rotated
      INTEGER :: l=-1        ! the l quantum number to which this U parameter belongs
      INTEGER :: atomType=0 ! The atom type to which this U parameter belongs
      LOGICAL :: l_amf=.FALSE. ! logical switch to choose the "around mean field" LDA+U limit
  END TYPE t_utype
  TYPE t_opctype
      SEQUENCE
      INTEGER :: l=-1,n=-1
      INTEGER :: atomType=0 ! The atom type to which this U parameter belongs
  END TYPE t_opctype
  TYPE,EXTENDS(t_fleurinput_base):: t_atoms
     !<no of types
  INTEGER :: ntype=-1
  !<total-no of atoms
  INTEGER :: nat=-1
  !<dimensions of LO's
  INTEGER ::nlod=0
  INTEGER ::llod=0
  INTEGER ::nlotot=0
  !lmaxd=maxval(lmax)
  INTEGER:: lmaxd=-1
  ! no of density matrices to calculate
  INTEGER :: n_denmat=0
  ! no of lda+us
  INTEGER ::n_u=0
  ! no of lda+hubbard1s
  INTEGER :: n_hia=0
  ! no of lda+orbital polarization corrections
  INTEGER :: n_opc=0
  ! dimensions
  INTEGER :: jmtd=-1
  INTEGER :: msh=0 !core state mesh was in dimension
  !No of element
  INTEGER, ALLOCATABLE ::nz(:)
  !atoms per type
  INTEGER, ALLOCATABLE::neq(:)
  ! type of each atom itype(atoms%nat) used for OMP unrolling
  INTEGER, ALLOCATABLE::itype(:)
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
  !no of sphhar for atom type(ntype
  INTEGER, ALLOCATABLE ::nlhtyp(:)
  !Calculate forces for this atom?
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
  !lda_u information(4*ntype)
  !lda+hubbard1 information is attached behind lda+u
  !so the dimension actually used is atoms%n_u+atoms%n_hia
  TYPE(t_utype), ALLOCATABLE::lda_u(:)
  TYPE(t_opctype), ALLOCATABLE::lda_opc(:)
  INTEGER, ALLOCATABLE :: relax(:, :) !<(3,ntype)
  !flipSpinTheta and flipSpinPhi are the angles which are given
  !in the input to rotate the charge den by these polar angles.
  !Typical one needs ntype angles.
  REAL, ALLOCATABLE :: flipSpinPhi(:)
  REAL, ALLOCATABLE :: flipSpinTheta(:)
  !Logical switch which decides if the rotated cdn should be scaled.
  !Yet untested feature.
  LOGICAL, ALLOCATABLE :: flipSpinScale(:)
  !Switches for the output for the calculation of crystal field coefficients
  LOGICAL, ALLOCATABLE :: l_outputCFcdn(:)
  LOGICAL, ALLOCATABLE :: l_outputCFpot(:)
  LOGICAL, ALLOCATABLE :: l_outputCFremove4f(:)
  !special
  LOGICAL, ALLOCATABLE :: lda_atom(:)
CONTAINS
  PROCEDURE :: init=>init_atoms
  PROCEDURE :: nsp => calc_nsp_atom
  PROCEDURE :: same_species
  PROCEDURE :: read_xml => read_xml_atoms
  PROCEDURE :: mpi_bc=>mpi_bc_atoms
END TYPE t_atoms

PUBLIC :: t_atoms,t_utype, readAtomAttribute

CONTAINS
SUBROUTINE mpi_bc_atoms(this,mpi_comm,irank)
 USE m_mpi_bc_tool
 CLASS(t_atoms),INTENT(INOUT)::this
 INTEGER,INTENT(IN):: mpi_comm
 INTEGER,INTENT(IN),OPTIONAL::irank
 INTEGER ::rank,myrank,ierr,n
 IF (PRESENT(irank)) THEN
    rank=irank
 ELSE
    rank=0
 END IF
 CALL mpi_bc(this%ntype,rank,mpi_comm)
 CALL mpi_bc(this%nat,rank,mpi_comm)
 CALL mpi_bc(this%nlod,rank,mpi_comm)
 CALL mpi_bc(this%llod,rank,mpi_comm)
 CALL mpi_bc(this%nlotot,rank,mpi_comm)
 CALL mpi_bc(this%lmaxd,rank,mpi_comm)
 CALL mpi_bc(this%n_u,rank,mpi_comm)
 CALL mpi_bc(this%n_hia,rank,mpi_comm)
 CALL mpi_bc(this%n_opc, rank, mpi_comm)
 CALL mpi_bc(this%n_denmat, rank, mpi_comm)
 CALL mpi_bc(this%jmtd,rank,mpi_comm)
 CALL mpi_bc(this%msh,rank,mpi_comm)
 CALL mpi_bc(this%nz,rank,mpi_comm)
 CALL mpi_bc(this%neq,rank,mpi_comm)
 CALL mpi_bc(this%jri,rank,mpi_comm)
 CALL mpi_bc(this%lmax,rank,mpi_comm)
 CALL mpi_bc(this%lnonsph,rank,mpi_comm)
 CALL mpi_bc(this%ncv,rank,mpi_comm)
 CALL mpi_bc(this%nlo,rank,mpi_comm)
 CALL mpi_bc(this%llo,rank,mpi_comm)
 CALL mpi_bc(this%lapw_l,rank,mpi_comm)
 CALL mpi_bc(this%lo1l,rank,mpi_comm)
 CALL mpi_bc(this%ulo_der,rank,mpi_comm)
 CALL mpi_bc(this%nlol,rank,mpi_comm)
 CALL mpi_bc(this%l_dulo,rank,mpi_comm)
 CALL mpi_bc(this%nlhtyp,rank,mpi_comm)
 CALL mpi_bc(this%l_geo,rank,mpi_comm)
 CALL mpi_bc(this%rmt,rank,mpi_comm)
 CALL mpi_bc(this%dx,rank,mpi_comm)
 CALL mpi_bc(this%volmts,rank,mpi_comm)
 CALL mpi_bc(this%rmsh,rank,mpi_comm)
 CALL mpi_bc(this%zatom,rank,mpi_comm)
 CALL mpi_bc(this%bmu,rank,mpi_comm)
 CALL mpi_bc(this%pos,rank,mpi_comm)
 CALL mpi_bc(this%taual,rank,mpi_comm)
 CALL mpi_bc(this%label,rank,mpi_comm)
 CALL mpi_bc(this%relax,rank,mpi_comm)
 CALL mpi_bc(this%flipSpinPhi,rank,mpi_comm)
 CALL mpi_bc(this%flipSpinTheta,rank,mpi_comm)
 CALL mpi_bc(this%flipSpinScale,rank,mpi_comm)
 CALL mpi_bc(this%l_outputCFcdn,rank,mpi_comm)
 CALL mpi_bc(this%l_outputCFpot,rank,mpi_comm)
 CALL mpi_bc(this%l_outputCFremove4f,rank,mpi_comm)
 call mpi_bc(this%itype,rank,mpi_comm)
 CALL mpi_bc(this%lda_atom, rank, mpi_comm)

#ifdef CPP_MPI
 CALL mpi_COMM_RANK(mpi_comm,myrank,ierr)
 IF (myrank.NE.rank) THEN
    IF (ALLOCATED(this%econf)) DEALLOCATE(this%econf)
    IF (ALLOCATED(this%lda_u)) DEALLOCATE(this%lda_u)
    IF (ALLOCATED(this%lda_opc)) DEALLOCATE(this%lda_opc)
    ALLOCATE(this%econf(this%ntype))
    ALLOCATE(this%lda_u(4*this%ntype))
    ALLOCATE(this%lda_opc(4*this%ntype))
 ENDIF
 DO n=1,this%ntype
    CALL this%econf(n)%broadcast(rank,mpi_comm)
 ENDDO
 DO n=1,this%n_u+this%n_hia
    CALL mpi_bc(this%lda_u(n)%j,rank,mpi_comm)
    CALL mpi_bc(this%lda_u(n)%u,rank,mpi_comm)
    CALL mpi_bc(this%lda_u(n)%theta,rank,mpi_comm)
    CALL mpi_bc(this%lda_u(n)%phi,rank,mpi_comm)
    CALL mpi_bc(this%lda_u(n)%l,rank,mpi_comm)
    CALL mpi_bc(this%lda_u(n)%atomType,rank,mpi_comm)
    CALL mpi_bc(this%lda_u(n)%l_amf,rank,mpi_comm)
 ENDDO
 DO n=1,this%n_opc
   CALL mpi_bc(this%lda_opc(n)%n,rank,mpi_comm)
   CALL mpi_bc(this%lda_opc(n)%l,rank,mpi_comm)
   CALL mpi_bc(this%lda_opc(n)%atomType,rank,mpi_comm)
ENDDO
#endif
END SUBROUTINE mpi_bc_atoms

LOGICAL FUNCTION same_species(atoms,n,nn)
 USE m_judft
 IMPLICIT NONE
 CLASS(t_atoms),INTENT(IN)::atoms
 INTEGER,INTENT(in)::n,nn
 LOGICAL same
 IF (n>atoms%ntype.OR.nn>atoms%ntype) CALL judft_error("Same species checked for non-existing atom")

 same=atoms%nz(n)==atoms%nz(nn)
 same=same.AND.atoms%jri(n)==atoms%jri(nn)
 same=same.AND.atoms%dx(n)==atoms%dx(nn)
 same=same.AND.atoms%rmt(n)==atoms%rmt(nn)
 same=same.AND.atoms%lmax(n)==atoms%lmax(nn)
 same=same.AND.atoms%lnonsph(n)==atoms%lnonsph(nn)
 same=same.AND.atoms%nlo(n)==atoms%nlo(nn)
 IF (atoms%nlo(n)==atoms%nlo(nn)) same=same.AND.ALL(atoms%llo(:,n)==atoms%llo(:,nn))
 same=same.AND.atoms%lapw_l(n)==atoms%lapw_l(nn)
 same=same.AND.(atoms%l_geo(n).EQV.atoms%l_geo(nn))
 same=same.AND.TRIM(atoms%econf(n)%coreconfig)==TRIM(atoms%econf(nn)%coreconfig)
 same=same.AND.TRIM(atoms%econf(n)%valenceconfig)==TRIM(atoms%econf(nn)%valenceconfig)
 same=same.AND.TRIM(atoms%econf(n)%valenceconfig)==TRIM(atoms%econf(nn)%valenceconfig)
 if (same) same=same.AND.ALL(abs(atoms%econf(n)%occupation-atoms%econf(nn)%occupation)<1E-8)
 same_species=same
END FUNCTION same_species
PURE FUNCTION calc_nsp_atom(self) RESULT(nsp)
 IMPLICIT NONE
 CLASS(t_atoms),INTENT(IN)      :: self
 INTEGER                        :: nsp

 nsp = (self%lmaxd+1+MOD(self%lmaxd+1,2))*(2*self%lmaxd+1)
END FUNCTION calc_nsp_atom

SUBROUTINE read_xml_atoms(this,xml)
 USE m_types_xml
 USE m_constants
 IMPLICIT NONE
 CLASS(t_atoms),INTENT(INOUT):: this
 TYPE(t_xml),INTENT(INOUT)    :: xml

 CHARACTER(len=200):: xpaths,xpathg,xpath,valueString,lstring,nstring,core,valence
 INTEGER           :: i,j,numberNodes,ilo,lNumCount,nNumCount,l,n,itype,na,jrc,numU
 INTEGER,ALLOCATABLE::lNumbers(:),nNumbers(:)
 LOGICAL           :: relaxx,relaxy,relaxz,l_flipElectronConfigSpins
 INTEGER,ALLOCATABLE :: itmp(:,:)
 REAL                :: down,up,dr,radius,vca_charge
 CHARACTER(len=20)   :: state
 this%ntype= xml%get_ntype()
 this%nat =  xml%get_nat()
 ALLOCATE(this%nz(this%ntype))     !nz and zatom have the same content!
 ALLOCATE(this%zatom(this%ntype))  !nz and zatom have the same content!
 ALLOCATE(this%jri(this%ntype))
 ALLOCATE(this%dx(this%ntype))
 ALLOCATE(this%lmax(this%ntype))
 ALLOCATE(this%nlo(this%ntype))
 ALLOCATE(this%lnonsph(this%ntype))
 ALLOCATE(this%flipSpinPhi(this%ntype))
 ALLOCATE(this%flipSpinTheta(this%ntype))
 ALLOCATE(this%flipSpinScale(this%ntype))
 ALLOCATE(this%l_geo(this%ntype))
 ALLOCATE(this%lda_u(4*this%ntype))
 ALLOCATE(this%lda_opc(4*this%ntype))
 ALLOCATE(this%bmu(this%ntype))
 ALLOCATE(this%relax(3,this%ntype))
 ALLOCATE(this%neq(this%ntype));this%neq=0
 ALLOCATE(this%taual(3,this%nat))
 ALLOCATE(this%label(this%nat))
 ALLOCATE(this%pos(3,this%nat))
 ALLOCATE(this%rmt(this%ntype))
 ALLOCATE(this%econf(this%ntype))
 ALLOCATE(this%ncv(this%ntype)) ! For what is this?
 ALLOCATE(this%lapw_l(this%ntype)) ! Where do I put this?
 ALLOCATE(this%l_outputCFcdn(this%ntype),source=.FALSE.)
 ALLOCATE(this%l_outputCFpot(this%ntype),source=.FALSE.)
 ALLOCATE(this%l_outputCFremove4f(this%ntype),source=.FALSE.)
 this%nlod=MAXVAL(xml%get_nlo())
 ALLOCATE(this%llo(this%nlod,this%ntype))
 this%llo=0
 ALLOCATE(this%ulo_der(this%nlod,this%ntype))
 ALLOCATE(this%speciesname(this%ntype))
 ALLOCATE(this%lda_atom(this%ntype))
 this%lda_atom = .FALSE.
 this%lapw_l(:) = -1
 this%n_u = 0
 this%bmu = 0.0
 na=0
 this%flipSpinPhi = 0.0
 this%flipSpinTheta = 0.0
 this%flipSpinScale = .FALSE.
 DO n = 1, this%ntype
    !in Species:
    !@name,element,atomicNumber,coreStates
    !@optional: magMom,flipSpin,magField,vcaAddCharge
    !mtSphere,atomicCutoffs
    !optional: energyParameters,prodBasis,special,force,electronConfig,nocoParams,ldaU(up to 4),lo(as many as needed)
    xpathg=xml%groupPath(n)
    xpaths=xml%speciesPath(n)
    this%speciesname(n)=TRIM(ADJUSTL(xml%getAttributeValue(TRIM(ADJUSTL(xPathg))//'/@species')))
    this%nz(n)=evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPaths))//'/@atomicNumber'))
    vca_charge=0.0;call readAtomAttribute(xml,n,'/special/@vca_charge',vca_charge)
    this%zatom(n) = this%nz(n)+vca_charge
    IF (this%nz(n).EQ.0) THEN
       WRITE(*,*) 'Note: Replacing atomic number 0 by 1.0e-10 on atom type ', n
       this%zatom(n) = 1.0e-10
    END IF

    CALL readAtomAttribute(xml,n,'/modInitDen/@flipSpinPhi',this%flipSpinPhi(n))
    CALL readAtomAttribute(xml,n,'/modInitDen/@flipSpinTheta',this%flipSpinTheta(n))
    CALL readAtomAttribute(xml,n,'/modInitDen/@flipSpinScale',this%flipSpinScale(n))
    CALL readAtomAttribute(xml,n,'/modInitDen/@magMom',this%bmu(n))
    IF (xml%versionNumber<32) THEN
       this%bmu(n) = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPaths))//'/@magMom'))
    END IF

    !Now the xml elements
    !mtSphere
    CALL readAtomAttribute(xml,n,'/mtSphere/@radius',this%rmt(n))
    CALL readAtomAttribute(xml,n,'/mtSphere/@gridPoints',this%jri(n))
    CALL readAtomAttribute(xml,n,'/mtSphere/@logIncrement',this%dx(n))
    !atomicCuttoffs
    CALL readAtomAttribute(xml,n,'/atomicCutoffs/@lmax',this%lmax(n))
    CALL readAtomAttribute(xml,n,'/atomicCutoffs/@lnonsphr',this%lnonsph(n))
    IF (this%lmax(n)<this%lnonsph(n)) THEN
      this%lnonsph(n)=this%lmax(n)
      call judft_warn("lnonsph cannot be larger than lmax")
    ENDIF
    this%lapw_l(n) = -1
    CALL readAtomAttribute(xml,n,'/atomicCutoffs/@lmaxAPW',this%lapw_l(n))
    !special
    CALL readAtomAttribute(xml,n,'/special/@lda',this%lda_atom(n))
    !force type
    xpath=''
    IF(xml%getNumberOfNodes(TRIM(ADJUSTL(xPaths))//'/force')==1) xpath=xpaths
    IF(xml%getNumberOfNodes(TRIM(ADJUSTL(xPathg))//'/force')==1) xpath=xpathg
    this%relax(:,n) = 0
    IF (xpath.NE.'') THEN
       this%l_geo(n) = evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/force/@calculate'))
       valueString = xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/force/@relaxXYZ')
       READ(valueString,'(3l1)') relaxX, relaxY, relaxZ
       IF (relaxX) this%relax(1,n) = 1
       IF (relaxY) this%relax(2,n) = 1
       IF (relaxZ) this%relax(3,n) = 1
    ELSE
       this%l_geo(n) = .FALSE.
    END IF
    !LO stuff
    this%nlo(n) = 0
    DO ilo = 1,xml%getNumberOfNodes(TRIM(ADJUSTL(xpaths))//'/lo')+xml%getNumberOfNodes(TRIM(ADJUSTL(xpathg))//'/lo')
       IF (ilo>xml%getNumberOfNodes(TRIM(ADJUSTL(xpaths))//'/lo')) THEN
          WRITE(xpath,*) TRIM(ADJUSTL(xpathg))//'/lo[',ilo-xml%getNumberOfNodes(TRIM(ADJUSTL(xpaths))//'/lo'),']'
       ELSE
          WRITE(xpath,*) TRIM(ADJUSTL(xpaths))//'/lo[',ilo,']'
       END IF
       lString = xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@l')
       nString = xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@n')
       CALL xml%getIntegerSequenceFromString(TRIM(ADJUSTL(lString)), lNumbers, lNumCount)
       CALL xml%getIntegerSequenceFromString(TRIM(ADJUSTL(nString)), nNumbers, nNumCount)
       IF(lNumCount.NE.nNumCount) THEN
          CALL judft_error('Error in LO input: l quantum number count does not equal n quantum number count')
       END IF
       IF (.NOT. ALLOCATED(this%llo)) ALLOCATE(this%llo(1,this%ntype),this%ulo_der(1,this%ntype))
       DO i=1,lNumCount
          this%llo(this%nlo(n)+i,n) = lNumbers(i)
          this%ulo_der(this%nlo(n)+i,n) = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@eDeriv'))
       ENDDO
       this%nlo(n) = this%nlo(n) +lnumcount
       DEALLOCATE (lNumbers, nNumbers)
    END DO
    !LDA+U
    DO i = 1,xml%getNumberOfNodes(TRIM(ADJUSTL(xPaths))//'/ldaU')
       WRITE(xpath,*) TRIM(ADJUSTL(xPaths))//'/ldaU[',i,']'
       IF (i.GT.4) CALL juDFT_error("Too many U parameters provided for a certain species (maximum is 4).",calledby ="read_xml_atoms")
       this%n_u = this%n_u + 1
       this%lda_u(this%n_u)%l = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@l'))
       this%lda_u(this%n_u)%u =  evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@U'))
       this%lda_u(this%n_u)%j = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@J'))
       this%lda_u(this%n_u)%l_amf =  evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@l_amf'))
       this%lda_u(this%n_u)%phi =  evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@phi'))
       this%lda_u(this%n_u)%theta =  evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@theta'))
       this%lda_u(this%n_u)%atomType = n
    END DO
    DO i = 1,xml%getNumberOfNodes(TRIM(ADJUSTL(xPaths))//'/ldaOPC')
      WRITE(xpath,*) TRIM(ADJUSTL(xPaths))//'/ldaOPC[',i,']'
      IF (i.GT.4) CALL juDFT_error("Too many OP corretions provided for a certain species (maximum is 4).",calledby ="read_xml_atoms")
      this%n_opc = this%n_opc + 1
      this%lda_opc(this%n_opc)%l = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@l'))
      this%lda_opc(this%n_opc)%n =  evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@n'))
      this%lda_opc(this%n_opc)%atomType = n
   END DO
    !electron config
    IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPaths))//'/electronConfig')==1) THEN
       l_flipElectronConfigSpins = .FALSE.
       IF (xml%versionNumber>=34) THEN
           l_flipElectronConfigSpins = evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPaths))//'/electronConfig/@flipSpins'))
       ENDIF
       core=xml%getAttributeValue(TRIM(ADJUSTL(xPaths))//'/electronConfig/coreConfig')
       valence=xml%getAttributeValue(TRIM(ADJUSTL(xPaths))//'/electronConfig/valenceConfig')
       CALL this%econf(n)%init(core,valence)
       numberNodes = xml%getNumberOfNodes(TRIM(ADJUSTL(xPaths))//'/electronConfig/stateOccupation')
       IF (numberNodes.GE.1) THEN
          DO i = 1, numberNodes
             WRITE(xpath,"(a,a,i0,a)") TRIM(ADJUSTL(xPaths)),'/electronConfig/stateOccupation[',i,']'
             state=xml%getAttributeValue(TRIM(xpath)//'/@state')
             up=evaluateFirstOnly(xml%getAttributeValue(TRIM(xpath)//'/@spinUp'))
             down=evaluateFirstOnly(xml%getAttributeValue(TRIM(xpath)//'/@spinDown'))
             IF(.NOT.l_flipElectronConfigSpins) THEN
                CALL this%econf(n)%set_occupation(state,up,down)
             ELSE
                CALL this%econf(n)%set_occupation(state,down,up)
             END IF
          END DO
       END IF
    ELSE IF (xml%versionNumber<32) then
       CALL this%econf(n)%init_num(evaluateFirstIntOnly(xml%getAttributeValue(TRIM(xpaths)//'/@coreStates')),this%nz(n))
       CALL this%econf(n)%set_initial_moment(evaluateFirstOnly(xml%getAttributeValue(TRIM(xpaths)//'/@magMom')))
    END IF
    !crystalField output
    numberNodes = xml%getNumberOfNodes(TRIM(ADJUSTL(xPathg))//'/cFCoeffs')
    IF(numberNodes.EQ.1) THEN
       this%l_outputCFcdn(n) = evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathg))//'/cFCoeffs/@chargeDensity'))
       this%l_outputCFpot(n) = evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathg))//'/cFCoeffs/@potential'))
       this%l_outputCFremove4f(n) = evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathg))//'/cFCoeffs/@remove4f'))
    ENDIF
    ! Read in atom positions
    numberNodes = xml%getNumberOfNodes(TRIM(ADJUSTL(xPathg))//'/relPos')
    DO i = 1, numberNodes
       this%neq(n)=this%neq(n)+1
       na = na + 1
       WRITE(xPath,*) TRIM(ADJUSTL(xPathg)),'/relPos[',i,']'
       IF(xml%getNumberOfNodes(TRIM(ADJUSTL(xPath))//'/@label').NE.0) THEN
          this%label(na) = xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@label')
       ELSE
          WRITE(this%label(na),'(i0)') na
       END IF
       valueString = xml%getAttributeValue(TRIM(ADJUSTL(xPath)))
       this%taual(1,na) = evaluatefirst(valueString)
       this%taual(2,na) = evaluatefirst(valueString)
       this%taual(3,na) = evaluatefirst(valueString)
       this%pos(:,na)=this%taual(:,na) !Use as flag that no rescaling is needed in init
    END DO

    numberNodes = xml%getNumberOfNodes(TRIM(ADJUSTL(xPathg))//'/absPos')
    DO i = 1, numberNodes
       na = na + 1
       STOP 'absPos not yet implemented!'
    END DO

    numberNodes = xml%getNumberOfNodes(TRIM(ADJUSTL(xPathg))//'/filmPos')
    DO i = 1, numberNodes
       this%neq(n)=this%neq(n)+1
       na = na + 1
       WRITE(xPath,*) TRIM(ADJUSTL(xPathg)),'/filmPos[',i,']'
       IF(xml%getNumberOfNodes(TRIM(ADJUSTL(xPath))//'/@label').NE.0) THEN
          this%label(na) = xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@label')
       ELSE
          WRITE(this%label(na),'(i0)') na
       END IF
       valueString = xml%getAttributeValue(TRIM(ADJUSTL(xPath)))
       this%taual(1,na) = evaluatefirst(valueString)
       this%taual(2,na) = evaluatefirst(valueString)
       this%taual(3,na) = evaluatefirst(valueString)
       this%pos(1:2,na)=this%taual(1:2,na) !Use as flag that no rescaling is needed in init
       this%pos(3,na)=this%taual(3,na)+1
    END DO

 END DO

 !Read in DFT+Hubbard1 information (only LDA+U relevant information rest is stored in hub1inp)
 !Stored behind LDA+U in the same array
 DO n = 1, this%ntype
    !To avoid segmentation faults read in the number of LDA+U for this atom again
    xpathS=xml%speciesPath(n)
    numU = xml%getNumberOfNodes(TRIM(ADJUSTL(xPathS))//'/ldaU')
    DO i = 1,xml%getNumberOfNodes(TRIM(ADJUSTL(xPathS))//'/ldaHIA')
       WRITE(xpath,*) TRIM(ADJUSTL(xPathS))//'/ldaHIA[',i,']'
       IF(i+numU.GT.4) CALL juDFT_error("Too many U parameters provided for a certain species (maximum is 4).",calledby ="read_xml_atoms")
       this%n_hia = this%n_hia + 1
       this%lda_u(this%n_u+this%n_hia)%l = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@l'))

       this%lda_u(this%n_u+this%n_hia)%u =  evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@U'))
       this%lda_u(this%n_u+this%n_hia)%j = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@J'))
       this%lda_u(this%n_u+this%n_hia)%l_amf = evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@l_amf'))
       this%lda_u(this%n_u+this%n_hia)%phi   = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@phi'))
       this%lda_u(this%n_u+this%n_hia)%theta = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@theta'))
       this%lda_u(this%n_u+this%n_hia)%atomType = n
    END DO
 ENDDO

 this%nlotot = 0
 DO n = 1, this%ntype
    DO l = 1,this%nlo(n)
       this%nlotot = this%nlotot + this%neq(n) * ( 2*this%llo(l,n) + 1 )
    ENDDO
 ENDDO

 ! Check the LO stuff and call setlomap (from inped):
 IF (SIZE(this%llo,1)>0) this%llod=MAXVAL(this%llo)
 ALLOCATE(this%lo1l(0:this%llod,this%ntype))
 ALLOCATE(this%nlol(0:this%llod,this%ntype))
 ALLOCATE(this%l_dulo(this%nlod,this%ntype))

 DO n = 1, this%ntype
    IF (this%nlo(n).GE.1) THEN
       IF (this%nlo(n).GT.this%nlod) THEN
          WRITE (oUnit,*) 'nlo(n) =',this%nlo(n),' > nlod =',this%nlod
          CALL juDFT_error("nlo(n)>nlod",calledby ="postprocessInput")
       END IF
       DO j=1,this%nlo(n)
          IF ( (this%llo(j,n).GT.this%llod).OR.(MOD(-this%llod,10)-1).GT.this%llod ) THEN
             WRITE (oUnit,*) 'llo(j,n) =',this%llo(j,n),' > llod =',this%llod
             CALL juDFT_error("llo(j,n)>llod",calledby ="postprocessInput")
          END IF
       END DO

       ! Replace call to setlomap with the following 3 loops (preliminary).
       ! this%nlol and this%lo1l arrays are strange. This should be solved differently.
       DO l = 0,this%llod
          this%nlol(l,n) = 0
          this%lo1l(l,n) = 0
       END DO

       DO ilo = 1,this%nlod
          this%l_dulo(ilo,n) = .FALSE.
       END DO

       DO ilo = 1,this%nlo(n)
          WRITE(oUnit,'(A,I2,A,I2)') 'I use',this%ulo_der(ilo,n),'. derivative of l =',this%llo(ilo,n)
          IF (this%llo(ilo,n)>this%llod) CALL juDFT_error(" l > llod!!!",calledby="postprocessInput")
          l = this%llo(ilo,n)
          IF (ilo.EQ.1) THEN
             this%lo1l(l,n) = ilo
          ELSE
             IF (l.NE.this%llo(ilo-1,n)) THEN
                this%lo1l(l,n) = ilo
             END IF
          END IF
          this%nlol(l,n) = this%nlol(l,n) + 1
       END DO
       WRITE (oUnit,*) 'atoms%lapw_l(n) = ',this%lapw_l(n)
    END IF

 END DO

 ! Check lda+u stuff (from inped)

 DO i = 1, this%n_u+this%n_hia
    n = this%lda_u(i)%atomType
    IF (this%nlo(n).GE.1) THEN
       DO j = 1, this%nlo(n)
          IF ((ABS(this%llo(j,n)).EQ.this%lda_u(i)%l) .AND. (.NOT.this%l_dulo(j,n)) ) &
               CALL juDFT_warn("LO and LDA+U for same l not implemented",calledby="read_xml_atoms")
       END DO
    END IF
 END DO
 IF(this%n_u.GE.1) THEN
    DO i = this%n_u+1, this%n_u+this%n_hia
       n = this%lda_u(i)%atomType
       l = this%lda_u(i)%l
       DO j = 1, this%n_u
          IF(this%lda_u(j)%atomType.EQ.n.AND.this%lda_u(j)%l.EQ.l) &
             CALL juDFT_error("LDA+U and LDA+Hubbard1 should not be used on the same orbital",calledby="read_xml_atoms")
       END DO
    END DO
 END IF
 DO i = 1, this%n_opc
   IF (this%lda_opc(i)%l/=2.OR.this%lda_opc(i)%n/=3) THEN
      CALL juDFT_warn("LDA+OP only implemented/tested for 3d orbitals",calledby="read_xml_atoms")
   END IF
END DO


 this%jmtd = MAXVAL(this%jri(:))
 ALLOCATE(this%rmsh(this%jmtd,this%ntype))
 ALLOCATE(this%volmts(this%ntype))
 this%rmsh = 0.0
 na = 0
 DO iType = 1, this%ntype
    ! Calculate mesh for valence states
    radius = this%rmt(iType)*EXP(this%dx(iType)*(1-this%jri(iType)))
    dr = EXP(this%dx(iType))
    DO i = 1, this%jri(iType)
       this%rmsh(i,iType) = radius
       radius = radius*dr
    END DO
    ! Calculate mesh dimension for core states
    radius = this%rmt(iType)
    jrc = this%jri(iType)
    DO WHILE (radius < this%rmt(iType) + 20.0)
       jrc = jrc + 1
       radius = radius*dr
    END DO
    this%msh=MAX(this%msh,jrc)

    this%volmts(iType) = (fpi_const/3.0)*this%rmt(iType)**3
 END DO
 this%nlotot = 0
 DO n = 1, this%ntype
    DO l = 1,this%nlo(n)
       this%nlotot = this%nlotot + this%neq(n) * ( 2*this%llo(l,n) + 1 )
    END DO
 END DO

 this%lmaxd=MAXVAL(this%lmax)

 ALLOCATE(this%nlhtyp(xml%get_ntype()))
END SUBROUTINE read_xml_atoms


SUBROUTINE readAtomAttributeString(xml, atomType, relAttPath, outString, l_error)

   USE m_types_xml

   IMPLICIT NONE

   TYPE(t_xml), INTENT(IN)           :: xml
   INTEGER, INTENT(IN)               :: atomType
   CHARACTER(LEN=*), INTENT(IN)      :: relAttPath
   CHARACTER(LEN=200), INTENT(INOUT) :: outString
   LOGICAL, INTENT(INOUT)            :: l_error

   CHARACTER(LEN=200) :: path, groupPath, speciesPath

   l_error = .FALSE.
   groupPath = xml%groupPath(atomType)
   speciesPath = xml%speciesPath(atomType)

   path = TRIM(groupPath)//TRIM(relAttPath)
   IF (xml%getNumberOfNodes(TRIM(PATH)).NE.1) THEN
      path = TRIM(speciesPath)//TRIM(relAttPath)
   END IF
   IF (xml%getNumberOfNodes(TRIM(PATH)).NE.1) THEN
      l_error = .TRUE.
      RETURN
   END IF

   outString = xml%getAttributeValue(TRIM(path))

END SUBROUTINE readAtomAttributeString


SUBROUTINE readAtomAttribute(xml, atomType, relAttPath, outValue)

   USE m_types_xml
   USE m_juDFT

   IMPLICIT NONE

   TYPE(t_xml), INTENT(IN)           :: xml
   INTEGER, INTENT(IN)               :: atomType
   CHARACTER(LEN=*), INTENT(IN)      :: relAttPath
   CLASS(*), INTENT(INOUT)           :: outValue

   LOGICAL            :: l_error
   CHARACTER(LEN=200) :: valueString

   l_error = .FALSE.

   CALL readAtomAttributeString(xml, atomType, relAttPath, valueString, l_error)
   IF (l_error) RETURN
   SELECT TYPE(outValue)
      TYPE IS(INTEGER)
         outValue = evaluateFirstIntOnly(valueString)
      TYPE IS(REAL)
         outValue = evaluatefirstOnly(valueString)
      TYPE IS(LOGICAL)
         outValue = evaluateFirstBoolOnly(valueString)
      CLASS DEFAULT
         CALL juDFT_error("unknown type passed to readAttribute", calledby = "types_atomes - readAttribute")
   END SELECT

END SUBROUTINE readAtomAttribute


SUBROUTINE init_atoms(this,cell)

   USE m_types_cell

   IMPLICIT NONE

   CLASS(t_atoms),INTENT(inout):: this
   TYPE(t_cell),INTENT(IN)   :: cell
   integer :: it, ineq, ic

   WHERE (ABS(this%pos(3,:)-this%taual(3,:))>0.5) this%taual(3,:) = this%taual(3,:) / cell%amat(3,3)
   this%pos(:,:) = MATMUL(cell%amat,this%taual(:,:))

   allocate(this%itype(this%nat))
   ic=0
   DO it = 1, this%ntype
      DO ineq = 1, this%neq(it)
         ic = ic + 1
         this%itype(ic) = it
      END DO
   END DO

   this%n_denmat = this%n_u + this%n_hia + this%n_opc

END SUBROUTINE init_atoms

END MODULE m_types_atoms
