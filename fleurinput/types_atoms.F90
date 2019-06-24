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
      REAL :: u, j         ! the actual U and J parameters
      REAL :: theta,phi   !the rotation angles by which the density metrics is rotated
      INTEGER :: l        ! the l quantum number to which this U parameter belongs
      INTEGER :: atomType ! The atom type to which this U parameter belongs
      LOGICAL :: l_amf ! logical switch to choose the "around mean field" LDA+U limit
   END TYPE t_utype
   TYPE,EXTENDS(t_fleurinput_base):: t_atoms
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
      PROCEDURE :: nsp => calc_nsp_atom
      PROCEDURE :: same_species
      PROCEDURE :: read_xml => read_xml_atoms
      procedure :: mpi_bc=>mpi_bc_atoms
   END TYPE t_atoms

   PUBLIC :: t_atoms

 CONTAINS
   subroutine mpi_bc_atoms(this,mpi_comm,irank)
     use m_mpi_bc_tool
     class(t_atoms),INTENT(INOUT)::this
     integer,INTENT(IN):: mpi_comm
     INTEGER,INTENT(IN),OPTIONAL::irank
     INTEGER ::rank
     if (present(irank)) THEN
        rank=0
     else
        rank=irank
     end if

     call mpi_bc(this% ntype,rank,mpi_comm)
     call mpi_bc(this% nat,rank,mpi_comm)
     call mpi_bc(this%nlod,rank,mpi_comm)
     call mpi_bc(this%llod,rank,mpi_comm)
     call mpi_bc(this%nlotot,rank,mpi_comm)
     call mpi_bc(this% lmaxd,rank,mpi_comm)
     call mpi_bc(this%n_u,rank,mpi_comm)
     call mpi_bc(this% jmtd,rank,mpi_comm)
     call mpi_bc(this%nz,rank,mpi_comm)
     call mpi_bc(this%neq,rank,mpi_comm)
     call mpi_bc(this%jri,rank,mpi_comm)
     call mpi_bc(this%lmax,rank,mpi_comm)
     call mpi_bc(this%lnonsph,rank,mpi_comm)
     call mpi_bc(this%ncv,rank,mpi_comm)
     call mpi_bc(this%nlo,rank,mpi_comm)
     call mpi_bc(this%llo,rank,mpi_comm)
     call mpi_bc(this%lapw_l,rank,mpi_comm)
     call mpi_bc(this%lo1l,rank,mpi_comm)
     call mpi_bc(this%ulo_der,rank,mpi_comm)
     call mpi_bc(this%nlol,rank,mpi_comm)
     call mpi_bc(this%l_dulo,rank,mpi_comm)
     call mpi_bc(this%ngopr,rank,mpi_comm)
     call mpi_bc(this%ntypsy,rank,mpi_comm)
     call mpi_bc(this%nlhtyp,rank,mpi_comm)
     call mpi_bc(this%invsat,rank,mpi_comm)
     call mpi_bc(this% l_geo,rank,mpi_comm)
     call mpi_bc(this%rmt,rank,mpi_comm)
     call mpi_bc(this%dx,rank,mpi_comm)
     call mpi_bc(this%volmts,rank,mpi_comm)
     call mpi_bc(this%rmsh,rank,mpi_comm)
     call mpi_bc(this%zatom,rank,mpi_comm)
     call mpi_bc(this%bmu,rank,mpi_comm)
     call mpi_bc(this%pos,rank,mpi_comm)
     call mpi_bc(this%taual,rank,mpi_comm)
     call mpi_bc(this% namex,rank,mpi_comm)
     call mpi_bc(this% icorr,rank,mpi_comm)
     call mpi_bc(this% igrd,rank,mpi_comm)
     call mpi_bc(this% krla,rank,mpi_comm)
     call mpi_bc(this% relcor,rank,mpi_comm)
     call mpi_bc(this% relax,rank,mpi_comm)
     call mpi_bc(this% nflip,rank,mpi_comm)



     !this needs work
     !TYPE(t_econfig),ALLOCATABLE::econf(:)
     !TYPE(t_utype), ALLOCATABLE::lda_u(:)
   end subroutine mpi_bc_atoms

   LOGICAL FUNCTION same_species(atoms,n,nn)
     USE m_judft
     IMPLICIT NONE
     CLASS(t_atoms),INTENT(IN)::atoms
     INTEGER,INTENT(in)::n,nn

     IF (n>atoms%ntype.OR.nn>atoms%ntype) CALL judft_error("Same species checked for non-existing atom")

     same_species=atoms%nz(n)==atoms%nz(nn)
     same_species=same_species.AND.atoms%jri(n)==atoms%jri(nn)
     same_species=same_species.AND.atoms%dx(n)==atoms%dx(nn)
     same_species=same_species.AND.atoms%rmt(n)==atoms%rmt(nn)
     same_species=same_species.AND.atoms%lmax(n)==atoms%lmax(nn)
     same_species=same_species.AND.atoms%lnonsph(n)==atoms%lnonsph(nn)
     same_species=same_species.AND.atoms%nlo(n)==atoms%nlo(nn)
     IF (atoms%nlo(n)==atoms%nlo(nn)) same_species=same_species.AND.ALL(atoms%llo(:,n)==atoms%llo(:,nn))
     same_species=same_species.AND.atoms%lapw_l(n)==atoms%lapw_l(nn)
     same_species=same_species.AND.atoms%l_geo(n)==atoms%l_geo(nn)
     same_species=same_species.AND.TRIM(atoms%econf(n)%coreconfig)==TRIM(atoms%econf(nn)%coreconfig)
     same_species=same_species.AND.TRIM(atoms%econf(n)%valenceconfig)==TRIM(atoms%econf(nn)%valenceconfig)
     same_species=same_species.AND.TRIM(atoms%econf(n)%valenceconfig)==TRIM(atoms%econf(nn)%valenceconfig)
  END FUNCTION
  PURE FUNCTION calc_nsp_atom(self) RESULT(nsp)
    IMPLICIT NONE
    CLASS(t_atoms),INTENT(IN)      :: self
    INTEGER                        :: nsp
      
    nsp = (self%lmaxd+1+MOD(self%lmaxd+1,2))*(2*self%lmaxd+1)
   END FUNCTION

   SUBROUTINE read_xml_atoms(this,xml)
    USE m_types_xml
    IMPLICIT NONE
    CLASS(t_atoms),INTENT(OUT):: this
    TYPE(t_xml),INTENT(IN)    :: xml

    CHARACTER(len=200):: xpaths,xpathg,xpath,valueString,lstring,nstring,core,valence
    INTEGER           :: i,j,numberNodes,ilo,lNumCount,nNumCount,l,n
    INTEGER,ALLOCATABLE::lNumbers(:),nNumbers(:)
    LOGICAL           :: relaxx,relaxy,relaxz
    INTEGER,ALLOCATABLE :: itmp(:,:)
    REAL                :: down,up
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
    ALLOCATE(this%nflip(this%ntype))
    ALLOCATE(this%l_geo(this%ntype))
    ALLOCATE(this%lda_u(4*this%ntype))
    ALLOCATE(this%bmu(this%ntype))
    ALLOCATE(this%relax(3,this%ntype))
    ALLOCATE(this%neq(this%ntype))
    ALLOCATE(this%taual(3,this%nat))
    ALLOCATE(this%label(this%nat))
    ALLOCATE(this%pos(3,this%nat))
    ALLOCATE(this%rmt(this%ntype))
    ALLOCATE(this%econf(this%ntype))
    ALLOCATE(this%ncv(this%ntype)) ! For what is this?
    ALLOCATE(this%ngopr(this%nat)) ! For what is this?
    ALLOCATE(this%lapw_l(this%ntype)) ! Where do I put this?
    ALLOCATE(this%invsat(this%nat)) ! Where do I put this?
    ALLOCATE(this%llo(MAXVAL(xml%get_nlo()),this%ntype))
    this%lapw_l(:) = -1
    this%n_u = 0
    DO n = 1, this%ntype
       !in Species:
       !@name,element,atomicNumber,coreStates
       !@optional: magMom,flipSpin,magField,vcaAddCharge
       !mtSphere,atomicCutoffs
       !optional: energyParameters,prodBasis,special,force,electronConfig,nocoParams,ldaU(up to 4),lo(as many as needed)
       xpathg=xml%groupPath(n)
       xpaths=xml%speciesPath(n)
       this%nz(n)=evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPaths))//'/@atomicNumber'))
       IF (this%nz(n).EQ.0) THEN
          WRITE(*,*) 'Note: Replacing atomic number 0 by 1.0e-10 on atom type ', n
          this%zatom(n) = 1.0e-10
       END IF
       this%zatom(n) = this%nz(n)
       
       IF (evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPaths))//'/@flipSpin'))) THEN
          this%nflip(n) = 1
       ELSE
          this%nflip(n) = 0
       ENDIF
       this%bmu(n) = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPaths))//'/@magMom'))
       !Now the xml elements
       !mtSphere
       xpath=xpaths
       this%rmt(n) =  evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/mtSphere/@radius'))
       this%jri(n) = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/mtSphere/@gridPoints'))
       this%dx(n) = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/mtSphere/@logIncrement'))
       !atomicCuttoffs
       xpath=xpaths
       this%lmax(n) = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/atomicCutoffs/@lmax'))
       this%lnonsph(n) = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/atomicCutoffs/@lnonsphr'))
       this%lapw_l(n) = -1
       IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPath))//'/atomicCutoffs/@lmaxAPW').EQ.1) &
            this%lapw_l(n) = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/atomicCutoffs/@lmaxAPW'))
       !force type
       xpath=''
       IF(xml%getNumberOfNodes(TRIM(ADJUSTL(xPaths))//'/force')==1) xpath=xpaths
       IF (xpath.NE.'') THEN
          this%l_geo(n) = evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPathg))//'/force/@calculate'))
          valueString = xml%getAttributeValue(TRIM(ADJUSTL(xPathg))//'force/@relaxXYZ')
          READ(valueString,'(3l1)') relaxX, relaxY, relaxZ
          IF (relaxX) this%relax(1,n) = 1
          IF (relaxY) this%relax(2,n) = 1
          IF (relaxZ) this%relax(3,n) = 1
       ELSE
          this%l_geo(n) = .FALSE.
          this%relax(:,n) = 0
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
          CALL getIntegerSequenceFromString(TRIM(ADJUSTL(lString)), lNumbers, lNumCount)
          CALL getIntegerSequenceFromString(TRIM(ADJUSTL(nString)), nNumbers, nNumCount)
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
          IF (i.GT.4) CALL juDFT_error("Too many U parameters provided for a certain species (maximum is 4).",calledby ="types_this")
          this%n_u = this%n_u + 1
          this%lda_u(this%n_u)%l = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@l'))
          
          this%lda_u(this%n_u)%u =  evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@U'))
          this%lda_u(this%n_u)%j = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@J'))
          this%lda_u(this%n_u)%l_amf =  evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@l_amf'))
          this%lda_u(this%n_u)%atomType = n
       END DO
       !electron config
       IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xml%getAttributeValue(xPaths)))//'/electronConfig')==1) THEN
          core=xml%getAttributeValue(TRIM(ADJUSTL(xml%getAttributeValue(xPaths)))//'/electronConfig/coreConfig')
          valence=xml%getAttributeValue(TRIM(ADJUSTL(xml%getAttributeValue(xPaths)))//'/electronConfig/valenceConfig')
          CALL this%econf(n)%init(core,valence)
          numberNodes = xml%getNumberOfNodes(TRIM(ADJUSTL(xml%getAttributeValue(xPaths)))//'/electronConfig/stateOccupation')
          IF (numberNodes.GE.1) THEN
             DO i = 1, numberNodes
                WRITE(xpath,"(a,a,i0,a)") TRIM(ADJUSTL(xPath)),'/electronConfig/stateOccupation[',i,']'
                state=xml%getAttributeValue(TRIM(xpath)//'/@state')
                up=evaluateFirstOnly(xml%getAttributeValue(TRIM(xpath)//'/@spinUp'))
                down=evaluateFirstOnly(xml%getAttributeValue(TRIM(xpath)//'/@spinDown'))
                CALL this%econf(n)%set_occupation(state,up,down)
             END DO
          END IF
       END IF
     
    END DO

    this%nlotot = 0
    DO n = 1, this%ntype
       DO l = 1,this%nlo(n)
          this%nlotot = this%nlotot + this%neq(n) * ( 2*this%llo(l,n) + 1 )
       ENDDO
    ENDDO
    
    ! Check the LO stuff and call setlomap (from inped):
    
    ALLOCATE(this%lo1l(0:this%llod,this%ntype))
    ALLOCATE(this%nlol(0:this%llod,this%ntype))
    
    
    DO n = 1, this%ntype
       IF (this%nlo(n).GE.1) THEN
          IF (this%nlo(n).GT.this%nlod) THEN
             WRITE (6,*) 'nlo(n) =',this%nlo(n),' > nlod =',this%nlod
             CALL juDFT_error("nlo(n)>nlod",calledby ="postprocessInput")
          END IF
          DO j=1,this%nlo(n)
             IF ( (this%llo(j,n).GT.this%llod).OR.(MOD(-this%llod,10)-1).GT.this%llod ) THEN
                WRITE (6,*) 'llo(j,n) =',this%llo(j,n),' > llod =',this%llod
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
             WRITE(6,'(A,I2,A,I2)') 'I use',this%ulo_der(ilo,n),'. derivative of l =',this%llo(ilo,n)
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
          WRITE (6,*) 'this%lapw_l(n) = ',this%lapw_l(n)
       END IF
       
    END DO
    
    ! Check lda+u stuff (from inped)
    
    DO i = 1, this%n_u
       n = this%lda_u(i)%atomType
       IF (this%nlo(n).GE.1) THEN
          DO j = 1, this%nlo(n)
             IF ((ABS(this%llo(j,n)).EQ.this%lda_u(i)%l) .AND. (.NOT.this%l_dulo(j,n)) ) &
                  WRITE (*,*) 'LO and LDA+U for same l not implemented'
          END DO
       END IF
    END DO
    
    
  END SUBROUTINE read_xml_atoms
  

 END MODULE m_types_atoms
