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

   TYPE t_gfelementtype
      SEQUENCE
      !defines the l and atomType elements for given greens function element (used for mapping index in types_greensf)
      INTEGER l
      INTEGER lp
      INTEGER atomType
      INTEGER atomTypep
   END TYPE t_gfelementtype

   TYPE t_j0calctype
      INTEGER atomType  !atom Type for which to calculate J0
      INTEGER l_min     !Minimum l considered
      INTEGER l_max     !Maximum l considered
      LOGICAL l_avgexc  !Determines wether we average over the exchange splittings for all l
      LOGICAL l_eDependence  !Switch to output J0 with variating fermi energy (only with contourDOS)
   END TYPE

  TYPE t_utype
      SEQUENCE
      REAL :: u=0.0, j=0.0         ! the actual U and J parameters
      REAL :: theta=0.0,phi=0.0   !the rotation angles by which the density metrics is rotated
      INTEGER :: l=-1        ! the l quantum number to which this U parameter belongs
      INTEGER :: atomType=0 ! The atom type to which this U parameter belongs
      LOGICAL :: l_amf=.false. ! logical switch to choose the "around mean field" LDA+U limit
   END TYPE t_utype
   TYPE,EXTENDS(t_fleurinput_base):: t_atoms
      !<no of types
      INTEGER :: ntype
      !<total-no of atoms
      INTEGER :: nat
      !<dimensions of LO's
      INTEGER ::nlod
      INTEGER ::llod=0
      INTEGER ::nlotot
      !lmaxd=maxval(lmax)
      INTEGER:: lmaxd
      ! no of lda+us
      INTEGER ::n_u=0
      INTEGER :: n_hia=0
      INTEGER :: n_j0=0
      ! no of greens function calculations (in total)
      INTEGER :: n_gf=0
      ! dimensions
      INTEGER :: jmtd
      INTEGER :: msh=0 !core state mesh was in dimension
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
      !no of sphhar for atom type(ntype
      INTEGER, ALLOCATABLE ::nlhtyp(:)
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
      !lda_u information(ntype)
      TYPE(t_utype), ALLOCATABLE::lda_u(:)
      TYPE(t_utype),ALLOCATABLE::lda_hia(:)
      !j0 calc information
      TYPE(t_gfelementtype), ALLOCATABLE::gfelem(:)
      TYPE(t_j0calctype), ALLOCATABLE::j0(:)

      INTEGER, ALLOCATABLE :: relax(:, :) !<(3,ntype)
      !flipSpinTheta and flipSpinPhi are the angles which are given
      !in the input to rotate the charge den by these polar angles.
      !Typical one needs ntype angles.
      REAL, ALLOCATABLE :: flipSpinPhi(:)
      REAL, ALLOCATABLE :: flipSpinTheta(:)
      !Logical switch which decides if the rotated cdn should be scaled.
      !Yet untested feature.
      LOGICAL, ALLOCATABLE :: flipSpinScale(:)
   CONTAINS
      PROCEDURE :: init=>init_atoms
      PROCEDURE :: nsp => calc_nsp_atom
      PROCEDURE :: same_species
      PROCEDURE :: read_xml => read_xml_atoms
      procedure :: mpi_bc=>mpi_bc_atoms
   END TYPE t_atoms

   PUBLIC :: t_atoms,t_utype

 CONTAINS
   subroutine mpi_bc_atoms(this,mpi_comm,irank)
     use m_mpi_bc_tool
     class(t_atoms),INTENT(INOUT)::this
     integer,INTENT(IN):: mpi_comm
     INTEGER,INTENT(IN),OPTIONAL::irank
     INTEGER ::rank,myrank,ierr,n
     if (present(irank)) THEN
        rank=irank
     else
        rank=0
     end if
     print *,"Attention, HUB1 parameters not in BC"
     call mpi_bc(this%ntype,rank,mpi_comm)
     call mpi_bc(this%nat,rank,mpi_comm)
     call mpi_bc(this%nlod,rank,mpi_comm)
     call mpi_bc(this%llod,rank,mpi_comm)
     call mpi_bc(this%nlotot,rank,mpi_comm)
     call mpi_bc(this%lmaxd,rank,mpi_comm)
     call mpi_bc(this%n_u,rank,mpi_comm)
     call mpi_bc(this%jmtd,rank,mpi_comm)
     call mpi_bc(this%msh,rank,mpi_comm)
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
     call mpi_bc(this%nlhtyp,rank,mpi_comm)
     call mpi_bc(this%l_geo,rank,mpi_comm)
     call mpi_bc(this%rmt,rank,mpi_comm)
     call mpi_bc(this%dx,rank,mpi_comm)
     call mpi_bc(this%volmts,rank,mpi_comm)
     call mpi_bc(this%rmsh,rank,mpi_comm)
     call mpi_bc(this%zatom,rank,mpi_comm)
     call mpi_bc(this%bmu,rank,mpi_comm)
     call mpi_bc(this%pos,rank,mpi_comm)
     call mpi_bc(this%taual,rank,mpi_comm)
     call mpi_bc(this%relax,rank,mpi_comm)
     call mpi_bc(this%flipSpinPhi,rank,mpi_comm)
     call mpi_bc(this%flipSpinTheta,rank,mpi_comm)
     call mpi_bc(this%flipSpinScale,rank,mpi_comm)

#ifdef CPP_MPI
     CALL mpi_COMM_RANK(mpi_comm,myrank,ierr)
     IF (myrank.ne.rank) Then
       if (allocated(this%econf)) DEALLOCATE(this%econf)
       if (allocated(this%lda_u)) DEALLOCATE(this%lda_u)
       ALLOCATE(this%econf(this%ntype))
       ALLOCATE(this%lda_u(4*this%ntype))
     ENDIF
     DO n=1,this%ntype
       call this%econf(n)%broadcast(rank,mpi_comm)
     endDO
     DO n=1,this%n_u
       call mpi_bc(this%lda_u(n)%j,rank,mpi_comm)
       call mpi_bc(this%lda_u(n)%u,rank,mpi_comm)
       call mpi_bc(this%lda_u(n)%theta,rank,mpi_comm)
       call mpi_bc(this%lda_u(n)%phi,rank,mpi_comm)
       call mpi_bc(this%lda_u(n)%l,rank,mpi_comm)
       call mpi_bc(this%lda_u(n)%atomType,rank,mpi_comm)
       call mpi_bc(this%lda_u(n)%l_amf,rank,mpi_comm)
     ENDDO
#endif
   end subroutine mpi_bc_atoms

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
     same=same.AND.(atoms%l_geo(n).eqv.atoms%l_geo(nn))
     same=same.AND.TRIM(atoms%econf(n)%coreconfig)==TRIM(atoms%econf(nn)%coreconfig)
     same=same.AND.TRIM(atoms%econf(n)%valenceconfig)==TRIM(atoms%econf(nn)%valenceconfig)
     same=same.AND.TRIM(atoms%econf(n)%valenceconfig)==TRIM(atoms%econf(nn)%valenceconfig)
     same_species=same
  END FUNCTION
  PURE FUNCTION calc_nsp_atom(self) RESULT(nsp)
    IMPLICIT NONE
    CLASS(t_atoms),INTENT(IN)      :: self
    INTEGER                        :: nsp

    nsp = (self%lmaxd+1+MOD(self%lmaxd+1,2))*(2*self%lmaxd+1)
   END FUNCTION

   SUBROUTINE read_xml_atoms(this,xml)
    USE m_types_xml
    use m_constants
    IMPLICIT NONE
    CLASS(t_atoms),INTENT(INOUT):: this
    TYPE(t_xml),INTENT(IN)    :: xml

    CHARACTER(len=200):: xpaths,xpathg,xpath,valueString,lstring,nstring,core,valence
    INTEGER           :: i,j,numberNodes,ilo,lNumCount,nNumCount,l,n,itype,na,jrc
    INTEGER,ALLOCATABLE::lNumbers(:),nNumbers(:)
    LOGICAL           :: relaxx,relaxy,relaxz
    INTEGER,ALLOCATABLE :: itmp(:,:)
    REAL                :: down,up,dr,radius
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
    ALLOCATE(this%bmu(this%ntype))
    ALLOCATE(this%relax(3,this%ntype))
    ALLOCATE(this%neq(this%ntype));this%neq=0
    ALLOCATE(this%taual(3,this%nat))
    ALLOCATE(this%label(this%nat))
    ALLOCATE(this%pos(3,this%nat))
    ALLOCATE(this%rmt(this%ntype))
    ALLOCATE(this%j0(this%ntype))
    ALLOCATE(this%gfelem(4*This%ntype))
    ALLOCATE(this%econf(this%ntype))
    ALLOCATE(this%ncv(this%ntype)) ! For what is this?
    ALLOCATE(this%lapw_l(this%ntype)) ! Where do I put this?
    this%nlod=MAXVAL(xml%get_nlo())
    ALLOCATE(this%llo(this%nlod,this%ntype))
    this%llo=0
    ALLOCATE(this%ulo_der(this%nlod,this%ntype))
    ALLOCATE(this%speciesname(this%ntype))
    this%lapw_l(:) = -1
    this%n_u = 0
    na=0
    DO n = 1, this%ntype
       !in Species:
       !@name,element,atomicNumber,coreStates
       !@optional: magMom,flipSpin,magField,vcaAddCharge
       !mtSphere,atomicCutoffs
       !optional: energyParameters,prodBasis,special,force,electronConfig,nocoParams,ldaU(up to 4),lo(as many as needed)
       xpathg=xml%groupPath(n)
       xpaths=xml%speciesPath(n)
       this%speciesname(n)=trim(adjustl(xml%getAttributeValue(TRIM(ADJUSTL(xPathg))//'/@species')))
       this%nz(n)=evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPaths))//'/@atomicNumber'))
       IF (this%nz(n).EQ.0) THEN
          WRITE(*,*) 'Note: Replacing atomic number 0 by 1.0e-10 on atom type ', n
          this%zatom(n) = 1.0e-10
       END IF
       this%zatom(n) = this%nz(n)
       this%flipSpinPhi(n) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPaths))//'/@flipSpinPhi'))
       this%flipSpinTheta(n) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xpaths))//'/@flipSpinTheta'))
       this%flipSpinScale(n) = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xpaths))//'/@flipSpinScale'))

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
       IF(xml%getNumberOfNodes(TRIM(ADJUSTL(xPathg))//'/force')==1) xpath=xpathg
       IF (xpath.NE.'') THEN
          this%l_geo(n) = evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/force/@calculate'))
          valueString = xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/force/@relaxXYZ')
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
          IF (i.GT.4) CALL juDFT_error("Too many U parameters provided for a certain species (maximum is 4).",calledby ="types_this")
          this%n_u = this%n_u + 1
          this%lda_u(this%n_u)%l = evaluateFirstIntOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@l'))

          this%lda_u(this%n_u)%u =  evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@U'))
          this%lda_u(this%n_u)%j = evaluateFirstOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@J'))
          this%lda_u(this%n_u)%l_amf =  evaluateFirstBoolOnly(xml%getAttributeValue(TRIM(ADJUSTL(xPath))//'/@l_amf'))
          this%lda_u(this%n_u)%atomType = n
       END DO
       !electron config
       IF (xml%getNumberOfNodes(TRIM(ADJUSTL(xPaths))//'/electronConfig')==1) THEN
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
                CALL this%econf(n)%set_occupation(state,up,down)
             END DO
          END IF
       END IF
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

    this%jmtd = maxval(this%jri(:))
    ALLOCATE(this%rmsh(this%jmtd,this%ntype))
    ALLOCATE(this%volmts(this%ntype))
    na = 0
    DO iType = 1, this%ntype
       ! Calculate mesh for valence states
       radius = this%rmt(iType)*exp(this%dx(iType)*(1-this%jri(iType)))
       dr = exp(this%dx(iType))
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
       this%msh=max(this%msh,jrc)

       this%volmts(iType) = (fpi_const/3.0)*this%rmt(iType)**3
    END DO
    this%nlotot = 0
    DO n = 1, this%ntype
       DO l = 1,this%nlo(n)
          this%nlotot = this%nlotot + this%neq(n) * ( 2*this%llo(l,n) + 1 )
       END DO
    END DO

    this%lmaxd=maxval(this%lmax)


  END SUBROUTINE read_xml_atoms

  subroutine init_atoms(this,cell)
    use m_types_cell
    class(t_atoms),intent(inout):: this
    type(t_cell),INTENT(IN)   :: cell

    where (abs(this%pos(3,:)-this%taual(3,:))>0.5) this%taual(3,:) = this%taual(3,:) / cell%amat(3,3)
    this%pos(:,:) = matmul(cell%amat,this%taual(:,:))
  end subroutine init_atoms
   SUBROUTINE add_gfjob(nType,lmin,lmax,atoms,l_off,l_inter,l_nn)

         USE m_juDFT

         INTEGER,          INTENT(IN)     :: nType
         INTEGER,          INTENT(IN)     :: lmin
         INTEGER,          INTENT(IN)     :: lmax
         TYPE(t_atoms),    INTENT(INOUT)  :: atoms
         LOGICAL,          INTENT(IN)     :: l_off !l!=lp
         LOGICAL,          INTENT(IN)     :: l_inter
         LOGICAL,          INTENT(IN)     :: l_nn

         INTEGER l,lp,i_gf
         LOGICAL l_found

         IF(l_inter) CALL juDFT_error("Intersite greens function not yet implemented",calledby="add_gfjob")

         !TODO: add the nearest neighbours jobs

         DO l = lmin, lmax
            DO lp = MERGE(lmin,l,l_off), MERGE(lmax,l,l_off)
               !Check if this job has already been added
               l_found = .FALSE.
               DO i_gf = 1, atoms%n_gf
                  IF(atoms%gfelem(i_gf)%l.NE.l) CYCLE
                  IF(atoms%gfelem(i_gf)%lp.NE.lp) CYCLE
                  IF(atoms%gfelem(i_gf)%atomType.NE.nType) CYCLE
                  IF(atoms%gfelem(i_gf)%atomTypep.NE.nType) CYCLE
                  l_found = .TRUE.
               ENDDO
               IF(l_found) CYCLE !This job is already in the array

               atoms%n_gf = atoms%n_gf + 1
               atoms%gfelem(atoms%n_gf)%l = l
               atoms%gfelem(atoms%n_gf)%atomType = nType
               atoms%gfelem(atoms%n_gf)%lp = lp
               atoms%gfelem(atoms%n_gf)%atomTypep = nType !For now

            ENDDO
         ENDDO

      END SUBROUTINE add_gfjob
 END MODULE m_types_atoms
