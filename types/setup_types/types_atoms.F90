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
     CHARACTER(LEN=20), ALLOCATABLE :: label(:)
     CHARACTER(len=20), ALLOCATABLE :: speciesName(:)
     !name and other data of explicitely provided xc functional
     CHARACTER(len=4), ALLOCATABLE :: namex(:)
     INTEGER,          ALLOCATABLE :: icorr(:)
     INTEGER,          ALLOCATABLE :: igrd(:)
     

     INTEGER,          ALLOCATABLE :: krla(:)
     LOGICAL,          ALLOCATABLE :: relcor(:)
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
    IMPLICIT NONE
    CLASS(t_atoms),INTENT(INOUT):: tt
    INTEGER,INTENT(IN)               :: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL      :: origin

#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER :: pe,ierr,irank,n(2)
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
    !Integer arrays of size ntype
    IF (pe.NE.irank) THEN
       IF (ALLOCATE(atoms%nz)) THEN
          DEALLOCATE(atoms%nz,atoms%neq,atoms%jri,atoms%jri,atoms%ncst,atoms%numStatesProvided)
          DEALLOCATE(atoms%lmax,atoms%lnonsph,atoms%ncv,atoms%nlo,atoms%lapw_l)
          DEALLOCATE(atoms%nlhtyp,atoms%l_geo,atoms%rmt,atoms%dx,atoms%volmts,atoms%zatom)
          DEALLOCATE(atoms%bmu,atoms%label
    
    
    CALL MPI_BCAST(tt%,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%,1,MPI_INTEGER,pe,mpi_comm,ierr)












    CALL MPI_BCAST(tt%nmzxy,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%layers,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%nvac,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%nstars,1,MPI_INTEGER,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%nstm,1,MPI_INTEGER,pe,mpi_comm,ierr)
    
    CALL MPI_BCAST(tt%delz,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%dvac,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%tworkf,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%locx,2,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
    CALL MPI_BCAST(tt%locy,2,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)

    CALL MPI_BCAST(tt%starcoeff,1,MPI_LOGICAL,pe,mpi_comm,ierr)

    IF (irank==pe) THEN
       n=SHAPE(tt%izlay)
       CALL MPI_BCAST(n,2,MPI_INTEGER,pe,mpi_comm,ierr)
    ELSE
       CALL MPI_BCAST(n,2,MPI_INTEGER,pe,mpi_comm,ierr)
       IF (ALLOCATED(tt%izlay)) DEALLOCATE(tt%izlay)
       ALLOCATE(tt%izlay(n(1),n(2)))
    ENDIF
    CALL MPI_BCAST(tt%izlay,SIZE(tt%izlay),MPI_INTEGER,pe,mpi_comm,ierr)
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

    WRITE(unit,*,IOSTAT=iostat) '"atoms":{'


    CALL json_print(unit,"nmz",tt%nmz)
    CALL JSON_PRINT(unit,"nmzxy",tt%nmzxy)
    CALL JSON_PRINT(unit,"layers",tt%layers) 
    CALL JSON_PRINT(unit, "nvac",tt%nvac)
    CALL JSON_PRINT(unit, "delz",tt%delz)  
    CALL JSON_PRINT(unit, "dvac",tt%dvac)  
    CALL JSON_PRINT(unit, "nstars",tt%nstars)  
    CALL JSON_PRINT(unit, "nstm",tt%nstm)  
    CALL JSON_PRINT(unit, "tworkf",tt%tworkf)  
    CALL JSON_PRINT(unit, "locx",tt%locx)  
    CALL JSON_PRINT(unit, "locy",tt%locy)  
    CALL JSON_PRINT(unit, "starcoeff",tt%starcoeff)
    CALL JSON_PRINT(unit, "izlay_shape",SHAPE(tt%izlay))
    CALL JSON_PRINT(unit,"izlay",RESHAPE(tt%izlay,(/SIZE(tt%izlay)/)))
    
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

    CHARACTER(len=40)::string
    INTEGER,ALLOCATABLE:: itemp(:),n(2)
    CALL json_open_class("atoms",unit,iostat)
    IF (iostat.NE.0)   RETURN

    CALL json_read(unit,"nmz",tt%nmz)
    CALL JSON_READ(unit,"nmzxy",tt%nmzxy)
    CALL JSON_READ(unit,"layers",tt%layers) 
    CALL JSON_READ(unit, "nvac",tt%nvac)
    CALL JSON_READ(unit, "delz",tt%delz)  
    CALL JSON_READ(unit, "dvac",tt%dvac)  
    CALL JSON_READ(unit, "nstars",tt%nstars)  
    CALL JSON_READ(unit, "nstm",tt%nstm)  
    CALL JSON_READ(unit, "tworkf",tt%tworkf)  
    CALL JSON_READ(unit, "locx",tt%locx)  
    CALL JSON_READ(unit, "locy",tt%locy)  
    CALL JSON_READ(unit, "starcoeff",tt%starcoeff)
    CALL JSON_READ(unit, "izlay_shape",n)

    IF (ALLOCATED(tt%izlay)) DEALLOCATE(tt%izlay)
    ALLOCATE(tt%izlay(n(1),n(2)),itemp(n(1)*n(2)))
    CALL JSON_READ(unit,"izlay",itemp)
    tt%izlay=RESHAPE(itemp,n)

    CALL json_close_class(unit)
    
  END SUBROUTINE read_atoms

 


  SUBROUTINE read_xml_atoms(atoms)
    USE m_xmlIntWrapFort
    USE m_constants
    USE m_calculator
    USE m_invs3
    IMPLICIT NONE
    CLASS(t_atoms),INTENT(OUT):: atoms


    LOGICAL::film,vacdos
    CHARACTER(len=200):: xpathA,valueString


    film=(xmlGetNumberOfNodes('/fleurInput/vacuum/filmLattice')==1)
    vacuum%nvac = 2

    vacuum%nmz = 250
    vacuum%delz = 25.0/vacuum%nmz
    vacuum%nmzxy = 100

    !DOS stuff for vacuum
    vacuum%layers = 1
    vacuum%starcoeff = .FALSE.
    vacuum%nstars = 0
    vacuum%locx = 0.0
    vacuum%locy = 0.0
    vacuum%nstm = 0
    vacuum%tworkf = 0.0
    !Check if vacdos=t
    vacdos=.FALSE.
    IF (xmlGetNumberOfNodes('/fleurInput/output')==1) &
         vacdos = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/output/@vacdos'))
    
    xPathA = '/fleurInput/output/vacuumDOS'
    IF (film.AND.vacdos.AND.xmlGetNumberOfNodes(xPathA).EQ.1) THEN
       vacuum%layers = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@layers'))
       vacuum%starcoeff = evaluateFirstBoolOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@star'))
       vacuum%nstars = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nstars'))
       vacuum%locx(1) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locx1'))
       vacuum%locx(2) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locx2'))
       vacuum%locy(1) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locy1'))
       vacuum%locy(2) = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locy2'))
       vacuum%nstm = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nstm'))
       vacuum%tworkf = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@tworkf'))
    END IF

    IF (vacuum%layers>1) THEN
       ALLOCATE(vacuum%izlay(vacuum%layerd,2))
       CALL judft_error("Layer mode not supported in types_vacuum")
    ELSE
       ALLOCATE(vacuum%izlay(1,1))
       vacuum%izlay=1.0
    END IF
    
  END SUBROUTINE read_xml_atoms


  SUBROUTINE init_atoms(atoms)
    USE types_sym
    IMPLICIT NONE
    CLASS(t_vacuum),INTENT(OUT):: atoms
    
  END SUBROUTINE init_atoms

    
END MODULE m_types_atoms
