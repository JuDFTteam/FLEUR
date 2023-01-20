!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_banddos
  USE m_juDFT
  USE m_types_fleurinput_base
  IMPLICIT NONE
  PRIVATE
  PUBLIC:: t_banddos
  TYPE,EXTENDS(t_fleurinput_base):: t_banddos
     LOGICAL :: dos =.FALSE.
     REAL    :: e1_dos=0.5
     REAL    :: e2_dos=-0.5
     REAL    :: sig_dos=0.015
     INTEGER :: ndos_points=1301
     LOGICAL :: l_storeEVData = .TRUE.


     LOGICAL :: vacdos =.FALSE.
     INTEGER :: layers=0
     INTEGER :: nstars=0
     INTEGER :: nstm=0
     REAL    :: tworkf=0.0
     REAL    :: locx(2)=[0.,0.]
     REAL    :: locy(2)=[0.,0.]
     LOGICAL :: starcoeff=.FALSE.
     INTEGER, ALLOCATABLE :: izlay(:, :)


     LOGICAL :: band =.FALSE.
     LOGICAL :: unfoldband =.FALSE.
     LOGICAL :: unfoldUseOlap = .TRUE.
     INTEGER :: s_cell_x=1
     INTEGER :: s_cell_y=1
     INTEGER :: s_cell_z=1
     INTEGER, ALLOCATABLE :: unfoldTransMat(:,:)


     LOGICAL :: l_mcd =.FALSE.
     REAL    :: e_mcd_lo =-10.0
     REAL    :: e_mcd_up= 0.0

     LOGICAL :: l_orb =.FALSE.

     LOGICAL :: l_jDOS = .FALSE.

     LOGICAL :: l_slab=.false.

     LOGICAL,ALLOCATABLE :: dos_atom(:) ! for each atom (not type) switch on DOS
     INTEGER,ALLOCATABLE :: dos_typelist(:) !list of types for which DOS is calculated
     INTEGER,ALLOCATABLE :: dos_atomlist(:) !list of atoms for which DOS is calculated
     INTEGER,ALLOCATABLE  :: map_atomtype(:) !map an atomtype to corresponding entry in DOS
     real,allocatable    :: alpha(:),beta(:),gamma(:)
     !INTEGER :: ndir =0
     !INTEGER :: orbCompAtom=0
     !INTEGER :: jDOSAtom=0
     !INTEGER :: projdos !selects one atomtype and prints the projected dos if there are to many atoms
   CONTAINS
     PROCEDURE :: read_xml=>read_xml_banddos
     PROCEDURE :: mpi_bc=>mpi_bc_banddos
  END TYPE t_banddos
CONTAINS
  SUBROUTINE mpi_bc_banddos(this,mpi_comm,irank)
    USE m_mpi_bc_tool
    CLASS(t_banddos),INTENT(INOUT)::this
    integer,INTENT(IN):: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL::irank
    INTEGER ::rank
    if (present(irank)) THEN
       rank=irank
    else
       rank=0
    end if
    CALL mpi_bc(this%dos ,rank,mpi_comm)
    CALL mpi_bc(this%band ,rank,mpi_comm)
    CALL mpi_bc(this%l_mcd ,rank,mpi_comm)
    CALL mpi_bc(this%l_orb ,rank,mpi_comm)
    CALL mpi_bc(this%l_jDOS,rank,mpi_comm)
    CALL mpi_bc(this%vacdos ,rank,mpi_comm)
    CALL mpi_bc(this%e1_dos,rank,mpi_comm)
    CALL mpi_bc(this%e2_dos,rank,mpi_comm)
    CALL mpi_bc(this%sig_dos,rank,mpi_comm)
    CALL mpi_bc(this%e_mcd_lo ,rank,mpi_comm)
    CALL mpi_bc(this%e_mcd_up,rank,mpi_comm)
    CALL mpi_bc(this%unfoldband ,rank,mpi_comm)
    CALL mpi_bc(this%unfoldUseOlap ,rank,mpi_comm)
    CALL mpi_bc(this%s_cell_x,rank,mpi_comm)
    CALL mpi_bc(this%s_cell_y,rank,mpi_comm)
    CALL mpi_bc(this%s_cell_z,rank,mpi_comm)
    CALL mpi_bc(this%unfoldTransMat, rank, mpi_comm)
    CALL mpi_bc(this%alpha,rank,mpi_comm)
    CALL mpi_bc(this%beta,rank,mpi_comm)
    CALL mpi_bc(this%gamma,rank,mpi_comm)
    CALL mpi_bc(this%dos_atom,rank,mpi_comm)
    CALL mpi_bc(this%l_slab,rank,mpi_comm)
    CALL mpi_bc(this%ndos_points,rank,mpi_comm)
    CALL mpi_bc(this%layers,rank,mpi_comm)
    CALL mpi_bc(this%nstars,rank,mpi_comm)
    CALL mpi_bc(this%nstm,rank,mpi_comm)
    CALL mpi_bc(this%tworkf,rank,mpi_comm)
    CALL mpi_bc(this%locx(1),rank,mpi_comm)
    CALL mpi_bc(this%locy(1),rank,mpi_comm)
    CALL mpi_bc(this%locx(2),rank,mpi_comm)
    CALL mpi_bc(this%locy(2),rank,mpi_comm)
    CALL mpi_bc(this%starcoeff,rank,mpi_comm)
    CALL mpi_bc(this%izlay,rank,mpi_comm)
    CALL mpi_bc(this%dos_atom,rank,mpi_comm)
    CALL mpi_bc(this%dos_atomlist,rank,mpi_comm)
    CALL mpi_bc(this%dos_typelist,rank,mpi_comm)
    CALL mpi_bc(this%map_atomtype,rank,mpi_comm)

  END SUBROUTINE mpi_bc_banddos
  SUBROUTINE read_xml_banddos(this,xml)
    USE m_types_xml
    CLASS(t_banddos),INTENT(INOUT)::this
    TYPE(t_xml),INTENT(INOUT)::xml

    CHARACTER(len=200) :: xPathA, xPathB, valueString
    INTEGER::numberNodes,iType,i,na,n,n_dos_atom,n_dos_type
    LOGICAL::all_atoms,dos_atom_found
    integer,allocatable:: dos_atomlist(:),dos_typelist(:),neq(:)

    this%band = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/@band'))
    this%dos = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/@dos'))
    !this%l_slab = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/@slab'))
    all_atoms=.true.
    numberNodes = xml%GetNumberOfNodes('/fleurInput/output/bandDOS')
    IF (numberNodes.EQ.1) THEN
       IF(xml%versionNumber>=34) THEN
          this%l_storeEVData=evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/bandDOS/@storeEVData'))
       ENDIF
       this%l_orb=evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/bandDOS/@orbcomp'))
       this%l_jDOS=evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/bandDOS/@jDOS'))
       all_atoms=evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/bandDOS/@all_atoms'))
       this%e2_dos = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/bandDOS/@minEnergy'))
       this%e1_dos = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/bandDOS/@maxEnergy'))
       this%sig_dos = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/bandDOS/@sigma'))
       this%ndos_points=evaluateFirstIntOnly(xml%GetAttributeValue('/fleurInput/output/bandDOS/@numberPoints'))
    END IF

    ! Read in optional magnetic circular dichroism parameters
    numberNodes = xml%GetNumberOfNodes('/fleurInput/output/magneticCircularDichroism')
    IF (numberNodes.EQ.1) THEN
      this%l_mcd=evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output//magneticCircularDichroism/@mcd'))
      this%e_mcd_lo = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/magneticCircularDichroism/@energyLo'))
      this%e_mcd_up = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/magneticCircularDichroism/@energyUp'))
    END IF

    allocate(this%dos_atom(xml%get_nat()))
    allocate(this%alpha(size(this%dos_atom)));this%alpha=0.0
    allocate(this%beta(size(this%dos_atom)));this%beta=0.0
    allocate(this%gamma(size(this%dos_atom)));this%gamma=0.0
    allocate(neq(xml%get_ntype()), source=0)
    ALLOCATE(this%unfoldTransMat(3,3))

    !check if cdinf is given
    IF (xml%GetNumberOfNodes('/fleurInput/output/checks').EQ.1) THEN
      this%dos_atom= evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/checks/@cdinf'))
    else
      this%dos_atom=.false.
    endif
    this%dos_atom=(all_atoms.and.(this%dos.or.this%band.or.this%dos_atom))
    na = 0
    IF (xml%versionNumber > 31) then
    DO iType = 1, xml%GetNumberOfNodes('/fleurInput/atomGroups/atomGroup')
       WRITE(xPathA,*) '/fleurInput/atomGroups/atomGroup[',iType,']'
       if (xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/relPos')>0) THEN
         neq(itype)= xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/relPos')
         xPathA=TRIM(ADJUSTL(xPathA))//'/relPos'
       elseif(xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/absPos')>0) THEN
         neq(itype)= xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/absPos')
         xPathA=TRIM(ADJUSTL(xPathA))//'/absPos'
       else
         neq(itype)= xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathA))//'/filmPos')
         xPathA=TRIM(ADJUSTL(xPathA))//'/filmPos'
      endif
       DO i = 1, neq(itype)
          na = na + 1
          WRITE(xPathB,*) TRIM(ADJUSTL(xPathA))//'[',i,']'
          if (xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathB))//'/@banddos')==1) &
          this%dos_atom(na) = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@banddos'))
          if (xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathB))//'/@alpha')==1) &
          this%alpha(na) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@alpha'))
          if (xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathB))//'/@beta')==1) &
          this%beta(na) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@beta'))
          if (xml%GetNumberOfNodes(TRIM(ADJUSTL(xPathB))//'/@gamma')==1) &
          this%gamma(na) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/@gamma'))
        ENDDO
    ENDDO
    endif
    ! Read in optional parameter for unfolding bandstructure of supercell
    numberNodes = xml%GetNumberOfNodes('/fleurInput/output/unfoldingBand')
    IF (numberNodes.EQ.1) THEN
       this%unfoldband = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/unfoldingBand/@unfoldBand'))
       this%s_cell_x = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/unfoldingBand/@supercellX'))
       this%s_cell_y = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/unfoldingBand/@supercellY'))
       this%s_cell_z = evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/output/unfoldingBand/@supercellZ'))
       numberNodes = xml%GetNumberOfNodes('/fleurInput/output/unfoldingBand/@useOlap')
       this%unfoldUseOlap = .TRUE.
       IF (numberNodes.EQ.1) this%unfoldUseOlap = evaluateFirstBoolOnly(xml%GetAttributeValue('/fleurInput/output/unfoldingBand/@useOlap'))
       this%unfoldTransMat(:,:) = 0
       this%unfoldTransMat(1,1) = 1
       this%unfoldTransMat(2,2) = 1
       this%unfoldTransMat(3,3) = 1
       numberNodes = xml%GetNumberOfNodes('/fleurInput/output/unfoldingBand/transMat')
       IF (numberNodes.EQ.1) THEN
          xPathB = '/fleurInput/output/unfoldingBand/transMat'
          valueString = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/row-1')))
          this%unfoldTransMat(1,1) = NINT(evaluateFirst(valueString))
          this%unfoldTransMat(2,1) = NINT(evaluateFirst(valueString))
          this%unfoldTransMat(3,1) = NINT(evaluateFirst(valueString))
          valueString = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/row-2')))
          this%unfoldTransMat(1,2) = NINT(evaluateFirst(valueString))
          this%unfoldTransMat(2,2) = NINT(evaluateFirst(valueString))
          this%unfoldTransMat(3,2) = NINT(evaluateFirst(valueString))
          valueString = TRIM(ADJUSTL(xml%GetAttributeValue(TRIM(ADJUSTL(xPathB))//'/row-3')))
          this%unfoldTransMat(1,3) = NINT(evaluateFirst(valueString))
          this%unfoldTransMat(2,3) = NINT(evaluateFirst(valueString))
          this%unfoldTransMat(3,3) = NINT(evaluateFirst(valueString))
       END IF
    END IF
    if (xml%versionNumber <32) return
    xPathA = '/fleurInput/output/vacuumDOS'
    IF (xml%GetNumberOfNodes(xpathA)==1) THEN
       this%vacdos = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@vacdos'))
       this%starcoeff = evaluateFirstBoolOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@star'))
       this%nstars = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nstars'))
       this%locx(1) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locx1'))
       this%locx(2) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locx2'))
       this%locy(1) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locy1'))
       this%locy(2) = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@locy2'))
       this%nstm = evaluateFirstIntOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@nstm'))
       this%tworkf = evaluateFirstOnly(xml%GetAttributeValue(TRIM(ADJUSTL(xPathA))//'/@tworkf'))

       this%layers = xml%GetNumberOfNodes('/fleurInput/output/vacuumDOS/layer')
       ALLOCATE(this%izlay(this%layers,2))
       DO n=1,this%layers
         write(xPathA,'(a,i0,a)') '/fleurInput/output/vacuumDOS/layer[',n,']'
         valueString = xml%GetAttributeValue(TRIM(ADJUSTL(xPathA)))
         this%izlay(n,1)=evaluateFirst(valueString)
         this%izlay(n,2)=evaluateFirst(valueString)
       ENDDO
     ELSE
       allocate(this%izlay(0,2))
    END IF

    !Create a list of all atoms and all types for which the DOS is calculated
    ALLOCATE(dos_atomlist(xml%get_nat()),source=0)
    allocate(dos_typelist(xml%get_ntype()),source=0)

    na=0
    n_dos_atom=0
    n_dos_type=0
    DO itype=1,xml%get_ntype()
      dos_atom_found=.false.
      DO n=1,neq(itype)
        na=na+1
        if (this%dos_atom(na)) then
          dos_atom_found=.true.
          n_dos_atom=n_dos_atom+1
          dos_atomlist(n_dos_atom)=na
        endif
      enddo
      if (dos_atom_found)then
        n_dos_type=n_dos_type+1
        dos_typelist(n_dos_type)=iType
      ENDIF
    ENDDO
    this%dos_atomlist=dos_atomlist(:n_dos_atom)
    this%dos_typelist=dos_typelist(:n_dos_type)
    !create map
    allocate(this%map_atomtype(size(neq)))
    this%map_atomtype=0
    DO itype=1,size(neq)
      DO n=1,size(dos_typelist)
        if (dos_typelist(n)==itype) this%map_atomtype(itype)=n
      ENDDO
    ENDDO



  END SUBROUTINE read_xml_banddos

END MODULE m_types_banddos
