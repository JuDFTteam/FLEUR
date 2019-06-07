MODULE m_make_crystal
  USE m_juDFT
  !********************************************************************
  !      generate space group operations from lattice information
  !********************************************************************
CONTAINS
  SUBROUTINE make_crystal(film, symor,atomid,atompos,atomlabel,amat,dvac,noco,&
       cell,sym,atoms)
    USE m_types
    USE m_spggen
    USE m_generator
    IMPLICIT NONE
    !===> Arguments
    LOGICAL, INTENT(IN)     :: film
    LOGICAL, INTENT(IN)     :: symor       ! on input: if true, reduce symmetry if oldfleur
    REAL,    INTENT(IN)     :: atomid(:)
    REAL,    INTENT(INOUT)  :: atompos(:,;)!might be shifted
    CHARACTER(len=*),INTENT(IN)::atomlabel(:)
    REAL,    INTENT(IN)     :: amat,dvac
    TYPE(t_noco),INTENT(in) :: noco

    TYPE(t_cell),INTENT(out)  ::cell
    TYPE(t_sym),INTENT(out)   ::sym
    TYPE(t_atoms),INTENT(out) ::atoms
  
 
    !===> Local Variables
    INTEGER :: i,j,k,n,m,na,nt,inversionOp
    REAL,PARAMETER :: eps7 = 1.0e-7  
  

    !Generate cell-datatype
    cell%amat=amat
    IF (film) THEN
       CALL cell%init(dvac)
    ELSE
       CALL cell%init()
    ENDIF

    
    !generate basic atom type
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !---> atomic positions:
    !--->
    !---> atomic positions input must be in
    !---> lattice vector units, as determined by logical cartesian.
    !
    !---> read in number of atoms or types

    !---> read in atomic positions and shift to (-1/2,1/2] in lattice
    !---> coords. also read in identification (atomic) number (atomid)
    !---> to distinguish different atom types (need not be atomic number)
    DO n=1,atoms%nat
       atompos(:,n) = atompos(:,n) - ANINT( atompos(:,n) - eps7 )
    ENDDO

    !--->    calculate space group symmetry
    sym%symor=symor
    CALL make_spacegroup(film,noco,cell,atompos,atomid,sym)
    ! Check whether there is an inversion center that is not at the
    ! origin and if one is found shift the crystal such that the
    ! inversion is with respect to the origin. Then recalculate
    ! symmetry operations.
    inversionOp = -1
    symOpLoop: DO k = 1, sym%nops
       IF (ALL(sym%mrot(:,:,k)==invs_matrix)) THEN
          inversionOp = k
          EXIT symOpLoop
       END IF
    END DO symOpLoop
    IF (inversionOp.GT.0) THEN
       IF(ANY(ABS(sym%tau(:,inversionOp)).GT.eps7)) THEN
          WRITE(*,*) 'Found inversion center at finite position.'
          WRITE(*,*) 'Shifting crystal by:'
          WRITE(*,'(3f15.10)') -0.5*sym%tau(:,inversionOp)            ! changed to minus
          WRITE(*,*) ''
          DO k = 1, size(atomid)
             atompos(:,k) = atompos(:,k) - 0.5*sym%tau(:,inversionOp) ! GB`18
          END DO
          CALL make_spacegroup(film,noco,cell,atompos,atomid,sym)
       END IF
    END IF

    !finish generation of symmetry....
    CALL sym%init(cell%amat,film)

    !Generate basic atom type
    call make_atom_groups(sym,atompos,atomid,atomlabel,atoms)
    
    !--->    determine a set of generators for this group
    CALL generator(sym%nops,sym%mrot,sym%tau,6,0)

    !---> output: the atomic positions, etc.

    WRITE (6,'(//," Atomic positions:",/,1x,17("-"))')
    WRITE (6,'(" atom types =",i5/,"      total =",i5)') ntype,nat
    WRITE (6,'(/,7x,"lattice coordinates",15x,"(scaled) Cartesian coordinates   atom")')

    na = 0
    DO nt=1,atoms%ntype
       DO n=1,atoms%neq(nt)
          WRITE (6,'(3f10.6,10x,3f10.6,i7)') atoms%taual(:,na+n),atoms%pos(:,na+n),na+n
       ENDDO
       na = na + atoms%neq(nt)
    ENDDO

    RETURN

  END SUBROUTINE crystal
END MODULE m_crystal
