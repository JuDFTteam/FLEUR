MODULE m_make_crystal
  USE m_juDFT
  !********************************************************************
  !      generate space group operations from lattice information
  !********************************************************************
CONTAINS
  SUBROUTINE make_crystal(film, atomid,atompos,mag_mom,atomlabel,dvac,noco,&
       cell,sym,atoms)

    USE m_types_cell
    USE m_types_sym
    USE m_types_atoms
    USE m_types_noco
    USE m_constants
    USE m_make_spacegroup
    USE m_make_atom_groups
    USE m_inv3
    !USE m_generator
    IMPLICIT NONE
    !===> Arguments
    LOGICAL, INTENT(IN)     :: film
    REAL,    INTENT(IN)     :: atomid(:)
    REAL,    INTENT(INOUT)  :: atompos(:,:)!might be shifted
    CHARACTER(len=*),INTENT(IN)::atomlabel(:)
    REAL,    INTENT(IN)     :: dvac
    TYPE(t_noco),INTENT(in) :: noco

    TYPE(t_cell),INTENT(inout)   ::cell
    TYPE(t_sym),INTENT(inout) ::sym !symor is checked
    TYPE(t_atoms),INTENT(out) ::atoms
    REAL,ALLOCATABLE,INTENT(INOUT)::mag_mom(:,:)

    !===> Local Variables
    INTEGER :: i,j,k,n,m,na,nt,inversionOp
    LOGICAL :: l_posCorrected(3)
    REAL    :: rest, det
    REAL,PARAMETER :: eps7 = 1.0e-7
    INTEGER,PARAMETER :: invs_matrix(3,3)=RESHAPE([-1,0,0,0,-1,0,0,0,-1],[3,3])
    INTEGER,ALLOCATABLE :: inpgen_atom_for_type(:)
    REAL,ALLOCATABLE    :: mag_mom_tmp(:,:)


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
    IF (film) THEN
       cell%amat(3,3)=max(cell%amat(3,3),8.0+2*maxval(abs(atompos(3,:))))
       CALL inv3(cell%amat,cell%bmat,det)
       cell%bmat = tpi_const*cell%bmat
       atompos(3,:)=atompos(3,:)/cell%amat(3,3)
       DO n=1,SIZE(atompos,2)
         atompos(:2,n) = atompos(:2,n) - ANINT( atompos(:2,n) - eps7 )
       ENDDO
    ELSE
      DO n=1,SIZE(atompos,2)
        atompos(:,n) = atompos(:,n) - ANINT( atompos(:,n) - eps7 )
      ENDDO
    ENDIF

    DO na = 1, SIZE(atompos,2)
       l_posCorrected(:) = .FALSE.
       DO i = 2, 40
          rest = ABS(i*atompos(1, na) - NINT(i*atompos(1, na)))
          IF (.NOT.l_posCorrected(1) .AND. (rest .LT. (i*0.000001))) THEN
             l_posCorrected(1) = .TRUE.
             atompos(1, na) = NINT(i*atompos(1, na)) / REAL(i)
          END IF
          rest = ABS(i*atompos(2, na) - NINT(i*atompos(2, na)))
          IF (.NOT.l_posCorrected(2) .AND. (rest .LT. (i*0.000001))) THEN
             l_posCorrected(2) = .TRUE.
             atompos(2, na) = NINT(i*atompos(2, na)) / REAL(i)
          END IF
          IF (.NOT.film) THEN
             rest = ABS(i*atompos(3, na) - NINT(i*atompos(3, na)))
             IF (.NOT.l_posCorrected(3) .AND. (rest .LT. (i*0.000001))) THEN
                l_posCorrected(3) = .TRUE.
                atompos(3, na) = NINT(i*atompos(3, na)) / REAL(i)
             END IF
          END IF
       END DO
    END DO

    !--->    calculate space group symmetry
    CALL make_spacegroup(film,noco,cell,atompos,atomid,mag_mom,sym)
    ! Check whether there is an inversion center that is not at the
    ! origin and if one is found shift the crystal such that the
    ! inversion is with respect to the origin. Then recalculate
    ! symmetry operations.
    inversionOp = -1
    symOpLoop: DO k = 1, sym%nop
       IF (ALL(sym%mrot(:,:,k)==invs_matrix)) THEN
          inversionOp = k
          EXIT symOpLoop
       END IF
    END DO symOpLoop
    IF (inversionOp.GT.0) THEN
       IF(ANY(ABS(sym%tau(:,inversionOp)).GT.eps7).and..not.(film.and.ABS(sym%tau(3,inversionOp))>eps7)) THEN
          WRITE(*,*) 'Found inversion center at finite position.'
          WRITE(*,*) 'Shifting crystal by:'
          WRITE(*,'(3f15.10)') -0.5*sym%tau(:,inversionOp)            ! changed to minus
          WRITE(*,*) ''
          DO k = 1, SIZE(atomid)
             atompos(:,k) = atompos(:,k) - 0.5*sym%tau(:,inversionOp) ! GB`18
             atompos(:,k) = atompos(:,k) - ANINT( atompos(:,k) - eps7 )
          END DO
          CALL make_spacegroup(film,noco,cell,atompos,atomid,mag_mom,sym)
       END IF
    END IF

    !finish generation of symmetry....
    CALL sym%init(cell,film)

    !Generate basic atom type
    ALLOCATE(inpgen_atom_for_type(size(atomid)))
    CALL make_atom_groups(sym,cell,atompos,atomid,atomlabel,atoms,inpgen_atom_for_type)

    !assign magnetic moments to types
    allocate(mag_mom_tmp(3,atoms%ntype))
    DO nt=1,atoms%ntype
      mag_mom_tmp(:,nt)=mag_mom(:,inpgen_atom_for_type(nt))
    enddo
    call move_alloc(mag_mom_tmp,mag_mom)

   

    !--->    determine a set of generators for this group
    !CALL generator(sym%nop,sym%mrot,sym%tau,6,0)

    !---> output: the atomic positions, etc.

    WRITE (oUnit,'(//," Atomic positions:",/,1x,17("-"))')
    WRITE (oUnit,'(" atom types =",i5/,"      total =",i5)') atoms%ntype,atoms%nat
    WRITE (oUnit,'(/,7x,"lattice coordinates",15x,"(scaled) Cartesian coordinates   atom")')

    na = 0
    DO nt=1,atoms%ntype
       DO n=1,atoms%neq(nt)
          WRITE (oUnit,'(3f10.6,10x,3f10.6,i7)') atoms%taual(:,na+n),atoms%pos(:,na+n),na+n
       ENDDO
       na = na + atoms%neq(nt)
    ENDDO

    RETURN

  END SUBROUTINE make_crystal
END MODULE m_make_crystal
