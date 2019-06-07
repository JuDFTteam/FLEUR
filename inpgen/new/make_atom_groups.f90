MODULE m_make_atom_groups
  USE m_juDFT
  !********************************************************************
  !     calculates the space group operations given the lattice vectors
  !     and the atomic positions.                              mw 12-99
  !     Modified by GM (2016)
  !********************************************************************
CONTAINS
  SUBROUTINE make_atom_groups(sym,atompos,atomid,atomlabel,atoms)
    !Use the symmetry to generate correct mapping of atoms into types
    IMPLICIT NONE
    TYPE(t_sym),INTENT(in)     :: sym
    REAL,INTENT(in)            :: atompos(:,:)
    REAL,INTENT(in)            :: atomid(:,:)
    CHARACTER(len=*),INTENT(in):: atomlabel(:)
    TYPE(t_atoms),INTENT(out)  :: atoms

    INTEGER, ALLOCATABLE :: natype(:),natrep(:),ity(:)  ! or  'nat'
    INTEGER              :: ntypm,n,i,j
    LOGICAL              :: lnew
    
    ALLOCATE(natype(SIZE(atomid)),natrep(SIZE(atomid)),ity(SIZE(atomid)))
  
    ntypm = 1
    ity(1) = 1
    DO n=2, SIZE(atomid)
       lnew = .TRUE.
       DO i=1,n-1
          IF ( ABS( atomid(i)-atomid(n) ) < eps7 ) THEN
             ity(n) = ity(i)
             lnew = .FALSE.
             EXIT
          ENDIF
       ENDDO
       IF (lnew) THEN
          ntypm = ntypm + 1
          ity(n) = ntypm
       ENDIF
    ENDDO

    
    natype(1:nat) = 0
    ntype = 0
    DO i =1,SIZE(atomid)
       IF ( natype(i) .NE. 0 ) CYCLE
       ntype = ntype + 1   ! new atom type
       natype(i) = ntype   ! atom type
       natrep(i) = i       ! this is the representative
       !--->    rotate representative and get symmetry equavalent atoms
       DO n=1,sym%nop
          tr(:) = MATMUL( sym%mrot(:,:,n) , atompos(:,i) ) + sym%tau(:,n)
          tr(:) = tr(:) - ANINT( tr(:) -eps7 )
          !--->       this (rotated) atom done already? (correct symmetry assumed)
          DO j=i+1,SIZE(atomid)
             IF ( natype(j) .NE. 0 ) CYCLE
             IF ( ity(j) .NE. ity(i) ) CYCLE
             IF ( ANY( ABS( tr(:) - atompos(:,j) ) > eps7 ) ) CYCLE
             natrep(j) = i      ! representative atom
             natype(j) = ntype  ! atom type
             EXIT
          ENDDO
       ENDDO
    ENDDO

    
    atoms%ntype=ntype
    atoms%nat=SIZE(atomid)
    
    ALLOCATE(atoms%neq(ntype),atoms%taual(3,atoms%nat),atoms%pos(3,atoms%nat),atoms%zatom(ntype),atoms%label(atoms%nat))

    atoms%neq(1:ntype) = 0
    DO n=1,nat
       atoms%neq( natype(n) ) = atoms%neq( natype(n) ) + 1
       atoms%zatom( natype(n) ) = atom_id(n)
    ENDDO
    atoms%taual=atompos(:,:atoms%nat)
    atoms%label=atomlabel
    WHERE ( ABS( atoms%taual ) < eps12 ) atoms%taual = 0.00
    
    !Generate postions in cartesian coordinates
    ALLOCATE(atoms%pos(3,atoms%nat)
    atoms%pos(:,n) = MATMUL( cell%amat , atoms%taual(:,n) )
    
  END SUBROUTINE make_atom_groups
END MODULE m_make_atom_groups
