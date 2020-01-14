MODULE m_make_atom_groups
  USE m_juDFT
  !********************************************************************
  !     calculates the space group operations given the lattice vectors
  !     and the atomic positions.                              mw 12-99
  !     Modified by GM (2016)
  !********************************************************************
CONTAINS
  SUBROUTINE make_atom_groups(sym,cell,atompos,atomid,atomlabel,atoms)
    !Use the symmetry to generate correct mapping of atoms into types
    USE m_types_sym
    USE m_types_cell
    USE m_types_atoms

    IMPLICIT NONE
    TYPE(t_sym),INTENT(in)     :: sym
    TYPE(t_cell),INTENT(IN)    :: cell
    REAL,INTENT(in)            :: atompos(:,:)
    REAL,INTENT(in)            :: atomid(:)
    CHARACTER(len=*),INTENT(in):: atomlabel(:)
    TYPE(t_atoms),INTENT(out)  :: atoms

    INTEGER, ALLOCATABLE :: natype(:),natrep(:),ity(:)  ! or  'nat'
    INTEGER              :: ntypm,n,i,j,ntype
    LOGICAL              :: lnew
    REAL                 :: tr(3)
    REAL,PARAMETER       :: eps7=1.e-7,eps12=1e-12

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


    natype(1:size(atomid)) = 0
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
    DO n=1,atoms%nat
       atoms%neq( natype(n) ) = atoms%neq( natype(n) ) + 1
       atoms%zatom( natype(n) ) = atomid(n)
    ENDDO
    DO n=1,atoms%ntype
      atoms%taual(1,sum(atoms%neq(:n-1))+1:sum(atoms%neq(:n)))=pack(atompos(1,:),natype==n)
      atoms%taual(2,sum(atoms%neq(:n-1))+1:sum(atoms%neq(:n)))=pack(atompos(2,:),natype==n)
      atoms%taual(3,sum(atoms%neq(:n-1))+1:sum(atoms%neq(:n)))=pack(atompos(3,:),natype==n)
      atoms%label(sum(atoms%neq(:n-1))+1:sum(atoms%neq(:n)))=pack(atomlabel(:),natype==n)
    enddo
    WHERE ( ABS( atoms%taual ) < eps12 ) atoms%taual = 0.00

    !Generate postions in cartesian coordinates
    atoms%pos(:,:) = MATMUL( cell%amat , atoms%taual(:,:) )

  END SUBROUTINE make_atom_groups
END MODULE m_make_atom_groups
