!--------------------------------------------------------------------------------
! Copyright (c) 2017 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_make_sphhar
CONTAINS
  SUBROUTINE make_sphhar(atoms,sphhar,sym,cell,oneD)
    USE m_types_setup
    USE m_localsym
    USE m_od_chisym

    TYPE(t_atoms),INTENT(inout)::atoms
    TYPE(t_sphhar),INTENT(inout)::sphhar
    TYPE(t_cell),INTENT(in)::cell
    TYPE(t_sym),INTENT(inout)::sym
    TYPE(t_oneD),INTENT(in)::oned



    INTEGER :: ii,i,j
    INTEGER, ALLOCATABLE :: lmx1(:), nq1(:), nlhtp1(:)

    ! Dimensioning of lattice harmonics

    ALLOCATE(sphhar%clnu(1,1,1),sphhar%nlh(1),sphhar%llh(1,1),sphhar%nmem(1,1),sphhar%mlh(1,1,1))
    sphhar%ntypsd = 0
    IF (.NOT.oneD%odd%d1) THEN
       CALL local_sym(atoms%lmaxd,atoms%lmax,sym%nop,sym%mrot,sym%tau,&
            atoms%nat,atoms%ntype,atoms%neq,cell%amat,cell%bmat,&
            atoms%taual,sphhar%nlhd,sphhar%memd,sphhar%ntypsd,.TRUE.,&
            atoms%nlhtyp,sphhar%nlh,sphhar%llh,&
            sphhar%nmem,sphhar%mlh,sphhar%clnu)
    ELSE IF (oneD%odd%d1) THEN
       WRITE(*,*) 'Note: I would be surprised if lattice harmonics generation works'
       WRITE(*,*) 'Dimensioning of local arrays seems to be inconsistent with routine local_sym'
       ALLOCATE (nq1(atoms%nat),lmx1(atoms%nat),nlhtp1(atoms%nat))
       ii = 1
       nq1=1
       DO i = 1,atoms%ntype
          DO j = 1,atoms%neq(i)
             lmx1(ii) = atoms%lmax(i)
             ii = ii + 1
          END DO
       END DO
       CALL local_sym(atoms%lmaxd,lmx1,sym%nop,sym%mrot,sym%tau,&
            atoms%nat,atoms%nat,nq1,cell%amat,cell%bmat,atoms%taual,&
            sphhar%nlhd,sphhar%memd,sphhar%ntypsd,.TRUE.,nlhtp1,&
            sphhar%nlh,sphhar%llh,sphhar%nmem,&
            sphhar%mlh,sphhar%clnu)
       ii = 1
       DO i = 1,atoms%ntype
          atoms%nlhtyp(i) = nlhtp1(ii)
          ii = ii + atoms%neq(i)
       END DO
       DEALLOCATE (nq1,lmx1,nlhtp1)
    END IF
    DEALLOCATE(sphhar%clnu,sphhar%nlh,sphhar%llh,sphhar%nmem,sphhar%mlh)

    ALLOCATE(sphhar%clnu(sphhar%memd,0:sphhar%nlhd,sphhar%ntypsd))
    ALLOCATE(sphhar%llh(0:sphhar%nlhd,sphhar%ntypsd))
    ALLOCATE(sphhar%mlh(sphhar%memd,0:sphhar%nlhd,sphhar%ntypsd))
    ALLOCATE(sphhar%nlh(sphhar%ntypsd),sphhar%nmem(0:sphhar%nlhd,sphhar%ntypsd))
    ! Generate lattice harmonics

    IF (.NOT.oneD%odd%d1) THEN
       CALL local_sym(atoms%lmaxd,atoms%lmax,sym%nop,sym%mrot,sym%tau,&
            atoms%nat,atoms%ntype,atoms%neq,cell%amat,cell%bmat,atoms%taual,&
            sphhar%nlhd,sphhar%memd,sphhar%ntypsd,.FALSE.,&
            atoms%nlhtyp,sphhar%nlh,sphhar%llh,sphhar%nmem,sphhar%mlh,sphhar%clnu)
       sym%nsymt = sphhar%ntypsd
       !     oneD%mrot1(:,:,:) = sym%mrot(:,:,:)
       !     oneD%tau1(:,:) = sym%tau(:,:)
    ELSE IF (oneD%odd%d1) THEN
       WRITE(*,*) 'Note: I would be surprised if lattice harmonics generation works'
       WRITE(*,*) 'Dimensioning of local arrays seems to be inconsistent with routine local_sym'
       CALL od_chisym(oneD%odd,oneD%mrot1,oneD%tau1,sym%zrfs,sym%invs,sym%invs2,cell%amat)
       ALLOCATE (nq1(atoms%nat),lmx1(atoms%nat),nlhtp1(atoms%nat))
       ii = 1
       DO i = 1,atoms%ntype
          DO j = 1,atoms%neq(i)
             nq1(ii) = 1
             lmx1(ii) = atoms%lmax(i)
             ii = ii + 1
          END DO
       END DO
       CALL local_sym(atoms%lmaxd,lmx1,sym%nop,sym%mrot,sym%tau,&
            atoms%nat,atoms%nat,nq1,cell%amat,cell%bmat,atoms%taual,&
            sphhar%nlhd,sphhar%memd,sphhar%ntypsd,.FALSE.,&
            nlhtp1,sphhar%nlh,sphhar%llh,sphhar%nmem,sphhar%mlh,sphhar%clnu)
       sym%nsymt = sphhar%ntypsd
       ii = 1
       DO i = 1,atoms%ntype
          atoms%nlhtyp(i) = nlhtp1(ii)
          ii = ii + atoms%neq(i)
       END DO
       DEALLOCATE (lmx1,nlhtp1)
    END IF
  END SUBROUTINE make_sphhar
END MODULE m_make_sphhar
