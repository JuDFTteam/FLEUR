!--------------------------------------------------------------------------------
! Copyright (c) 2017 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_make_sphhar
CONTAINS
  SUBROUTINE make_sphhar(l_write,atoms,sphhar,sym,cell )
    USE m_types_setup
    USE m_localsym
     
    LOGICAL,INTENT(IN) :: l_write
    TYPE(t_atoms),INTENT(inout)::atoms
    TYPE(t_sphhar),INTENT(inout)::sphhar
    TYPE(t_cell),INTENT(in)::cell
    TYPE(t_sym),INTENT(inout)::sym
     



    INTEGER :: ii,i,j
    INTEGER, ALLOCATABLE :: lmx1(:), nq1(:), nlhtp1(:)

    ! Dimensioning of lattice harmonics

    ALLOCATE(sphhar%clnu(1,1,1),sphhar%nlh(1),sphhar%llh(1,1),sphhar%nmem(1,1),sphhar%mlh(1,1,1))
    sphhar%ntypsd = 0
       CALL local_sym(l_write,atoms%lmaxd,atoms%lmax,sym%nop,sym%mrot,sym%tau,&
            atoms%nat,atoms%ntype,atoms%neq,cell%amat,cell%bmat,&
            atoms%taual,sphhar%nlhd,sphhar%memd,sphhar%ntypsd,.TRUE.,&
            atoms%nlhtyp,sphhar%nlh,sphhar%llh,&
            sphhar%nmem,sphhar%mlh,sphhar%clnu)
    
    DEALLOCATE(sphhar%clnu,sphhar%nlh,sphhar%llh,sphhar%nmem,sphhar%mlh)

    ALLOCATE(sphhar%clnu(sphhar%memd,0:sphhar%nlhd,sphhar%ntypsd))
    ALLOCATE(sphhar%llh(0:sphhar%nlhd,sphhar%ntypsd))
    ALLOCATE(sphhar%mlh(sphhar%memd,0:sphhar%nlhd,sphhar%ntypsd))
    ALLOCATE(sphhar%nlh(sphhar%ntypsd),sphhar%nmem(0:sphhar%nlhd,sphhar%ntypsd))

    sphhar%clnu = CMPLX(0.0,0.0)
    sphhar%llh = 0
    sphhar%mlh = 0
    sphhar%nlh = 0
    sphhar%nmem = 0

    ! Generate lattice harmonics

       CALL local_sym(l_write,atoms%lmaxd,atoms%lmax,sym%nop,sym%mrot,sym%tau,&
            atoms%nat,atoms%ntype,atoms%neq,cell%amat,cell%bmat,atoms%taual,&
            sphhar%nlhd,sphhar%memd,sphhar%ntypsd,.FALSE.,&
            atoms%nlhtyp,sphhar%nlh,sphhar%llh,sphhar%nmem,sphhar%mlh,sphhar%clnu)
       sym%nsymt = sphhar%ntypsd
       !     oneD%mrot1(:,:,:) = sym%mrot(:,:,:)
       !     oneD%tau1(:,:) = sym%tau(:,:)
  END SUBROUTINE make_sphhar
END MODULE m_make_sphhar
