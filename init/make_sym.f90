 !--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_make_sym
  USE m_judft
  IMPLICIT NONE
  PRIVATE
  PUBLIC make_sym
CONTAINS
  SUBROUTINE make_sym(sym,cell,atoms,noco,oneD,input)
    !Generates missing symmetry info.
    !tau,mrot and nop have to be specified alread
    USE m_dwigner
    USE m_mapatom
    USE m_od_mapatom
    use m_ptsym
    USE m_types_sym
    USE m_types_cell
    USE m_types_atoms
    USE m_types_noco
    USE m_types_oneD
    use m_types_input
    TYPE(t_sym),INTENT(INOUT) :: sym
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_noco),INTENT(IN)   :: noco
    TYPE(t_oneD),INTENT(INOUT):: oneD
    TYPE(t_input),INTENT(IN)  :: input

    integer :: nsymt
    integer,allocatable::nrot(:),locops(:,:)

    !Check for additional time-reversal symmetry
    IF( sym%invs .OR. noco%l_soc ) THEN
       sym%nsym = sym%nop
    ELSE
       ! combine time reversal symmetry with the spatial symmetry opera
       ! thus the symmetry operations are doubled
       sym%nsym = 2*sym%nop
    END IF

    !Generated wigner symbols for LDA+U (includes DFT+HubbardI)
    IF (ALLOCATED(sym%d_wgn)) DEALLOCATE(sym%d_wgn)
    ALLOCATE(sym%d_wgn(-3:3,-3:3,3,sym%nop))
    IF (atoms%n_u+atoms%n_hia.GT.0) THEN
       CALL d_wigner(sym%nop,sym%mrot,cell%bmat,3,sym%d_wgn)
    END IF

    !Atom specific symmetries

    allocate(locops(sym%nop,atoms%nat),nrot(atoms%nat))
    allocate(sym%ntypsy(atoms%nat))
    call ptsym(atoms%ntype,atoms%nat,atoms%neq,atoms%taual,sym%nop,sym%mrot,sym%tau,atoms%lmax,&
         nsymt,sym%ntypsy,nrot,locops)


    IF (.NOT.oneD%odd%d1) THEN
     CALL mapatom(sym,atoms,cell,input,noco)
     allocate(oneD%ngopr1(atoms%nat))
     oneD%ngopr1 = sym%ngopr
  ELSE
     CALL juDFT_error("The oneD version is broken here. Compare call to mapatom with old version")
     CALL mapatom(sym,atoms,cell,input,noco)
     !CALL od_mapatom(oneD,atoms,sym,cell)
  END IF


END SUBROUTINE make_sym
END MODULE m_make_sym
