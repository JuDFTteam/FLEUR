!--------------------------------------------------------------------------------
! Copyright (c) 2024 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_vgen_constraint
    IMPLICIT NONE
    CONTAINS
    subroutine vgen_constraint(atoms,noco,nococonv,vtot)
        use m_types
        TYPE(t_atoms),INTENT(in)    :: atoms
        TYPE(t_noco),INTENT(IN)     :: noco
        TYPE(t_nococonv),INTENT(IN) :: nococonv
        TYPE(t_potden),intent(INOUT):: vtot

        integer:: n
        DO n=1,atoms%ntype
            if (.not.noco%l_constrained(n)) cycle
            vtot%mt(:,0,n,3)=-0.5*nococonv%b_con(1,n)
            vtot%mt(:,0,n,4)=0.5*nococonv%b_con(2,n)
        ENDDO
    END subroutine
end module            