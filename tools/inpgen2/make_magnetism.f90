 !--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_make_magnetism
    implicit none
    contains
    subroutine make_magnetism(input,noco,atoms,mag_mom)
        use m_types_input
        use m_types_noco
        use m_types_atoms
        use m_polangle
        TYPE(t_input),INTENT(INOUT):: input
        TYPE(t_noco),INTENT(INOUT) :: noco 
        TYPE(t_atoms),INTENT(INOUT):: atoms 
        REAL,INTENT(IN)            :: mag_mom(0:,:)

        INTEGER:: n 

        !Magnetic defaults for noco
        ALLOCATE ( noco%alph_inp(atoms%ntype), noco%beta_inp(atoms%ntype))
        noco%alph_inp(:) = 0.0
        noco%beta_inp(:) = 0.0
            
        noco%qss_inp = MERGE(noco%qss_inp, [0.0, 0.0, 0.0], noco%l_ss)
        ALLOCATE(noco%l_constrained(atoms%ntype))
        noco%l_constrained(:) = .FALSE.
        ALLOCATE(noco%l_unrestrictMT(atoms%ntype))
        noco%l_unrestrictMT(:) = .FALSE.
        ALLOCATE(noco%l_alignMT(atoms%ntype))
        noco%l_alignMT(:)=.false.
        ALLOCATE(noco%mix_RelaxWeightOffD(atoms%ntype))
        noco%mix_RelaxWeightOffD(:)=1.0
        noco%mag_mixing_scheme=0
        Allocate(atoms%bmu(atoms%ntype))
        atoms%bmu=0.0

        if (all(abs(mag_mom)<1E-5)) RETURN !No magnetic moments given
        input%jspins=2
        !check if we have a collinear setup
        if (all(abs(mag_mom(2:,:))<1E-5).and.all(mag_mom(0,:)==0)) THEN
            atoms%bmu(:)=mag_mom(1,:)
            RETURN
        endif

        !Now do the noco-setup
        noco%l_noco=.true. 
        
        DO n=1,atoms%ntype
            if (mag_mom(0,n)==0) THEN
                atoms%bmu(n)=sqrt(dot_product(mag_mom(:,n),mag_mom(:,n))) !set mag moment to absolute value
                call pol_angle(mag_mom(1,n),mag_mom(2,n),mag_mom(3,n), noco%beta_inp(n), noco%alph_inp(n),.true.)
            else !angle are given
                noco%alph_inp(n)=mag_mom(1,n)
                noco%beta_inp(n)=mag_mom(2,n)
                atoms%bmu(n)=mag_mom(3,n)
            endif      
        ENDDO
    END subroutine
END module            