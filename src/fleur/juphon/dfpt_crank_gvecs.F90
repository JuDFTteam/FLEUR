!--------------------------------------------------------------------------------
! Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------


MODULE m_dfpt_crank_gvecs


    IMPLICIT NONE 


CONTAINS
    SUBROUTINE crank_gvecs(fi,fmpi,sym,cell,input,sphhar,vacuum,noco,local_stars,local_potden,local_atoms,qvec,iDir,iDtype)
        USE m_types_fleurinput
        USE m_types
        use m_constants
        USE m_types_stars
        USE m_types_atoms
        USE m_types_sym
        USE m_types_vacuum
        USE m_types_input
        USE m_types_cell
        USE m_convn
        USE m_npy
        USE m_types_mpi
        USE m_make_stars
        TYPE(t_fleurinput), INTENT(IN) :: fi
        TYPE(t_mpi), INTENT(IN) :: fmpi
        TYPE(t_sym), INTENT(IN) :: sym
        TYPE(t_cell), INTENT(IN) :: cell
        TYPE(t_input), INTENT(IN) :: input
        TYPE(t_sphhar), INTENT(IN) :: sphhar
        TYPE(t_vacuum), INTENT(IN) :: vacuum
        TYPE(t_noco), INTENT(IN)   :: noco
        TYPE(t_stars), INTENT(INOUT) :: local_stars
        TYPE(t_potden), INTENT(INOUT) :: local_potden
        TYPE(t_atoms),INTENT(OUT)   :: local_atoms
        REAL, OPTIONAL, INTENT(IN)    :: qvec(3)
        INTEGER, OPTIONAL, INTENT(IN) :: iDir,iDtype
        LOGICAL :: l_starsq



        !local_stars%gmax=12
        l_starsq=.FALSE.

        IF (PRESENT(qvec)) l_starsq=.TRUE.

        IF (l_starsq) THEN
            call make_stars(local_stars,sym,fi%atoms,vacuum,sphhar,input,cell,noco,fmpi,qvec,iDtype,iDir,gfactor=fi%juphon%gfactor)
        ELSE
            call make_stars(local_stars,sym,fi%atoms,vacuum,sphhar,input,cell,noco,fmpi,gfactor=fi%juphon%gfactor)
        END IF 
    

        local_atoms = fi%atoms

        call convn(fmpi%irank == 0, local_atoms, local_stars)
        
        call local_potden%init(local_stars,local_atoms,sphhar,vacuum,noco,input%jspins,POTDEN_TYPE_POTTOT)
    END SUBROUTINE crank_gvecs

    SUBROUTINE cast_smaller_grid(pot1,pot2,global_stars,input)
        USE m_types_input
        USE m_types

        TYPE(t_potden), INTENT(INOUT) :: pot1
        TYPE(t_potden), INTENT(IN)    :: pot2
        TYPE(t_stars), INTENT(IN)     :: global_stars
        TYPE(t_input), INTENT(IN)     :: input 

        
        pot1%pw(:,:) = pot2%pw(:global_stars%ng3,:)
        pot1%mt(:,0:,:,:) = pot2%mt(:,0:,:,:)
    
        IF ( input%film) THEN
            pot1%vac(:,:,:,:)=pot2%vac(:,:global_stars%ng2,:,:)
        END IF 

    

    END SUBROUTINE cast_smaller_grid

    SUBROUTINE cast_onto_larger_grid(pot1,pot2,global_stars,input)
        USE m_types_input
        USE m_types

        TYPE(t_potden), INTENT(INOUT) :: pot1
        TYPE(t_potden), INTENT(IN)    :: pot2
        TYPE(t_stars), INTENT(IN)     :: global_stars
        TYPE(t_input), INTENT(IN)     :: input 


        pot1%pw(:global_stars%ng3,:) = pot2%pw(:,:) 
        pot1%mt(:,0:,:,:) =  pot2%mt(:,0:,:,:) 
    
        IF ( input%film) THEN
            pot1%vac(:,:global_stars%ng2,:,:)=pot2%vac(:,:,:,:)
        END IF 

    END SUBROUTINE cast_onto_larger_grid

    SUBROUTINE copy_stars(stars1,stars2)
        USE m_types_stars

        CLASS(t_stars),INTENT(INOUT) :: stars1 
        CLASS(t_stars),INTENT(IN) :: stars2 

        CALL stars1%reset_stars()


        IF ( .NOT. ALLOCATED(stars1%kv3)) allocate(stars1%kv3,mold=stars2%kv3)
        IF ( .NOT. ALLOCATED(stars1%sk3)) allocate(stars1%sk3,mold=stars2%sk3)
        IF ( .NOT. ALLOCATED(stars1%ig)) allocate(stars1%ig,mold=stars2%ig)
        IF ( .NOT. ALLOCATED(stars1%nstr)) allocate(stars1%nstr,mold=stars2%nstr)
        IF ( .NOT. ALLOCATED(stars1%rgphs)) allocate(stars1%rgphs,mold=stars2%rgphs)
        IF ( .NOT. ALLOCATED(stars1%ustep)) allocate(stars1%ustep,mold=stars2%ustep)
        IF ( .NOT. ALLOCATED(stars1%ufft)) allocate(stars1%ufft,mold=stars2%ufft)
        IF ( .NOT. ALLOCATED(stars1%gq)) allocate(stars1%gq,mold=stars2%gq)
        IF ( .NOT. ALLOCATED(stars1%gq2)) allocate(stars1%gq2,mold=stars2%gq2)
        IF ( .NOT. ALLOCATED(stars1%ufft1)) allocate(stars1%ufft1,mold=stars2%ufft1)
        IF ( .NOT. ALLOCATED(stars1%kv2)) allocate(stars1%kv2,mold=stars2%kv2)
        IF ( .NOT. ALLOCATED(stars1%sk2)) allocate(stars1%sk2,mold=stars2%sk2)
        IF ( .NOT. ALLOCATED(stars1%nstr2)) allocate(stars1%nstr2,mold=stars2%nstr2)
        IF ( .NOT. ALLOCATED(stars1%i2g)) allocate(stars1%i2g,mold=stars2%i2g)
        IF ( .NOT. ALLOCATED(stars1%ig2)) allocate(stars1%ig2,mold=stars2%ig2)
        IF ( .NOT. ALLOCATED(stars1%igvac)) allocate(stars1%igvac,mold=stars2%igvac)
        IF ( .NOT. ALLOCATED(stars1%phi2)) allocate(stars1%phi2,mold=stars2%phi2)
        IF ( .NOT. ALLOCATED(stars1%r2gphs)) allocate(stars1%r2gphs,mold=stars2%r2gphs)
   


        stars1%gmax=stars2%gmax
        stars1%center = stars2%center


       stars1%kv3=stars2%kv3
       stars1%ng3=stars2%ng3
       stars1%sk3=stars2%sk3
       stars1%ig=stars2%ig
       stars1%nstr=stars2%nstr
       stars1%rgphs =stars2%rgphs
       stars1%ustep =stars2%ustep
       stars1%ufft =stars2%ufft
       stars1%gq =stars2%gq
       stars1%gq2=stars2%gq2
       stars1%ufft1 =stars2%ufft1
       stars1%kv2 =stars2%kv2
       stars1%ng2 =stars2%ng2
       stars1%sk2 =stars2%sk2
       stars1%nstr2 =stars2%nstr2
       stars1%i2g =stars2%i2g
       stars1%ig2 =stars2%ig2
       stars1%igvac =stars2%igvac
       stars1%phi2 =stars2%phi2
       stars1%r2gphs =stars2%r2gphs
   



    END SUBROUTINE copy_stars

END MODULE m_dfpt_crank_gvecs