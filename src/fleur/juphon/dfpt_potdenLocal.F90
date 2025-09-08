!--------------------------------------------------------------------------------
! Copyright (c) 2024 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------


MODULE m_dfpt_potdenLocal


    IMPLICIT NONE 


CONTAINS
    SUBROUTINE create_typesLocal(fi,fmpi,sym,cell,input,sphhar,vacuum,noco,stars,potdenLocal,atomsLocal,qvec,iDir,iDtype)
        
        ! This subroutine creates the types with a bigger Gmaxz Cutoff
        ! Nessesary for the Film-Mode Calcaultion 
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
        TYPE(t_stars), INTENT(INOUT) :: stars
        TYPE(t_potden), INTENT(INOUT) :: potdenLocal
        TYPE(t_atoms),INTENT(OUT)   :: atomsLocal
        REAL, OPTIONAL, INTENT(IN)    :: qvec(3)
        INTEGER, OPTIONAL, INTENT(IN) :: iDir,iDtype


        call make_stars(stars,sym,fi%atoms,vacuum,sphhar,input,cell,noco,fmpi,qvec,iDtype,iDir,gmaxzLocal=fi%juphon%gmaxzLocal)
        atomsLocal = fi%atoms
        call convn(fmpi%irank == 0, atomsLocal, stars)
        call potdenLocal%init(stars,atomsLocal,sphhar,vacuum,noco,input%jspins,POTDEN_TYPE_POTTOT)

    END SUBROUTINE create_typesLocal

    SUBROUTINE cast_smaller_grid(pot1,pot2,starsGlobal,input)
        
        USE m_types_input
        USE m_types

        TYPE(t_potden), INTENT(INOUT) :: pot1
        TYPE(t_potden), INTENT(IN)    :: pot2
        TYPE(t_stars), INTENT(IN)     :: starsGlobal
        TYPE(t_input), INTENT(IN)     :: input 

        
        pot1%pw(:,:) = pot2%pw(:starsGlobal%ng3,:)
        pot1%mt(:,0:,:,:) = pot2%mt(:,0:,:,:)
    
        IF ( input%film) THEN
            pot1%vac(:,:,:,:)=pot2%vac(:,:starsGlobal%ng2,:,:)
        END IF 

    END SUBROUTINE cast_smaller_grid

    SUBROUTINE cast_onto_larger_grid(pot1,pot2,starsGlobal,input)
        
        USE m_types_input
        USE m_types

        TYPE(t_potden), INTENT(INOUT) :: pot1
        TYPE(t_potden), INTENT(IN)    :: pot2
        TYPE(t_stars), INTENT(IN)     :: starsGlobal
        TYPE(t_input), INTENT(IN)     :: input 


        pot1%pw(:starsGlobal%ng3,:) = pot2%pw(:,:) 
        pot1%mt(:,0:,:,:) =  pot2%mt(:,0:,:,:) 
    
        IF ( input%film) THEN
            pot1%vac(:,:starsGlobal%ng2,:,:)=pot2%vac(:,:,:,:)
        END IF 

    END SUBROUTINE cast_onto_larger_grid

END MODULE m_dfpt_potdenLocal
