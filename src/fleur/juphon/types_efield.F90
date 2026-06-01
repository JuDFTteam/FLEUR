!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_types_efield
    use m_juDFT
    use m_types_dfpt_scf

    implicit none 


       
    private
    public :: t_efield_pert

    type, extends(t_dfpt_scf) :: t_efield_pert
        private
        complex , allocatable :: dielecTen(:,:) ! (3, 3)

    contains
        procedure :: set_dielecTen, get_dielecTen

        ! maybe add an init that allocates the private variables --> should do
        procedure :: init_child => init_child_efield 
        procedure :: q_indepent_properties => q_indepent_properties_efield
        procedure :: postprocessing_scf => postprocessing_scf_efield
        procedure :: postprocessing_qpoint => postprocessing_qpoint_efield
        procedure :: write_outfiles => write_outfiles_efield 
    end type t_efield_pert

    contains

    subroutine set_dielecTen(this,inDielecTen)

        class(t_efield_pert), intent(inout) :: this
        complex, intent(in)           :: inDielecTen(:,:)
        ! maybe add a check if sizes match --> after init the size should be relatively fixxed
        this%dielecTen(:,:) = inDielecTen(:,:)
    end subroutine set_dielecTen


    subroutine get_dielecTen(this,outDielecTen)
        class(t_efield_pert), intent(in) :: this
        complex, allocatable, intent(out) :: outDielecTen(:,:)

        if (allocated(outDielecTen)) deallocate(outDielecTen)
        allocate(outDielecTen,mold=this%dielecTen)

        outDielecTen = this%dielecTen        
    end subroutine get_dielecTen


    subroutine init_child_efield(this,fi,nqpts,dynMatNac)
        use m_types
        class(t_efield_pert), intent(inout) :: this
        type(t_fleurinput), intent(in) :: fi 
        integer, intent(in)  :: nqpts
        complex, optional, intent(in)  :: dynMatNac(:,:)
        

        allocate(this%dielecTen(3,3))

        this%dielecTen = cmplx(0.0,0.0)
        
    end subroutine init_child_efield

    subroutine q_indepent_properties_efield(this,sternheimerJob,fi,fmpi,sphhar,hybdat,xcpot,nococonv,stars,rho,vTot,grRho3,grVtot3,grVC3,grVext3,grgrVext3x3)
        
        use m_types
        

        class(t_efield_pert), intent(inout) :: this
        type(t_sternheimerjob),intent(in) :: sternheimerJob
        type(t_fleurinput), intent(in)  :: fi 
        type(t_mpi), intent(in)         :: fmpi
        type(t_stars),intent(in)      :: stars
        type(t_sphhar),intent(in)     :: sphhar
        class(t_xcpot),intent(in)     :: xcpot
        type(t_nococonv),intent(in)   :: nococonv
        type(t_hybdat),intent(inout)     :: hybdat
        type(t_potden),intent(in)     :: rho 
        type(t_potden),intent(in)     :: vTot
        type(t_potden), intent(inout) :: grRho3(3), grVtot3(3), grVC3(3), grVext3(3),grgrVext3x3(3,3)

    end subroutine q_indepent_properties_efield

    subroutine postprocessing_scf_efield(this,sternheimerJob,fi,stars,starsq,sphhar,xcpot,nococonv,hybdat,fmpi,qpts,q_list,iQ,iDtype,iDir,eig_id,dfpt_eig_id, &
                                          dfpt_eig_id2,enpara,results,results1,l_real,juphon,rho,vTot,grRho3,grVext3,grVc3,den1,vTot1,den1Im,vTot1Im,vC1,vC1Im)
        
        
        use m_types
        use m_dfpt_dielecten
        

        class(t_efield_pert),intent(inout) :: this
        type(t_sternheimerjob),intent(in) :: sternheimerJob 
        type(t_fleurinput), intent(in)  :: fi 
        type(t_stars),intent(in)      :: stars
        type(t_stars),intent(in)      :: starsq
        type(t_sphhar),intent(in)     :: sphhar
        class(t_xcpot),intent(in)     :: xcpot
        type(t_nococonv),intent(in)   :: nococonv
        type(t_hybdat),intent(inout)     :: hybdat
        type(t_mpi), intent(in)         :: fmpi
        type(t_kpts), intent(in)       :: qpts
        integer, allocatable,intent(in) :: q_list(:)
        integer, intent(in)            :: iQ,iDtype,iDir,eig_id,dfpt_eig_id,dfpt_eig_id2
        type(t_enpara),intent(inout)     :: enpara
        type(t_results),intent(inout)    :: results, results1
        logical,intent(in)               :: l_real
        type(t_juphon),intent(in)        :: juphon
        type(t_potden),intent(in)        :: rho,vTot,grRho3(3),grVext3(3),grVC3(3),den1,vTot1,den1Im,vTot1Im
        type(t_potden),intent(inout)     :: vC1,vC1Im

        
        complex,allocatable :: diel_tensor(:,:)

        
        call this%get_dielecTen(diel_tensor)


        call timestart("Dielecten row")
        call dfpt_dielecten_HF_int(fi,stars,starsq,sphhar,fmpi,den1,den1Im,diel_tensor(iDir,:),rho,iDir,1)
        call timestop("Dielecten row")
        
        call this%set_dielecTen(diel_tensor)


    end subroutine postprocessing_scf_efield


    subroutine postprocessing_qpoint_efield(this,fi,fmpi,juphon,qpts,iQ,q_list)

        use m_types
        use m_dfpt_dielecten

        class(t_efield_pert), intent(inout) :: this
        type(t_fleurinput),intent(in)  :: fi 
        type(t_mpi),intent(in)         :: fmpi
        type(t_juphon),intent(in)      :: juphon
        type(t_kpts),intent(in)        :: qpts
        integer,intent(in)             :: iQ
        integer, intent(in)            :: q_list(:)

    end subroutine postprocessing_qpoint_efield


    subroutine write_outfiles_efield(this,fi,fmpi,juPhon)

        use m_types 
        use m_dfpt_dielecten

        class(t_efield_pert),intent(inout)   :: this         
        type(t_fleurinput),intent(in) :: fi
        type(t_mpi),intent(in)        :: fmpi
        type(t_juPhon),intent(in)     :: juPhon


        complex, allocatable :: diel_tensor(:,:)


        call this%get_dielecTen(diel_tensor)


        ! IO of results 
        if (fmpi%irank==0) then
            call timestart("diel_tensor")
            if (juPhon%l_efield_scr) then
                write(*,*) "Scf calculation for screened electric field perturbation finished"
                call dfpt_dielecten_final_old(fi,diel_tensor(:,:))
            else
                write(*,*) "Scf calculation for bare electric field perturbation finished"
                call dfpt_dielecten_final_new(fi,diel_tensor(:,:))
            end if 
            call timestop("diel_tensor")
        end if 


    end subroutine write_outfiles_efield




end module m_types_efield