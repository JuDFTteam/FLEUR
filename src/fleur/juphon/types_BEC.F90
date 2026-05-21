!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_types_BEC
    use m_juDFT
    use m_types_dfpt

    implicit none 


       
    private
    public :: t_BEC

    type, extends(t_dfpt) :: t_BEC
        private
        complex , allocatable :: borneffcharge(:,:,:),borneffcharge_ctrb(:,:,:,:) 

    contains
        procedure :: set_BEC, get_BEC

        ! maybe add an init that allocates the private variables --> should do
        procedure :: init_child => init_child_BEC 
        procedure :: q_indepent_properties => q_indepent_properties_BEC
        procedure :: postprocessing_scf => postprocessing_scf_BEC
        procedure :: postprocessing_qpoint => postprocessing_qpoint_BEC
        procedure :: write_outfiles => write_outfiles_BEC 
    end type t_BEC

    contains

    subroutine set_BEC(this,inborneffcharge,inborneffcharge_ctrb)

        class(t_BEC), intent(inout) :: this
        complex, intent(in)         :: inborneffcharge(:,:,:),inborneffcharge_ctrb(:,:,:,:)
        ! maybe add a check if sizes match --> after init the size should be relatively fixxed
        this%borneffcharge(:,:,:) = inborneffcharge(:,:,:)
        this%borneffcharge_ctrb(:,:,:,:) = inborneffcharge_ctrb(:,:,:,:)
    end subroutine set_BEC


    subroutine get_BEC(this,outborneffcharge,outborneffcharge_ctrb)
        class(t_BEC), intent(in) :: this
        complex, allocatable, intent(out) :: outborneffcharge(:,:,:),outborneffcharge_ctrb(:,:,:,:)

        if (allocated(outborneffcharge)) deallocate(outborneffcharge)
        if (allocated(outborneffcharge_ctrb)) deallocate(outborneffcharge_ctrb)
        allocate(outborneffcharge,mold=this%borneffcharge)
        allocate(outborneffcharge_ctrb,mold=this%borneffcharge_ctrb)

        outborneffcharge = this%borneffcharge        
        outborneffcharge_ctrb = this%borneffcharge_ctrb     
    end subroutine get_BEC


subroutine init_child_BEC(this,fi,nqpts,dynMatNac)
        use m_types
        class(t_BEC), intent(inout) :: this
        type(t_fleurinput), intent(in) :: fi 
        integer, intent(in)  :: nqpts
        complex, optional, intent(in)  :: dynMatNac(:,:)
        

        allocate(this%borneffcharge(fi%atoms%ntype,3,3))
        allocate(this%borneffcharge_ctrb(fi%atoms%ntype,3,3,8+fi%atoms%nat))

        this%borneffcharge = cmplx(0.0,0.0)
        this%borneffcharge_ctrb = cmplx(0.0,0.0)
        
    end subroutine init_child_BEC

    subroutine q_indepent_properties_BEC(this,sternheimerJob,fi,fmpi,sphhar,hybdat,xcpot,nococonv,stars,rho,vTot,grRho3,grVtot3,grVC3,grVext3,grgrVext3x3)
        
        use m_types
        
        use m_dfpt_generate_gradient

        class(t_BEC), intent(inout) :: this
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


        ! calculate the gradients as this is currently a necessary input for the sternheimer in any perturbation 

        ! Generate the gradients of the density and the various potentials, that will be used at different points in the programm.
        ! The density gradient is calculated by numerical differentiation, while the potential gradients are constructed (from the
        ! density gradient) by a Weinert construction, just like the potentials are from the density.
        ! This is done to ensure good continuity.
        call timestart("Gradient generation")
        call dfpt_generate_gradient(sternheimerJob,fi,fmpi,sphhar,hybdat,xcpot,nococonv,stars,rho,vTot,grRho3,grVtot3,grVC3,grVext3,grgrVext3x3)
        call timestop("Gradient generation")


    end subroutine q_indepent_properties_BEC

    subroutine postprocessing_scf_BEC(this,sternheimerJob,fi,stars,starsq,sphhar,xcpot,nococonv,hybdat,fmpi,qpts,q_list,iQ,iDtype,iDir,eig_id,dfpt_eig_id, &
                                          dfpt_eig_id2,enpara,results,results1,l_real,juphon,rho,vTot,grRho3,grVext3,grVc3,den1,vTot1,den1Im,vTot1Im,vC1,vC1Im)
        
        
        use m_types
        use m_dfpt_born_effcharge
        

        class(t_BEC),intent(inout) :: this
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

        
        complex,allocatable :: borneffcharge(:,:,:),borneffcharge_ctrb(:,:,:,:)

        
        call this%get_BEC(borneffcharge,borneffcharge_ctrb)
        !call this%get_BEC(borneffcharge_ctrb)


        call timestart("BEC element")
        CALL dfpt_born_eff_charge_element(fi,stars,starsq,sphhar,fmpi,rho,den1,den1Im,grRho3(iDir),borneffcharge(iDtype,iDir,iQ),&
                                                          borneffcharge_ctrb(iDtype,iDir,iQ,:),iDir,iDtype,iQ,1)
        call timestop("BEC element")
        
        call this%set_BEC(borneffcharge,borneffcharge_ctrb)
        !call this%set_BEC(borneffcharge_ctrb)


    end subroutine postprocessing_scf_BEC


    subroutine postprocessing_qpoint_BEC(this,fi,fmpi,juphon,qpts,iQ,q_list)

        use m_types
        use m_dfpt_dielecten

        class(t_BEC), intent(inout) :: this
        type(t_fleurinput),intent(in)  :: fi 
        type(t_mpi),intent(in)         :: fmpi
        type(t_juphon),intent(in)      :: juphon
        type(t_kpts),intent(in)        :: qpts
        integer,intent(in)             :: iQ
        integer, intent(in)            :: q_list(:)

    end subroutine postprocessing_qpoint_BEC


    subroutine write_outfiles_BEC(this,fi,fmpi,juPhon)

        use m_types 
        use m_dfpt_born_effcharge

        class(t_BEC),intent(inout)   :: this         
        type(t_fleurinput),intent(in) :: fi
        type(t_mpi),intent(in)        :: fmpi
        type(t_juPhon),intent(in)     :: juPhon


        complex, allocatable :: borneffcharge(:,:,:),borneffcharge_ctrb(:,:,:,:)


        call this%get_BEC(borneffcharge,borneffcharge_ctrb)
        !call this%get_BEC(borneffcharge_ctrb)


        ! IO of results 
        if (fmpi%irank==0) then
            call timestart("borneffcharge")
            call dfpt_born_eff_charge_final(fi,borneffcharge,borneffcharge_ctrb)
            call timestop("borneffcharge")
        end if 


    end subroutine write_outfiles_BEC




end module m_types_BEC