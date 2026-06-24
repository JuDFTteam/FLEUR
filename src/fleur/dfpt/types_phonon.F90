!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------


module m_types_phonon

    use m_juDFT
    use m_types_dfpt_scf

    implicit none 


       
    private
    public :: t_phonon

    type, extends(t_dfpt_scf) :: t_phonon
        private
        complex , allocatable :: dynMat(:,:,:) ! (q-Point , 3*natoms, 3*natoms)
        complex , allocatable :: dynMatNac(:,:) ! (q-Point , 3*natoms, 3*natoms) NAC correction LO-TO 
        complex, allocatable  :: Eii2(:,:) ! q = 0 of ion-ion interaction (3*natoms, 3*natoms)

    contains
        procedure :: set_Eii2, get_Eii2
        procedure :: set_dynMat, get_dynMat
        procedure :: get_dynNAC

        ! maybe add an init that allocates the private variables --> should do
        procedure :: init_child => init_child_phonon 
        procedure :: q_indepent_properties => q_indepent_properties_phonon
        procedure :: postprocessing_scf => postprocessing_scf_phonon
        procedure :: postprocessing_qpoint => postprocessing_qpoint_phonon
        procedure :: write_outfiles => write_outfiles_phonon
    end type t_phonon


    contains

    subroutine set_dynMat(this,inDynMat)

        class(t_phonon), intent(inout) :: this
        complex, intent(in)           :: inDynMat(:,:,:)
        ! maybe add a check if sizes match --> after init the size should be relatively fixxed
        this%dynMat(:,:,:) = inDynMat(:,:,:)
    end subroutine set_dynMat


    subroutine get_dynMat(this,outDynMat)
        class(t_phonon), intent(in) :: this
        complex, allocatable, intent(out) :: outDynmat(:,:,:)

        if (allocated(outDynmat)) deallocate(outDynmat)
        allocate(outDynmat,mold=this%dynMat)

        outDynmat = this%dynMat        
    end subroutine get_dynMat


    subroutine set_Eii2(this,inEii2)
        class(t_phonon), intent(inout) :: this
        complex, intent(in)           :: inEii2(:,:)
        ! maybe add a check if sizes match --> after init the size should be relatively fixxed
        this%Eii2(:,:) = inEii2(:,:)
    end subroutine set_Eii2

    
    subroutine get_Eii2(this,outEii2)
        class(t_phonon), intent(in) :: this
        complex, allocatable, intent(out) :: outEii2(:,:)

        if (allocated(outEii2)) deallocate(outEii2)
        allocate(outEii2,mold=this%Eii2)

        outEii2 = this%Eii2        
    end subroutine get_Eii2
    
    subroutine get_dynNAC(this,outDynNAC)
        class(t_phonon), intent(in) :: this
        complex, allocatable, intent(out) :: outDynNAC(:,:)

        if (allocated(outDynNAC)) deallocate(outDynNAC)
        allocate(outDynNAC,mold=this%dynMatNac)

        outDynNAC = this%dynMatNac        
    end subroutine get_dynNAC


    subroutine init_child_phonon(this,fi,nqpts,dynMatNac)
        use m_types
        class(t_phonon), intent(inout) :: this
        type(t_fleurinput), intent(in) :: fi 
        integer, intent(in)  :: nqpts
        complex, optional, intent(in)  :: dynMatNac(:,:)
        

        allocate(this%dynMat(nqpts,3*fi%atoms%nat,3*fi%atoms%nat))
        allocate(this%Eii2(3*fi%atoms%nat,3*fi%atoms%nat))
        allocate(this%dynMatNac(3*fi%atoms%nat,3*fi%atoms%nat))

        this%dynMat = cmplx(0.0,0.0)
        this%Eii2 = cmplx(0.0,0.0)
        this%dynMatNac = cmplx(0.0,0.0)

    end subroutine init_child_phonon



    subroutine q_indepent_properties_phonon(this,sternheimerJob,fi,fmpi,sphhar,hybdat,xcpot,nococonv,stars,rho,vTot,grRho3,grVtot3,grVC3,grVext3,grgrVext3x3)
        
        use m_types
        use m_dfpt_eii2    
        
        use m_dfpt_generate_gradient


        class(t_phonon), intent(inout) :: this
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

        ! local variables 
        type(t_potden) :: potdummy
        real :: e2_vm(fi%atoms%nat,3,3)
        complex :: E2ndOrdII(3*fi%atoms%nat,3*fi%atoms%nat) 
        integer :: iDir, iDir2, iDtype


        E2ndOrdII = cmplx(0.0,0.0)
        e2_vm = 0.0

        ! Generate the gradients of the density and the various potentials, that will be used at different points in the programm.
        ! The density gradient is calculated by numerical differentiation, while the potential gradients are constructed (from the
        ! density gradient) by a Weinert construction, just like the potentials are from the density.
        ! This is done to ensure good continuity.
        call timestart("Gradient generation")
        call dfpt_generate_gradient(sternheimerJob,fi,fmpi,sphhar,hybdat,xcpot,nococonv,stars,rho,vTot,grRho3,grVtot3,grVC3,grVext3,grgrVext3x3)
        call timestop("Gradient generation")


        
        call potdummy%copyPotDen(rho)
        call potdummy%resetPotDen()

        call timestart("Eii2 q=0")
        do iDir = 1, 3 
            do iDir2 = 1, 3 
                call dfpt_e2_madelung(fi%atoms,fi%input%jspins,potdummy%mt(:,0,:,:),grgrVext3x3(iDir2,iDir)%mt(:,0,:,1),e2_vm(:,iDir2,iDir))
            end do !iDir2
        end do !iDir 
        call timestop("Eii2 q=0")

        ! map the q=0 part of Eii(2) to the corresponding element in the dynMat shape 
        do iDtype = 1, fi%atoms%ntype
            do iDir2 = 1, 3
                do iDir = 1, 3
                    E2ndOrdII(3*(iDtype-1)+iDir2,3*(iDtype-1)+iDir) = e2_vm(iDtype,iDir2,iDir)
                end do 
            end do 
        end do 

        ! write to type 
        call this%set_Eii2(E2ndOrdII)
    
    end subroutine q_indepent_properties_phonon


    subroutine postprocessing_scf_phonon(this,sternheimerJob,fi,stars,starsq,sphhar,xcpot,nococonv,hybdat,fmpi,qpts,q_list,iQ,iDtype,iDir,eig_id,dfpt_eig_id, &
                                          dfpt_eig_id2,enpara,results,results1,l_real,dfpt,rho,vTot,grRho3,grVext3,grVc3,den1,vTot1,den1Im,vTot1Im,vC1,vC1Im)
        
        
        use m_types
        use m_dfpt_dynmat
        

        class(t_phonon),intent(inout) :: this
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
        type(t_dfpt),intent(in)        :: dfpt
        type(t_potden),intent(in)        :: rho,vTot,grRho3(3),grVext3(3),grVC3(3),den1,vTot1,den1Im,vTot1Im
        type(t_potden),intent(inout)     :: vC1,vC1Im

        
        complex,allocatable :: E2ndOrdII(:,:)
        complex,allocatable :: dyn_mat(:,:,:)

        
        call this%get_Eii2(E2ndOrdII)
        call this%get_dynMat(dyn_mat)

        call timestart("Dynmat row")
        call dfpt_dynmat_row(sternheimerJob, fi, stars, starsq, sphhar, xcpot, nococonv, hybdat, fmpi, qpts, q_list(iQ), iDtype, iDir, &
                                eig_id, dfpt_eig_id, dfpt_eig_id2, enpara, results, results1, l_real, dfpt, &
                                rho, vTot, grRho3, grVext3, grVC3, &
                                den1, vTot1, den1Im, vTot1Im, vC1, vC1Im, dyn_mat(iQ,3 *(iDtype-1)+iDir,:), E2ndOrdII)

        call timestop("Dynmat row")
        
        dyn_mat(iQ,3 *(iDtype-1)+iDir,:) = dyn_mat(iQ,3 *(iDtype-1)+iDir,:) + conjg(E2ndOrdII(3 *(iDtype-1)+iDir,:))

        call this%set_dynMat(dyn_mat)


    end subroutine postprocessing_scf_phonon

    subroutine postprocessing_qpoint_phonon(this,fi,fmpi,dfpt,qpts,iQ,q_list)

        use m_types
        use m_dfpt_dynmat_eig

        class(t_phonon), intent(inout) :: this
        type(t_fleurinput),intent(in)  :: fi 
        type(t_mpi),intent(in)         :: fmpi
        type(t_dfpt),intent(in)        :: dfpt
        type(t_kpts),intent(in)        :: qpts
        integer,intent(in)             :: iQ
        integer, intent(in)            :: q_list(:)

        complex,allocatable :: dyn_mat(:,:,:)

        real, allocatable    :: eigenVals(:)
        complex, allocatable :: eigenFreqs(:), eigenVecs(:,:),dynMatNac(:,:)
        logical :: l_gamma

        l_gamma = (norm2(qpts%bk(:,q_list(iQ))) .lt. 1e-8)

        call this%get_dynMat(dyn_mat)

        if (fmpi%irank==0) then 
            ! Add NAC contribution
            if (l_gamma .and. dfpt%l_polar) then 
                call this%get_dynNAC(dynMatNac)
                dyn_mat(iQ,:,:) = dyn_mat(iQ,:,:) + dynMatNac(:,:) 
            end if 

            call timestart("Dynmat diagonalization")
            call DiagonalizeDynMat(fi%atoms, qpts%bk(:,q_list(iQ)), dfpt%calcEigenVec, dyn_mat(iQ,:,:), eigenVals, eigenVecs, q_list(iQ),.TRUE.,"raw",dfpt%l_sumrule_scf,l_writeOutput=.true.)
            call timestop("Dynmat diagonalization")

            call timestart("Frequency calculation")
            call CalculateFrequencies(fi%atoms, q_list(iQ), eigenVals, eigenFreqs,"raw",qpts%bk(:,q_list(iQ)))
            call timestop("Frequency calculation")
        end if 

    end subroutine postprocessing_qpoint_phonon

    subroutine write_outfiles_phonon(this,fi,fmpi,dfpt)

        use m_types 

        class(t_phonon),intent(inout)   :: this         
        type(t_fleurinput),intent(in) :: fi
        type(t_mpi),intent(in)        :: fmpi
        type(t_dfpt),intent(in)        :: dfpt

    end subroutine write_outfiles_phonon




end module m_types_phonon