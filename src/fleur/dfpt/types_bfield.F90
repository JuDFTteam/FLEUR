!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_types_bfield
    use m_juDFT
    use m_types_dfpt_scf
    use m_types
    use m_dfpt_magsusc

    implicit none

    private
    public :: t_bfield

    type, extends(t_dfpt_scf) :: t_bfield
        private
        complex :: magnetic_susc

    contains
        procedure :: set_magnetic_susc, get_magnetic_susc
        procedure :: init_child => init_child_bfield
        procedure :: q_indepent_properties => q_indepent_properties_bfield
        procedure :: postprocessing_scf => postprocessing_scf_bfield
        procedure :: postprocessing_qpoint => postprocessing_qpoint_bfield
        procedure :: write_outfiles => write_outfiles_bfield
    end type t_bfield

contains

    subroutine set_magnetic_susc(this, in_susc)
        class(t_bfield), intent(inout) :: this
        complex, intent(in)            :: in_susc

        this%magnetic_susc = in_susc
    end subroutine set_magnetic_susc

    subroutine get_magnetic_susc(this, out_susc)
        class(t_bfield), intent(in) :: this
        complex, intent(out)        :: out_susc

        out_susc = this%magnetic_susc
    end subroutine get_magnetic_susc

    subroutine init_child_bfield(this,fi,nqpts,dynMatNac)
        use m_types
        class(t_bfield), intent(inout) :: this
        type(t_fleurinput), intent(in) :: fi
        integer, intent(in)            :: nqpts
        complex, optional, intent(in)  :: dynMatNac(:,:)

        this%magnetic_susc = cmplx(0.0,0.0)

    end subroutine init_child_bfield

    subroutine q_indepent_properties_bfield(this,sternheimerJob,fi,fmpi,sphhar,hybdat,xcpot,nococonv,stars,rho,vTot,grRho3,grVtot3,grVC3,grVext3,grgrVext3x3)
        use m_types
        class(t_bfield), intent(inout) :: this
        type(t_sternheimerjob),intent(in) :: sternheimerJob
        type(t_fleurinput), intent(in)  :: fi
        type(t_mpi), intent(in)         :: fmpi
        type(t_stars),intent(in)        :: stars
        type(t_sphhar),intent(in)       :: sphhar
        class(t_xcpot),intent(in)       :: xcpot
        type(t_nococonv),intent(in)     :: nococonv
        type(t_hybdat),intent(inout)    :: hybdat
        type(t_potden),intent(in)       :: rho
        type(t_potden),intent(in)       :: vTot
        type(t_potden), intent(inout)      :: grRho3(3), grVtot3(3), grVC3(3), grVext3(3), grgrVext3x3(3,3)
    end subroutine q_indepent_properties_bfield

    subroutine postprocessing_scf_bfield(this,sternheimerJob,fi,stars,starsq,sphhar,xcpot,nococonv,hybdat,fmpi,qpts,q_list,iQ,iDtype,iDir,eig_id,dfpt_eig_id, &
                                          dfpt_eig_id2,enpara,results,results1,l_real,dfpt,rho,vTot,grRho3,grVext3,grVC3,den1,vTot1,den1Im,vTot1Im,vC1,vC1Im)
        use m_types
        class(t_bfield),intent(inout) :: this
        type(t_sternheimerjob),intent(in) :: sternheimerJob
        type(t_fleurinput), intent(in)  :: fi
        type(t_stars),intent(in)       :: stars, starsq
        type(t_sphhar),intent(in)      :: sphhar
        class(t_xcpot),intent(in)      :: xcpot
        type(t_nococonv),intent(in)    :: nococonv
        type(t_hybdat),intent(inout)   :: hybdat
        type(t_mpi), intent(in)        :: fmpi
        type(t_kpts), intent(in)       :: qpts
        integer, allocatable,intent(in) :: q_list(:)
        integer, intent(in)            :: iQ, iDtype, iDir, eig_id, dfpt_eig_id, dfpt_eig_id2
        type(t_enpara),intent(inout)   :: enpara
        type(t_results),intent(inout)  :: results, results1
        logical, intent(in)            :: l_real
        type(t_dfpt),intent(in)      :: dfpt
        type(t_potden),intent(in)      :: rho, vTot, grRho3(3), grVext3(3), grVC3(3), den1, vTot1, den1Im, vTot1Im
        type(t_potden),intent(inout)   :: vC1, vC1Im

        complex :: magnetic_susc_local

        call dfpt_magnetic_susc(fi,stars,starsq,sphhar,fmpi,den1,den1Im,magnetic_susc_local)
        this%magnetic_susc = this%magnetic_susc + magnetic_susc_local
    end subroutine postprocessing_scf_bfield

    subroutine postprocessing_qpoint_bfield(this,fi,fmpi,dfpt,qpts,iQ,q_list)
        use m_types
        class(t_bfield), intent(inout) :: this
        type(t_fleurinput),intent(in) :: fi
        type(t_mpi),intent(in)        :: fmpi
        type(t_dfpt),intent(in)     :: dfpt
        type(t_kpts),intent(in)       :: qpts
        integer,intent(in)            :: iQ
        integer, intent(in)           :: q_list(:)
    end subroutine postprocessing_qpoint_bfield

    subroutine write_outfiles_bfield(this,fi,fmpi,dfpt)
        use m_types
        class(t_bfield),intent(inout)   :: this
        type(t_fleurinput),intent(in)  :: fi
        type(t_mpi),intent(in)         :: fmpi
        type(t_dfpt),intent(in)      :: dfpt

        complex :: magnetic_susc

        call this%get_magnetic_susc(magnetic_susc)
        if (fmpi%irank == 0) then
            write(*,*) 'Scf calculation for Zeeman field perturbation finished'
            write(*,*) 'Magnetic susceptibility in a.u. per unit cell: ', magnetic_susc
        end if
    end subroutine write_outfiles_bfield

end module m_types_bfield
