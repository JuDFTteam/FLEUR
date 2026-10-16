!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------



module m_dfpt_interpolation


    use m_juDFT
    use m_constants
    use m_dfpt_NAC
    use m_fleur_init
    use m_dfpt_dynmat_fourier
    use m_dfpt_dynmat_eig
    use m_dfpt_dynmat_io
    use m_make_dos
    use m_types_eigdos
    use m_types_banddos
    use m_types_dos
    use m_types_enpara
    use m_types_fleurinput
    use m_types_forcetheo
    use m_types_hybdat
    use m_types_kpts
    use m_types_mpdata
    use m_types_mpi
    use m_types_nococonv
    use m_types_misc
    use m_types_sphhar
    use m_types_stars
    use m_types_xcpot

    implicit none 
    private
    public :: dfpt_interpolation

contains 

    subroutine dfpt_interpolation(fi,fmpi,nococonv,results)


        type(t_fleurinput), intent(in) :: fi 
        type(t_mpi), intent(in)        :: fmpi
        type(t_nococonv), intent(in)   :: nococonv 
        type(t_results), intent(in)    :: results

        ! Full symmetrized type variables:
        type(t_mpi)        :: fmpi_fullsym
        type(t_fleurinput) :: fi_fullsym
        type(t_sphhar)     :: sphhar_fullsym
        type(t_stars)      :: stars_fullsym
        type(t_nococonv)   :: nococonv_fullsym
        type(t_enpara)     :: enpara_fullsym
        type(t_results)    :: results_fullsym
        type(t_hybdat)     :: hybdat_fullsym
        type(t_mpdata)     :: mpdata_fullsym

        class(t_xcpot),     allocatable :: xcpot_fullsym
        class(t_forcetheo), allocatable :: forcetheo_fullsym

        type(t_kpts)  :: qpts
        integer, allocatable :: q_list(:)
        
        integer :: iQ, iSpecial

        character(len=100) :: inp_pref
        character(len=50)  :: label
        character(len=:), allocatable :: dynfilename
        logical :: l_sumrule
        real    :: weight

        ! dynMat properties
        complex, allocatable   :: eigenFreqs(:), eigenVecs(:,:)
        complex, allocatable   :: dyn_mat(:,:,:), dyn_mat_interp(:,:,:), dyn_mat_gamma(:,:,:)
        real,    allocatable   :: eigenVals(:), eigenValsFull(:,:,:)


        ! helper types
        type(t_banddos)                 :: banddosLocal
        type(t_eigdos_list),allocatable :: eigdos(:)
        type(t_dos), target             :: dos


        ! If the Dynmats-Files were already created, we can read them in and do postprocessing.
        ! a) Transform the q-Mesh onto real space.
        ! b) Transform it back onto a dense q-path.
        ! c) Transform it back to a denser grid
        ! d) Perform a DOS calculation for the denser grid.


        ! determine the Job 
        if (fi%dfpt%l_band) then 
            dynfilename = "dynmat_band.hdf"
        else if (fi%dfpt%l_dos) then 
            dynfilename = "dynmat_dos.hdf"
        else
            call juDFT_error("Invalid option for post-procession of dynMatfiles.", calledby="dfpt_interpolation.F90")
        end if 

        ! Read in the dynMats in the IBZ
        ! We currently rely on a seperate set of files fullsym_* that contain the symmetry
        ! Maybe we should change this in the future.
        ! Read qpoints from fullsym_inp.xml and fullsym_kpts.xml
        inp_pref = ADJUSTL("fullsym_")
        fmpi_fullsym%l_mpi_multithreaded = fmpi%l_mpi_multithreaded
        fmpi_fullsym%mpi_comm = fmpi%mpi_comm
        ! Skip setupMPI: this fleur_init only re-reads the fullsym_ input set; the
        ! parallel-solver setup it would otherwise do is never used here and fails
        ! when the fullsym q-mesh does not factor evenly onto the MPI ranks.
        call fleur_init(fmpi_fullsym, fi_fullsym, sphhar_fullsym, stars_fullsym, nococonv_fullsym, forcetheo_fullsym, &
                        enpara_fullsym, xcpot_fullsym, results_fullsym, hybdat_fullsym, mpdata_fullsym, &
                        inp_pref, l_skip_setupmpi=.true.)
        qpts = fi_fullsym%kpts


        IF (fi%dfpt%l_dos) allocate(eigenValsFull(3*fi%atoms%nat,fi%dfpt%qpts_interpol%nkpt,fi%input%jspins))


        ALLOCATE(q_list(SIZE(qpts%bk,2)))
        q_list = (/(iQ, iQ=1,SIZE(qpts%bk,2), 1)/)

        ALLOCATE(dyn_mat(3*fi%atoms%ntype,3*fi%atoms%ntype,SIZE(q_list)))

        dyn_mat = cmplx(0.0,0.0)


        if (fmpi%irank==0) then 
            do iQ = 1, qpts%nkpt
                call read_dynmat_q("dynmat.hdf", iQ, dyn_mat(:,:,iQ))
            end do

            !subtract Long range part
            !if (fi%dfpt%l_polar) then
            !    do iQ = 1, qpts%nkpt
            !        dyn_mat_NAC_q = cmplx(0.0,0.0)
            !        call get_NAC_ewald(fi,qpts,stars_fullsym,dyn_mat_NAC_q,qpts%bk(:,iQ),iQ)
            !        dyn_mat(:,:,iQ) = dyn_mat(:,:,iQ) - dyn_mat_NAC_q
            !    end do
            !end if

            ! Fourier-interpolate the coarse dynamical matrices onto the fine q-mesh
            call interpolate_dynmat(fi_fullsym%atoms, fi_fullsym%sym, fi%cell, qpts, dyn_mat, &
                                    fi%dfpt%l_WSinterpol, fi%dfpt%qpts_interpol%bk, dyn_mat_interp, &
                                    l_bornhuang=fi%dfpt%l_bornhuang)

            ! The projection already imposes the translational rule on the force constants.
            l_sumrule = fi%dfpt%l_sumrule_intp .and. .not. fi%dfpt%l_bornhuang
            if (l_sumrule) then
                ! The Gamma matrix is needed before the first fine q-point
                call interpolate_dynmat(fi_fullsym%atoms, fi_fullsym%sym, fi%cell, qpts, dyn_mat, &
                                        fi%dfpt%l_WSinterpol, reshape([0.0, 0.0, 0.0], [3, 1]), dyn_mat_gamma, &
                                        l_bornhuang=fi%dfpt%l_bornhuang)
                call hermitize_scale_dynmat(fi%atoms, dyn_mat_gamma(:,:,1))
            end if

            call init_dynmat_file(dynfilename, fi%atoms, fi%cell, fi%dfpt%qpts_interpol, l_truncate=.true.)

            ! diagonalize the interpolated dynamical matrix at each fine q-point
            do iQ = 1, fi%dfpt%qpts_interpol%nkpt
                call timestart("Dynmat diagonalization")
                call hermitize_scale_dynmat(fi%atoms, dyn_mat_interp(:,:,iQ))
                if (l_sumrule) call apply_sumrule(dyn_mat_interp(:,:,iQ), dyn_mat_gamma(:,:,1))
                call DiagonalizeDynMat(fi%atoms, fi%dfpt%qpts_interpol%bk(:,iQ), dyn_mat_interp(:,:,iQ), eigenVals, eigenVecs, .not.fi%dfpt%l_dos)
                call timestop("Dynmat diagonalization")

                call timestart("Frequency calculation")
                call CalculateFrequencies(fi%atoms, iQ, eigenVals, eigenFreqs, fi%dfpt%qpts_interpol%bk(:,iQ), l_writeOutput=(.not.fi%dfpt%l_dos))
                call timestop("Frequency calculation")

                call timestart("dynMat IO")
                weight = 0.0
                if (allocated(fi%dfpt%qpts_interpol%wtkpt)) weight = fi%dfpt%qpts_interpol%wtkpt(iQ)
                if (fi%dfpt%l_dos) then
                    call write_dynmat_q(dynfilename, iQ, fi%dfpt%qpts_interpol%bk(:,iQ), weight, eigenVals, eigenFreqs)
                else
                    label = ''
                    if (allocated(fi%dfpt%qpts_interpol%specialPointIndices)) then
                        do iSpecial = 1, fi%dfpt%qpts_interpol%numSpecialPoints
                            if (fi%dfpt%qpts_interpol%specialPointIndices(iSpecial) == iQ) label = fi%dfpt%qpts_interpol%specialPointNames(iSpecial)
                        end do
                    end if
                    call write_dynmat_q(dynfilename, iQ, fi%dfpt%qpts_interpol%bk(:,iQ), weight, eigenVals, eigenFreqs, dynMat=dyn_mat_interp(:,:,iQ), eigVecs=eigenVecs, label=label)
                end if
                call timestop("dynMat IO")

                if (fi%dfpt%l_dos) eigenValsFull(:,iQ,1) = sqrt(eigenVals(:)) * 1e3 ! factor so banddos saves in meV 

                deallocate(eigenVals, eigenVecs, eigenFreqs)
            end do ! iQ

            if (fi%dfpt%l_dos) then 
                banddosLocal = fi%banddos 
                banddosLocal%dos = .true.
                banddosLocal%e1_dos = minval(eigenValsFull)
                banddosLocal%e2_dos = maxval(eigenValsFull)
                call dos%init(fi%input,fi%atoms,fi%dfpt%qpts_interpol,banddosLocal,.false.,eigenValsFull)
                allocate(eigdos(1))
                eigdos(1)%p=>dos 
                call make_dos(fi%dfpt%qpts_interpol,fi%atoms,fi%vacuum,fi%input,banddosLocal,fi%sliceplot,fi%noco,nococonv,fi%sym,fi%cell,results,eigdos,fi%dfpt)
            end if 
        end if !irank==0 


    end subroutine dfpt_interpolation 


end module m_dfpt_interpolation