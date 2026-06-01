!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_dfpt_bfield

    use m_juDFT
    use m_constants
    use m_types

    implicit none 

contains 

    subroutine dfpt_bfield(fi_ext,fmpi,stars,sphhar,xcpot,forcetheo,enpara,nococonv,hybdat,rho,vTot,vxc,&
                            grRho3,grVtot3,grVext3,results,resultsq,results1, eig_id,q_eig_id, dfpt_eig_id,dfpt_eig_id2)



        use m_eigen
        use m_fermie
        use m_dfpt_sternheimer
        use m_dfpt_magsusc


        type(t_fleurinput),intent(in) :: fi_ext 
        type(t_mpi),intent(in)        :: fmpi
        type(t_stars),intent(in)      :: stars
        type(t_sphhar),intent(in)     :: sphhar
        class(t_xcpot),intent(in)     :: xcpot
        type(t_forcetheo),intent(in)  :: forcetheo
        type(t_enpara),intent(inout)     :: enpara
        type(t_nococonv),intent(in)   :: nococonv
        type(t_hybdat),intent(inout)     :: hybdat
        type(t_potden),intent(in)     :: rho 
        type(t_potden),intent(in)     :: vTot
        type(t_potden),intent(in)     :: vxc
        type(t_potden), intent(in)    :: grRho3(3),grVtot3(3),grVext3(3)
        type(t_results),intent(inout)       :: results
        type(t_results),intent(inout)       :: resultsq
        type(t_results),intent(inout)       :: results1

        integer, intent(in) :: eig_id, q_eig_id, dfpt_eig_id,dfpt_eig_id2

        ! scf variables 
        type(t_potden) :: den1, den1Im, vTot1, vTot1Im, vC1,vC1im , vTot1m, vTot1mIm
        type(t_kpts) :: kqpts, qpts 
        type(t_stars)  :: starsq

        ! dielectric tensor 
        type(t_fleurinput) :: fi
        complex :: magnetic_susc

        ! helper types
        type(t_hub1data) :: hub1data

        integer :: iDir, ikpt
        logical :: l_real
        character(len=20) :: dfpt_tag

        !complex :: sigma_coul(2) , sigma_ext(2) ! this was a previous idea for discontinuities at the vac-ir boundary --> not used anylonger

#ifdef CPP_MPI
        integer :: ierr
#endif 
    
        !set l_phonon to true, not nice WIP
        fi= fi_ext
        fi%dfpt%l_phonon=.FALSE.


        !sigma_coul = cmplx(0.0,0.0)
        !sigma_ext  = cmplx(0.0,0.0)

        l_real = fi%sym%invs.and.(.not.fi%noco%l_soc).and.(.not.fi%noco%l_noco).and.fi%atoms%n_hia==0 !change tomorra 


        magnetic_susc = CMPLX(0,0)

        if (fmpi%irank==0) write(*,*) "Scf calculation for Zeeman field perturbation"

        do iDir = 1,1 !for all cartesian directions
            call timestart("Dirloop Bfield")

            dfpt_tag = ''
            write(dfpt_tag,'(a1,i0,a2,i0)') 'q', 1, '_j', iDir

            kqpts = fi%kpts
            do ikpt = 1, fi%kpts%nkpt
                kqpts%bk(:, ikpt) = kqpts%bk(:, ikpt) + fi%dfpt%qvec(:,1)
            end do 
            !stop
            qpts = fi%kpts
            deallocate(qpts%bk)
            allocate(qpts%bk(3,1)) ! this is not nice. Maybe change the expected type in dfpt_sternheimer/make_stars
            qpts%bk(:,1) = fi%dfpt%qvec(:,1)


            call timestart("Eigenstuff at k+q")
            ! Get the eigenstuff at k+q
            call resultsq%reset_results(fi%input)

            call eigen(fi, fmpi, stars, sphhar, xcpot, forcetheo, enpara, nococonv,  &
                    hybdat, 1, q_eig_id, resultsq, rho, vTot, vxc, hub1data, &
                    qpts%bk(:,1))

            ! Fermi level and occupancies
            call timestart("determination of fermi energy")
            call fermie(q_eig_id, fmpi, kqpts, fi%input, fi%noco, enpara%epara_min, fi%cell, resultsq)
            call timestop("determination of fermi energy")

#ifdef CPP_MPI
            call MPI_BCAST(resultsq%ef, 1, MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
            call MPI_BCAST(resultsq%w_iks, size(resultsq%w_iks), MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
#endif
            call timestop("Eigenstuff at k+q")

            if (fmpi%irank==0) then
                write(*,*) 'Starting calculation for:'
                write(*,*) ' direction = ', iDir
                write(*,*) ' q         = ', qpts%bk(:,1)
            end if 

            ! reset arrays 
            call starsq%reset_stars()
            call den1%reset_dfpt()
            call den1Im%reset_dfpt()
            call vTot1%reset_dfpt()
            call vTot1Im%reset_dfpt()
            call results1%reset_results(fi%input)

            if (fmpi%irank==0) write(*,*) '-------------------------'
            call timestart("Sternheimer")
            call dfpt_sternheimer(fi, xcpot, sphhar, stars, starsq, nococonv, qpts, fmpi, results, resultsq, enpara, hybdat, fi%dfpt, &
                                  rho, vTot, grRho3(iDir), grVtot3(iDir), grVext3(iDir), 1, 1, iDir, &
                                  dfpt_tag, eig_id, l_real, results1, dfpt_eig_id, dfpt_eig_id2, q_eig_id, &
                                  den1, vTot1, den1Im, vTot1Im, vC1, vC1Im)
            call timestop("Sternheimer")
            if (fmpi%irank==0) write(*,*) '-------------------------'  
            call dfpt_magnetic_susc(fi,stars,starsq,sphhar,fmpi,den1,den1Im,magnetic_susc)



            call timestop("Dirloop Bfield")
        end do !iDir 


! #if defined(CPP_MPI)
!    CALL MPI_BARRIER(fmpi%MPI_COMM,ierr)
! #endif  

        ! IO of results 
        !if (fmpi%irank==0) then
        !    call timestart("magnetic susceptibilty")
        !    write(*,*) "Scf calculation for Zeeman field perturbation finished"
        !    call dfpt_dielecten_final_old(fi,magnetic_susc)
        !    call timestop("magnetic susceptibilty")
        !end if 


    end subroutine dfpt_bfield


end module m_dfpt_bfield