!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------


module m_dfpt_borncharges

    use m_juDFT
    use m_constants
    use m_types

    implicit none 

contains 

    subroutine dfpt_borncharges(fi,fmpi,stars,sphhar,xcpot,forcetheo,enpara,nococonv,hybdat,rho,vTot,vxc,&
                            grRho3,grVtot3,grVc3,grVext3, results,resultsq, results1, eig_id,q_eig_id, & 
                            dfpt_eig_id,dfpt_eig_id2,l_minusq,resultsqm,results1m,qm_eig_id, dfpt_eigm_id, dfpt_eigm_id2)

        use m_eigen
        use m_fermie
        use m_dfpt_sternheimer
        use m_types_BEC
        USE m_dfpt_born_effcharge



        type(t_fleurinput),intent(in) :: fi 
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
        type(t_potden), intent(in)    :: grRho3(3),grVtot3(3),grVc3(3),grVext3(3)
        type(t_results),intent(inout)       :: results
        type(t_results),intent(inout)       :: resultsq
        type(t_results),intent(inout)       :: resultsqm
        type(t_results),intent(inout)       :: results1 , results1m



        integer, intent(in) :: eig_id, q_eig_id, dfpt_eig_id,dfpt_eig_id2,qm_eig_id,dfpt_eigm_id,dfpt_eigm_id2
        logical, intent(in) :: l_minusq

        ! scf variables 
        type(t_potden) :: den1, den1Im, vTot1, vTot1Im, vC1,vC1Im , vTot1m, vTot1mIm
        type(t_kpts) :: kqpts, kmqpts , qpts 
        type(t_stars)  :: starsq, starsmq

        ! BEC properties
        type(t_BEC) :: bec
        complex,allocatable :: born_eff_charge_contributions(:,:,:,:)
        ! helper types
        type(t_hub1data) :: hub1data

        integer, allocatable :: q_list(:)
        integer :: iQ,iArray, q_start, q_stop , ikpt ,iDir, iDir2, iDtype

        character(len=20) :: dfpt_tag
        logical :: l_real 

        !complex :: sigma_coul(2) , sigma_ext(2) ! this was a previous idea for discontinuities at the vac-ir boundary --> not used anylonger

#ifdef CPP_MPI
        integer :: ierr
#endif 

        !sigma_coul = cmplx(0.0,0.0)
        !sigma_ext = cmplx(0.0,0.0)

        l_real = fi%sym%invs.and.(.not.fi%noco%l_soc).and.(.not.fi%noco%l_noco).and.fi%atoms%n_hia==0

        ! create a kpts type that contains the necessary q-vectors 
        qpts = fi%kpts
        deallocate(qpts%bk)
        allocate(qpts%bk,mold=fi%juphon%qvec_efield) ! this is not nice. Maybe change the expected type in dfpt_sternheimer/make_stars
        qpts%bk(:, :size(fi%juPhon%qvec,2)) = fi%juPhon%qvec_efield
        allocate(q_list(size(fi%juPhon%qvec_efield,2)))  
        q_list = (/(iArray, iArray=1,SIZE(fi%juPhon%qvec_efield,2), 1)/)


        allocate(born_eff_charge_contributions(fi%atoms%nat,3,3,8+fi%atoms%nat))
        born_eff_charge_contributions = cmplx(0.0)
        call bec%init(fi)

     

        do iQ = 1, 3 ! all cartesian directions
            call timestart("q-Point")
            
            kqpts = fi%kpts
            do ikpt = 1, fi%kpts%nkpt
                kqpts%bk(:, ikpt) = kqpts%bk(:, ikpt) + qpts%bk(:,q_list(iQ))
            end do 
   
            if (l_minusq) then 
                kmqpts = fi%kpts
                do ikpt = 1, fi%kpts%nkpt
                    kmqpts%bk(:, ikpt) = kmqpts%bk(:, ikpt) - qpts%bk(:,q_list(iQ))
                end do 
            end if 

            call timestart("Eigenstuff at k+q")
            ! Get the eigenstuff at k+q
            call resultsq%reset_results(fi%input)

            call eigen(fi, fmpi, stars, sphhar, xcpot, forcetheo, enpara, nococonv,  &
                    hybdat, 1, q_eig_id, resultsq, rho, vTot, vxc, hub1data, &
                    qpts%bk(:,q_list(iQ)))

            ! Fermi level and occupancies
            call timestart("determination of fermi energy")
            call fermie(q_eig_id, fmpi, kqpts, fi%input, fi%noco, enpara%epara_min, fi%cell, resultsq)
            call timestop("determination of fermi energy")

#ifdef CPP_MPI
            call MPI_BCAST(resultsq%ef, 1, MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
            call MPI_BCAST(resultsq%w_iks, size(resultsq%w_iks), MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
#endif
            call timestop("Eigenstuff at k+q")

   
            if (l_minusq) then
                call timestart("Eigenstuff at k-q")
                ! Solve eigenvalue equation on k-q grid
                call resultsqm%reset_results(fi%input)

                call eigen(fi, fmpi, stars, sphhar, xcpot, forcetheo, enpara, nococonv,  &
                        hybdat, 1, qm_eig_id, resultsqm, rho, vTot, vxc, hub1data, &
                        -qpts%bk(:,q_list(iQ)))

                ! Fermi level and occupancies
                call timestart("determination of fermi energy")
                call fermie(qm_eig_id, fmpi, kmqpts, fi%input, fi%noco, enpara%epara_min, fi%cell, resultsqm)
                call timestop("determination of fermi energy")

#ifdef CPP_MPI
                call MPI_BCAST(resultsqm%ef, 1, MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
                call MPI_BCAST(resultsqm%w_iks, SIZE(resultsqm%w_iks), MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
#endif

                call timestop("Eigenstuff at k-q")
            end if


            do iDtype = 1 , fi%atoms%ntype
                call timestart("Typeloop Phonon")
                    do iDir = 1 , 3 
                        call timestart("Dirloop Phonon")
                        dfpt_tag = ''
                        write(dfpt_tag,'(a1,i0,a2,i0,a2,i0)') 'q', q_list(iQ), '_b', iDtype, '_j', iDir

                        if (fmpi%irank==0) then 
                            write(*,*) 'Starting calculation for:'
                            write(*,*) ' q         = ', fi%juPhon%qvec(:,q_list(iQ))
                            write(*,*) ' atom      = ', iDtype
                            write(*,*) ' direction = ', iDir
                        end if 

                        ! reset arrays 
                        call starsq%reset_stars()
                        call den1%reset_dfpt()
                        call den1Im%reset_dfpt()
                        call vTot1%reset_dfpt()
                        call vTot1Im%reset_dfpt()
                        call vC1%reset_dfpt()
                        call vC1Im%reset_dfpt
                        call results1%reset_results(fi%input)
                        if (l_minusq) then 
                            call starsmq%reset_stars()
                            call vTot1m%reset_dfpt()
                            call vTot1mIm%reset_dfpt()
                            call results1m%reset_results(fi%input)
                            ! I am unsure why there is no vC1m 
                        end if 
                            
                        if (fmpi%irank==0) write(*,*) '-------------------------'
                        ! Legacy Comment: Dont dare to delete it
                        ! This is where the magic happens. The Sternheimer equation is solved
                        ! iteratively, providing the scf part of dfpt calculations.
                        if (l_minusq) then 
                            call timestart("Sternheimer with -q")
                            call dfpt_sternheimer(fi, xcpot, sphhar, stars, starsq, nococonv, qpts, fmpi, results, resultsq, enpara, hybdat, &
                                                rho, vTot, grRho3(iDir), grVtot3(iDir), grVext3(iDir), q_list(iQ), iDtype, iDir, &
                                                dfpt_tag, eig_id, l_real, results1, dfpt_eig_id, dfpt_eig_id2, q_eig_id, &
                                                den1, vTot1, den1Im, vTot1Im, vC1, vC1Im,&
                                                starsmq, resultsqm, dfpt_eigm_id, dfpt_eigm_id2, qm_eig_id, results1m, vTot1m, vTot1mIm)
                            call timestop("Sternheimer with -q")
                        else
                            call timestart("Sternheimer")
                            call dfpt_sternheimer(fi, xcpot, sphhar, stars, starsq, nococonv, qpts, fmpi, results, resultsq, enpara, hybdat, &
                                                rho, vTot, grRho3(iDir), grVtot3(iDir), grVext3(iDir), q_list(iQ), iDtype, iDir, &
                                                dfpt_tag, eig_id, l_real, results1, dfpt_eig_id, dfpt_eig_id2, q_eig_id, &
                                                den1, vTot1, den1Im, vTot1Im, vC1, vC1Im)
                            call timestop("Sternheimer")
                        end if 

        
                        CALL dfpt_born_eff_charge_element(fi,stars,starsq,sphhar,fmpi,rho,den1,den1Im,grRho3(iDir),BEC%BEC_element(iDtype,iDir,iQ),&
                                                          born_eff_charge_contributions(iDtype,iDir,iQ,:),iDir,iDtype,iQ,1)

                        
                        if (fmpi%irank == 0 .and. fi%juphon%l_rm_qhdf) call system("rm "//trim(dfpt_tag)//".hdf")
                        call timestop("Dirloop Phonon")
                    end do ! iDir
                call timestop("Typeloop Phonon")
            end do ! iDtype
            call timestop("q-Point")
        end do ! iQ

#if defined(CPP_MPI)
        call MPI_BARRIER(fmpi%MPI_COMM,ierr)
#endif 
        ! IO of BEC properties
        if (fmpi%irank==0) call dfpt_born_eff_charge_final(fi,BEC%BEC_element,born_eff_charge_contributions(:,:,:,:))


    end subroutine dfpt_borncharges

end module m_dfpt_borncharges