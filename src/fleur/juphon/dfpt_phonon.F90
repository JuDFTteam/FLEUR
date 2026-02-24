!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_dfpt_phonon

    use m_juDFT
    use m_constants
    use m_types

    implicit none 

contains 

    subroutine dfpt_phonon(fi,fmpi,stars,sphhar,xcpot,forcetheo,enpara,nococonv,hybdat,rho,vTot,vxc,&
                            grRho3,grVtot3,grVc3,grVext3,grgrVext3x3, results,resultsq, results1, eig_id,q_eig_id, & 
                            dfpt_eig_id,dfpt_eig_id2,l_minusq,resultsqm,results1m,qm_eig_id, dfpt_eigm_id, dfpt_eigm_id2)

        use m_dfpt_eii2
        use m_eigen
        use m_fermie
        use m_dfpt_sternheimer
        use m_dfpt_dynmat
        use m_dfpt_dynmat_eig


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
        type(t_potden), intent(in)    :: grRho3(3),grVtot3(3),grVc3(3),grVext3(3),grgrVext3x3(3,3)
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

        ! dynMat properties
        complex, allocatable :: dyn_mat(:,:,:)
        real, allocatable    :: eigenVals(:)
        complex, allocatable :: eigenFreqs(:), eigenVecs(:,:)

        ! helper types
        type(t_potden) :: imagrhodummy 
        type(t_hub1data) :: hub1data

        integer, allocatable :: q_list(:)
        integer :: iQ,iArray, q_start, q_stop , ikpt ,iDir, iDir2, iDtype

        character(len=20) :: dfpt_tag
        logical :: l_real 

        real,allocatable :: e2_vm(:,:,:) ! q=0 part of Eii(2) term
        complex,allocatable :: E2ndOrdII(:,:) ! Eii(2) potential

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
        allocate(qpts%bk,mold=fi%juphon%qvec) ! this is not nice. Maybe change the expected type in dfpt_sternheimer/make_stars
        qpts%bk(:, :size(fi%juPhon%qvec,2)) = fi%juPhon%qvec
        allocate(q_list(size(fi%juPhon%qvec,2)))  
        q_list = (/(iArray, iArray=1,SIZE(fi%juPhon%qvec,2), 1)/)

        allocate(dyn_mat(size(q_list),3*fi%atoms%ntype,3*fi%atoms%ntype))
        dyn_mat = cmplx(0.0,0.0)

        ! calculate the q=0 part of the Eii(2) term 
        allocate(e2_vm(fi%atoms%nat,3,3))
        allocate( E2ndOrdII(3 * fi%atoms%nat, 3 * fi%atoms%nat) )

        call imagrhodummy%copyPotDen(rho)
        call imagrhodummy%resetpotden()
        
        call timestart("Eii2 q=0")
        do iDir = 1 ,3
            do iDir2 = 1 ,3  
                call dfpt_e2_madelung(fi%atoms,fi%input%jspins,imagrhodummy%mt(:,0,:,:),grgrVext3x3(iDir2,iDir)%mt(:,0,:,1),e2_vm(:,iDir2,iDir))
            end do 
        end do 
        call timestop("Eii2 q=0")

        q_start = fi%juPhon%startq
        q_stop = merge(fi%juPhon%stopq,size(q_list),fi%juPhon%stopq/=0)

        do iQ = q_start, q_stop
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

            
            E2ndOrdII = CMPLX(0.0,0.0)
            ! map the q=0 part of Eii(2) to the corresponding element in the dynMat shape 
            do iDtype = 1, fi%atoms%ntype
                do iDir2 = 1, 3
                    do iDir = 1, 3
                        E2ndOrdII(3*(iDtype-1)+iDir2,3*(iDtype-1)+iDir) = e2_vm(iDtype,iDir2,iDir)
                    end do 
                end do 
            end do 

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
                        call vC1Im%reset_dfpt()
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
                                                den1, vTot1, den1Im, vTot1Im, vC1, vC1Im, &
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

                        call timestart("Dynmat")
                        ! Once the first order quantities are converged, we can construct all
                        ! additional necessary quantities and from that the dynamical matrix.
                        call dfpt_dynmat_row(fi, stars, starsq, sphhar, xcpot, nococonv, hybdat, fmpi, qpts, q_list(iQ), iDtype, iDir, &
                                            eig_id, dfpt_eig_id, dfpt_eig_id2, enpara, results, results1, l_real,&
                                            rho, vTot, grRho3, grVext3, grVC3, &
                                            den1, vTot1, den1Im, vTot1Im, vC1, vC1Im, dyn_mat(iQ,3 *(iDtype-1)+iDir,:), E2ndOrdII)

                        call timestop("Dynmat")
                        dyn_mat(iQ,3 *(iDtype-1)+iDir,:) = dyn_mat(iQ,3 *(iDtype-1)+iDir,:) + conjg(E2ndOrdII(3 *(iDtype-1)+iDir,:))
                        if (fmpi%irank==0) write(*,*) "dynmat row for ", dfpt_tag
                        if (fmpi%irank==0) write(*,*) dyn_mat(iQ,3 *(iDtype-1)+iDir,:)
                        if (fmpi%irank == 0 .and. fi%juphon%l_rm_qhdf) call system("rm "//trim(dfpt_tag)//".hdf")
                        call timestop("Dirloop Phonon")
                    end do ! iDir
                call timestop("Typeloop Phonon")
            end do ! iDtype

            
            if (fmpi%irank==0) then
                ! Diagonalize the dynMat 
                write(*,*) '-------------------------'
                call timestart("Dynmat diagonalization")
                call DiagonalizeDynMat(fi%atoms, qpts%bk(:,q_list(iQ)), fi%juPhon%calcEigenVec, dyn_mat(iQ,:,:), eigenVals, eigenVecs, q_list(iQ),.TRUE.,"raw",fi%juphon%l_sumrule)
                call timestop("Dynmat diagonalization")

                call timestart("Frequency calculation")
                call CalculateFrequencies(fi%atoms, q_list(iQ), eigenVals, eigenFreqs,"raw",qpts%bk(:,q_list(iQ)))
                call timestop("Frequency calculation")
            end if 

            IF (fmpi%irank==0) DEALLOCATE(eigenVals, eigenVecs, eigenFreqs)

            call timestop("q-Point")
        end do ! iQ



    end subroutine dfpt_phonon

end module m_dfpt_phonon