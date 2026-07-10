!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_dfpt_postprocess_pot

#ifdef CPP_MPI
    USE mpi
#endif
    USE m_juDFT


    USE m_types 
    USE m_constants
    
    implicit none


contains 

    subroutine dfpt_postprocess_elph(fmpi,fi,stars,sphhar,xcpot,forcetheo,enpara,nococonv,hybdat, &
                                  rho, vTot, vxc,results,eig_id,resultsq,q_eig_id,l_real)

        use m_types 
        use m_cdn_io
        use m_make_stars
        use m_dfpt_dynmat_eig
        use m_eigen 
        use m_dfpt_vgen
        use m_dfpt_elph_mat
        use m_fermie
        use m_dfpt_generate_gradient
        use m_dfpt_vgen

        type(t_mpi), intent(in)       :: fmpi
        type(t_fleurinput),intent(in) :: fi 
        type(t_stars), intent(in)     :: stars
        type(t_sphhar), intent(in)    :: sphhar
        class(t_xcpot), intent(in)       :: xcpot
        type(t_forcetheo),intent(inout) :: forcetheo 
        type(t_enpara), intent(inout) :: enpara
        type(t_nococonv), intent(in)  :: nococonv
        type(t_hybdat), intent(inout) :: hybdat
        type(t_potden), intent(in)    :: rho
        type(t_potden), intent(in)    :: vTot
        type(t_potden), intent(in)    :: vxc
        type(t_results), intent(in)   :: results
        type(t_results),intent(inout) :: resultsq

        integer, intent(in) :: eig_id, q_eig_id
        logical, intent(in) :: l_real

        type(t_hub1data) :: hub1data
        type(t_stars) :: starsq
        type(t_kpts)  :: kqpts , qpts
        type(t_sternheimerJob) :: sternheimerJob
        type(t_potden) :: vTot1, vTot1Im, den1, den1Im, rho_local
        type(t_results) :: dummy_results


        type(t_potden) :: grRho3(3), grVtot3(3), grVext3(3), grVc3(3),grgrVext3x3(3,3)

        integer :: ikpt, iQ ,iDir, iDtype, iPerturb ,iArray, iMode, killcont(6), bandWindowSize
        integer :: bandWindow(2)
        logical :: l_dummy , l_exist
        complex :: pref 
        complex, allocatable :: dynMat(:,:) , eigenVecs(:,:)
        real, allocatable :: eigenVals(:)
        complex, allocatable :: gmatCart(:,:,:,:) ! (nu',nu,kpoints,jsp)
        complex, allocatable :: gmat(:,:,:,:,:,:) ! (nu',nu,kpoints,jsp,iMode,iQ)
        integer, allocatable :: q_list(:)

        character(len=20) :: dfpt_tag
        character(len=100)  :: filename

#ifdef CPP_MPI 
        integer :: ierr
#endif 
        real    :: atomic_mass_array(118)

        atomic_mass_array = atomicMasses_const * massInElectronMasses
        ! killcont can be used to blot out certain contricutions to the
        ! perturbed matrices.
        ! In this order: V1_pw_pw, T1_pw, S1_pw, V1_MT, ikGH0_MT, ikGS0_MT
        killcont = [1,1,1,1,1,1]

        call sternheimerJob%init(fi,l_phonon=.true.)
        call dummy_results%init(fi%input, fi%atoms, fi%kpts, fi%noco)

        call rho_local%copyPotDen(rho)
        
         ! create a kpts type that contains the necessary q-vectors 
        qpts = fi%kpts
        deallocate(qpts%bk)
        allocate(qpts%bk,mold=fi%dfpt%qvec) ! this is not nice. Maybe change the expected type in dfpt_sternheimer/make_stars
        qpts%bk(:, :size(fi%dfpt%qvec,2)) = fi%dfpt%qvec
        allocate(q_list(size(fi%dfpt%qvec,2)))  
        q_list = (/(iArray, iArray=1,SIZE(fi%dfpt%qvec,2), 1)/)

        ! Determine the band window for the electron-phonon matrix elements.
        if (fi%wannierlib%l_wannierize) then
           bandWindow = [fi%wannierlib%min_band, fi%wannierlib%max_band]
        else
           if (.not. allocated(fi%dfpt%bandWindow)) &
              call juDFT_error("dfpt bandWindow is required when not wannierizing",calledby="construct_elph_mat")
           bandWindow = fi%dfpt%bandWindow
        end if

        if (bandWindow(1) < 1 .or. bandWindow(2) < bandWindow(1)) &
           call juDFT_error("Invalid band window in construct_elph_mat",calledby="construct_elph_mat")

        bandWindowSize = bandWindow(2) - bandWindow(1) + 1

        allocate(dynMat(3*fi%atoms%nat,3*fi%atoms%nat))
        allocate(gmat(bandWindowSize,bandWindowSize,fi%kpts%nkpt,fi%input%jspins,3*fi%atoms%nat,size(q_list)))
        allocate(gmatCart(bandWindowSize,bandWindowSize,fi%kpts%nkpt,fi%input%jspins))
        dynMat = cmplx(0,0)
        gmat = cmplx(0,0)
        gmatCart= cmplx(0.0,0.0)

        call timestart("Gradient generation")
        call dfpt_generate_gradient(sternheimerJob,fi,fmpi,sphhar,hybdat,xcpot,nococonv,stars,rho,vTot,grRho3,grVtot3,grVc3,grVext3,grgrVext3x3)
        call timestop("Gradient generation")

        do iQ = 1 , size(q_list)
            call timestart("q-point elph")

            if (fmpi%irank == 0 ) then 
                call timestart("dynMat IO")
                ! Read in eigenvectors and eigenvalues for given q-point
                ! Be careful only irank 0 has eigenVals and eigenVecs allocated 
                call read_dynmats(fi%atoms%nat,iQ,dynMat)
                CALL DiagonalizeDynMat(fi%atoms, qpts%bk(:,q_list(iQ)), fi%dfpt%calcEigenVec, dynMat(:,:), eigenVals, eigenVecs, q_list(iQ),.false.,"raw",.false.,l_writeOutput=.false.)
                call timestop("dynMat IO")
            end if 

#ifdef CPP_MPI
            call MPI_BARRIER(fmpi%mpi_comm, ierr)
#endif
            do iDtype = 1 , fi%atoms%nat
                call timestart("Typeloop")
                do iDir = 1 , 3
                    call timestart("Dirloop")
                    
                    write(dfpt_tag,'(a1,i0,a2,i0,a2,i0)') 'q', q_list(iQ), '_b', iDtype, '_j', iDir

                    iPerturb = iDir+3*(iDtype-1)


                    call make_stars(starsq, fi%sym, fi%atoms, fi%vacuum, sphhar, fi%input, fi%cell, fi%noco, fmpi, qpts%bk(:,iQ), iDtype, iDir,sternheimerJob%l_efield)
                    starsq%ufft = stars%ufft

                    call den1%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
                    call den1Im%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)

                    call vTot1%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
                    call vTot1Im%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)
                    
                    ! allocate the pw_w part
                    allocate( vTot1%pw_w(size(vTot1%pw,1),size(vTot1%pw,2)))

                    if (fmpi%irank==0) then 
                        call timestart("den1 IO")
                        ! We write out the density response in the sternheimer iteration 
                        ! read in the densities
                        filename = trim(dfpt_tag)
                        inquire(file=trim(filename)//".hdf",EXIST=l_exist)
                        if (l_exist) call readDensity(starsq, fi%noco, fi%vacuum, fi%atoms, fi%cell, sphhar, &
                                                    fi%input, fi%sym, CDN_ARCHIVE_TYPE_CDN1_const, CDN_INPUT_DEN_const, 0, &
                                                    dummy_results%ef, dummy_results%last_distance, l_dummy, den1,  &
                                                    inFilename=trim(filename),denIm=den1Im)
                        call timestop("den1 IO")
                        
                        ! add the gradient to the density that we store
                        ! as we need den1 = z^(1) - grad as the input of dfpt_vgen
                        den1%mt(:,0:,iDtype,:) = den1%mt(:,0:,iDtype,:) - grRho3(iDir)%mt(:,0:,iDtype,:)
                    end if                     

#ifdef CPP_MPI
                    call den1%distribute(fmpi%mpi_comm)
                    call den1Im%distribute(fmpi%mpi_comm)
#endif 
                    call timestart("Generating Potential Perturbation")
                    call dfpt_vgen(sternheimerJob,hybdat,fi%field,fi%input,xcpot,fi%atoms,sphhar,stars,fi%vacuum,fi%sym,&
                                    fi%dfpt,fi%cell,fmpi,fi%noco,nococonv,rho_local,vTot,&
                                    starsq,den1Im,vTot1,.TRUE.,vTot1Im,den1,iDtype,iDir,[1,1])
                    call timestop("Generating Potential Perturbation")
                    ! The matrix element needs the gradient correction in the MT 
                    vTot1%mt(:,0:,iDtype,:) = vTot1%mt(:,0:,iDtype,:) + grVtot3(iDir)%mt(:,0:,iDtype,:)

                    ! construct the electron-phonon element in cartesian basis 
                    call timestart("construction elph element")
                    call construct_elph_element(sternheimerJob,fi,sphhar,results,fmpi,enpara,nococonv,starsq,vTot1,vTot1Im,vTot,rho, qpts%bk(:, iQ),&
                                        eig_id,q_eig_id,iDir,iDtype,killcont,l_real,gmatCart,bandWindow)

                    call timestop("construction elph element")

                    ! construct the electron-phonon element in the normal basis  
                    if (fmpi%irank == 0 ) then 
                        do iMode = 1 , 3*fi%atoms%nat    
                            pref = 1.0 
                            if (eigenVals(iMode) .lt. 0.0) pref = -1*ImagUnit
                            gmat(:,:,:,:,iMode,iQ) =  gmat(:,:,:,:,iMode,iQ) + eigenVecs(iPerturb,iMode) * & 
                                                pref / sqrt(2* atomic_mass_array(fi%atoms%nz(ceiling(iPerturb/3.0))) * sqrt(abs(eigenVals(iMode))) ) * gmatCart(:,:,:,:)     
                        end do 
                    end if 
                    
                    ! reset some variables 
                    call starsq%reset_stars()
                    call vTot1%reset_dfpt()
                    call vTot1Im%reset_dfpt()
                    call den1%reset_dfpt()
                    call den1Im%reset_dfpt()
                    
                call timestop("Dirloop")
                end do ! iDir
                call timestop("Typeloop")
            end do ! iDtype

            call timestop("q-point elph")
  
     
        end do !iQ
        ! free some memory
        deallocate(gmatCart)

        if (fmpi%irank==0) print * , "Starting thhe construction of the interpolation"

        ! Perform Wannier interpolation
        if (fi%wannierlib%l_wannierize) call el_ph_wannier_interpolate(fmpi,fi,results,gmat)


    end subroutine dfpt_postprocess_elph


    subroutine read_dynmats(natoms,iQ,dynMat)

        integer, intent(in) :: natoms
        integer, intent(in)  :: iQ
        complex, intent(out) :: dynMat(:,:) 

        integer :: iread 
        character(len=100) :: trash 
        real,allocatable    :: numbers(:,:)

        allocate(numbers(3*natoms,6*natoms))
        numbers = 0.0
        if (iQ<=9) then
            open( 3001, file="dynMatq=000"//int2str(iQ), status="old")
        else if (iQ<=99) then  
            open( 3001, file="dynMatq=00"//int2str(iQ), status="old")
        else if (iQ<=999) then
            open( 3001, file="dynMatq=0"//int2str(iQ), status="old")
        else 
            open( 3001, file="dynMatq="//int2str(iQ), status="old")
        end if 

        
        do iread = 1, 3 + 3*natoms !Loop over dynmat rows
            if (iread<4) then
                read( 3001,*) trash
                write(*,*) iread, trash
            else
                read( 3001,*) numbers(iread-3,:)
                write(*,*) iread, numbers(iread-3,:)
                dynMat(iread-3,:) = cmplx(numbers(iread-3,::2),numbers(iread-3,2::2))
            end if
        end do ! iread
        close(3001)



    end subroutine read_dynmats







end module m_dfpt_postprocess_pot
