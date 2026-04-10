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

    subroutine construct_elph_mat(fmpi,fi,stars,sphhar,xcpot,forcetheo,enpara,nococonv,hybdat, &
                                  rho, vTot, vxc,results,eig_id,resultsq,q_eig_id,l_real)

        use m_types 
        use m_cdn_io
        use m_make_stars
        use m_dfpt_dynmat_eig
        use m_eigen 
        use m_dfpt_vgen
        use m_dfpt_elph_mat
        use m_fermie

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
        type(t_potden) :: vTot1, vTot1Im

        integer :: ikpt, iQ ,iDir, iDtype, iPerturb ,iArray, iMode, killcont(6), bandWindowSize
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
        
         ! create a kpts type that contains the necessary q-vectors 
        qpts = fi%kpts
        deallocate(qpts%bk)
        allocate(qpts%bk,mold=fi%juPhon%qvec) ! this is not nice. Maybe change the expected type in dfpt_sternheimer/make_stars
        qpts%bk(:, :size(fi%juPhon%qvec,2)) = fi%juPhon%qvec
        allocate(q_list(size(fi%juPhon%qvec,2)))  
        q_list = (/(iArray, iArray=1,SIZE(fi%juPhon%qvec,2), 1)/)

        bandWindowSize = fi%juPhon%bandWindow(2) - fi%juPhon%bandWindow(1) + 1 

        allocate(dynMat(3*fi%atoms%nat,3*fi%atoms%nat))
        allocate(gmat(bandWindowSize,bandWindowSize,fi%kpts%nkpt,fi%input%jspins,3*fi%atoms%nat,size(q_list)))        
        allocate(gmatCart(bandWindowSize,bandWindowSize,fi%kpts%nkpt,fi%input%jspins))
        dynMat = cmplx(0,0)
        gmat = cmplx(0,0)
        gmatCart= cmplx(0.0,0.0)




        do iQ = 1 , size(q_list)
            call timestart("q-point elph")

            kqpts = fi%kpts
            ! Construct the shifted k-grid
            do ikpt = 1, fi%kpts%nkpt
               kqpts%bk(:, ikpt) = kqpts%bk(:, ikpt) + qpts%bk(:,q_list(iQ))
            end do 

            ! Note:
            ! I currently think we do not need the minus solution for the el-ph matrix element 
            ! since we construct < \psi_{k+q} | V(q) | \psi_{k} }. The density never enters.  

            call timestart("Eigenstuff at k+q")
            ! construct eigenfunction on the k+q grid 
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
            call MPI_BCAST(resultsq%w_iks, SIZE(resultsq%w_iks), MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
#endif

            call timestop("Eigenstuff at k+q")

            if (fmpi%irank == 0 ) then 
                call timestart("dynMat IO")
                ! Read in eigenvectors and eigenvalues for given q-point
                ! Be careful only irank 0 has eigenVals and eigenVecs allocated 
                call read_dynmats(fi%atoms%nat,iQ,dynMat)
                CALL DiagonalizeDynMat(fi%atoms, qpts%bk(:,q_list(iQ)), fi%juPhon%calcEigenVec, dynMat(:,:), eigenVals, eigenVecs, q_list(iQ),.TRUE.,"raw",.FALSE.)
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


                    CALL make_stars(starsq, fi%sym, fi%atoms, fi%vacuum, sphhar, fi%input, fi%cell, fi%noco, fmpi, qpts%bk(:,iQ), iDtype, iDir,sternheimerJob%l_efield)
                    starsq%ufft = stars%ufft

                    CALL vTot1%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
                    CALL vTot1Im%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)

                    CALL timestart("pot1 IO")
                    ! read in the potential perturbation 
                    filename = TRIM("pot1_"//dfpt_tag)
                    INQUIRE(FILE=TRIM(filename)//".hdf",EXIST=l_exist)
                    if (l_exist) CALL readDensity(starsq, fi%noco, fi%vacuum, fi%atoms, fi%cell, sphhar, &
                                                fi%input, fi%sym, CDN_ARCHIVE_TYPE_CDN1_const, CDN_INPUT_DEN_const, 0, &
                                                resultsq%ef, resultsq%last_distance, l_dummy, vTot1,  &
                                                inFilename=TRIM(filename),denIm=vTot1Im)
                    CALL timestop("pot1 IO")                    

                    ! construct electron-phonon element in cartesian basis 
                    call timestart("elph element")
                    CALL matrix_element(sternheimerJob,fi,sphhar,results,resultsq,fmpi,enpara,nococonv,starsq,vTot1,vTot1Im,vTot,rho, qpts%bk(:, iQ),&
                                        eig_id,q_eig_id,iDir,iDtype,killcont,l_real,gmatCart,fi%juPhon%bandWindow) 

                    call timestop("elph element")

                    if (fmpi%irank == 0 ) then 
                        do iMode = 1 , 3*fi%atoms%nat    
                            pref = 1.0 
                            if (eigenVals(iMode) .lt. 0.0) pref = -1*ImagUnit
                            gmat(:,:,:,:,iMode,iQ) =  gmat(:,:,:,:,iMode,iQ) + eigenVecs(iPerturb,iMode) & 
                                                / sqrt(2* atomic_mass_array(fi%atoms%nz(ceiling(iPerturb/3.0))) * sqrt(eigenVals(iMode))) * gmatCart(:,:,:,:)          
                        end do 
                    end if 
                call timestop("Dirloop")
                end do ! iDir
                call timestop("Typeloop")
            end do ! iDtype

            ! reset some variables 
            call starsq%reset_stars()
            call vTot1%reset_dfpt()
            call vTot1Im%reset_dfpt()
            call timestop("q-point elph")
  
            ! Now do IO with el-ph matrix element 
  
        end do !iQ



    end subroutine construct_elph_mat







    subroutine read_dynmats(natoms,iQ,dynMat)

        integer, intent(in) :: natoms
        integer, intent(in)  :: iQ
        complex, allocatable, intent(out) :: dynMat(:,:) 

        integer :: iread 
        character(len=100) :: trash 
        real,allocatable    :: numbers(:,:)

        allocate(numbers(3*natoms,3*natoms))
        numbers = cmplx(0.0,0.0)
        if (iQ<=9) then
            open( 3001, file="dynMatq=000"//int2str(iQ), status="old")
        else if (iQ<=99) then  
            open( 3001, file="dynMatq=00"//int2str(iQ), status="old")
        else 
            open( 3001, file="dynMatq=0"//int2str(iQ), status="old")
        end if 

        
        do iread = 1, 3 + 3*natoms !Loop over dynmat rows
            if (iread<4) then
                read( 3001,*) trash
                write(*,*) iread, trash
            else
                read( 3001,*) numbers(iread-3,:)
                write(*,*) iread, numbers(iread-3,:)
                dynMat(iread-3,:) = CMPLX(numbers(iread-3,::2),numbers(iread-3,2::2))
            end if
        end do ! iread
        close(3001)



    end subroutine read_dynmats







end module m_dfpt_postprocess_pot
