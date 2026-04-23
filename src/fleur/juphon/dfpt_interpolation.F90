!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------



module m_dfpt_interpolation


    use m_juDFT
    use m_constants
    use m_types
    use m_dfpt_NAC

    implicit none 

contains 

    subroutine dfpt_interpolation(fi,fmpi,nococonv,results)

        use m_fleur_init
        use m_dfpt_dynmat_sym
        use m_dfpt_dynmat_eig
        use m_make_dos
        use m_types_eigdos

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
        type(t_wann)       :: wann_fullsym
        type(t_hybdat)     :: hybdat_fullsym
        type(t_mpdata)     :: mpdata_fullsym

        class(t_xcpot),     allocatable :: xcpot_fullsym
        class(t_forcetheo), allocatable :: forcetheo_fullsym

        type(t_kpts)  :: qpts
        integer, allocatable :: q_list(:)
        
        integer :: iQ, iread, iDir, iDir2

        ! for IO of dynMats
        real    :: numbers(3*fi%atoms%nat,6*fi%atoms%nat)
        character(len=100) :: trash, inp_pref 
        character(len=4)   :: dynfiletag 

        ! dynMat properties
        complex, allocatable   :: eigenFreqs(:), eigenVecs(:,:)
        complex, allocatable   :: dyn_mat(:,:,:), dyn_mat_r(:,:,:), dyn_mat_q_full(:,:,:), dyn_mat_pathq(:,:),dyn_mat_NAC_q(:,:),dyn_mat_NAC_r(:,:,:)
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
        if (fi%juPhon%l_band) then 
            dynfiletag = "band"
        else if (fi%juPhon%l_intp .or. fi%juPhon%l_dos) then 
            dynfiletag = "full"
        else
            call juDFT_error("Invalid option for post-procession of dynMatfiles.", calledby="dfpt_interpolation.F90")
        end if 

        IF (fi%juPhon%l_dos) allocate(eigenValsFull(3*fi%atoms%nat,fi%kpts%nkpt,fi%input%jspins))

        ! Read in the dynMats in the IBZ
        ! We currently rely on a seperate set of files fullsym_* that contain the symmetry
        ! Maybe we should change this in the future.
        ! Read qpoints from fullsym_inp.xml and fullsym_kpts.xml
        if (fi%juPhon%qmode .ne. 1 )   call juDFT_error("Invalid option for post-procession of dynMatfiles. qmode not set to 1", calledby="dfpt_interpolation.F90")
        inp_pref = ADJUSTL("fullsym_")
        fmpi_fullsym%l_mpi_multithreaded = fmpi%l_mpi_multithreaded
        fmpi_fullsym%mpi_comm = fmpi%mpi_comm
        call fleur_init(fmpi_fullsym, fi_fullsym, sphhar_fullsym, stars_fullsym, nococonv_fullsym, forcetheo_fullsym, &
                        enpara_fullsym, xcpot_fullsym, results_fullsym, wann_fullsym, hybdat_fullsym, mpdata_fullsym, &
                        inp_pref)
        qpts = fi_fullsym%kpts

        ALLOCATE(q_list(SIZE(qpts%bk,2)))
        q_list = (/(iQ, iQ=1,SIZE(qpts%bk,2), 1)/)

        ALLOCATE(dyn_mat(SIZE(q_list),3*fi%atoms%ntype,3*fi%atoms%ntype))
        ALLOCATE(dyn_mat_NAC_q(3*fi%atoms%ntype,3*fi%atoms%ntype))
        ALLOCATE(dyn_mat_NAC_r(qpts%nkptf,3*fi%atoms%ntype,3*fi%atoms%ntype))

        dyn_mat = cmplx(0.0,0.0)


        if (fmpi%irank==0) then 
            ! this was copied from dfpt.F90 durin refactor --> think about a more clever way
            do iQ = 1, qpts%nkpt ! Loop over dynmat files to read
                if (iQ<=9) THEN
                    open( 3001, file="dynMatq=000"//int2str(iQ), status="old")
                else if(iQ<=99) THEN 
                    open( 3001, file="dynMatq=00"//int2str(iQ), status="old")
                else 
                    open( 3001, file="dynMatq=0"//int2str(iQ), status="old")
                end if 
                do iread = 1, 3 + 3*fi%atoms%nat ! Loop over dynmat rows
                    if (iread<4) then
                        read( 3001,*) trash
                        write(*,*) iread, trash
                    else
                        read( 3001,*) numbers(iread-3,:)
                        write(*,*) iread, numbers(iread-3,:)
                        dyn_mat(iQ,iread-3,:) = cmplx(numbers(iread-3,::2),numbers(iread-3,2::2))
                    end if
                end do  ! iread
                close(3001)
            end do  ! iQ

            !subtract Long range part
            !if (fi%juPhon%l_polar) then
            !    do iQ = 1, qpts%nkpt
            !        dyn_mat_NAC_q = cmplx(0.0,0.0)
            !        call get_NAC_ewald(fi,qpts,stars_fullsym,dyn_mat_NAC_q,qpts%bk(:,iQ),iQ)
            !        dyn_mat(iQ,:,:) = dyn_mat(iQ,:,:) - dyn_mat_NAC_q
            !    end do
            !end if
            !stop


            ! Real space transformation
            ALLOCATE(dyn_mat_r(qpts%nkptf,3*fi%atoms%nat,3*fi%atoms%nat))
            call ft_dyn(fi_fullsym%atoms, qpts, fi_fullsym%sym, fi%cell%amat, dyn_mat, dyn_mat_r, dyn_mat_q_full)
            
            ! In order to call the normal diagonalisation routines
            ! The FCM must be not-normalized --> otherwise we find the wrong unit
            ! Either change here or in dfpt_dynmat_eig.F90 if tag != raw 
            do iDir = 1, 3*fi%atoms%nat
                do iDir2 = 1, 3*fi%atoms%nat
                    dyn_mat_r(:,iDir, iDir2) = dyn_mat_r(:,iDir, iDir2) * massInElectronMasses*  &
                                               SQRT(atomicMasses_const(fi%atoms%nz(CEILING(iDir/3.0)))*atomicMasses_const(fi%atoms%nz(CEILING(iDir2/3.0))))
                end do
            end do
            !subtract Long range part
            if (fi%juPhon%l_polar) then
                dyn_mat_NAC_r = cmplx(0.0,0.0)
                call get_NAC_ewald_r(fi,qpts,stars_fullsym,dyn_mat_NAC_r,[0.,0.,0.],iQ)
                dyn_mat_r(:,:,:) = dyn_mat_r(:,:,:) - dyn_mat_NAC_r(:,:,:)

            end if


            ! interpolate to dense grid on a arbitrary q-point
            ! specified in the inp.xml kpts.xml

            !print*,"fi%kpts%bk(:,iQ)",fi%kpts%bk(:,:)
            !stop
            !print*,"qpts%nkptf",qpts%nkptf
            !print*,"fi%kpts%nkptf",fi%kpts%nkptf
            !stop
            !print*,"fi%kpts%nkpt",fi%kpts%nkpt
            !stop
            do iQ = 1, fi%kpts%nkpt
                call ift_dyn(fi_fullsym%atoms,fi_fullsym%kpts,fi_fullsym%sym,fi%cell%amat,fi%kpts%bk(:,iQ),dyn_mat_r,dyn_mat_pathq)
                !print*,"sum(dyn_mat_pathq) in inter 1",sum(dyn_mat_pathq(:,:))
                WRITE(*,*) '-------------------------'
                !print*,"fi%kpts%bk(:,iQ)",fi%kpts%bk(:,iQ)
                !if ((fi%juPhon%l_polar) .and. (norm2(fi%kpts%bk(:,iQ)) .lt. 1e-8)) then
                !    print*,"at gamma"
                !    call get_NAC(fi,fi%kpts,dyn_mat_pathq,iQ)
                !end if
                if (fi%juPhon%l_polar) then
                    dyn_mat_NAC_q = cmplx(0.0,0.0)
                    print*,"before adding it again"
                    call get_NAC_ewald(fi,qpts,stars_fullsym,dyn_mat_NAC_q,fi%kpts%bk(:,iQ),iQ)
                    print*,"sum(dyn_mat_NAC_q)",sum(dyn_mat_NAC_q)
                    dyn_mat_pathq(:,:) = dyn_mat_pathq(:,:)+ dyn_mat_NAC_q(:,:)
                end if

                print*,"sum(dyn_mat_pathq) in inter 2",sum(dyn_mat_pathq(:,:))
                call timestart("Dynmat diagonalization")
                call DiagonalizeDynMat(fi%atoms, fi%kpts%bk(:,iQ), fi%juPhon%calcEigenVec, dyn_mat_pathq, eigenVals, eigenVecs, iQ,.TRUE.,TRIM(dynfiletag),fi%juPhon%l_sumrule)
                call timestop("Dynmat diagonalization")
                !stop
                call timestart("Frequency calculation")
                call CalculateFrequencies(fi%atoms, iQ, eigenVals, eigenFreqs,TRIM(dynfiletag),fi%kpts%bk(:,iQ))
                call timestop("Frequency calculation")
                !stop
                if (fi%juPhon%l_dos) eigenValsFull(:,iQ,1) = eigenFreqs(:) ! save eigenfrequencies in case of dos

                deallocate(eigenVals, eigenVecs, eigenFreqs, dyn_mat_pathq)
                !stop
            end do ! iQ
            stop

            if (fi%juPhon%l_dos) then 
                banddosLocal = fi%banddos 
                banddosLocal%dos = .true.
                call dos%init(fi%input,fi%atoms,fi%kpts,banddosLocal,.false.,eigenValsFull)
                eigdos(1)%p=>dos 
                call make_dos(fi%kpts,fi%atoms,fi%vacuum,fi%input,fi%banddos,fi%sliceplot,fi%noco,nococonv,fi%sym,fi%cell,results,eigdos,fi%juPhon)
            end if 
        end if !irank==0 


    end subroutine dfpt_interpolation 


end module m_dfpt_interpolation