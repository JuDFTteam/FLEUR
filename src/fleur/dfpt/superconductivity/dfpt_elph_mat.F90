!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_dfpt_elph_mat

#ifdef CPP_MPI
    USE mpi
#endif
    USE m_juDFT
    
#ifdef _OPENACC_later
    USE cublas
#define CPP_zgemv cublaszgemv
#else
#define CPP_zgemv zgemv
#endif

    USE m_types 
    USE m_constants
    USE m_npy 

    IMPLICIT NONE

CONTAINS
    SUBROUTINE construct_elph_element(sternheimerJob,fi,sphhar,results,fmpi,enpara,nococonv,starsq,v1real,v1imag,vTot,inden,bqpt,eig_id,q_eig_id,iDir,iDtype,killcont,l_real,gmatBuffer,nuWindow)
        ! This routine is very similar to dfpt_eigen
        ! However, we do not need the gmat which is slightly different to z1
        ! Output needs to be different 
        USE m_dfpt_eigen_hssetup
        USE m_util !this needed?
        USE m_types_mpimat
        USE m_dfpt_tlmplm
        USE m_local_hamiltonian
        
        USE m_eig66_io, ONLY : write_eig, read_eig

        IMPLICIT NONE 

        TYPE(t_sternheimerJob),INTENT(IN) :: sternheimerJob
        TYPE(t_fleurinput), INTENT(IN) :: fi 
        TYPE(t_sphhar), INTENT(IN) :: sphhar
        TYPE(t_results), INTENT(IN) :: results
        TYPE(t_mpi),INTENT(IN) :: fmpi
        TYPE(t_enpara), INTENT(IN) :: enpara
        TYPE(t_nococonv), INTENT(IN) :: nococonv
        TYPE(t_stars),INTENT(IN) :: starsq
        TYPE(t_potden), INTENT(IN) :: v1real,v1imag,vtot,inden
        REAL,  INTENT(IN) :: bqpt(3)
        INTEGER, INTENT(IN) :: eig_id, q_eig_id,iDir, iDtype ,killcont(6) 
        LOGICAL, INTENT(IN) :: l_real
        COMPLEX,INTENT(INOUT) :: gmatBuffer(:,:,:,:) !(nu',nu,kpoints,jsp)
        INTEGER, INTENT(IN) :: nuWindow(2) ! Eigenstate value that we want to consider 


        TYPE(t_tlmplm)  :: td, tdV1
        TYPE(t_potden) :: vx
        TYPE(t_usdus)  :: ud
        TYPE(t_lapw)   :: lapw,lapwq
        CLASS(t_mat), ALLOCATABLE :: zMatk, zMatq, gmat ! this we propably rename to something better
        CLASS(t_mat), ALLOCATABLE :: hmat,smat
        TYPE(t_hub1data) :: hub1data

#ifdef CPP_MPI
        INTEGER :: ierr
#endif 

        INTEGER :: jsp, nk_i, nk ,noccbd,noccbdq,nbasfcn,nbasfcnq , i , neigk,neigq, nu , nkqpt
        INTEGER :: tmp1 ,  dealloc_stat
        character(len=300)        :: errmsg
        REAL :: bkpt(3) , bkqpt(3)
        REAL, ALLOCATABLE :: eigk(:),eigq(:)
        COMPLEX, ALLOCATABLE :: gmatH(:,:),gmatS(:,:),tempVec(:),tempMat1(:)

        INTEGER, ALLOCATABLE      :: ev_list(:), q_ev_list(:)

        CALL vx%copyPotDen(vTot)
        ALLOCATE(vx%pw_w, mold=vx%pw)
        vx%pw_w = vTot%pw_w

        gmatBuffer=0.0 

        ! Get the (lm) matrix elements for V1 and H0
        CALL ud%init(fi%atoms,fi%input%jspins)
        CALL dfpt_tlmplm(fi%atoms,fi%sym,sphhar,fi%input,fi%noco,enpara,fi%hub1inp,hub1data,vTot,fmpi,tdV1,v1real,v1imag,.FALSE.)
        CALL local_ham(sphhar,fi%atoms,fi%sym,fi%noco,nococonv,enpara,fmpi,vTot,vx,inden,fi%input,fi%hub1inp,hub1data,td,ud,0.0,.TRUE.)

#ifndef _OPENACC  
!newer nvhpc versions fail here with ICE        

        DO jsp = 1, merge(1,fi%input%jspins,fi%noco%l_noco)
            DO nk_i = 1, size(fmpi%k_list)
                nk=fmpi%k_list(nk_i)

                ! Get the required eigenvectors and values at k for occupied bands:
                bkpt = fi%kpts%bk(:, nk)
                
                ! find the corresponding kqpt point
                ! folding to first BZ is done in %get_nk
                bkqpt = bkpt + bqpt 
                nkqpt = fi%kpts%get_nk(bkqpt)


                CALL lapw%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, nk, fi%cell, fmpi)
                CALL lapwq%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, nkqpt, fi%cell, fmpi)
 

                noccbd  = COUNT(results%w_iks(:,nk,jsp)*2.0/fi%input%jspins>1.e-8)
                noccbdq = COUNT(results%w_iks(:,nkqpt,jsp)*2.0/fi%input%jspins>1.e-8)
                
                
                nbasfcn  = MERGE(lapw%nv(1)+lapw%nv(2)+2*fi%atoms%nlotot,lapw%nv(1)+fi%atoms%nlotot,fi%noco%l_noco)
                nbasfcnq = MERGE(lapwq%nv(1)+lapwq%nv(2)+2*fi%atoms%nlotot,lapwq%nv(1)+fi%atoms%nlotot,fi%noco%l_noco)
                
                IF (fmpi%n_size == 1) THEN
                    ALLOCATE (t_mat::zMatk)
                    ALLOCATE (t_mat::zMatq)
                ELSE
                    ALLOCATE (t_mpimat::zMatk)
                    ALLOCATE (t_mpimat::zMatq)
                END IF

                ! Initialize the expansion coefficient matrices at k and k+q
                ! Then read all the stuff into it
                CALL zMatk%init(l_real,nbasfcn,nuWindow(2))
                CALL zMatq%init(l_real,nbasfcnq,nuWindow(2))


                ! We only need this for read in 
                ! Data inside we do not care about
                ALLOCATE(ev_list(nuWindow(2)))
                ev_list = (/(i, i=1,nuWindow(2), 1)/)
                ALLOCATE(q_ev_list(nuWindow(2)))
                q_ev_list = (/(i, i=1,nuWindow(2), 1)/)

                ALLOCATE(eigk(nuWindow(2)))
                ALLOCATE(eigq(nuWindow(2)))

                CALL timestart("Read eigenstuff at k/k+q")
                CALL read_eig(eig_id, nk, jsp, list=ev_list, neig=neigk, eig=eigk, zmat=zMatk)
                CALL read_eig(eig_id, nkqpt, jsp, list=q_ev_list, neig=neigq, eig=eigq, zmat=zMatq)
                CALL timestop("Read eigenstuff at k/k+q")


                ! Construct the perturbed Hamiltonian and Overlap matrix perturbations:
                CALL timestart("Setup of matrix perturbations")
                CALL dfpt_eigen_hssetup(sternheimerJob,jsp,fmpi,fi,enpara,nococonv,starsq,ud,td,tdV1,vTot,v1real,lapw,lapwq,iDir,iDtype,hmat,smat,nk,killcont)
                CALL timestop("Setup of matrix perturbations")
    
                IF (fmpi%n_size == 1) THEN
                   ALLOCATE (t_mat::gmat)
                ELSE
                   ALLOCATE (t_mpimat::gmat)
                END IF
    
                ! Initialize the electron-phonon matrix
                ! CALL gmat%init(.FALSE.,nbasfcnq,nuWindow(2))
                ! gmat%data_c(:,:) = cmplx(0.0,0.0)
                ALLOCATE(gmatH(nuWindow(2),nuWindow(2)))
                ALLOCATE(gmatS,mold=gmatH)
                
                gmatH = CMPLX(0.0,0.0)
                gmatS = CMPLX(0.0,0.0)

                ! Allocate auxiliary 
                ! tempVec  size G vectors 
                ! tempMat1 size bandWindow we are interested in  
                ALLOCATE(tempVec(nbasfcnq))
                ALLOCATE(tempMat1(nuWindow(2)))
                tempVec = cmplx(0.0,0.0)
                tempMat1 = cmplx(0.0,0.0)
                CALL timestart("Matrix multiplication")
                DO nu = 1, nuWindow(2) 
                    IF (l_real) THEN ! l_real for zMatk
                        tempVec(:nbasfcnq) = MATMUL(hmat%data_c,zMatk%data_r(:nbasfcn,nu))
                    ELSE
                        CALL CPP_zgemv('N',nbasfcnq,nbasfcn,CMPLX(1.0,0.0),hmat%data_c,nbasfcnq,zMatk%data_c(:nbasfcn,nu),1,CMPLX(0.0,0.0),tempVec,1)
                    END IF
                    
                    IF (zMatq%l_real) THEN ! l_real for zMatq
                        tempMat1(:nuWindow(2)) = MATMUL(TRANSPOSE(zMatq%data_r),tempvec)
                    ELSE
                        CALL CPP_zgemv('C',nbasfcnq,nuWindow(2),CMPLX(1.0,0.0),zmatq%data_c,nbasfcnq,tempvec,1,CMPLX(0.0,0.0),tempMat1,1)
                    END IF
                    ! tempMat1 = H^{(1}_{\nu'\nu}
                    ! We have to think about what happens if nu=nu' at q=0 ---> Is there any rule?
                    ! gmatH exists for testing 
                    ! gmat%data_c(:nbasfcnq,nu) = tempMat1(:nbasfcnq)
                    gmatH(:nuWindow(2),nu) = tempMat1(:nuWindow(2))

                    IF (l_real) THEN ! l_real for zMatk
                        tempVec(:nbasfcnq) = MATMUL(smat%data_c,zMatk%data_r(:nbasfcn,nu))
                    ELSE
                        CALL CPP_zgemv('N',nbasfcnq,nbasfcn,CMPLX(1.0,0.0),smat%data_c,nbasfcnq,zMatk%data_c(:nbasfcn,nu),1,CMPLX(0.0,0.0),tempVec,1)
                    END IF

                    IF (zMatq%l_real) THEN ! l_real for zMatq
                        tempMat1(:nuwindow(2)) = MATMUL(TRANSPOSE(zMatq%data_r),tempvec)
                    ELSE
                        CALL CPP_zgemv('C',nbasfcnq,nuWindow(2),CMPLX(1.0,0.0),zmatq%data_c,nbasfcnq,tempvec,1,CMPLX(0.0,0.0),tempMat1,1)
                    END IF
                     
                    ! tempMat1 = S^{(1}_{\nu'\nu}
                    ! We have to think about what happens if nu=nu' at q=0 ---> Is there any rule?
                    ! gmatS exists for testing 
                    !gmat%data_c(:nbasfcnq,nu) = -eigk(nu)*tempMat1(:nbasfcnq)
                    gmatS(:nuWindow(2),nu) = -eigk(nu)*tempMat1(:nuWindow(2))

                END DO !nu
                CALL timestop("Matrix multiplication")

                ! gmat%data_c(:nuWindow(2),:nuWindow(2)) = gmatH(:nuWindow(2),:nuWindow(2)) + gmatS(:nuWindow(2),:nuWindow(2))

                ! gmatBuffer(:,:,nk,jsp) = gmat%data_c(nuWindow(1):nuWindow(2), nuWindow(1):nuWindow(2))

                gmatBuffer(:,:,nk,jsp) = gmatH(nuWindow(1):nuWindow(2),nuWindow(1):nuWindow(2)) + gmatS(nuWindow(1):nuWindow(2),nuWindow(1):nuWindow(2))
                
                CALL smat%free()
                CALL hmat%free()
                DEALLOCATE(hmat,smat, stat=dealloc_stat, errmsg=errmsg)
                IF(dealloc_stat /= 0) CALL juDFT_error("Deallocation failed for hmat or smat", hint=errmsg, calledby="dfpt_elph_mat.F90")

                
                IF (ALLOCATED(ev_list)) DEALLOCATE(ev_list)
                IF (ALLOCATED(q_ev_list)) DEALLOCATE(q_ev_list)
                IF (ALLOCATED(eigk)) DEALLOCATE(eigk)
                IF (ALLOCATED(eigq)) DEALLOCATE(eigq)
                IF (ALLOCATED(gmatH)) DEALLOCATE(gmatH)
                IF (ALLOCATED(gmatS)) DEALLOCATE(gmatS)
                IF (ALLOCATED(tempVec)) DEALLOCATE(tempVec)
                IF (ALLOCATED(tempMat1)) DEALLOCATE(tempMat1)
                IF (ALLOCATED(zMatk)) THEN
                    CALL zMatk%free()
                    DEALLOCATE(zMatk)
                 END IF
                 IF (ALLOCATED(zmatq)) THEN
                    CALL zMatq%free()
                    DEALLOCATE(zMatq)
                 END IF
                 IF (ALLOCATED(gmat)) THEN
                    CALL gmat%free()
                    DEALLOCATE(gmat)
                 END IF

            END DO !nk_i
        END DO !jsp
#endif 

#ifdef CPP_MPI
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,gmatBuffer,size(gmatBuffer),MPI_DOUBLE_COMPLEX,MPI_SUM,fmpi%mpi_comm,ierr)
#endif


    END SUBROUTINE construct_elph_element

    subroutine el_ph_wannier_interpolate(fmpi,fi,results,dynMats,gmatCart)

        use m_eig66_io, only : open_eig,close_eig,read_eig
        use m_wannier_interpolate
        use m_matrix_interpolation
        use m_dfpt_dynmat_fourier
        use m_dfpt_dynmat_eig
        use m_dfpt_elph_linewidth 

        type(t_mpi), intent(in)         :: fmpi
        type(t_fleurinput) , intent(in) :: fi 
        type(t_results),     intent(in) :: results
        complex,intent(in)              :: dynMats(:,:,:) !(3*nat,3*nat,nqcoarse)
        complex,intent(in)              :: gmatCart(:,:,:,:,:,:) !(nu,nu',kpoints,spin,iPerturb,iQ)

        integer :: eig_id_interpol,q_eig_id_interpol, num_wann , num_bands , ikpt, iMode, ispin, iPerturb
        integer :: iQ, nKept, nSurv, iKept, ne, iWann, i, nLocK, nInterpol
        real    :: qvec(3), pref
        real,allocatable :: kqpts_interpol(:,:)
        integer, allocatable :: myIdx(:)       ! this rank's global fine-k indices
        real,    allocatable :: myKpts(:,:)    ! (3, nLocK) this rank's slice of the fine k-mesh
        integer, allocatable :: fermiKidx(:)   ! local fine-mesh k-point indices near E_F
        real,    allocatable :: fermiKpts(:,:) ! (3, nKept) coords of those k-points
        integer, allocatable :: kqKeptIdx(:)   ! indices (into fermiK*) where k+q also near E_F
        real,    allocatable :: kqKeptKpts(:,:), survKpts(:,:), qsingle(:,:)
        class(t_mat), allocatable :: zMatk, zMatkq ! interpolated eigenvectors (num_wann x num_wann)
        real,    allocatable :: eigBuff(:), eigBuffq(:)
        integer, allocatable :: wann_list(:)
        real,    allocatable :: eigk(:,:,:), eigkq(:,:,:)

        real,    allocatable :: wtkpt_line(:), ph_linewidth(:)
        real,    allocatable :: eigenValsQ(:)
        complex, allocatable :: eigenVecs(:,:), eigenFreqs(:)
        integer              :: lw_unit
        real                  :: sqrtOmegaMax

#ifdef CPP_MPI
        integer :: ierr
#endif
        complex, allocatable :: U_full(:,:,:,:) !(num_bands,num_wann,ikpt,jspin)
        complex, allocatable :: gmatInterpol_q(:,:,:,:,:,:) !(nwann,nwann,nSurv,1,jspin,mode)
        complex, allocatable :: gmatEig(:,:,:,:,:) !(num_wann,num_wann,nSurv,jspin,mode) eigenbasis
        complex, allocatable :: dynMat_interpol(:,:,:)
        type(t_wann_ft), allocatable :: ftRealspace(:,:) !(3*nat, jspins) q-independent forward transforms, built once
        real    :: atomic_mass_array(118)

        num_bands = fi%wannierlib%max_band - fi%wannierlib%min_band + 1
        num_wann  = fi%wannierlib%num_wann

        atomic_mass_array = atomicMasses_const * massInElectronMasses

        if ( (.not. allocated(fi%dfpt%qpts_interpol%bk)) .or. (.not. allocated(fi%dfpt%kpts_interpol%bk))) call juDFT_error("No meshes to interpoalte on.",calledby="dfpt_elph_mat.F90")

        ! distribute the interpolated k-mesh across ranks (round-robin / strided)
        nInterpol = fi%dfpt%kpts_interpol%nkpt
        myIdx = (/(i, i = fmpi%irank+1, nInterpol, fmpi%isize)/)
        nLocK  = size(myIdx)
        allocate(myKpts(3, nLocK))
        do ikpt = 1, nLocK
            myKpts(:, ikpt) = fi%dfpt%kpts_interpol%bk(:, myIdx(ikpt))
        end do

        ! load the Umats from the Wannier step
        call timestart("Load Wannier U matrices")
        allocate(U_full(num_bands,num_wann,fi%kpts%nkptf,fi%input%jspins))
        if (num_bands > num_wann) then
            do ispin = 1 , fi%input%jspins
                do ikpt = 1 , fi%kpts%nkptf
                    U_full(:,:,ikpt,ispin) = matmul(results%U_dis(:,:,ikpt,ispin),results%U_mat(:,:,ikpt,ispin))
                end do
            end do
        else
            U_full(:,:,:,:) = results%U_mat(:,:,:,:)
        end if
        call timestop("Load Wannier U matrices")

        ! diagnostic: verify the matrix-element interpolation
        if (fmpi%irank==0) call check_matrixq_roundtrip(fi, gmatCart, U_full)

        sqrtOmegaMax = global_phonon_energy_bound(fi, dynMats)

        if (nLocK > 0) then
            ! Size of the interpolated hamiltonian is only num_wann x num_wann.
            ! MPI_COMM_SELF: each rank keeps a private (in-memory) eig store for its own k-slice.
#ifdef CPP_MPI
            eig_id_interpol = open_eig(MPI_COMM_SELF, num_wann, num_wann, nLocK, fi%input%jspins, fi%noco%l_noco, &
                                    .true., .false., fi%noco%l_soc, .false., .FALSE., 1)
#else
            eig_id_interpol = open_eig(fmpi%mpi_comm, num_wann, num_wann, nLocK, fi%input%jspins, fi%noco%l_noco, &
                                    .true., .false., fi%noco%l_soc, .false., .FALSE., 1)
#endif
            call timestart("Bandstructure interpolation (k-mesh)")
            call interpolate_bandstructure(fi,results,myKpts,eig_id_interpol,.false.)
            call timestop("Bandstructure interpolation (k-mesh)")

            ! only keep k points where the eigenvalues are close to the fermi energy
            call timestart("Fermi k-point selection")
            call select_fermi_kpoints(fi, results, eig_id_interpol, nLocK, myKpts, fermiKidx, fermiKpts, omegaMax=sqrtOmegaMax)
            call timestop("Fermi k-point selection")
        else
            allocate(fermiKidx(0), fermiKpts(3,0))
        end if
        nKept = size(fermiKidx)

        ! per-rank scratch for the k+q interpolation over the surviving k-points
        if (nKept > 0) then
#ifdef CPP_MPI
            q_eig_id_interpol = open_eig(MPI_COMM_SELF, num_wann, num_wann, nKept, fi%input%jspins, fi%noco%l_noco, &
                                    .true., .false., fi%noco%l_soc, .false., .FALSE., 1)
#else
            q_eig_id_interpol = open_eig(fmpi%mpi_comm, num_wann, num_wann, nKept, fi%input%jspins, fi%noco%l_noco, &
                                    .true., .false., fi%noco%l_soc, .false., .FALSE., 1)
#endif
            allocate(kqpts_interpol(3, nKept))
            allocate(eigBuff(num_wann))
            allocate(eigBuffq(num_wann))
            allocate(wann_list(num_wann))
            wann_list = (/(iWann, iWann=1, num_wann)/)
            allocate(t_mat::zMatk)
            allocate(t_mat::zMatkq)
            call zMatk%init(.false., num_wann, num_wann)
            call zMatkq%init(.false., num_wann, num_wann)
        end if

        allocate(qsingle(3, 1))
        allocate(ph_linewidth(3*fi%atoms%nat))

        if (fmpi%irank==0) &
            open(newunit=lw_unit, file="linewidth", status='replace', action='write', form='formatted')

        ! Build the real-space Wannier-gauge element once.
        call timestart("Forward FT elph-element")
        allocate(ftRealspace(3*fi%atoms%nat, fi%input%jspins))
        do iMode = 1 , 3*fi%atoms%nat
            do ispin = 1 , fi%input%jspins
                call wannier_matrixq_forward(fi,gmatCart(:,:,:,ispin,iMode,:),U_full(:,:,:,ispin),fi%kpts,fi%kpts,ftRealspace(iMode,ispin))
            end do !ispin
        end do !iMode
        call timestop("Forward FT elph-element")

        do iQ = 1 , fi%dfpt%qpts_interpol%nkpt
            call timestart("q-point (el-ph interpolation)")
            ! do one q-point at a time
            qvec = fi%dfpt%qpts_interpol%bk(:, iQ)
            qsingle(:, 1) = qvec
            ph_linewidth = 0.0

            ! interpolate the dynMat and find the eigenvalue for this q
            call timestart("Dynmat Interpolation")
            call interpolate_dynmat(fi%atoms,fi%sym,fi%cell,fi%dfpt%qvec,dynMats,fi%dfpt%l_WSinterpol,qsingle,dynMat_interpol)
            call timestop("Dynmat Interpolation")
            call timestart("Dynmat diagonalization")
            call DiagonalizeDynMat(fi%atoms, qvec, fi%dfpt%calcEigenVec, dynMat_interpol(:,:,1), eigenValsQ, eigenVecs, iQ, .true., &
                                   'band', .false., l_writeOutput=(fmpi%irank==0))
            call timestop("Dynmat diagonalization")
            if (fmpi%irank==0) then
                call timestart("Frequency calculation")
                call CalculateFrequencies(fi%atoms,iQ,eigenValsQ,eigenFreqs,"band",qvec)
                call timestop("Frequency calculation")
            end if

            nSurv = 0
            if (nKept > 0) then
                do ikpt = 1, nKept
                    kqpts_interpol(:, ikpt) = fermiKpts(:, ikpt) + qvec
                end do
                ! interpolate the bands at k+q and keep those also near E_F
                call timestart("k+q bandstructure + Fermi selection")
                call interpolate_bandstructure(fi,results,kqpts_interpol,q_eig_id_interpol,.false.)
                call select_fermi_kpoints(fi, results, q_eig_id_interpol, nKept, kqpts_interpol, kqKeptIdx, kqKeptKpts, omegaMax=maxval(sqrt(abs(eigenValsQ))))
                call timestop("k+q bandstructure + Fermi selection")
                nSurv = size(kqKeptIdx)
            end if

            if (nSurv > 0) then
                ! surviving k-points: k and k+q at E_F
                allocate(survKpts(3, nSurv))
                do ikpt = 1, nSurv
                    survKpts(:, ikpt) = fermiKpts(:, kqKeptIdx(ikpt))
                end do

                ! interpolate the matrix element only for the surviving k-points, single q
                call timestart("Matrix element interpolation")
                allocate(gmatInterpol_q(num_wann, num_wann, nSurv, 1, fi%input%jspins, 3*fi%atoms%nat))
                gmatInterpol_q = cmplx(0.0, 0.0)
                do iMode = 1 , 3*fi%atoms%nat
                    do ispin = 1 , fi%input%jspins
                        call wannier_matrixq_backward(ftRealspace(iMode,ispin),survKpts,qsingle,gmatInterpol_q(:,:,:,:,ispin,iMode))
                    end do !ispin
                end do !iMode
                call timestop("Matrix element interpolation")

                ! rotate the Wannier-gauge matrix element into the eigenbasis:
                !   g_eig(k+q,k) = zmat^dagger(k+q) . g_wann(k+q,k) . zmat(k)
                call timestart("Eigenbasis rotation")
                allocate(gmatEig(num_wann, num_wann, nSurv, fi%input%jspins, 3*fi%atoms%nat))
                allocate(eigk(num_wann, nSurv, fi%input%jspins))
                allocate(eigkq(num_wann, nSurv, fi%input%jspins))
                gmatEig = cmplx(0.0, 0.0)
                do ispin = 1, fi%input%jspins
                    do ikpt = 1, nSurv
                        iKept = kqKeptIdx(ikpt)
                        call read_eig(eig_id_interpol,  fermiKidx(iKept), ispin, list=wann_list,neig=ne, eig=eigBuff, zmat=zMatk)
                        call read_eig(q_eig_id_interpol, iKept, ispin, list=wann_list, neig=ne, eig=eigBuffq, zmat=zMatkq)
                        eigk(:,ikpt,ispin)  = eigBuff
                        eigkq(:,ikpt,ispin) = eigBuffq
                        do iMode = 1, 3*fi%atoms%nat
                            pref = 1.0
                            if (eigenValsQ(iMode) .lt. 0.0) pref = -1*ImagUnit
                            do iPerturb = 1, 3*fi%atoms%nat
                                gmatEig(:,:,ikpt,ispin,iMode) =  gmatEig(:,:,ikpt,ispin,iMode)  +  eigenVecs(iPerturb,iMode) * &
                                                                pref / sqrt(2* atomic_mass_array(fi%atoms%nz(ceiling(iPerturb/3.0))) * sqrt(abs(eigenValsQ(iMode))) ) *&
                                                                matmul(conjg(transpose(zMatkq%data_c)),matmul( gmatInterpol_q(:,:,ikpt,1,ispin,iPerturb),zMatk%data_c ))
                            end do ! iPerturb
                        end do !iMode
                    end do !ikpt
                end do !ispin
                deallocate(gmatInterpol_q)
                call timestop("Eigenbasis rotation")

                allocate(wtkpt_line(nSurv))
                wtkpt_line = 1.0 / real(nInterpol)  

                call timestart("Phonon linewidth construction")
                call dfpt_ph_linewidth(fi, wtkpt_line, eigk, eigkq, gmatEig, eigenValsQ, results%ef, ph_linewidth)
                call timestop("Phonon linewidth construction")

                deallocate(wtkpt_line)
                deallocate(gmatEig, eigk, eigkq, survKpts)
            end if
            if (allocated(kqKeptIdx))  deallocate(kqKeptIdx)
            if (allocated(kqKeptKpts)) deallocate(kqKeptKpts)

            ! sum the per-rank partial linewidths into the full linewidth for this q
#ifdef CPP_MPI
            CALL mpi_allreduce(mpi_in_place, ph_linewidth, size(ph_linewidth), mpi_double_precision, mpi_sum, fmpi%mpi_comm, ierr)
#endif
            if (fmpi%irank==0) then
                write(lw_unit,*) "q-Point", qvec
                write(lw_unit,*) ph_linewidth(:)
            end if
            call timestop("q-point (el-ph interpolation)")
        end do !iQ

        if (fmpi%irank==0) close(lw_unit)

        if (nKept > 0) then
            call zMatk%free()
            deallocate(zMatk)
            call zMatkq%free() 
            deallocate(zMatkq)
            deallocate(eigBuff, eigBuffq, wann_list, kqpts_interpol)
            call close_eig(q_eig_id_interpol)
        end if
        if (nLocK > 0) call close_eig(eig_id_interpol)
        deallocate(ph_linewidth, qsingle)
        deallocate(ftRealspace)

#ifdef CPP_MPI
        CALL MPI_BARRIER(fmpi%mpi_comm, ierr)
#endif

    end subroutine el_ph_wannier_interpolate

    function global_phonon_energy_bound(fi, dynMats) result(omegaMax)
        ! compute maxvalue of the eigenValues of dynMat
        use m_dfpt_dynmat_eig, only : DiagonalizeDynMat

        type(t_fleurinput), intent(in) :: fi
        complex,             intent(in) :: dynMats(:,:,:) !(3*nat,3*nat,nqcoarse)
        real :: omegaMax

        integer :: iQc
        real,    allocatable :: wTmp(:)
        complex, allocatable :: aTmp(:,:)

        omegaMax = 0.0
        if (fi%dfpt%i_integration /= 1) return

        do iQc = 1, fi%dfpt%qvec%nkpt
            call DiagonalizeDynMat(fi%atoms, fi%dfpt%qvec%bk(:,iQc), .false., dynMats(:,:,iQc), &
                                   wTmp, aTmp, iQc, .true., "prescan", .false., .false.)
            omegaMax = max(omegaMax, maxval(sqrt(abs(wTmp))))
            deallocate(wTmp, aTmp)
        end do

    end function global_phonon_energy_bound

    subroutine select_fermi_kpoints(fi, results, eig_id, npoints, coords, keptIdx, keptKpts, omegaMax)
        ! Read the interpolated eigenvalues stored in eig_id for the point set given by
        ! coords(:,1:npoints) and keep only those points that have at least one eigenvalue
        ! within window of the Fermi energy.
        use m_eig66_io, only : read_eig

        type(t_fleurinput), intent(in)  :: fi
        type(t_results),    intent(in)  :: results
        integer,            intent(in)  :: eig_id, npoints
        real,               intent(in)  :: coords(:,:)      ! (3, npoints) coords in eig_id order
        integer, allocatable, intent(out) :: keptIdx(:)     ! indices into 1..npoints near E_F
        real,    allocatable, intent(out) :: keptKpts(:,:)  ! (3, nKept) packed coords
        real,    optional,    intent(in)  :: omegaMax        ! phonon energy scale (Hartree), i_integration=1 only

        integer :: num_wann, jspin, ik, ne, nKept
        real    :: window
        real,    allocatable :: eigvals(:)
        logical, allocatable :: mask(:)

        num_wann = fi%wannierlib%num_wann
        window   = 6.0 * fi%input%tkb
        if (fi%dfpt%i_integration == 1 .and. present(omegaMax)) window = 6.0 * max(fi%input%tkb, omegaMax)

        allocate(eigvals(num_wann))
        allocate(mask(npoints))
        mask = .false.

        ! A point is kept if any band at any spin lies within the window of E_F.
        do jspin = 1, fi%input%jspins
            do ik = 1, npoints
                call read_eig(eig_id, ik, jspin, neig=ne, eig=eigvals)
                if (any(abs(eigvals(:ne) - results%ef) < window)) mask(ik) = .true.
            end do
        end do

        nKept = count(mask)
        allocate(keptIdx(nKept))
        allocate(keptKpts(3, nKept))
        nKept = 0
        do ik = 1, npoints
            if (mask(ik)) then
                nKept = nKept + 1
                keptIdx(nKept)    = ik
                keptKpts(:,nKept) = coords(:, ik)
            end if
        end do

        deallocate(eigvals, mask)

    end subroutine select_fermi_kpoints

    subroutine check_matrixq_roundtrip(fi, gmat, U_full)
        ! Round-trip test of the e-ph matrix-element interpolation:
        ! interpolate gmat onto the SAME coarse (k,q) grid and compare (in the
        ! Wannier gauge) to the directly rotated coarse element. On the coarse
        ! grid the Fourier interpolation is exact, so this must match to ~1e-10.
        use m_matrix_interpolation
        type(t_fleurinput), intent(in) :: fi
        complex,            intent(in) :: gmat(:,:,:,:,:,:)   ! (nb,nb,k,jsp,mode,q) Bloch, coarse
        complex,            intent(in) :: U_full(:,:,:,:)     ! (nb,nwann,nkptf,jsp)

        integer :: nwann, nkpt, nq, ispin, iMode, ik, iq, ikq
        complex, allocatable :: gInterp(:,:,:,:), gRef(:,:,:,:)
        real :: dmax, dref,dint

        nwann = fi%wannierlib%num_wann
        nkpt  = fi%kpts%nkpt
        nq    = size(gmat,6)                 ! coarse q-mesh (== nkpt in the current setup)
        allocate(gInterp(nwann,nwann,nkpt,nq), gRef(nwann,nwann,nkpt,nq))
        dmax = 0.0
        dref = 0.0
        dint = 0.0 

        do ispin = 1, fi%input%jspins
            do iMode = 1, size(gmat,5)
                ! "after": interpolate the coarse gmat back onto the coarse (k,q) grid
                call wannier_matrixq_interpolate(fi, gmat(:,:,:,ispin,iMode,:), U_full(:,:,:,ispin), &
                         fi%kpts, fi%kpts%bk(:,1:nkpt), gInterp, fi%kpts, fi%kpts%bk(:,1:nq))
                ! "before": coarse Bloch element rotated into the Wannier gauge  (U^dag(k+q) g U(k))
                do iq = 1, nq
                    do ik = 1, nkpt
                        ikq = fi%kpts%get_nk(fi%kpts%bk(:,ik) + fi%kpts%bk(:,iq))
                        gRef(:,:,ik,iq) = matmul(conjg(transpose(U_full(:,:,ikq,ispin))), &
                                          matmul(gmat(:,:,ik,ispin,iMode,iq), U_full(:,:,ik,ispin)))
                    end do
                end do
                dmax = max(dmax, maxval(abs(gInterp - gRef)))
                dref = max(dref, maxval(abs(gRef)))
                dint = max(dint, maxval(abs(gInterp)))
            end do
        end do

        write(44,'(a)')        '--- e-ph matrix-element interpolation round-trip (coarse grid) ---'
        write(44,'(a,es12.4)') '   max |g_after - g_before| = ', dmax
        write(44,'(a,es12.4)') '   max |g_before| (scale)   = ', dref
        write(44,'(a,es12.4)') '   max |g_after| (scale)   = ', dint
        write(44,'(a,es12.4)') '   relative deviation        = ', dmax/max(dref,1e-30)
        write(44,'(a)')        '   (should be ~1e-10; large => interpolation loses/rescales gmat)'
        deallocate(gInterp, gRef)
    end subroutine check_matrixq_roundtrip

END MODULE  m_dfpt_elph_mat