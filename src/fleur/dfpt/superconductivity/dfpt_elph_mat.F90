!--------------------------------------------------------------------------------
! Copyright (c) 2024 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
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
        COMPLEX,ALLOCATABLE,INTENT(INOUT) :: gmatBuffer(:,:,:,:) !(nu',nu,kpoints,jsp)
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

    subroutine el_ph_wannier_interpolate(fmpi,fi,results,gmat)

        use m_eig66_io, only : open_eig,close_eig,read_eig
        use m_wannier_interpolate
        use m_matrix_interpolation
        use m_dosbin


        type(t_mpi), intent(in)         :: fmpi
        type(t_fleurinput) , intent(in) :: fi 
        type(t_results),     intent(in) :: results
        complex,intent(in)              :: gmat(:,:,:,:,:,:) !(nu,nu',kpoints,spin,iMode,iQ)

        integer :: eig_id_interpol,q_eig_id_interpol, num_wann , num_bands , ikpt, iMode, ispin
        integer :: iQ, nKept, nSurv, iKept, ne, iWann
        real    :: qvec(3)
        real,allocatable :: kqpts_interpol(:,:)
        integer, allocatable :: fermiKidx(:)   ! fine-mesh k-point indices near E_F
        real,    allocatable :: fermiKpts(:,:) ! (3, nKept) coords of those k-points
        integer, allocatable :: kqKeptIdx(:)   ! indices (into fermiK*) where k+q also near E_F
        real,    allocatable :: kqKeptKpts(:,:), survKpts(:,:), qsingle(:,:)
        class(t_mat), allocatable :: zMatk, zMatkq ! interpolated eigenvectors (num_wann x num_wann)
        real,    allocatable :: eigBuff(:), eigBuffq(:)
        integer, allocatable :: wann_list(:)
        real,    allocatable :: eigk(:,:,:), eigkq(:,:,:) 

        type(t_mpi)          :: fmpi_line
        real,    allocatable :: wtkpt_line(:), eGrid(:), linewidth(:,:), ph_linewidth(:)
        real,    allocatable :: eigenValsQ(:)
        integer              :: gridPoint, nFermi
        logical              :: l_firstWrite

#ifdef CPP_MPI
        integer :: ierr
#endif
        complex, allocatable :: U_full(:,:,:,:) !(num_bands,num_wann,ikpt,jspin)
        complex, allocatable :: gmatInterpol_q(:,:,:,:,:,:) !(nwann,nwann,nSurv,1,jspin,mode)
        complex, allocatable :: gmatEig(:,:,:,:,:) !(num_wann,num_wann,nSurv,jspin,mode) eigenbasis

        num_bands = fi%wannierlib%max_band - fi%wannierlib%min_band + 1
        num_wann  = fi%wannierlib%num_wann
        
        if ( (.not. allocated(fi%dfpt%qvec_interpolate)) .or. (.not. allocated(fi%dfpt%kMesh_interpol))) call juDFT_error("No meshes to interpoalte on.",calledby="dfpt_elph_mat.F90") 

        if (fmpi%irank==0) then

            ! Size of the interpolated hamiltonian is only num_wann x num_wann
            ! MPI_COMM_SELF used, as it would lead to wrong RMA problems 
#ifdef CPP_MPI
            eig_id_interpol = open_eig(MPI_COMM_SELF, num_wann, num_wann, size(fi%dfpt%kMesh_interpol,2), fi%input%jspins, fi%noco%l_noco, &
                                    .true., .false., fi%noco%l_soc, .false., .FALSE., 1)
            q_eig_id_interpol = open_eig(MPI_COMM_SELF, num_wann, num_wann, size(fi%dfpt%kMesh_interpol,2), fi%input%jspins, fi%noco%l_noco, &
                                    .true., .false., fi%noco%l_soc, .false., .FALSE., 1)
#else
            eig_id_interpol = open_eig(fmpi%mpi_comm, num_wann, num_wann, size(fi%dfpt%kMesh_interpol,2), fi%input%jspins, fi%noco%l_noco, &
                                    .true., .false., fi%noco%l_soc, .false., .FALSE., 1)
            q_eig_id_interpol = open_eig(fmpi%mpi_comm, num_wann, num_wann, size(fi%dfpt%kMesh_interpol,2), fi%input%jspins, fi%noco%l_noco, &
                                    .true., .false., fi%noco%l_soc, .false., .FALSE., 1)
#endif

            call interpolate_bandstructure(fi,results,fi%dfpt%kMesh_interpol,eig_id_interpol,.false.)

            ! only keep k points where the eigenvalues are close to the fermi energy
            call select_fermi_kpoints(fi, results, eig_id_interpol, size(fi%dfpt%kMesh_interpol,2), fi%dfpt%kMesh_interpol,fermiKidx, fermiKpts)


            ! load the Umats from the Wannier step
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

            nKept = size(fermiKidx)
            allocate(kqpts_interpol(3, nKept))
            allocate(qsingle(3, 1))
            allocate(eigBuff(num_wann))
            allocate(eigBuffq(num_wann))
            allocate(wann_list(num_wann))
            wann_list = (/(iWann, iWann=1, num_wann)/)

            allocate(t_mat::zMatk)
            allocate(t_mat::zMatkq)
            call zMatk%init(.false., num_wann, num_wann)
            call zMatkq%init(.false., num_wann, num_wann)

            allocate(eGrid(fi%banddos%ndos_points))
            allocate(linewidth(fi%banddos%ndos_points, fi%input%jspins))
            allocate(ph_linewidth(3*fi%atoms%nat))
            allocate(eigenValsQ(3*fi%atoms%nat))
            do gridPoint = 1, fi%banddos%ndos_points
                eGrid(gridPoint) = -8*fi%input%tkb + (16*fi%input%tkb)/(fi%banddos%ndos_points-1.0)*(gridPoint-1.0)
            end do
            nFermi = minloc(abs(eGrid), dim=1)   ! grid point closest to E_F
            ! local serial MPI descriptor: dos_bin_double indexes via fmpi%k_list only
            fmpi_line%irank    = 0
            fmpi_line%isize    = 1
            fmpi_line%mpi_comm = fmpi%mpi_comm
            l_firstWrite = .true.

            do iQ = 1 , size(fi%dfpt%qvec_interpolate,2)
                ! do one q-point at a time
                qvec = fi%dfpt%qvec_interpolate(:, iQ)

                do ikpt = 1, nKept
                    kqpts_interpol(:, ikpt) = fermiKpts(:, ikpt) + qvec
                end do

                ! interpolate the bands at k+q and keep those also near E_F
                call interpolate_bandstructure(fi,results,kqpts_interpol,q_eig_id_interpol,.false.)
                call select_fermi_kpoints(fi, results, q_eig_id_interpol, nKept, kqpts_interpol, kqKeptIdx, kqKeptKpts)

                nSurv = size(kqKeptIdx)
                if (nSurv == 0) then
                    ! no k-point has both k and k+q at E_F for this q
                    deallocate(kqKeptIdx, kqKeptKpts)
                    cycle
                end if

                ! surviving k-points: k and k+q at E_F
                allocate(survKpts(3, nSurv))
                do ikpt = 1, nSurv
                    survKpts(:, ikpt) = fermiKpts(:, kqKeptIdx(ikpt))
                end do

                ! interpolate the matrix element only for the surviving k-points, single q
                qsingle(:, 1) = qvec
                allocate(gmatInterpol_q(num_wann, num_wann, nSurv, 1, fi%input%jspins, 3*fi%atoms%nat))
                gmatInterpol_q = cmplx(0.0, 0.0)
                do iMode = 1 , 3*fi%atoms%nat
                    do ispin = 1 , fi%input%jspins
                        call wannier_matrixq_interpolate(fi,gmat(:,:,:,ispin,iMode,:),U_full(:,:,:,ispin),fi%kpts,survKpts,gmatInterpol_q(:,:,:,:,ispin,iMode),fi%kpts,qsingle)
                    end do !ispin
                end do !iMode

                ! rotate the Wannier-gauge matrix element into the eigenbasis:
                !   g_eig(k+q,k) = zmat^dagger(k+q) . g_wann(k+q,k) . zmat(k)
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
                            gmatEig(:,:,ikpt,ispin,iMode) = matmul(conjg(transpose(zMatkq%data_c)),matmul( gmatInterpol_q(:,:,ikpt,1,ispin,iMode),zMatk%data_c ))
                        end do !iMode
                    end do !ikpt
                end do !ispin
                deallocate(gmatInterpol_q)

                allocate(fmpi_line%k_list(nSurv))
                fmpi_line%k_list = (/(ikpt, ikpt=1,nSurv)/)
                allocate(wtkpt_line(nSurv))
                wtkpt_line = 1.0 / real(size(fi%dfpt%kMesh_interpol,2))

                eigenValsQ = 1.0   ! TODO(PartA): replace with interpolated eigenVals from dfpt_interpolation

                ph_linewidth = 0.0
                do iMode = 1, 3*fi%atoms%nat
                    if (eigenValsQ(iMode) < 0.0) cycle   ! imaginary mode -> linewidth 0
                    linewidth = 0.0
                    call dos_bin_double(fmpi_line, fi%input%jspins, wtkpt_line, eGrid, eigk, eigkq, abs(gmatEig(:,:,:,:,iMode))**2, fi%dfpt%smearingGauss, linewidth, results%ef)
                    do ispin = 1, fi%input%jspins
                        ph_linewidth(iMode) = ph_linewidth(iMode) + tpi_const*sqrt(eigenValsQ(iMode))*linewidth(nFermi,ispin)
                    end do
                end do

                ! write the linewidth for this q
                if (l_firstWrite) then
                    open(110, file="linewidth", status='replace', action='write', form='formatted')
                    l_firstWrite = .false.
                else
                    open(110, file="linewidth", status='old', position='append', action='write')
                end if
                write(110,*) "q-Point", qvec
                write(110,*) ph_linewidth(:)
                close(110)

                deallocate(wtkpt_line, fmpi_line%k_list)
                deallocate(gmatEig, eigk, eigkq, survKpts, kqKeptIdx, kqKeptKpts)
            end do !iQ

            call zMatk%free()
            deallocate(zMatk)
            call zMatkq%free()
            deallocate(zMatkq)
            
            deallocate(eigBuff,eigBuffq, wann_list)
            deallocate(eGrid, linewidth, ph_linewidth, eigenValsQ)
            call close_eig(eig_id_interpol)
            call close_eig(q_eig_id_interpol)

        end if

#ifdef CPP_MPI
        CALL MPI_BARRIER(fmpi%mpi_comm, ierr)
#endif

    end subroutine el_ph_wannier_interpolate

    subroutine select_fermi_kpoints(fi, results, eig_id, npoints, coords, keptIdx, keptKpts)
        ! Read the interpolated eigenvalues stored in eig_id for the point set given by
        ! coords(:,1:npoints) and keep only those points that have at least one eigenvalue
        ! within 6*fermiSmearingEnergy of the Fermi energy.
        use m_eig66_io, only : read_eig

        type(t_fleurinput), intent(in)  :: fi
        type(t_results),    intent(in)  :: results
        integer,            intent(in)  :: eig_id, npoints
        real,               intent(in)  :: coords(:,:)      ! (3, npoints) coords in eig_id order
        integer, allocatable, intent(out) :: keptIdx(:)     ! indices into 1..npoints near E_F
        real,    allocatable, intent(out) :: keptKpts(:,:)  ! (3, nKept) packed coords

        integer :: num_wann, jspin, ik, ne, nKept
        real    :: window
        real,    allocatable :: eigvals(:)
        logical, allocatable :: mask(:)

        num_wann = fi%wannierlib%num_wann
        window   = 6.0 * fi%input%tkb

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

END MODULE  m_dfpt_elph_mat