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
    SUBROUTINE dfpt_elph_mat(fi,xcpot,sphhar,stars,nococonv,qpts,fmpi,results, resultsq, results1, enpara,hybdat, rho,vTot,grRho,grVtot,iQ,eig_id,q_eig_id,l_real,denIn1,denIn1Im,eigenVecs,eigenVals)

        USE m_vgen
        USE m_make_stars
        USE m_dfpt_vgen
        USE m_eig66_io, ONLY : write_eig, read_eig
        USE m_dosbin
        USE m_smooth
        USE m_dfpt_fermie, ONLY : sfermi
        USE m_dfpt_elph_linewidth



        IMPLICIT NONE 


        TYPE(t_fleurinput), INTENT(IN) :: fi
        CLASS(t_xcpot),     INTENT(IN)    :: xcpot
        TYPE(t_stars),INTENT(IN) :: stars
        TYPE(t_nococonv), INTENT(IN) :: nococonv
        TYPE(t_mpi), INTENT(IN) :: fmpi
        TYPE(t_results), INTENT(IN) :: results,resultsq,results1
        TYPE(t_enpara), INTENT(IN) :: enpara
        TYPE(t_sphhar), INTENT(IN)  :: sphhar
        TYPE(t_kpts), INTENT(IN) :: qpts 
        TYPE(t_hybdat),     INTENT(INOUT) :: hybdat
        TYPE(t_potden), INTENT(IN) :: rho,vTot,grRho(3),grVtot(3)
        INTEGER,INTENT(IN)         :: iQ,eig_id,q_eig_id
        LOGICAL,INTENT(IN)         :: l_real
        TYPE(t_potden), ALLOCATABLE,  INTENT(IN)     :: denIn1(:) , denIn1Im(:)
        COMPLEX, ALLOCATABLE, INTENT(INOUT) :: eigenVecs(:,:) ! Only allocated on irank 0
        REAL,ALLOCATABLE, INTENT(INOUT) :: eigenVals(:) ! Only allocated on irank 0

        
        TYPE(t_potden) :: vTot1,vTot1Im,denIn1_loc, denIn1Im_loc, rho_loc


        TYPE(t_stars) :: starsq
        INTEGER :: iDtype, iDir, killcont(6) ,iMode , iPerturb
        REAL :: bqpt(3)
        COMPLEX :: sigma_loc(2)
        COMPLEX,ALLOCATABLE:: gmatCart(:,:,:,:) !(nu',nu,kpoints,jsp)
        COMPLEX,ALLOCATABLE:: gmat(:,:,:,:,:) !(nu',nu,kpoints,jsp,normal_mode)
        REAL, ALLOCATABLE :: ph_linewidth(:) !(normal_mode)
        INTEGER ::  nuWindow(2,2)  ! ,nbasfcnq_min


#ifdef CPP_MPI
        INTEGER :: ierr
#endif 

        REAL                                      :: atomic_mass_array(118)

        atomic_mass_array = atomicMasses_const * massInElectronMasses
        ! killcont can be used to blot out certain contricutions to the
        ! perturbed matrices.
        ! In this order: V1_pw_pw, T1_pw, S1_pw, V1_MT, ikGH0_MT, ikGS0_MT
        killcont = [1,1,1,1,1,1]
        
        !Up to now only irank == 0 knows the eigenvecs plus eigenvals 
        IF (.NOT. ALLOCATED(eigenVecs)) ALLOCATE(eigenVecs(3*fi%atoms%nat,3*fi%atoms%nat))
        IF (.NOT. ALLOCATED(eigenVals)) ALLOCATE(eigenVals(3*fi%atoms%nat))

#ifdef CPP_MPI
        CALL MPI_BCAST(eigenVecs, size(eigenVecs), MPI_DOUBLE_COMPLEX, 0, fmpi%mpi_comm, ierr)
        CALL MPI_BCAST(eigenVals, size(eigenVals), MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
#endif 

        CALL rho_loc%copyPotDen(rho)
        IF (fmpi%irank==0) WRITE(*,*) 'Generating Potentials for Electron-Phonon Matrix Elements'
        
        ! Introduce Energy Window for states
        ! Reduce memory and computational effort
        bqpt = qpts%bk(:, iQ)
        CALL timestart("Generating the energy window")
        CALL energy_window(fi,fmpi,results,resultsq,nococonv,bqpt,eigenVals,nuWindow)
        CALL timestop("Generating the energy window")
        
        DO iDtype=1,fi%atoms%nat
            DO iDir=1,3
                CALL denIn1_loc%copyPotDen(denIn1(iDir+3*(iDtype-1)))
                CALL denIn1Im_loc%copyPotDen(denIn1Im(iDir+3*(iDtype-1)))

                denIn1_loc%mt(:,0:,iDtype,:) = denIn1_loc%mt(:,0:,iDtype,:) - grRho(iDir)%mt(:,0:,iDtype,:)
                
                CALL make_stars(starsq, fi%sym, fi%atoms, fi%vacuum, sphhar, fi%input, fi%cell, fi%noco, fmpi, qpts%bk(:,iQ), iDtype, iDir)
                starsq%ufft = stars%ufft

                CALL vTot1%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
                CALL vTot1Im%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)


                iPerturb = iDir+3*(iDtype-1)
                
                CALL timestart("Generating Potential Perturbation")
                IF (fmpi%irank==0) WRITE(oUnit, *) "vEff1", iDir
                sigma_loc = cmplx(0.0,0.0)
                CALL dfpt_vgen(hybdat,fi%field,fi%input,xcpot,fi%atoms,sphhar,stars,fi%vacuum,fi%sym,&
                           fi%juphon,fi%cell,fmpi,fi%noco,nococonv,rho_loc,vTot,&
                           starsq,denIn1Im_loc,vTot1,.TRUE.,vTot1Im,denIn1_loc,iDtype,iDir,[1,1],sigma_loc)
                    
                CALL timestop("Generating Potential Perturbation")

                ! Add the gradient to potential 
                vTot1%mt(:,0:,iDtype,:) = vTot1%mt(:,0:,iDtype,:) + grVtot(iDir)%mt(:,0:,iDtype,:)

                CALL timestart("Generate electron-phonon matrix element")
                CALL matrix_element(fi,sphhar,results,resultsq,fmpi,enpara,nococonv,starsq,vTot1,vTot1Im,vTot,rho_loc,bqpt,eig_id,q_eig_id,iDir,iDtype,killcont,l_real,gmatCart,nuWindow) ! nbasfcnq_min
                CALL timestop("Generate electron-phonon matrix element")

                IF (.NOT. ALLOCATED(gmat)) THEN
                    ALLOCATE(gmat(size(gmatCart,1),size(gmatCart,2),size(gmatCart,3),size(gmatCart,4),3*fi%atoms%nat))
                    gmat=CMPLX(0.0,0.0)
                END IF 
                !TODO Read in the eigenvecotrs from Dynmats, here we can take them from memory
                !IF (fmpi%irank==0) THEN 
                    ! Numerics saves the day 
                    ! Think about Gamma if Frequencies are approximately zero
                    DO iMode = 1 , 3*fi%atoms%nat
                        IF (eigenVals(iMode) .LT. 0.0 ) THEN 
                            gmat(:,:,:,:,iMode) = gmat(:,:,:,:,iMode) + eigenVecs(iPerturb,iMode) * &
                            &                   ( (-1*ImagUnit) / SQRT(2* atomic_mass_array(fi%atoms%nz(CEILING(iPerturb/3.0))) * SQRT(ABS(eigenVals(iMode)))) ) * gmatCart(:,:,:,:) 
                        ELSE
                            gmat(:,:,:,:,iMode) = gmat(:,:,:,:,iMode) + eigenVecs(iPerturb,iMode) / SQRT(2* atomic_mass_array(fi%atoms%nz(CEILING(iPerturb/3.0))) * SQRT(eigenVals(iMode))) * gmatCart(:,:,:,:) 
                        END IF 
                    END DO  
                !END IF 
                
                CALL starsq%reset_stars()
                CALL denIn1_loc%reset_dfpt()
                CALL denIn1Im_loc%reset_dfpt()
                CALL vTot1%reset_dfpt()
                CALL vTot1Im%reset_dfpt()
                DEALLOCATE(gmatCart)
            END DO !iDir 
        END DO !iDtype 


        ! Construct the Superconduction temperature 
#ifdef CPP_MPI
        CALL MPI_BARRIER(fmpi%MPI_COMM,ierr)
#endif
        CALL timestart("Linewidth elph")
        !Set this code block behind a logical in the future 
        CALL dfpt_ph_linewidth(fi,fmpi,qpts,results,resultsq,results1,eigenVals,gmat,iQ,nuWindow, ph_linewidth) !nbasfcnq_min 
        CALL timestop("Linewidth elph")

#ifdef CPP_MPI
        CALL MPI_BARRIER(fmpi%MPI_COMM,ierr)
#endif

        DEALLOCATE(ph_linewidth)
        DEALLOCATE(gmat)
        IF (.NOT. fmpi%irank==0) THEN 
            DEALLOCATE(eigenVals)
            DEALLOCATE(eigenVecs)
        END IF 

    END SUBROUTINE dfpt_elph_mat


    SUBROUTINE matrix_element(fi,sphhar,results, resultsq,fmpi,enpara,nococonv,starsq,v1real,v1imag,vTot,inden,bqpt,eig_id,q_eig_id,iDir,iDtype,killcont,l_real,gmatBuffer,nuWindow)
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

        TYPE(t_fleurinput), INTENT(IN) :: fi 
        TYPE(t_sphhar), INTENT(IN) :: sphhar
        TYPE(t_results), INTENT(IN) :: results,resultsq
        TYPE(t_mpi),INTENT(IN) :: fmpi
        TYPE(t_enpara), INTENT(IN) :: enpara
        TYPE(t_nococonv), INTENT(IN) :: nococonv
        TYPE(t_stars),INTENT(IN) :: starsq
        TYPE(t_potden), INTENT(IN) :: v1real,v1imag,vtot,inden
        REAL,  INTENT(IN) :: bqpt(3)
        INTEGER, INTENT(IN) :: eig_id, q_eig_id,iDir, iDtype ,killcont(6) 
        LOGICAL, INTENT(IN) :: l_real
        COMPLEX,ALLOCATABLE,INTENT(INOUT) :: gmatBuffer(:,:,:,:) !(nu',nu,kpoints,jsp)
        !INTEGER,INTENT(INOUT):: nbasfcnq_min
        INTEGER, INTENT(IN) :: nuWindow(2,2)


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

        INTEGER :: jsp, nk_i, nk ,noccbd,noccbdq,nbasfcn,nbasfcnq , i , neigk,neigq, nu , noccbd_max !, nbasfcnq_max
        INTEGER :: tmp1 ,  dealloc_stat
        character(len=300)        :: errmsg
        REAL :: bkpt(3)
        REAL, ALLOCATABLE :: eigk(:),eigq(:)
        COMPLEX, ALLOCATABLE :: gmatH(:,:),gmatS(:,:),tempVec(:),tempMat1(:)

        INTEGER, ALLOCATABLE      :: ev_list(:), q_ev_list(:)


        ! lapw and lapwq can be different if k and k+q do not align
        !nbasfcnq_max = lapwq%dim_nvd()
        !noccbd_max = 0 
        !nbasfcnq_min = 1000000000! ridiculous but should ensure 
        !IF (fmpi%irank==0) THEN 
        !    DO jsp = 1, MERGE(1,fi%input%jspins,fi%noco%l_noco)
        !        DO nk_i = 1, fi%kpts%nkpt
        !            CALL lapw%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, nk_i, fi%cell, fmpi)
        !            CALL lapwq%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, nk_i, fi%cell, fmpi, bqpt)
        !            noccbd  = COUNT(results%w_iks(:,nk_i,jsp)*2.0/fi%input%jspins>1.e-8)
        !            nbasfcnq = MERGE(lapwq%nv(1)+lapwq%nv(2)+2*fi%atoms%nlotot,lapwq%nv(1)+fi%atoms%nlotot,fi%noco%l_noco)
        !            IF (noccbd.GT.noccbd_max) noccbd_max = noccbd
        !            ! maybe we can reduce the calculation to nbasfcnq_min --> saves memory and time 
        !            IF (nbasfcnq .LT. nbasfcnq_min) nbasfcnq_min =nbasfcnq
        !        END DO !nk_i
        !    END DO !jsp 
        !END IF 

        noccbd_max = nuWindow(1,2)

!#ifdef CPP_MPI
!        CALL MPI_BCAST(noccbd_max, 1, MPI_INTEGER, 0, fmpi%mpi_comm, ierr)
!        CALL MPI_BCAST(nbasfcnq_min, 1, MPI_INTEGER, 0, fmpi%mpi_comm, ierr)
!#endif 

        CALL vx%copyPotDen(vTot)
        ALLOCATE(vx%pw_w, mold=vx%pw)
        vx%pw_w = vTot%pw_w
        ALLOCATE(gmatBuffer(nuWindow(2,2)-nuWindow(2,1)+1,nuWindow(1,2)-nuWindow(1,1)+1,size(fmpi%k_list),fi%input%jspins))
        gmatBuffer=0.0 

        ! Get the (lm) matrix elements for V1 and H0
        CALL ud%init(fi%atoms,fi%input%jspins)
        CALL dfpt_tlmplm(fi%atoms,fi%sym,sphhar,fi%input,fi%noco,enpara,fi%hub1inp,hub1data,vTot,fmpi,tdV1,v1real,v1imag,.FALSE.)
        CALL local_ham(sphhar,fi%atoms,fi%sym,fi%noco,nococonv,enpara,fmpi,vTot,vx,inden,fi%input,fi%hub1inp,hub1data,td,ud,0.0,.TRUE.)

#ifndef _OPENACC  
!newer nvhpc versions fail here with ICE        

        DO jsp = 1, MERGE(1,fi%input%jspins,fi%noco%l_noco)
            DO nk_i = 1, size(fmpi%k_list)
                nk=fmpi%k_list(nk_i)

                ! Get the required eigenvectors and values at k for occupied bands:
                bkpt = fi%kpts%bk(:, nk)

                CALL lapw%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, nk, fi%cell, fmpi)
                CALL lapwq%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, nk, fi%cell, fmpi, bqpt)

                noccbd  = COUNT(results%w_iks(:,nk,jsp)*2.0/fi%input%jspins>1.e-8)
                noccbdq = COUNT(resultsq%w_iks(:,nk,jsp)*2.0/fi%input%jspins>1.e-8)
                
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
                CALL zMatk%init(l_real,nbasfcn,noccbd_max)
                CALL zMatq%init(l_real,nbasfcnq,nbasfcnq)


                ! We only need this for read in 
                ! Data inside we do not care about
                ALLOCATE(ev_list(noccbd_max))
                ev_list = (/(i, i=1,noccbd_max, 1)/)
                ALLOCATE(q_ev_list(nbasfcnq))
                q_ev_list = (/(i, i=1,nbasfcnq, 1)/)

                ALLOCATE(eigk(noccbd_max))
                ALLOCATE(eigq(nbasfcnq))

                CALL timestart("Read eigenstuff at k/k+q")
                CALL read_eig(eig_id, nk, jsp, list=ev_list, neig=neigk, eig=eigk, zmat=zMatk)
                CALL read_eig(q_eig_id, nk, jsp, list=q_ev_list, neig=neigq, eig=eigq, zmat=zMatq)
                CALL timestop("Read eigenstuff at k/k+q")


                ! Construct the perturbed Hamiltonian and Overlap matrix perturbations:
                CALL timestart("Setup of matrix perturbations")
                CALL dfpt_eigen_hssetup(jsp,fmpi,fi,enpara,nococonv,starsq,ud,td,tdV1,vTot,v1real,lapw,lapwq,iDir,iDtype,hmat,smat,nk,killcont)
                CALL timestop("Setup of matrix perturbations")
    
                IF (fmpi%n_size == 1) THEN
                   ALLOCATE (t_mat::gmat)
                ELSE
                   ALLOCATE (t_mpimat::gmat)
                END IF
    
                ! Initialize the electron-phonon matrix
                CALL gmat%init(.FALSE.,nbasfcnq,noccbd_max)
                gmat%data_c(:,:) = cmplx(0.0,0.0)
                ALLOCATE(gmatH,mold=gmat%data_c)
                ALLOCATE(gmatS,mold=gmat%data_c)
                
                gmatH = CMPLX(0.0,0.0)
                gmatS = CMPLX(0.0,0.0)

                ! Allocate auxiliary 
                ! tempVec and tempMat1 could be made same variable --> do not need tempVec after usage 
                ALLOCATE(tempVec(nbasfcnq))
                ALLOCATE(tempMat1(nbasfcnq))
                tempVec = cmplx(0.0,0.0)
                tempMat1 = cmplx(0.0,0.0)
                CALL timestart("Matrix multiplication")
                DO nu = 1, noccbd_max ! this might need to be adjusted due to k-integration --> more states needed
                    IF (l_real) THEN ! l_real for zMatk
                        tempVec(:nbasfcnq) = MATMUL(hmat%data_c,zMatk%data_r(:nbasfcn,nu))
                    ELSE
                        CALL CPP_zgemv('N',nbasfcnq,nbasfcn,CMPLX(1.0,0.0),hmat%data_c,nbasfcnq,zMatk%data_c(:nbasfcn,nu),1,CMPLX(0.0,0.0),tempVec,1)
                    END IF
                    
                    IF (zMatq%l_real) THEN ! l_real for zMatq
                        tempMat1(:nbasfcnq) = MATMUL(TRANSPOSE(zMatq%data_r),tempvec)
                    ELSE
                        CALL CPP_zgemv('C',nbasfcnq,nbasfcnq,CMPLX(1.0,0.0),zmatq%data_c,nbasfcnq,tempvec,1,CMPLX(0.0,0.0),tempMat1,1)
                    END IF
                    ! tempMat1 = H^{(1}_{\nu'\nu}
                    ! We have to think about what happens if nu=nu' at q=0 ---> Is there any rule?
                    ! gmatH exists for testing 
                    ! gmat%data_c(:nbasfcnq,nu) = tempMat1(:nbasfcnq)
                    gmatH(:nbasfcnq,nu) = tempMat1(:nbasfcnq)

                    IF (l_real) THEN ! l_real for zMatk
                        tempVec(:nbasfcnq) = MATMUL(smat%data_c,zMatk%data_r(:nbasfcn,nu))
                    ELSE
                        CALL CPP_zgemv('N',nbasfcnq,nbasfcn,CMPLX(1.0,0.0),smat%data_c,nbasfcnq,zMatk%data_c(:nbasfcn,nu),1,CMPLX(0.0,0.0),tempVec,1)
                    END IF

                    IF (zMatq%l_real) THEN ! l_real for zMatq
                        tempMat1(:nbasfcnq) = MATMUL(TRANSPOSE(zMatq%data_r),tempvec)
                    ELSE
                        CALL CPP_zgemv('C',nbasfcnq,nbasfcnq,CMPLX(1.0,0.0),zmatq%data_c,nbasfcnq,tempvec,1,CMPLX(0.0,0.0),tempMat1,1)
                    END IF
                     
                    ! tempMat1 = S^{(1}_{\nu'\nu}
                    ! We have to think about what happens if nu=nu' at q=0 ---> Is there any rule?
                    ! gmatS exists for testing 
                    !gmat%data_c(:nbasfcnq,nu) = -eigk(nu)*tempMat1(:nbasfcnq)
                    gmatS(:nbasfcnq,nu) = -eigk(nu)*tempMat1(:nbasfcnq)

                END DO !nu
                CALL timestop("Matrix multiplication")

                gmat%data_c(:nbasfcnq,:noccbd_max) = gmatH(:nbasfcnq,:noccbd_max) + gmatS(:nbasfcnq,:noccbd_max)

                gmatBuffer(:,:,nk_i,jsp) = gmat%data_c(nuWindow(2,1):nuWindow(2,2), nuWindow(1,1):nuWindow(1,2))

                
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
!#ifdef CPP_MPI
!        CALL MPI_ALLREDUCE(MPI_IN_PLACE,gmatBuffer,size(gmatBuffer),MPI_DOUBLE_COMPLEX,MPI_SUM,fmpi%mpi_comm,ierr)
!        CALL MPI_BARRIER(fmpi%MPI_COMM,ierr)
!#endif

    END SUBROUTINE matrix_element


    SUBROUTINE energy_window(fi,fmpi,results,resultsq,nococonv,bqpt,eigenVals,nuWindow)
        ! This subroutine calculates the state window for the el-ph calculation 
        ! Since we only states close to the fermi level contribute, we will introduce 
        ! For the k and k+q states an effective window that is needed for the calculation

        IMPLICIT NONE 

        TYPE(t_fleurinput), INTENT(IN) :: fi 
        TYPE(t_mpi), INTENT(IN) :: fmpi
        TYPE(t_results), INTENT(IN) :: results,resultsq
        TYPE(t_nococonv), INTENT(IN) :: nococonv
        REAL, INTENT(IN) :: bqpt(3)
        REAL, INTENT(IN) :: eigenVals(:)
        INTEGER,INTENT(OUT) :: nuWindow(2,2)

        REAL :: eig,eigq
        REAL :: efermi, smearing, borderWindow 

        INTEGER :: nbasfcn,nbasfcnq , nu , iNupr, jsp, nk_i
        TYPE(t_lapw)   :: lapw,lapwq


#ifdef CPP_MPI
        INTEGER :: ierr
#endif 

        nuWindow = 0 
        efermi = results%ef 
        !smearing = fi%juPhon%smearingGauss ! lets look into a 5 smearing window around efermi

        
        ! Preset the Windows 
        IF (fmpi%irank == 0 ) THEN 
            CALL lapw%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, 1, fi%cell, fmpi)
            CALL lapwq%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, 1, fi%cell, fmpi, bqpt)

            nbasfcn = MERGE(lapw%nv(1)+lapw%nv(2)+2*fi%atoms%nlotot,lapw%nv(1)+fi%atoms%nlotot,fi%noco%l_noco)
            nbasfcnq = MERGE(lapwq%nv(1)+lapwq%nv(2)+2*fi%atoms%nlotot,lapwq%nv(1)+fi%atoms%nlotot,fi%noco%l_noco)
            
            DO nu = 1, nbasfcn ! If 1 is already above fermi we have no electrons 
                eig = results%eig(nu,1,1)
                IF ( (eig-efermi) .GT. 1e-8 ) THEN 
                    write(2200,*)  "eig that is above" , eig-efermi
                    nuWindow(1,1) = nu - 1 
                    nuWindow(1,2) = nu 
                    EXIT
                END IF
            END DO 

            DO iNupr = 1, nbasfcnq ! If 1 is already above fermi we have no electrons 
                eigq = resultsq%eig(iNupr,1,1)
                IF ( (eigq-efermi) .GT. 1e-8 ) THEN 
                    write(2200,*) "eigq that is above" , eigq-efermi
                    nuWindow(2,1) = iNupr - 1 
                    nuWindow(2,2) = iNupr 
                    EXIT
                END IF 
            END DO     

            write(2200,*) "preguess nu window", nuWindow(1,:)
            write(2200,*) "preguess iNupr window", nuWindow(2,:)
        
        END IF 
        
        ! If the Phonon freq. become larger, then we should set our window to that 
        borderWindow = fi%input%tkb
        IF ( MAXVAL(SQRT(ABS(eigenVals))) .GT. fi%input%tkb ) borderWindow =  MAXVAL(SQRT(ABS(eigenVals)))

        IF (fmpi%irank==0) THEN 
            DO jsp = 1, MERGE(1,fi%input%jspins,fi%noco%l_noco)
                DO nk_i = 1, fi%kpts%nkpt
                    CALL lapw%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, nk_i, fi%cell, fmpi)
                    CALL lapwq%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, nk_i, fi%cell, fmpi, bqpt)
                    nbasfcn = MERGE(lapw%nv(1)+lapw%nv(2)+2*fi%atoms%nlotot,lapw%nv(1)+fi%atoms%nlotot,fi%noco%l_noco)
                    nbasfcnq = MERGE(lapwq%nv(1)+lapwq%nv(2)+2*fi%atoms%nlotot,lapwq%nv(1)+fi%atoms%nlotot,fi%noco%l_noco)

                    DO nu = 1 , nbasfcn 
                        eig = results%eig(nu,nk_i,jsp)  
                        IF ( ABS(eig - efermi)  .LT.   8*borderWindow ) THEN  
                            ! eigenvalue is within energy window 
                            IF ( (eig - efermi) .LT. 1e-8 ) THEN 
                                ! This is the lower end of the window 
                                IF( nu .LT. nuWindow(1,1)) nuWindow(1,1) = nu 
                            ELSE
                                ! This is the upper end of the window 
                                IF( nu .GT. nuWindow(1,2)) nuWindow(1,2) = nu
                            END IF 
                        END IF
                    END DO 

                    DO iNupr= 1, nbasfcnq
                        eigq = resultsq%eig(iNupr,nk_i,jsp)
                        IF ( ABS(eigq - efermi)  .LT.   8*borderWindow ) THEN 
                            ! eigenvalue is within energy window 
                            IF ( (eigq - efermi) .LT. 1e-8 ) THEN 
                                ! This is the lower end of the window 
                                IF( iNupr .LT. nuWindow(2,1)) nuWindow(2,1) = iNupr 
                            ELSE
                                ! This is the upper end of the window 
                                IF( iNupr .GT. nuWindow(2,2)) nuWindow(2,2) = iNupr
                            END IF 
                        END IF
                    END DO 
                END DO !nk_i
            END DO !jsp 

            write(2200,*) "final nu window", nuWindow(1,:)
            write(2200,*) "final iNupr window", nuWindow(2,:)

             ! This is a first thing --> revert when talked with gustav 
            nuWindow(:,1) = MIN(nuWindow(1,1),nuWindow(2,1))
            nuWindow(:,2) = MAX(nuWindow(1,2),nuWindow(2,2))

            write(2200,*) "minval nu window", nuWindow(1,:)
            write(2200,*) "maxval iNupr window", nuWindow(2,:)

        ! This is not the ideal check 
        IF ( (nuWindow(1,2) -nuWindow(1,1)+1) .EQ. size(results%eig,1)  ) CALL juDFT_warn("The phonon calculation appears to have not converged. All states contribute to electron-phonon coupling",calledby="dfpt_elph_mat")
        END IF 



#ifdef CPP_MPI
        CALL MPI_BCAST(nuWindow, size(nuWindow), MPI_INTEGER, 0, fmpi%mpi_comm, ierr)
#endif 


    END SUBROUTINE energy_window


END MODULE  m_dfpt_elph_mat