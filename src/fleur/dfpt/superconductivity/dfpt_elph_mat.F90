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
        type(t_sternheimerJob) :: sternheimerJob
        INTEGER :: iDtype, iDir, killcont(6) ,iMode , iPerturb
        REAL :: bqpt(3)
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
        

        call sternheimerJob%init(fi,l_phonon=.true.)

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
        ! CALL timestart("Generating the energy window")
        ! CALL energy_window(fi,fmpi,results,resultsq,nococonv,bqpt,eigenVals,nuWindow)
        ! CALL timestop("Generating the energy window")
        
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
                CALL dfpt_vgen(sternheimerJob,hybdat,fi%field,fi%input,xcpot,fi%atoms,sphhar,stars,fi%vacuum,fi%sym,&
                           fi%dfpt,fi%cell,fmpi,fi%noco,nococonv,rho_loc,vTot,&
                           starsq,denIn1Im_loc,vTot1,.TRUE.,vTot1Im,denIn1_loc,iDtype,iDir,[1,1])
                    
                CALL timestop("Generating Potential Perturbation")

                ! Add the gradient to potential 
                vTot1%mt(:,0:,iDtype,:) = vTot1%mt(:,0:,iDtype,:) + grVtot(iDir)%mt(:,0:,iDtype,:)

                CALL timestart("Generate electron-phonon matrix element")
                ! currenlty this call will lead to segfalse here! 
                ! reactivate in the future. 
                CALL matrix_element(sternheimerJob,fi,sphhar,results,fmpi,enpara,nococonv,starsq,vTot1,vTot1Im,vTot,rho_loc,bqpt,eig_id,q_eig_id,iDir,iDtype,killcont,l_real,gmatCart,nuWindow) ! nbasfcnq_min
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


    SUBROUTINE matrix_element(sternheimerJob,fi,sphhar,results,fmpi,enpara,nococonv,starsq,v1real,v1imag,vTot,inden,bqpt,eig_id,q_eig_id,iDir,iDtype,killcont,l_real,gmatBuffer,nuWindow)
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


    END SUBROUTINE matrix_element

END MODULE  m_dfpt_elph_mat