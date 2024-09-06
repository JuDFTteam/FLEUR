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
    SUBROUTINE dfpt_elph_mat(fi,xcpot,sphhar,stars,nococonv,qpts,fmpi,results, resultsq, enpara,hybdat, rho,vTot,grRho,grVtot,iQ,eig_id,q_eig_id,l_real,denIn1,denIn1Im)

        USE m_vgen
        USE m_make_stars
        USE m_dfpt_vgen


        IMPLICIT NONE 


        TYPE(t_fleurinput), INTENT(IN) :: fi
        CLASS(t_xcpot),     INTENT(IN)    :: xcpot
        TYPE(t_stars),INTENT(IN) :: stars
        TYPE(t_nococonv), INTENT(IN) :: nococonv
        TYPE(t_mpi), INTENT(IN) :: fmpi
        TYPE(t_results), INTENT(IN) :: results,resultsq
        TYPE(t_enpara), INTENT(IN) :: enpara
        TYPE(t_sphhar), INTENT(IN)  :: sphhar
        TYPE(t_kpts), INTENT(IN) :: qpts 
        TYPE(t_hybdat),     INTENT(INOUT) :: hybdat
        TYPE(t_potden), INTENT(IN) :: rho,vTot,grRho(3),grVtot(3)
        INTEGER,INTENT(IN)         :: iQ,eig_id,q_eig_id
        LOGICAL,INTENT(IN)         :: l_real
        TYPE(t_potden), ALLOCATABLE,  INTENT(IN)     :: denIn1(:) , denIn1Im(:)
        
        TYPE(t_potden) :: vTot1,vTot1Im,denIn1_loc, denIn1Im_loc, rho_loc


        TYPE(t_stars) :: starsq
        INTEGER :: iDtype, iDir, killcont(6) 
        REAL :: bqpt(3)
        COMPLEX :: sigma_loc(2)

        ! killcont can be used to blot out certain contricutions to the
        ! perturbed matrices.
        ! In this order: V1_pw_pw, T1_pw, S1_pw, V1_MT, ikGH0_MT, ikGS0_MT
        killcont = [1,1,1,1,1,1]

        CALL rho_loc%copyPotDen(rho)
        IF (fmpi%irank==0) WRITE(*,*) 'Generating Potentials for Electron-Phonon Matrix Elements'
        
        DO iDtype=1,fi%atoms%ntype
            DO iDir=1,3
                CALL denIn1_loc%copyPotDen(denIn1(iDir+3*(iDtype-1)))
                CALL denIn1Im_loc%copyPotDen(denIn1Im(iDir+3*(iDtype-1)))

                denIn1_loc%mt(:,0:,iDtype,:) = denIn1_loc%mt(:,0:,iDtype,:) - grRho(iDir)%mt(:,0:,iDtype,:)
                
                CALL make_stars(starsq, fi%sym, fi%atoms, fi%vacuum, sphhar, fi%input, fi%cell, fi%noco, fmpi, qpts%bk(:,iQ), iDtype, iDir)
                starsq%ufft = stars%ufft

                CALL vTot1%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
                CALL vTot1Im%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)

                bqpt = qpts%bk(:, iQ)
                
                CALL timestart("Generating Potential Perturbation")
                IF (fmpi%irank==0) WRITE(oUnit, *) "vEff1", iDir
                sigma_loc = cmplx(0.0,0.0)
                CALL dfpt_vgen(hybdat,fi%field,fi%input,xcpot,fi%atoms,sphhar,stars,fi%vacuum,fi%sym,&
                           fi%juphon,fi%cell,fmpi,fi%noco,nococonv,rho_loc,vTot,&
                           starsq,denIn1Im_loc,vTot1,.TRUE.,vTot1Im,denIn1_loc,iDtype,iDir,[1,1],sigma_loc)
                    
                CALL timestop("Generating Potential Perturbation")

                CALL timestart("Generate electron-phonon matrix element")
                CALL matrix_element(fi,sphhar,results,resultsq,fmpi,enpara,nococonv,starsq,vTot1,vTot1Im,vTot,rho_loc,bqpt,eig_id,q_eig_id,iDir,iDtype,killcont,l_real)
                CALL timestop("Generate electron-phonon matrix element")

                CALL starsq%reset_stars()
                CALL denIn1_loc%resetpotden()
                CALL denIn1Im_loc%resetpotden()
            END DO !iDir 
        END DO !iDtype 

    END SUBROUTINE dfpt_elph_mat


    SUBROUTINE matrix_element(fi,sphhar,results, resultsq,fmpi,enpara,nococonv,starsq,v1real,v1imag,vTot,inden,bqpt,eig_id,q_eig_id,iDir,iDtype,killcont,l_real)
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


        TYPE(t_tlmplm)  :: td, tdV1
        TYPE(t_potden) :: vx
        TYPE(t_usdus)  :: ud
        TYPE(t_lapw)   :: lapw,lapwq
        CLASS(t_mat), ALLOCATABLE :: zMatk, zMatq, gmat ! this we propably rename to something better
        CLASS(t_mat), ALLOCATABLE :: hmat,smat
        TYPE(t_hub1data) :: hub1data


        INTEGER :: jsp, nk_i, nk ,noccbd,noccbdq,nbasfcn,nbasfcnq , i , neigk,neigq, nu
        REAL :: bkpt(3)
        REAL, ALLOCATABLE :: eigk(:),eigq(:)
        COMPLEX, ALLOCATABLE :: gmatH(:,:),gmatS(:,:),tempVec(:),tempMat1(:)

        INTEGER, ALLOCATABLE      :: ev_list(:), q_ev_list(:)

        CALL vx%copyPotDen(vTot)
        ALLOCATE(vx%pw_w, mold=vx%pw)
        vx%pw_w = vTot%pw_w

        ! Get the (lm) matrix elements for V1 and H0
        CALL ud%init(fi%atoms,fi%input%jspins)
        CALL dfpt_tlmplm(fi%atoms,fi%sym,sphhar,fi%input,fi%noco,enpara,fi%hub1inp,hub1data,vTot,fmpi,tdV1,v1real,v1imag,.FALSE.)
        CALL local_ham(sphhar,fi%atoms,fi%sym,fi%noco,nococonv,enpara,fmpi,vTot,vx,inden,fi%input,fi%hub1inp,hub1data,td,ud,0.0,.TRUE.)


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
                CALL zMatk%init(l_real,nbasfcn,noccbd)
                CALL zMatq%init(l_real,nbasfcnq,nbasfcnq)


                ! We only need this for read in 
                ! Data inside we do not care about
                ALLOCATE(ev_list(noccbd))
                ev_list = (/(i, i=1,noccbd, 1)/)
                ALLOCATE(q_ev_list(nbasfcnq))
                q_ev_list = (/(i, i=1,nbasfcnq, 1)/)

                ALLOCATE(eigk(noccbd))
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
                CALL gmat%init(.FALSE.,nbasfcnq,noccbd)
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
                DO nu = 1, noccbd ! this might need to be adjusted due to k-integration --> more states needed
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

                END DO 

                gmat%data_c(:nbasfcnq,:noccbd) = gmatH(:nbasfcnq,:noccbd) + gmatS(:nbasfcnq,:noccbd)

                CALL timestop("Matrix multiplication")
                IF (ALLOCATED(ev_list)) DEALLOCATE(ev_list)
                IF (ALLOCATED(q_ev_list)) DEALLOCATE(q_ev_list)
                IF (ALLOCATED(eigk)) DEALLOCATE(eigk)
                IF (ALLOCATED(eigq)) DEALLOCATE(eigq)
                IF (ALLOCATED(gmatH)) DEALLOCATE(gmatH)
                IF (ALLOCATED(gmatS)) DEALLOCATE(gmatS)
                IF (ALLOCATED(tempVec)) DEALLOCATE(tempVec)
                IF (ALLOCATED(tempMat1)) DEALLOCATE(tempMat1)
                IF (ALLOCATED(zmatk)) THEN
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
    END SUBROUTINE matrix_element


END MODULE  m_dfpt_elph_mat