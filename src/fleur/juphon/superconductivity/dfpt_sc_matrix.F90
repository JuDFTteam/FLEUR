!--------------------------------------------------------------------------------
! Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------


MODULE m_dfpt_sc_matrix

#ifdef CPP_MPI
    USE mpi
#endif
    USE m_juDFT

    IMPLICIT NONE

CONTAINS
    SUBROUTINE dfpt_sc_matrix(fi,sphhar,results,resultsq,fmpi,enpara,nococonv,starsq,v1real,v1imag,vTot,inden,bqpt, &
                                eig_id, q_eig_id,iDir, iDtype,killcont,l_real,g_mat)
        !> Setup of the Electron Phonon Matrix Elements
        !! Adapted from dfpt_eigen()
        !! This calculates the matrix element and already sums over the states \nu' and \nu
        !! Output: Matrix elements for specific k point, q point and cartesian direction
        USE m_types
        USE m_constants
        USE m_dfpt_eigen_hssetup
        USE m_pot_io
        USE m_util
        USE m_eig66_io, ONLY : write_eig, read_eig
        USE m_xmlOutput
        USE m_types_mpimat
        USE m_dfpt_tlmplm
        USE m_local_hamiltonian
        USE m_npy

        USE m_dfpt_tetra


        TYPE(t_fleurinput), INTENT(IN) :: fi
        TYPE(t_sphhar),     INTENT(IN) :: sphhar
        TYPE(t_results),INTENT(INOUT)  :: results, resultsq 
        TYPE(t_mpi),INTENT(IN)         :: fmpi
        TYPE(t_enpara),INTENT(IN)      :: enpara
        TYPE(t_nococonv),INTENT(IN)    :: nococonv
        TYPE(t_stars),INTENT(IN)       :: starsq
        TYPE(t_potden),INTENT(IN)      :: inden, v1real, v1imag, vTot
        REAL,         INTENT(IN)       :: bqpt(3)
        INTEGER,      INTENT(IN)       :: eig_id, q_eig_id, iDir, iDtype, killcont(6)
        LOGICAL,      INTENT(IN)       :: l_real
        REAL,         INTENT(OUT)      :: g_mat  !! store the matrix element at each k point 
        TYPE(t_lapw)                   :: lapw, lapwq
        TYPE(t_tlmplm)                 :: td, tdV1
        TYPE(t_potden)                 :: vx
        TYPE(t_hub1data)               :: hub1data
        TYPE(t_usdus)                  :: ud
        CLASS(t_mat), ALLOCATABLE      :: zMatk, zMatq , tmp_zMatk
        CLASS(t_mat), ALLOCATABLE      :: hmat,smat !smat is not needed here

        INTEGER                        :: i, jsp, nk_i, nk, iNupr, nu
        REAL                           :: bkpt(3), q_loop(3)
        INTEGER                        :: dealloc_stat, nbasfcnq, nbasfcn, neigk, neigq, noccbd, noccbdq
        CHARACTER(len=300)             :: errmsg
        INTEGER, ALLOCATABLE           :: ev_list(:), q_ev_list(:)
        REAL,    ALLOCATABLE           :: eigk(:), eigq(:)

        !!!! this array stores all eigenvalues for all k points
        REAL, ALLOCATABLE              :: eig(:,:,:), eigqq(:,:,:)  ! nu, k , spin
        COMPLEX, ALLOCATABLE           :: el_ph_mat(:,:,:,:) ! nu , nu', k, spin
        INTEGER, ALLOCATABLE           :: noccbd_k(:,:), nbasfcnq_k(:,:) ! number of function, spin
        !!!! delete the upper if neded


        COMPLEX, ALLOCATABLE           :: tempVec(:), tempMat1(:), tempMat2(:)
        !REAL, ALLOCATABLE              :: local_mat(:)  !! store the matrix element at each k point 
        REAL :: temp1,temp2


        COMPLEX, allocatable :: testmat1(:,:,:)
        COMPLEX, allocatable    :: testmat2(:,:)
        COMPLEX, allocatable    :: testmat3(:)
        integer :: nbasfcn_max

        !!!
        !!! What we might have to do is calculate the maximum length of noccbd and nbasfcnq
        !!! as we need the matrix element for each k point and each nu' and nu 
        !!! all of this needs to be at a certain q point 

        !!! What we should store is the sum over each nu and nu' at each k point

        !ALLOCATE(local_mat(size(fmpi%k_list))) !!! THIS STILL HAS TO DEAL WITH SPIIIIIIIN
        !local_mat=0.0
        ! Get the (lm) matrix elements for V1 and H0
        !!! Understand if this is necessary in this portion of code
        CALL vx%copyPotDen(vTot)
        ALLOCATE(vx%pw_w, mold=vx%pw)
        vx%pw_w = vTot%pw_w
        print *,"Set up hamiltonian"
        CALL ud%init(fi%atoms,fi%input%jspins)
        CALL dfpt_tlmplm(fi%atoms,fi%sym,sphhar,fi%input,fi%noco,enpara,fi%hub1inp,hub1data,vTot,fmpi,tdV1,v1real,v1imag,.FALSE.)
        CALL local_ham(sphhar,fi%atoms,fi%sym,fi%noco,nococonv,enpara,fmpi,vTot,vx,inden,fi%input,fi%hub1inp,hub1data,td,ud,0.0,.TRUE.)
        print *,"DONE Set up hamiltonian"
        !!! TEST variables
        nbasfcn_max = max(lapw%dim_nvd(),lapwq%dim_nvd())
        print *, nbasfcn_max
        print *, "Before allocation "
        ALLOCATE(testmat1(size(fmpi%k_list),nbasfcn_max,nbasfcn_max))
        ALLOCATE(testmat2(size(fmpi%k_list),nbasfcn_max))
        ALLOCATE(testmat3(size(fmpi%k_list)))

        !!!
        !!! This is for the tetra tests
        !!!

        print *, "Fermi in sc elements"
        print *, results%ef
        print *, resultsq%ef

        ALLOCATE(noccbd_k(size(fmpi%k_list),MERGE(1,fi%input%jspins,fi%noco%l_noco))) ! occupied for each k point 
        ALLOCATE(nbasfcnq_k,mold=noccbd_k) ! number of states for each k point 
        ALLOCATE(eig(nbasfcn_max,size(fmpi%k_list),MERGE(1,fi%input%jspins,fi%noco%l_noco)))   ! eigenvalues at each k point
        ALLOCATE(eigqq,mold=eig) ! eigenvalues at each k+q point 

        !!!!
        !!!! Number of states for nu has to be adjuted 
        !!!!
        ALLOCATE(el_ph_mat(nbasfcn_max,nbasfcn_max,size(fmpi%k_list),MERGE(1,fi%input%jspins,fi%noco%l_noco))) !electron phonon matrix element for k+q v' , k v 
        eig=0
        eigqq=0
        el_ph_mat=0
        !!!
        !!!
        !!!
        !!!

        DO jsp = 1, MERGE(1,fi%input%jspins,fi%noco%l_noco)
            DO nk_i = 1,size(fmpi%k_list)
                nk=fmpi%k_list(nk_i)
                ! Get the required eigenvectors and values at k for occupied bands:
                bkpt = fi%kpts%bk(:, nk)
                
                q_loop = bqpt

                CALL lapw%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, nk, fi%cell, fmpi)
                CALL lapwq%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, nk, fi%cell, fmpi, q_loop)
                !number of occupied bands
                noccbd  = COUNT(results%w_iks(:,nk,jsp)*2.0/fi%input%jspins>1.e-8) 
                noccbdq = COUNT(resultsq%w_iks(:,nk,jsp)*2.0/fi%input%jspins>1.e-8)
                !number of basis functions
                nbasfcn  = MERGE(lapw%nv(1)+lapw%nv(2)+2*fi%atoms%nlotot,lapw%nv(1)+fi%atoms%nlotot,fi%noco%l_noco)
                nbasfcnq = MERGE(lapwq%nv(1)+lapwq%nv(2)+2*fi%atoms%nlotot,lapwq%nv(1)+fi%atoms%nlotot,fi%noco%l_noco)

                !!! TETRA TEST
                noccbd_k(nk_i,jsp) = noccbd
                nbasfcnq_k(nk_i,jsp) = nbasfcnq
                !!!

                IF (fmpi%n_size == 1) THEN
                    ALLOCATE (t_mat::zMatk)
                    ALLOCATE (t_mat::tmp_zMatk)
                    ALLOCATE (t_mat::zMatq)
                ELSE
                    ALLOCATE (t_mpimat::zMatk)
                    ALLOCATE (t_mpimat::tmp_zMatk)
                    ALLOCATE (t_mpimat::zMatq)
                END IF
        
                ! Initialize the expansion coefficient matrices at k and k+q
                ! Then read all the stuff into it
                CALL zMatk%init(l_real,nbasfcn,noccbd)
                CALL tmp_zMatk%init(l_real,nbasfcn,nbasfcn)
                CALL zMatq%init(l_real,nbasfcnq,nbasfcnq)
                
                !ALLOCATE(ev_list(noccbd))
                !ev_list = (/(i, i=1,noccbd, 1)/)
                ALLOCATE(ev_list(nbasfcn))
                ev_list = (/(i, i=1,nbasfcn, 1)/)
                ALLOCATE(q_ev_list(nbasfcnq))
                q_ev_list = (/(i, i=1,nbasfcnq, 1)/)
    
                !ALLOCATE(eigk(noccbd))
                ALLOCATE(eigk(nbasfcn))
                ALLOCATE(eigq(nbasfcnq))
                
                CALL timestart("Read eigenstuff at k/k+q")
                CALL read_eig(eig_id, nk, jsp, list=ev_list, neig=neigk, eig=eigk, zmat=tmp_zMatk)
                write(1233,*), ev_list ,"ev_list"
                write(1234,*), eigk   ,"eigk"
                write(1232,*), neigk ,"neigk"
                CALL read_eig(q_eig_id, nk, jsp, list=q_ev_list, neig=neigq, eig=eigq, zmat=zMatq)
                CALL timestop("Read eigenstuff at k/k+q")

                call save_npy("now_eigk.npy",eigk)
                call save_npy("now_eigq.npy",eigq)
                IF (l_real) THEN
                    zMatk%data_r(:nbasfcn,:noccbd) = tmp_zMatk%data_r(:nbasfcn,:noccbd)
                ELSE
                    zMatk%data_c(:nbasfcn,:noccbd) = tmp_zMatk%data_c(:nbasfcn,:noccbd)
                END IF
                IF (ALLOCATED(tmp_zMatk)) THEN
                    CALL tmp_zMatk%free()
                    DEALLOCATE(tmp_zMatk)
                END IF

                !!! tetra test
                eig(:nbasfcn,nk_i,jsp) = eigk
                eigqq(:nbasfcnq,nk_i,jsp) = eigq
                !!!!

                ! Construct the perturbed Hamiltonian and Overlap matrix perturbations:
                CALL timestart("Setup of matrix perturbations")
                !! Neends to change to the correct killcont to set smat=0  
                CALL dfpt_eigen_hssetup(jsp,fmpi,fi,enpara,nococonv,starsq,ud,td,tdV1,vTot,v1real,lapw,lapwq,iDir,iDtype,hmat,smat,nk,killcont)
                CALL timestop("Setup of matrix perturbations")

                ! Allocate auxiliary quantities
                ALLOCATE(tempVec(nbasfcnq))
                ALLOCATE(tempMat1(nbasfcnq))
                ALLOCATE(tempMat2(nbasfcnq))
                !TODO: Optimize this with (SCA)LAPACK CALLS
                ! start the H(1) setup
                CALL timestart("Matrix multiplications")
                DO nu = 1, noccbd
                    IF (l_real) THEN ! l_real for zMatk
                        tempVec(:nbasfcnq) = MATMUL(hmat%data_c,zMatk%data_r(:nbasfcn,nu))
                    ELSE
                        tempVec(:nbasfcnq) = MATMUL(hmat%data_c,zMatk%data_c(:nbasfcn,nu))
                    END IF  
                    IF (zMatq%l_real) THEN ! l_real for zMatq
                        tempMat1(:nbasfcnq) = MATMUL(TRANSPOSE(zMatq%data_r),tempvec)
                    ELSE
                        tempMat1(:nbasfcnq) = MATMUL(CONJG(TRANSPOSE(zMatq%data_c)),tempvec)
                    END IF
                    !! end the h(1) setup

                    !! start the s(1) setup 
                    IF (l_real) THEN ! l_real for zMatk
                        tempVec(:nbasfcnq) = MATMUL(smat%data_c,zMatk%data_r(:nbasfcn,nu))
                    ELSE
                        tempVec(:nbasfcnq) = MATMUL(smat%data_c,zMatk%data_c(:nbasfcn,nu))
                    END IF
    
                    IF (zMatq%l_real) THEN ! l_real for zMatq
                        tempMat2(:nbasfcnq) = -eigk(nu) * MATMUL(TRANSPOSE(zMatq%data_r),tempvec)
                    ELSE
                        tempMat2(:nbasfcnq) = -eigk(nu) * MATMUL(CONJG(TRANSPOSE(zMatq%data_c)),tempvec)
                    END IF
                    !! end the s(1) setup

                    !! Now construct the el-ph element $H^1 - e_{k,\nu} S^1$ for each k+q nu', k nu combination
                    DO iNupr = 1, nbasfcnq
                        el_ph_mat(nu,iNupr,nk_i,jsp)  = tempMat1(iNupr) + tempMat2(iNupr) 
                    END DO !iNupr
                END DO ! nu
                CALL timestop("Matrix multiplications")
                
                
                !!!
                !!! We now constructed the phonon matrix element for a specific k-point and q point 
                !!! Summation that are left over k,q, direction, atoms 
                !!! Summation over k within routine 
                !!!


                CALL smat%free()
                CALL hmat%free()
                DEALLOCATE(hmat,smat, stat=dealloc_stat, errmsg=errmsg)
                IF(dealloc_stat /= 0) CALL juDFT_error("Deallocation failed for hmat or smat", hint=errmsg, calledby="dfpt_sc_matrix.F90")


                IF (ALLOCATED(ev_list)) DEALLOCATE(ev_list)
                IF (ALLOCATED(q_ev_list)) DEALLOCATE(q_ev_list)
                IF (ALLOCATED(eigk)) DEALLOCATE(eigk)
                IF (ALLOCATED(eigq)) DEALLOCATE(eigq)
                IF (ALLOCATED(tempVec)) DEALLOCATE(tempVec)
                IF (ALLOCATED(tempMat1)) DEALLOCATE(tempMat1)
                IF (ALLOCATED(tempMat2)) DEALLOCATE(tempMat2)
                IF (ALLOCATED(zMatk)) THEN
                    CALL zMatk%free()
                    DEALLOCATE(zMatk)
                END IF
                IF (ALLOCATED(zMatq)) THEN
                    CALL zMatq%free()
                    DEALLOCATE(zMatq)
                END IF
            
                g_mat = g_mat + temp1
                testmat3(nk_i)=g_mat

            END DO  !k loop 
        END DO ! spin loop ends

        
        !call save_npy(int2str(iDir)//"_dir_k_nu_nu_prime.npy",testmat1)
        !call save_npy(int2str(iDir)//"_dir_k_nu.npy",testmat2)
        !call save_npy(int2str(iDir)//"_dir_k.npy",testmat3)
        DEALLOCATE(testmat1)
        DEALLOCATE(testmat2)
        DEALLOCATE(testmat3)
        
        call save_npy("eig.npy",eig)
        call save_npy("eigqq.npy",eigqq)
        print *, MINVAL(nbasfcnq_k)
        print *, MINVAL(noccbd_k) 
        CALL timestart("K integration")
        CALL dfpt_k_int(fi, fmpi, resultsq,  eig,eigqq,noccbd_k,nbasfcnq_k,el_ph_mat)
        CALL timestop("K integration")
    END SUBROUTINE dfpt_sc_matrix


END MODULE m_dfpt_sc_matrix