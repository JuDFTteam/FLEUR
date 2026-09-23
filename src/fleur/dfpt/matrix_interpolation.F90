!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_matrix_interpolation 
    use m_juDFT
    use m_types
    use m_constants
    use m_npy

    implicit none

    ! Carrier for the q/fine-k-independent part of the double-mesh Wannier
    ! interpolation (real-space tensor + Wigner-Seitz data).
    type t_wann_ft
        integer :: nwann
        integer :: nk1, nk2, nk3, nq1, nq2, nq3
        integer :: ft_lim_k(2,3), ft_lim_q(2,3)
        integer :: bigBox_lim_k(2,3), bigBox_lim_q(2,3)
        real,    allocatable :: FTweight_k(:), FTweight_q(:)
        complex, allocatable :: matWannier(:,:,:,:)   ! nwann,nwann, R_e flat, R_p flat
        ! packed nonzero Wigner-Seitz weights for the backward FT (built once):
        !   *_Rvecs(:,i) big-box coords (phase), *_indFlat(i) 0-based storage idx, *_weightNZ(i) weight
        integer :: ws_k_nNZ = 0, ws_q_nNZ = 0
        integer, allocatable :: ws_k_Rvecs(:,:), ws_k_indFlat(:)
        integer, allocatable :: ws_q_Rvecs(:,:), ws_q_indFlat(:)
        real,    allocatable :: ws_k_weightNZ(:), ws_q_weightNZ(:)
    end type t_wann_ft

contains
    subroutine wannier_matrix_interpolate(fi,matElement,U_mat,kpts_coarse,kpts_fine,matInterpol,qpts_coarse,qpts_fine)

        use m_dfpt_dynmat_fourier , only : ft_dyn_direct, build_ws_ft, ft_fcm_weight_packed

        type(t_fleurinput), intent(in) :: fi
        complex, intent(in) :: matElement(:,:,:)                 ! nu',nu, kpoints
        complex, intent(in) :: U_mat(:,:,:)                      ! bloch, wannier, kpoint
        type(t_kpts),intent(in) :: kpts_coarse                   ! on the coarse Wannier k-mesh
        real,intent(in) :: kpts_fine(:,:)                        ! to interpolate onto
        complex,intent(inout) :: matInterpol(:,:,:)    ! interpolate matrix element  (nwann,nwann,kpts)
        type(t_kpts),optional,intent(in) :: qpts_coarse          ! on the coarse Wannier k-mesh
        type(t_kpts),optional,intent(in) :: qpts_fine            ! to interpolate onto


        integer :: ikpt , nu , iNupr, boxSize, nx, ny, nz, iGrid, ix , iy , iz
        integer :: ft_lim(2,3), bigBox_lim(2,3) ! discrete fourier box limits
        type(t_cell) :: cell
        integer, allocatable  :: supercellR(:,:)
        real, allocatable  :: FTweight(:)

        complex,allocatable :: matWannier(:,:,:)                    ! matrix elements in realspace (wannier gauge)
        complex,allocatable :: matRot(:,:)                          ! rotated matrix elements in wannier gauge
        real :: bkpt(3)

        integer :: nNZ                              ! packed nonzero-weight list (built once)
        integer, allocatable :: Rvecs(:,:), indFlat(:)
        real,    allocatable :: weightNZ(:)

        integer :: nwann, nfine, nGrid

        nwann = size(U_mat, 2)
        nfine = size(kpts_fine,2)

        matInterpol = cmplx(0.0, 0.0)
        allocate(matRot(nwann,nwann))
        matRot = cmplx(0.0,0.0)


        ft_lim(2,:) = kpts_coarse%nkpt3(:)/2
        ft_lim(1,:) = ft_lim(2,:) - kpts_coarse%nkpt3(:) + 1

        nGrid = kpts_coarse%nkpt3(1) * kpts_coarse%nkpt3(2) * kpts_coarse%nkpt3(3)
        allocate(matWannier(nwann,nwann,0:nGrid-1))
        matWannier = cmplx(0.0,0.0)


        call timestart("Forward FT + Wannier gauge")
        do ikpt = 1 , kpts_coarse%nkpt

            bkpt = kpts_coarse%bk(:, ikpt)

            ! rotate the matrix elements into wannier gauge
            ! U^dagger M U
            matRot(:,:) = matmul(conjg(transpose(U_mat(:,:,ikpt))),matmul(matElement(:,:,ikpt),U_mat(:,:,ikpt)))

            call ft_dyn_direct(ft_lim,1,bkpt,matRot,matWannier)
        end do ! ikpt

        matWannier(:,:,:)=matWannier(:,:,:)/kpts_coarse%nkpt
        call timestop("Forward FT + Wannier gauge")

        ! interpolate to fine mesh with WS construction
        call timestart("Wigner-Seitz weights")
        bigBox_lim(2,:) =   2*kpts_coarse%nkpt3(:)
        bigBox_lim(1,:) = - 2*kpts_coarse%nkpt3(:) 
        
        boxSize =  (4*kpts_coarse%nkpt3(1)+1) * (4*kpts_coarse%nkpt3(2)+1) * (4*kpts_coarse%nkpt3(3)+1) 
        allocate(FTweight(boxSize))
        allocate(supercellR(3,boxSize))
        FTweight = 0.0
        supercellR = 0 
        iGrid = 1
        do nz=bigBox_lim(1,3),bigBox_lim(2,3)
            do ny=bigBox_lim(1,2),bigBox_lim(2,2)
                do nx=bigBox_lim(1,1),bigBox_lim(2,1)
                    supercellR(:,iGrid) = (/nx,ny,nz/)
                    iGrid = iGrid+1
                end do 
            end do 
        end do 
        ! compute WS weights for the WS cell
        cell = fi%cell
        call cell%calculate_WSweight(supercellR,FTweight,scaleSupercell=kpts_coarse%nkpt3(:))
        call build_ws_ft(ft_lim, bigBox_lim, FTweight, nNZ, Rvecs, indFlat, weightNZ)
        call timestop("Wigner-Seitz weights")

        call timestart("Backward FT")
        do ikpt = 1 , nfine
            call ft_fcm_weight_packed(-1,nNZ,Rvecs,indFlat,weightNZ,kpts_fine(:,ikpt),matInterpol(:,:,ikpt),matWannier)
        end do !ikpt
        call timestop("Backward FT")



    end subroutine wannier_matrix_interpolate
    
    subroutine wannier_matrixq_interpolate(fi,matElement,U_mat,kpts_coarse,kpts_fine,matInterpol,qpts_coarse,qpts_fine)

        ! Wannier interpolation of a matrix element in bloch basis of type M_{m n}(k,q)
        ! Rotation is U^dagger(k+q) M(k,q) U(k)
        ! The q points have to connect between k points
        ! double fourier transform (R,R) --> (k,q) is done

        type(t_fleurinput), intent(in) :: fi
        complex, intent(in) :: matElement(:,:,:,:)                ! nu',nu, kpoints, qpts (nu' at k+q, nu at k)
        complex, intent(in) :: U_mat(:,:,:)                       ! bloch, wannier, kpoint
        type(t_kpts),intent(in) :: kpts_coarse                    ! on the coarse Wannier k-mesh
        real,intent(in) :: kpts_fine(:,:)                         ! fine k-mesh to interpolate onto
        complex,intent(inout) :: matInterpol(:,:,:,:)   ! interpolated matrix element  (nwann,nwann,kpts,qpts)
        type(t_kpts),intent(in) :: qpts_coarse                    ! on the coarse Wannier q-mesh
        real,intent(in) :: qpts_fine(:,:)                         ! fine q-mesh to interpolate onto

        type(t_wann_ft) :: ft

        call wannier_matrixq_forward(fi, matElement, U_mat, kpts_coarse, qpts_coarse, ft)
        call wannier_matrixq_backward(ft, kpts_fine, qpts_fine, matInterpol)

    end subroutine wannier_matrixq_interpolate

    subroutine wannier_matrixq_forward(fi,matElement,U_mat,kpts_coarse,qpts_coarse,ft)

        ! Build the real-space Wannier-gauge matWannier and the Wigner-Seitz weights.
        
        use m_dfpt_dynmat_fourier , only : ft_dyn_direct, build_ws_ft

        type(t_fleurinput), intent(in)  :: fi
        complex,            intent(in)  :: matElement(:,:,:,:)    ! nu',nu, kpoints, qpts (nu' at k+q, nu at k)
        complex,            intent(in)  :: U_mat(:,:,:)           ! bloch, wannier, kpoint
        type(t_kpts),       intent(in)  :: kpts_coarse            ! coarse Wannier k-mesh
        type(t_kpts),       intent(in)  :: qpts_coarse            ! coarse Wannier q-mesh, summed over its full zone (nkptf/bkf)
        type(t_wann_ft),    intent(out) :: ft

        integer :: ikpt, iqpt, ikqpt, ix, iy, iz, iGrid, iqR
        integer :: boxSize_k, boxSize_q
        integer :: nwann, nk_c, nq_c, nk1, nk2, nk3, nq1, nq2, nq3, nGrid_k, nGrid_q
        type(t_cell) :: cell
        integer, allocatable :: supercellR_k(:,:), supercellR_q(:,:)
        complex, allocatable :: matRot(:,:)                       ! wannier-gauge block for one (k,q)
        complex, allocatable :: matRpK(:,:,:,:)                  ! after q->R_p FT (nwann,nwann, R_p flat, nk_c)
        real :: bkqpt(3)

        nwann = size(U_mat, 2)
        nk_c  = kpts_coarse%nkpt
        nq_c  = qpts_coarse%nkptf

        nk1 = kpts_coarse%nkpt3(1); nk2 = kpts_coarse%nkpt3(2); nk3 = kpts_coarse%nkpt3(3)
        nq1 = qpts_coarse%nkpt3(1); nq2 = qpts_coarse%nkpt3(2); nq3 = qpts_coarse%nkpt3(3)

        ft%nwann = nwann
        ft%nk1 = nk1; ft%nk2 = nk2; ft%nk3 = nk3
        ft%nq1 = nq1; ft%nq2 = nq2; ft%nq3 = nq3

        ! fourier box limits (Monkhorst box -N/2 .. N/2) for each mesh
        ft%ft_lim_k(2,:) = kpts_coarse%nkpt3(:)/2
        ft%ft_lim_k(1,:) = ft%ft_lim_k(2,:) - kpts_coarse%nkpt3(:) + 1
        ft%ft_lim_q(2,:) = qpts_coarse%nkpt3(:)/2
        ft%ft_lim_q(1,:) = ft%ft_lim_q(2,:) - qpts_coarse%nkpt3(:) + 1

        nGrid_k = nk1*nk2*nk3
        nGrid_q = nq1*nq2*nq3

        allocate(matRot(nwann,nwann))
        allocate(matRpK(nwann,nwann, 0:nGrid_q-1, nk_c))
        allocate(ft%matWannier(nwann,nwann, 0:nGrid_k-1, 0:nGrid_q-1))
        ft%matWannier = cmplx(0.0,0.0)

        ! Wigner-Seitz weights and big-box supercell R vectors, one set per mesh
        call timestart("Wigner-Seitz weights (k and q)")
        ft%bigBox_lim_k(2,:) =   2*kpts_coarse%nkpt3(:)
        ft%bigBox_lim_k(1,:) = - 2*kpts_coarse%nkpt3(:)
        boxSize_k = (4*nk1+1) * (4*nk2+1) * (4*nk3+1)
        allocate(ft%FTweight_k(boxSize_k), supercellR_k(3,boxSize_k))
        ft%FTweight_k = 0.0; supercellR_k = 0
        iGrid = 1
        do iz=ft%bigBox_lim_k(1,3),ft%bigBox_lim_k(2,3)
            do iy=ft%bigBox_lim_k(1,2),ft%bigBox_lim_k(2,2)
                do ix=ft%bigBox_lim_k(1,1),ft%bigBox_lim_k(2,1)
                    supercellR_k(:,iGrid) = (/ix,iy,iz/)
                    iGrid = iGrid+1
                end do
            end do
        end do

        ft%bigBox_lim_q(2,:) =   2*qpts_coarse%nkpt3(:)
        ft%bigBox_lim_q(1,:) = - 2*qpts_coarse%nkpt3(:)
        boxSize_q = (4*nq1+1) * (4*nq2+1) * (4*nq3+1)
        allocate(ft%FTweight_q(boxSize_q), supercellR_q(3,boxSize_q))
        ft%FTweight_q = 0.0; supercellR_q = 0
        iGrid = 1
        do iz=ft%bigBox_lim_q(1,3),ft%bigBox_lim_q(2,3)
            do iy=ft%bigBox_lim_q(1,2),ft%bigBox_lim_q(2,2)
                do ix=ft%bigBox_lim_q(1,1),ft%bigBox_lim_q(2,1)
                    supercellR_q(:,iGrid) = (/ix,iy,iz/)
                    iGrid = iGrid+1
                end do
            end do
        end do

        cell = fi%cell
        call cell%calculate_WSweight(supercellR_k,ft%FTweight_k,scaleSupercell=kpts_coarse%nkpt3(:))
        call cell%calculate_WSweight(supercellR_q,ft%FTweight_q,scaleSupercell=qpts_coarse%nkpt3(:))
        call build_ws_ft(ft%ft_lim_k, ft%bigBox_lim_k, ft%FTweight_k, ft%ws_k_nNZ, ft%ws_k_Rvecs, ft%ws_k_indFlat, ft%ws_k_weightNZ)
        call build_ws_ft(ft%ft_lim_q, ft%bigBox_lim_q, ft%FTweight_q, ft%ws_q_nNZ, ft%ws_q_Rvecs, ft%ws_q_indFlat, ft%ws_q_weightNZ)
        call timestop("Wigner-Seitz weights (k and q)")

        ! forward fourier transform ( k-space ---> Realspace )
        ! fourier transform of the q-mesh --> R_p
        call timestart("Forward FT + Wannier gauge")
        matRpK = cmplx(0.0,0.0)
        do ikpt = 1 , nk_c
            do iqpt = 1 , nq_c
                ! fold k+q back into the BZ to find the stored U(k+q)=U(tilde k)
                bkqpt  = kpts_coarse%bk(:,ikpt) + qpts_coarse%bkf(:,iqpt)
                ikqpt = kpts_coarse%get_nk(bkqpt)
                if (ikqpt < 1) call juDFT_error("k+q not on the coarse mesh (q must connect k-points)", &
                                              calledby="wannier_matrixq_forward")
                ! rotate the matrix element into wannier gauge: U^dagger(k+q) M(k,q) U(k)
                matRot(:,:) = matmul(conjg(transpose(U_mat(:,:,ikqpt))),matmul(matElement(:,:,ikpt,iqpt), U_mat(:,:,ikpt)))
                call ft_dyn_direct(ft%ft_lim_q, 1, qpts_coarse%bkf(:,iqpt), matRot, matRpK(:,:,:,ikpt))
            end do !iqpt
            matRpK(:,:,:,ikpt) = matRpK(:,:,:,ikpt) / nq_c
        end do !ikpt
        call timestop("Forward FT + Wannier gauge")

        ! fourier transfrom k -> R_e, for each R_p grid point ----
        call timestart("Forward FT")
        do iqR = 0, nGrid_q-1
            do ikpt = 1 , nk_c
                call ft_dyn_direct(ft%ft_lim_k, 1, kpts_coarse%bk(:,ikpt), matRpK(:,:,iqR,ikpt), ft%matWannier(:,:,:,iqR))
            end do !ikpt
            ft%matWannier(:,:,:,iqR) = ft%matWannier(:,:,:,iqR) / nk_c
        end do
        call timestop("Forward FT")

    end subroutine wannier_matrixq_forward

    subroutine wannier_matrixq_backward(ft,kpts_fine,qpts_fine,matInterpol)

        ! Backward Fourier transform of the precomputed real-space to fine k space
    
        use m_dfpt_dynmat_fourier , only : ft_fcm_weight_packed

        type(t_wann_ft), intent(in)    :: ft
        real,            intent(in)    :: kpts_fine(:,:)          ! fine k-mesh to interpolate onto
        real,            intent(in)    :: qpts_fine(:,:)          ! fine q-mesh to interpolate onto
        complex,         intent(inout) :: matInterpol(:,:,:,:)    ! interpolated matrix element (nwann,nwann,kpts,qpts)

        integer :: ik_fine, iq_fine, ikR, nk_fine, nq_fine, nGrid_k
        complex, allocatable :: matReQ(:,:,:)                  ! after R_p->q' FT (nwann,nwann, R_e flat)

        nk_fine = size(kpts_fine, 2)
        nq_fine = size(qpts_fine, 2)
        nGrid_k = ft%nk1 * ft%nk2 * ft%nk3

        allocate(matReQ(ft%nwann,ft%nwann, 0:nGrid_k-1))
        matInterpol = cmplx(0.0, 0.0)

        ! backwards fourier transform ( Realspace --> kspace )
        call timestart("Backward FT")
        do iq_fine = 1 , nq_fine
            ! ---- stage 1: R_p -> q', for each R_e grid point ----
            matReQ = cmplx(0.0,0.0)
            do ikR = 0, nGrid_k-1
                call ft_fcm_weight_packed(-1, ft%ws_q_nNZ, ft%ws_q_Rvecs, ft%ws_q_indFlat, ft%ws_q_weightNZ, qpts_fine(:,iq_fine), &
                                   matReQ(:,:,ikR), ft%matWannier(:,:,ikR,:))
            end do
            ! ---- stage 2: R_e -> k' ----
            do ik_fine = 1 , nk_fine
                call ft_fcm_weight_packed(-1, ft%ws_k_nNZ, ft%ws_k_Rvecs, ft%ws_k_indFlat, ft%ws_k_weightNZ, kpts_fine(:,ik_fine), &
                                   matInterpol(:,:,ik_fine,iq_fine), matReQ)
            end do !ik_fine
        end do !iq_fine
        call timestop("Backward FT")

    end subroutine wannier_matrixq_backward

end module m_matrix_interpolation 