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

contains 
    subroutine wannier_matrix_interpolate(fi,matElement,U_mat,kpts_coarse,kpts_fine,matInterpol,qpts_coarse,qpts_fine)

        use m_dfpt_dynmat_fourier , only : ft_dyn_direct, ft_fcm_weight, unfold_grid

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

        complex,allocatable :: fft_grid(:,:,:)                     ! matrix elements in realspace on grid (wannier gauge)
        complex,allocatable :: matWannier(:,:,:,:,:)               ! matrix elements in realspace (wannier gauge)
        complex,allocatable :: matRot(:,:)                          ! rotated matrix elements in wannier gauge
        real :: bkpt(3)

        integer :: nwann, nfine

        nwann = size(U_mat, 2)
        nfine = size(kpts_fine,2)

        allocate(fft_grid(nwann,nwann,size(matElement,3)))
        matInterpol = cmplx(0.0, 0.0)
        allocate(matRot(nwann,nwann))
        matRot = cmplx(0.0,0.0)


        ft_lim(2,:) = kpts_coarse%nkpt3(:)/2
        ft_lim(1,:) = ft_lim(2,:) - kpts_coarse%nkpt3(:) + 1

        allocate(matWannier(nwann,nwann,0:(kpts_coarse%nkpt3(1)-1),&
                            0:(kpts_coarse%nkpt3(2)-1),0:(kpts_coarse%nkpt3(3)-1)))
        matWannier = cmplx(0.0,0.0)


        call timestart("Forward FT + Wannier gauge")
        fft_grid = cmplx(0.0,0.0)
        do ikpt = 1 , kpts_coarse%nkpt

            bkpt = kpts_coarse%bk(:, ikpt)

            ! rotate the matrix elements into wannier gauge
            ! U^dagger M U
            matRot(:,:) = matmul(conjg(transpose(U_mat(:,:,ikpt))),matmul(matElement(:,:,ikpt),U_mat(:,:,ikpt)))

            call ft_dyn_direct(ft_lim,1,bkpt,matRot,fft_grid(:,:,:))
        end do ! ikpt

        ! unfold the grid to nx,ny,nz indexing

        fft_grid(:,:,:)=fft_grid(:,:,:)/kpts_coarse%nkpt
        call unfold_grid(ft_lim, fft_grid, matWannier(:,:,:,:,:))
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
        call timestop("Wigner-Seitz weights")

        call timestart("Backward FT")
        do ikpt = 1 , nfine
            call ft_fcm_weight(-1,ft_lim,bigBox_lim,FTweight,kpts_fine(:,ikpt),matInterpol(:,:,ikpt),matWannier(:,:,:,:,:))
        end do !ikpt
        call timestop("Backward FT")



    end subroutine wannier_matrix_interpolate
    
    subroutine wannier_matrixq_interpolate(fi,matElement,U_mat,kpts_coarse,kpts_fine,matInterpol,qpts_coarse,qpts_fine)

        ! Wannier interpolation of a matrix element in bloch basis of type M_{m n}(k,q)
        ! Rotationi is U^dagger(k+q) M(k,q) U(k)
        ! The q points have to connect between k points
        ! double fourier transform (R,R) --> (k,q) is done 

        use m_dfpt_dynmat_fourier , only : ft_dyn_direct, ft_fcm_weight, unfold_grid

        type(t_fleurinput), intent(in) :: fi
        complex, intent(in) :: matElement(:,:,:,:)                ! nu',nu, kpoints, qpts (nu' at k+q, nu at k)
        complex, intent(in) :: U_mat(:,:,:)                       ! bloch, wannier, kpoint
        type(t_kpts),intent(in) :: kpts_coarse                    ! on the coarse Wannier k-mesh
        real,intent(in) :: kpts_fine(:,:)                         ! fine k-mesh to interpolate onto
        complex,intent(inout) :: matInterpol(:,:,:,:)   ! interpolated matrix element  (nwann,nwann,kpts,qpts)
        type(t_kpts),intent(in) :: qpts_coarse                    ! on the coarse Wannier q-mesh
        real,intent(in) :: qpts_fine(:,:)                         ! fine q-mesh to interpolate onto


        integer :: ikpt, iqpt, ik_fine, iq_fine, ikqpt, ix, iy, iz
        integer :: ft_lim_k(2,3), ft_lim_q(2,3), bigBox_lim_k(2,3), bigBox_lim_q(2,3) ! fourier box limits
        integer :: boxSize_k, boxSize_q, iGrid
        integer :: nwann, nk_c, nq_c, nk_fine, nq_fine
        integer :: nk1, nk2, nk3, nq1, nq2, nq3
        type(t_cell) :: cell
        integer, allocatable :: supercellR_k(:,:), supercellR_q(:,:)
        real, allocatable    :: FTweight_k(:), FTweight_q(:)

        complex,allocatable :: matRot(:,:)                                  ! wannier-gauge block for one (k,q)
        complex,allocatable :: fftq(:,:,:), fftk(:,:,:)                 ! flat accumulators for ft_dyn_direct
        complex,allocatable :: tempMat(:,:,:,:,:,:)                        ! after q->R_p FT   (nwann,nwann,nk_c, R_p)
        complex,allocatable :: matWannier(:,:,:,:,:,:,:,:)                    ! real space        (nwann,nwann, R_e, R_p)
        complex,allocatable :: tempMat2(:,:,:,:,:)                          ! after R_p->q' FT  (nwann,nwann, R_e)
        real :: bkqpt(3)

        nwann   = size(U_mat, 2)
        nk_c    = kpts_coarse%nkpt
        nq_c    = qpts_coarse%nkpt
        nk_fine = size(kpts_fine, 2)
        nq_fine = size(qpts_fine, 2)

        nk1 = kpts_coarse%nkpt3(1)
        nk2 = kpts_coarse%nkpt3(2)
        nk3 = kpts_coarse%nkpt3(3)
        
        nq1 = qpts_coarse%nkpt3(1)
        nq2 = qpts_coarse%nkpt3(2)
        nq3 = qpts_coarse%nkpt3(3)

        ! fourier box limits (Monkhorst box -N/2 .. N/2) for each mesh
        ft_lim_k(2,:) = kpts_coarse%nkpt3(:)/2
        ft_lim_k(1,:) = ft_lim_k(2,:) - kpts_coarse%nkpt3(:) + 1
        ft_lim_q(2,:) = qpts_coarse%nkpt3(:)/2
        ft_lim_q(1,:) = ft_lim_q(2,:) - qpts_coarse%nkpt3(:) + 1

        allocate(matRot(nwann,nwann))
        allocate(fftq(nwann,nwann, nq1*nq2*nq3))
        allocate(fftk(nwann,nwann, nk1*nk2*nk3))
        allocate(tempMat(nwann,nwann, nk_c, 0:nq1-1,0:nq2-1,0:nq3-1))
        allocate(matWannier(nwann,nwann, 0:nk1-1,0:nk2-1,0:nk3-1, 0:nq1-1,0:nq2-1,0:nq3-1))
        allocate(tempMat2(nwann,nwann, 0:nk1-1,0:nk2-1,0:nk3-1))
        matInterpol = cmplx(0.0, 0.0)

        ! Wigner-Seitz weights and big-box supercell R vectors, one set per mesh
        call timestart("Wigner-Seitz weights (k and q)")
        bigBox_lim_k(2,:) =   2*kpts_coarse%nkpt3(:)
        bigBox_lim_k(1,:) = - 2*kpts_coarse%nkpt3(:)
        boxSize_k = (4*nk1+1) * (4*nk2+1) * (4*nk3+1)
        allocate(FTweight_k(boxSize_k), supercellR_k(3,boxSize_k))
        FTweight_k = 0.0; supercellR_k = 0
        iGrid = 1
        do iz=bigBox_lim_k(1,3),bigBox_lim_k(2,3)
            do iy=bigBox_lim_k(1,2),bigBox_lim_k(2,2)
                do ix=bigBox_lim_k(1,1),bigBox_lim_k(2,1)
                    supercellR_k(:,iGrid) = (/ix,iy,iz/)
                    iGrid = iGrid+1
                end do
            end do
        end do

        bigBox_lim_q(2,:) =   2*qpts_coarse%nkpt3(:)
        bigBox_lim_q(1,:) = - 2*qpts_coarse%nkpt3(:)
        boxSize_q = (4*nq1+1) * (4*nq2+1) * (4*nq3+1)
        allocate(FTweight_q(boxSize_q), supercellR_q(3,boxSize_q))
        FTweight_q = 0.0; supercellR_q = 0
        iGrid = 1
        do iz=bigBox_lim_q(1,3),bigBox_lim_q(2,3)
            do iy=bigBox_lim_q(1,2),bigBox_lim_q(2,2)
                do ix=bigBox_lim_q(1,1),bigBox_lim_q(2,1)
                    supercellR_q(:,iGrid) = (/ix,iy,iz/)
                    iGrid = iGrid+1
                end do
            end do
        end do

        cell = fi%cell
        call cell%calculate_WSweight(supercellR_k,FTweight_k,scaleSupercell=kpts_coarse%nkpt3(:))
        call cell%calculate_WSweight(supercellR_q,FTweight_q,scaleSupercell=qpts_coarse%nkpt3(:))
        call timestop("Wigner-Seitz weights (k and q)")

        ! forward fourier transform ( k-space ---> Realspace )
        ! fourier transform of the q-mesh --> R_p
        call timestart("Forward FT + Wannier gauge")
        tempMat = cmplx(0.0,0.0)
        do ikpt = 1 , nk_c
            fftq = cmplx(0.0,0.0)
            do iqpt = 1 , nq_c
                ! fold k+q back into the BZ to find the stored U(k+q)=U(tilde k)
                bkqpt  = kpts_coarse%bk(:,ikpt) + qpts_coarse%bk(:,iqpt)
                ikqpt = kpts_coarse%get_nk(bkqpt)
                if (ikqpt < 1) call juDFT_error("k+q not on the coarse mesh (q must connect k-points)", &
                                              calledby="wannier_matrixq_interpolate")
                ! rotate the matrix element into wannier gauge: U^dagger(k+q) M(k,q) U(k)
                matRot(:,:) = matmul(conjg(transpose(U_mat(:,:,ikqpt))),matmul(matElement(:,:,ikpt,iqpt), U_mat(:,:,ikpt)))
                call ft_dyn_direct(ft_lim_q, 1, qpts_coarse%bk(:,iqpt), matRot, fftq)
            end do !iqpt
            fftq = fftq / nq_c
            call unfold_grid(ft_lim_q, fftq, tempMat(:,:,ikpt,:,:,:))
        end do !ikpt
        call timestop("Forward FT + Wannier gauge")

        ! fourier transfrom k -> R_e, for each R_p grid point ----
        call timestart("Forward FT")
        do iz=0,nq3-1
            do iy=0,nq2-1
                do ix=0,nq1-1
                    fftk = cmplx(0.0,0.0)
                    do ikpt = 1 , nk_c
                        call ft_dyn_direct(ft_lim_k, 1, kpts_coarse%bk(:,ikpt), tempMat(:,:,ikpt,ix,iy,iz), fftk)
                    end do !ikpt
                    fftk = fftk / nk_c
                    call unfold_grid(ft_lim_k, fftk, matWannier(:,:,:,:,:,ix,iy,iz))
                end do
            end do
        end do
        call timestop("Forward FT")

        ! backwards fourier transform ( Realspace --> kspace )
        call timestart("Backward FT")
        do iq_fine = 1 , nq_fine
            ! ---- stage 1: R_p -> q', for each R_e grid point ----
            tempMat2 = cmplx(0.0,0.0)
            do iz=0,nk3-1
                do iy=0,nk2-1
                    do ix=0,nk1-1
                        call ft_fcm_weight(-1, ft_lim_q, bigBox_lim_q, FTweight_q, qpts_fine(:,iq_fine), &
                                           tempMat2(:,:,ix,iy,iz), matWannier(:,:,ix,iy,iz,:,:,:))
                    end do
                end do
            end do
            ! ---- stage 2: R_e -> k' ----
            do ik_fine = 1 , nk_fine
                call ft_fcm_weight(-1, ft_lim_k, bigBox_lim_k, FTweight_k, kpts_fine(:,ik_fine), &
                                   matInterpol(:,:,ik_fine,iq_fine), tempMat2(:,:,:,:,:))
            end do !ik_fine
        end do !iq_fine
        call timestop("Backward FT")

    end subroutine

end module m_matrix_interpolation 