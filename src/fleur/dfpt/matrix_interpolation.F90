!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_matrix_interpolation 
    use m_juDFT
    use m_types
    use m_constants

    implicit none 

contains 
    subroutine wannier_matrix_interpolate(fi,matElement,U_mat,kpts_coarse,kpts_fine,matInterpol,qpts_coarse,qpts_fine)

        use m_dfpt_dynmat_sym , only : ft_dyn_direct, ft_fcm_weight

        type(t_fleurinput), intent(in) :: fi  
        complex, intent(in) :: matElement(:,:,:,:)               ! nu',nu, kpoints, spin
        complex, intent(in) :: U_mat(:,:,:,:)                    ! bloch, wannier, kpoint, jspin 
        type(t_kpts),intent(in) :: kpts_coarse                   ! on the coarse Wannier k-mesh  
        real,intent(in) :: kpts_fine(:,:)                        ! to interpolate onto
        complex,allocatable,intent(out) :: matInterpol(:,:,:,:)  ! interpoalte matrix element  
        type(t_kpts),optional,intent(in) :: qpts_coarse          ! on the coarse Wannier k-mesh  
        type(t_kpts),optional,intent(in) :: qpts_fine            ! to interpolate onto 
        

        integer :: ikpt, iSpin , nu , iNupr, boxSize, nx, ny, nz, iGrid, ix , iy , iz
        integer :: ft_lim(2,3), bigBox_lim(2,3) ! discrete fourier box limits
        type(t_cell) :: cell
        integer, allocatable  :: supercellR(:,:)
        real, allocatable  :: FTweight(:) 

        complex,allocatable :: fft_grid(:,:,:)                     ! matrix elements in realspace on grid (wannier gauge)
        complex,allocatable :: matWannier(:,:,:,:,:,:)               ! matrix elements in realspace (wannier gauge)
        complex,allocatable :: matRot(:,:)                          ! rotated matrix elements in wannier gauge 
        real :: bkpt(3)

        integer :: nwann, nfine
        
        nwann = size(U_mat, 2)
        nfine = size(kpts_fine,2)

        allocate(fft_grid(nwann,nwann,size(matElement,3)))
        allocate(matInterpol(nwann,nwann,nfine,size(matElement,4)))
        allocate(matRot(nwann,nwann))

        fft_grid = cmplx(0.0,0.0)
        matRot = cmplx(0.0,0.0)


        ft_lim(2,:) = kpts_coarse%nkpt3(:)/2
        ft_lim(1,:) = ft_lim(2,:) - kpts_coarse%nkpt3(:) + 1       
        
        allocate(matWannier(nwann,nwann,0:(kpts_coarse%nkpt3(1)-1),&
                            0:(kpts_coarse%nkpt3(2)-1),0:(kpts_coarse%nkpt3(3)-1),size(matElement,4)))
        matWannier = cmplx(0.0,0.0)


        do iSpin = 1 , size(matElement,4)
            do ikpt = 1 , kpts_coarse%nkpt

                bkpt = kpts_coarse%bk(:, ikpt)

                ! rotate the matrix elements into wannier gauge
                ! U^dagger M U
                matRot(:,:) = matmul(conjg(transpose(U_mat(:,:,ikpt,ispin))),matmul(matElement(:,:,ikpt,iSpin),U_mat(:,:,ikpt,ispin)))
                
                call ft_dyn_direct(ft_lim,1,bkpt,matRot,fft_grid(:,:,:))
            end do ! ikpt 
            
            ! unfold the grid to nx,ny,nz indexing 
            fft_grid(:,:,:)=fft_grid(:,:,:)/kpts_coarse%nkpt
            iGrid = 1 
            do iz=ft_lim(1,3),ft_lim(2,3)
                do iy=ft_lim(1,2),ft_lim(2,2)
                    do ix=ft_lim(1,1),ft_lim(2,1)
                    ! shift to storage indices (0-based)
                    nx = ix - ft_lim(1,1)
                    ny = iy - ft_lim(1,2)
                    nz = iz - ft_lim(1,3)
                    matWannier(:,:,nx,ny,nz,iSpin)= fft_grid(:,:,iGrid)
                    iGrid = iGrid + 1 
                    end do !ix
                end do !iy
            end do !iz       
        end do ! iSpin 

        ! interpolate to fine mesh with WS construction 
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

        do ispin = 1 , size(matElement,4)
            do ikpt = 1 , nfine
                call ft_fcm_weight(-1,ft_lim,bigBox_lim,FTweight,kpts_fine(:,ikpt),matInterpol(:,:,ikpt,ispin),matWannier(:,:,:,:,:,iSpin))
            end do !ikpt 
        end do ! ispin 



    end subroutine wannier_matrix_interpolate

end module m_matrix_interpolation 