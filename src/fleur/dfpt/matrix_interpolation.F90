!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_matrix_interpolation
#ifdef CPP_MPI
    use mpi
#endif
    use m_juDFT
    use m_constants
    use m_npy
    use m_dfpt_dynmat_fourier, only: ft_ws_forward, build_ws_ft, ft_fcm_weight_packed, build_ws_shifted
    use m_types_cell
    use m_types_fleurinput
    use m_types_kpts

    implicit none
    private
    public :: wannier_matrix_interpolate, wannier_matrixq_interpolate, wannier_matrixq_forward, &
       wannier_matrixq_backward, t_wann_ft

    ! Carrier for the q/fine-k-independent part of the double-mesh Wannier
    ! interpolation (real-space tensor + Wigner-Seitz data).
    type t_wann_ft
        integer :: nwann
        complex, allocatable :: matWannier(:,:,:,:)   ! nwann,nwann, 0:ws_k_nNZ-1 (R_e), 0:ws_q_nNZ-1 (R_p)
        ! Wigner-Seitz vectors *_Rvecs(:,i), stored at index i-1
        integer :: ws_k_nNZ = 0, ws_q_nNZ = 0
        integer, allocatable :: ws_k_Rvecs(:,:)
        integer, allocatable :: ws_q_Rvecs(:,:)
        real,    allocatable :: ws_k_weight(:,:,:)    ! (nwann,nwann,ws_k_nNZ) for the bond R_e + c_a - c_b
        real,    allocatable :: ws_q_weight(:,:)      ! (nwann,ws_q_nNZ) for the bond R_p + c_a - tau of the displaced atom
    end type t_wann_ft

contains
    subroutine wannier_matrix_interpolate(fi,matElement,U_mat,kpts_coarse,kpts_fine,matInterpol,qpts_coarse,qpts_fine)


        type(t_fleurinput), intent(in) :: fi
        complex, intent(in) :: matElement(:,:,:)                 ! nu',nu, kpoints
        complex, intent(in) :: U_mat(:,:,:)                      ! bloch, wannier, kpoint
        type(t_kpts),intent(in) :: kpts_coarse                   ! on the coarse Wannier k-mesh
        real,intent(in) :: kpts_fine(:,:)                        ! to interpolate onto
        complex,intent(inout) :: matInterpol(:,:,:)    ! interpolate matrix element  (nwann,nwann,kpts)
        type(t_kpts),optional,intent(in) :: qpts_coarse          ! on the coarse Wannier k-mesh
        type(t_kpts),optional,intent(in) :: qpts_fine            ! to interpolate onto


        integer :: ikpt , boxSize, nx, ny, nz, iGrid
        integer :: bigBox_lim(2,3)
        type(t_cell) :: cell
        integer, allocatable  :: supercellR(:,:)
        real, allocatable  :: FTweight(:)

        complex,allocatable :: matWannier(:,:,:)                    ! matrix elements in realspace (wannier gauge)
        complex,allocatable :: matRot(:,:)                          ! rotated matrix elements in wannier gauge

        integer :: nNZ                              ! packed nonzero-weight list (built once)
        integer, allocatable :: Rvecs(:,:)
        real,    allocatable :: weightNZ(:)

        integer :: nwann, nfine

        nwann = size(U_mat, 2)
        nfine = size(kpts_fine,2)

        matInterpol = cmplx(0.0, 0.0)
        allocate(matRot(nwann,nwann))
        matRot = cmplx(0.0,0.0)

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
        call build_ws_ft(bigBox_lim, FTweight, nNZ, Rvecs, weightNZ)
        call timestop("Wigner-Seitz weights")

        allocate(matWannier(nwann,nwann,0:nNZ-1))
        matWannier = cmplx(0.0,0.0)

        call timestart("Forward FT + Wannier gauge")
        do ikpt = 1 , kpts_coarse%nkpt
            ! rotate the matrix elements into wannier gauge
            ! U^dagger M U
            matRot(:,:) = matmul(conjg(transpose(U_mat(:,:,ikpt))),matmul(matElement(:,:,ikpt),U_mat(:,:,ikpt)))

            call ft_ws_forward(nNZ,Rvecs,kpts_coarse%bk(:,ikpt),matRot,matWannier)
        end do ! ikpt

        matWannier(:,:,:)=matWannier(:,:,:)/kpts_coarse%nkpt
        call timestop("Forward FT + Wannier gauge")

        call timestart("Backward FT")
        do ikpt = 1 , nfine
            call ft_fcm_weight_packed(-1,nNZ,Rvecs,weightNZ,kpts_fine(:,ikpt),matInterpol(:,:,ikpt),matWannier)
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
        real, allocatable :: centres(:,:)

        allocate(centres(3,size(U_mat,2)), source=0.0)
        call wannier_matrixq_forward(fi, matElement, U_mat, kpts_coarse, qpts_coarse, centres, ft)
        call wannier_matrixq_backward(ft, kpts_fine, qpts_fine, matInterpol)

    end subroutine wannier_matrixq_interpolate

    subroutine wannier_matrixq_forward(fi,matElement,U_mat,kpts_coarse,qpts_coarse,centres,ft,k_list,mpi_comm,iAtom)

        ! Build the real-space Wannier-gauge matWannier and the Wigner-Seitz weights of each Wannier pair (R_e) and bra function (R_p).
        ! With k_list/mpi_comm, matElement holds only the k-points k_list and matWannier is summed over all ranks.
        ! With iAtom, R_p runs over the union of the per-atom lists and carries the weights of the displaced atom iAtom.

        type(t_fleurinput), intent(in)  :: fi
        complex,            intent(in)  :: matElement(:,:,:,:)    ! nu',nu, local kpoints, qpts (nu' at k+q, nu at k)
        complex,            intent(in)  :: U_mat(:,:,:)           ! bloch, wannier, kpoint (full coarse mesh)
        type(t_kpts),       intent(in)  :: kpts_coarse            ! coarse Wannier k-mesh
        type(t_kpts),       intent(in)  :: qpts_coarse            ! coarse Wannier q-mesh, summed over its full zone (nkptf/bkf)
        real,               intent(in)  :: centres(:,:)           ! (3,nwann) Wannier centres, cartesian bohr
        type(t_wann_ft),    intent(out) :: ft
        integer, optional,  intent(in)  :: k_list(:)              ! global k index of each local slot of matElement
        integer, optional,  intent(in)  :: mpi_comm
        integer, optional,  intent(in)  :: iAtom                  ! atom displaced by this matrix element

        integer :: ikpt, nk_i, iqpt, ikqpt, iqR, iw, jw, jAtom, nAtom, iAtomW
        integer :: nwann, nk_c, nq_c
        real    :: tau(3)
        real,    allocatable :: sites(:,:), shifts(:,:), weightShift(:,:)
        integer, allocatable :: kpt_list(:)
        complex, allocatable :: matRot(:,:)                       ! wannier-gauge block for one (k,q)
        complex, allocatable :: matRpK(:,:,:,:)                  ! after q->R_p FT (nwann,nwann, R_p, local k)
        real :: bkqpt(3)
#ifdef CPP_MPI
        integer :: ierr
#endif

        nwann = size(U_mat, 2)
        nk_c  = kpts_coarse%nkpt
        nq_c  = qpts_coarse%nkptf

        if (present(k_list)) then
            kpt_list = k_list
        else
            kpt_list = [(ikpt, ikpt = 1, nk_c)]
        end if
        if (size(matElement,3) /= size(kpt_list)) call juDFT_error("matElement k dimension does not match the k list", &
                                                              calledby="wannier_matrixq_forward")

        ft%nwann = nwann

        ! The FT conventions put the bra function a (at k+q) at the origin, the ket function b at -R_e + c_b
        ! and the displaced atom at -R_p + tau, so the bonds are R_e + c_a - c_b and R_p + c_a - tau.
        call timestart("Wigner-Seitz weights (k and q)")
        sites = wannier_sites(fi, centres)
        allocate(shifts(3,nwann*nwann))
        do jw = 1, nwann
            do iw = 1, nwann
                shifts(:,iw+nwann*(jw-1)) = sites(:,iw) - sites(:,jw)
            end do
        end do
        call build_ws_shifted(fi%cell%amat, kpts_coarse%nkpt3, fi%input%film, shifts, ft%ws_k_nNZ, ft%ws_k_Rvecs, weightShift)
        ft%ws_k_weight = reshape(weightShift, [nwann, nwann, ft%ws_k_nNZ])
        deallocate(shifts, weightShift)

        nAtom = 1
        if (present(iAtom)) nAtom = fi%atoms%nat
        allocate(shifts(3,nwann*nAtom))
        do jAtom = 1, nAtom
            tau = 0.0
            if (present(iAtom)) tau = matmul(fi%cell%amat, fi%atoms%taual(:,jAtom))
            do iw = 1, nwann
                shifts(:,iw+nwann*(jAtom-1)) = sites(:,iw) - tau
            end do
        end do
        call build_ws_shifted(fi%cell%amat, qpts_coarse%nkpt3, fi%input%film, shifts, ft%ws_q_nNZ, ft%ws_q_Rvecs, weightShift)
        iAtomW = 1
        if (present(iAtom)) iAtomW = iAtom
        ft%ws_q_weight = weightShift(nwann*(iAtomW-1)+1:nwann*iAtomW,:)
        call timestop("Wigner-Seitz weights (k and q)")

        allocate(matRot(nwann,nwann))
        allocate(matRpK(nwann,nwann, 0:ft%ws_q_nNZ-1, size(kpt_list)))
        allocate(ft%matWannier(nwann,nwann, 0:ft%ws_k_nNZ-1, 0:ft%ws_q_nNZ-1))
        ft%matWannier = cmplx(0.0,0.0)

        ! forward fourier transform ( k-space ---> Realspace ) onto the WS vectors
        ! fourier transform of the q-mesh --> R_p
        call timestart("Forward FT + Wannier gauge")
        matRpK = cmplx(0.0,0.0)
        do nk_i = 1 , size(kpt_list)
            ikpt = kpt_list(nk_i)
            do iqpt = 1 , nq_c
                ! fold k+q back into the BZ to find the stored U(k+q)=U(tilde k)
                bkqpt  = kpts_coarse%bk(:,ikpt) + qpts_coarse%bkf(:,iqpt)
                ikqpt = kpts_coarse%get_nk(bkqpt)
                if (ikqpt < 1) call juDFT_error("k+q not on the coarse mesh (q must connect k-points)", &
                                              calledby="wannier_matrixq_forward")
                ! rotate the matrix element into wannier gauge: U^dagger(k+q) M(k,q) U(k)
                matRot(:,:) = matmul(conjg(transpose(U_mat(:,:,ikqpt))),matmul(matElement(:,:,nk_i,iqpt), U_mat(:,:,ikpt)))
                call ft_ws_forward(ft%ws_q_nNZ, ft%ws_q_Rvecs, qpts_coarse%bkf(:,iqpt), matRot, matRpK(:,:,:,nk_i))
            end do !iqpt
            matRpK(:,:,:,nk_i) = matRpK(:,:,:,nk_i) / nq_c
        end do !nk_i
        call timestop("Forward FT + Wannier gauge")

        ! fourier transfrom k -> R_e, for each R_p grid point ----
        call timestart("Forward FT")
        do iqR = 0, ft%ws_q_nNZ-1
            do nk_i = 1 , size(kpt_list)
                call ft_ws_forward(ft%ws_k_nNZ, ft%ws_k_Rvecs, kpts_coarse%bk(:,kpt_list(nk_i)), matRpK(:,:,iqR,nk_i), ft%matWannier(:,:,:,iqR))
            end do !nk_i
            ft%matWannier(:,:,:,iqR) = ft%matWannier(:,:,:,iqR) / nk_c
        end do
        call timestop("Forward FT")

#ifdef CPP_MPI
        if (present(mpi_comm)) then
            call timestart("Forward FT reduction")
            call MPI_ALLREDUCE(MPI_IN_PLACE, ft%matWannier, size(ft%matWannier), MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm, ierr)
            call timestop("Forward FT reduction")
        end if
#endif

    end subroutine wannier_matrixq_forward

    subroutine wannier_matrixq_backward(ft,kpts_fine,qpts_fine,matInterpol)

        ! Backward Fourier transform of the precomputed real-space to fine k space
    

        type(t_wann_ft), intent(in)    :: ft
        real,            intent(in)    :: kpts_fine(:,:)          ! fine k-mesh to interpolate onto
        real,            intent(in)    :: qpts_fine(:,:)          ! fine q-mesh to interpolate onto
        complex,         intent(inout) :: matInterpol(:,:,:,:)    ! interpolated matrix element (nwann,nwann,kpts,qpts)

        integer :: ik_fine, iq_fine, ikR, nk_fine, nq_fine
        complex, allocatable :: matReQ(:,:,:)                  ! after R_p->q' FT (nwann,nwann, R_e)

        nk_fine = size(kpts_fine, 2)
        nq_fine = size(qpts_fine, 2)

        allocate(matReQ(ft%nwann,ft%nwann, 0:ft%ws_k_nNZ-1))
        matInterpol = cmplx(0.0, 0.0)

        ! backwards fourier transform ( Realspace --> kspace )
        call timestart("Backward FT")
        do iq_fine = 1 , nq_fine
            ! ---- stage 1: R_p -> q', for each R_e grid point ----
            matReQ = cmplx(0.0,0.0)
            do ikR = 0, ft%ws_k_nNZ-1
                call ft_backward_rows(ft%ws_q_nNZ, ft%ws_q_Rvecs, ft%ws_q_weight, qpts_fine(:,iq_fine), &
                                      matReQ(:,:,ikR), ft%matWannier(:,:,ikR,:))
            end do
            ! ---- stage 2: R_e -> k' ----
            do ik_fine = 1 , nk_fine
                call ft_backward_pairs(ft%ws_k_nNZ, ft%ws_k_Rvecs, ft%ws_k_weight, kpts_fine(:,ik_fine), &
                                       matInterpol(:,:,ik_fine,iq_fine), matReQ)
            end do !ik_fine
        end do !iq_fine
        call timestop("Backward FT")

    end subroutine wannier_matrixq_backward

    function wannier_sites(fi, centres) result(sites)
        !! Wannier centres moved onto the nearest atom image within snapRadius; centres farther from every atom are kept.
        type(t_fleurinput), intent(in) :: fi
        real,               intent(in) :: centres(:,:)    ! (3,nwann) cartesian bohr
        real, allocatable :: sites(:,:)

        real, parameter :: snapRadius = 1.0
        integer :: iw, iAtom, ix, iy, iz
        real    :: pos(3), dist, distMin

        sites = centres
        do iw = 1, size(centres,2)
            distMin = snapRadius
            do iAtom = 1, fi%atoms%nat
                do iz = -2, 2
                    do iy = -2, 2
                        do ix = -2, 2
                            pos = matmul(fi%cell%amat, fi%atoms%taual(:,iAtom) + [ix,iy,iz])
                            dist = norm2(pos - centres(:,iw))
                            if (dist < distMin) then
                                distMin = dist
                                sites(:,iw) = pos
                            end if
                        end do
                    end do
                end do
            end do
            if (distMin >= snapRadius) write(oUnit,'(a,i0,a)') "Wannier function ", iw, &
                                                              ": centre not on an atom, Wigner-Seitz weights use the Wannier90 centre"
        end do
    end function wannier_sites

    subroutine ft_backward_rows(nNZ, Rvecs, weight, bqpt, mat_q, mat_r)
        !! Backward FT r -> q with the weight weight(a,i) on row a of the block at Rvecs(:,i), stored at index i-1.
        integer, intent(in)    :: nNZ
        integer, intent(in)    :: Rvecs(:,:)
        real,    intent(in)    :: weight(:,:)
        real,    intent(in)    :: bqpt(3)
        complex, intent(inout) :: mat_q(:,:)
        complex, intent(in)    :: mat_r(:,:,0:)

        integer :: i, iw
        real    :: phas
        complex :: phase_fac

        do i = 1, nNZ
            phas = -tpi_const*dot_product(bqpt, real(Rvecs(:,i)))
            phase_fac = cmplx(cos(phas), sin(phas))
            do iw = 1, size(weight,1)
                if (weight(iw,i) == 0.0) cycle
                mat_q(iw,:) = mat_q(iw,:) + phase_fac*weight(iw,i)*mat_r(iw,:,i-1)
            end do
        end do
    end subroutine ft_backward_rows

    subroutine ft_backward_pairs(nNZ, Rvecs, weight, bqpt, mat_q, mat_r)
        !! Backward FT r -> q with the weight weight(a,b,i) on element (a,b) of the block at Rvecs(:,i), stored at index i-1.
        integer, intent(in)    :: nNZ
        integer, intent(in)    :: Rvecs(:,:)
        real,    intent(in)    :: weight(:,:,:)
        real,    intent(in)    :: bqpt(3)
        complex, intent(inout) :: mat_q(:,:)
        complex, intent(in)    :: mat_r(:,:,0:)

        integer :: i
        real    :: phas

        do i = 1, nNZ
            phas = -tpi_const*dot_product(bqpt, real(Rvecs(:,i)))
            mat_q = mat_q + cmplx(cos(phas), sin(phas))*weight(:,:,i)*mat_r(:,:,i-1)
        end do
    end subroutine ft_backward_pairs

end module m_matrix_interpolation 