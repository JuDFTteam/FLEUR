!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_wannier_interpolate

    use m_juDFT
    use m_types
    use m_constants
    use m_matrix_interpolation
    use m_npy
    implicit none

contains
    subroutine interpolate_bandstructure(fi, results, kpts_fine,eig_id_interpol,l_write_output)

        use m_available_solvers
        use m_types_solver
        use m_eig66_io
        use m_banddos_io
#ifdef CPP_HDF
        use hdf5
        use m_hdf_tools
#endif

        type(t_fleurinput), intent(in) :: fi
        type(t_results),    intent(in) :: results
        real,       intent(in) :: kpts_fine(:,:)
        integer,    intent(in) :: eig_id_interpol
        logical,    intent(in) :: l_write_output 

        integer  :: num_bands, num_wann, jspin, ib, ikpt, info, lwork
        real,    allocatable :: eig_interpol(:,:,:)
        real,    allocatable :: rwork(:)
        real,    allocatable :: band_weight(:,:,:)
        complex, allocatable :: H_bloch(:,:,:,:)
        complex, allocatable :: H_interpol(:,:,:,:)
        complex, allocatable :: work(:)
        complex, allocatable :: U_full(:,:,:)
        integer :: nfine

        type(t_kpts)    :: kpts_band
        type(t_banddos) :: banddosLocal
#ifdef CPP_HDF
        integer(HID_T) :: banddosFile_id
#endif

        class(t_mat),allocatable :: hmat_tmp
        class(t_mat),allocatable :: zMat_tmp
        class(t_solver), allocatable   :: solver, transform

        if (.not. fi%wannierlib%l_wannierize) return
        if (.not. allocated(results%U_mat)) &
            call juDFT_error('U_mat not allocated; run wannierization first', &
                             calledby='interpolate_bandstructure')

        num_bands = fi%wannierlib%max_band - fi%wannierlib%min_band + 1
        num_wann  = fi%wannierlib%num_wann
        nfine = size(kpts_fine,2)

        if (num_bands /= num_wann .and. .not. allocated(results%U_dis)) &
            call juDFT_error('U_dis not allocated; disentanglement data missing', &
                             calledby='interpolate_bandstructure')


        ! construct diagonal Bloch Hamiltonian H_k
        allocate(H_bloch(num_bands, num_bands, fi%kpts%nkptf, fi%input%jspins))
        H_bloch = cmplx(0.0, 0.0)
        call timestart("Setup Bloch Hamiltonain")
        do jspin = 1, fi%input%jspins
            do ikpt = 1, fi%kpts%nkptf
                do ib = 1, num_bands
                    H_bloch(ib, ib, ikpt, jspin) = &
                        cmplx(results%eig(fi%wannierlib%min_band + ib - 1, fi%kpts%bkp(ikpt), jspin), 0.0)
                end do
            end do
        end do
        call timestop("Setup Bloch Hamiltonain")


        ! Wannier interpolation
        allocate(H_interpol(num_wann, num_wann, nfine , fi%input%jspins))
        H_interpol = cmplx(0.0, 0.0)
        allocate(U_full(num_bands, num_wann, fi%kpts%nkptf))

        do jspin = 1, fi%input%jspins
            if (num_bands /= num_wann) then
                ! U_full = U_dis @ U_wann  shape: (num_bands x num_wann)
                do ikpt = 1, fi%kpts%nkptf
                    U_full(:,:,ikpt) = matmul(results%U_dis(:,:,ikpt,jspin), &
                                              results%U_mat(:,:,ikpt,jspin))
                end do
            else
                U_full(:,:,:) = results%U_mat(:,:,:,jspin)
            end if
            call timestart("Wannier Interpolation")
            call wannier_matrix_interpolate(fi, H_bloch(:,:,:,jspin), &
                                            U_full,                      &
                                            fi%kpts, kpts_fine,          &
                                            H_interpol(:,:,:,jspin))
            call timestop("Wannier Interpolation")
        end do
        
        deallocate(U_full)
        deallocate(H_bloch)

  
        ! diagonalize to get the eigenvalues 
        ! In wannier gauge, hamiltonian and zmat complex, even if its real 
        ! beforehand 
        allocate(eig_interpol(num_wann, nfine, fi%input%jspins))
        allocate (t_mat::hmat_tmp)
        call hmat_tmp%init(.false.,num_wann,num_wann)
        call select_solver(.false.,diag_solver=solver,diag_transform=transform)

        do jspin = 1 , fi%input%jspins
            do ikpt  = 1 , nfine 
                hmat_tmp%data_c = H_interpol(:,:,ikpt,jspin)
                call timestart("Diogonalization Wannier")
                call solver%solve_std_dp(hmat_tmp,num_wann,eig_interpol(:,ikpt,jspin),zMat_tmp)
                call timestop("Diogonalization Wannier")
                call write_eig(eig_id_interpol, ikpt, jspin, num_wann, num_wann, eig_interpol(:,ikpt,jspin),zMat=zMat_tmp)
                if (allocated(zMat_tmp)) then
                    call zMat_tmp%free()
                    deallocate(zMat_tmp)
                end if 
            end do !ikpt 
        end do !jspins


        if (l_write_output) then
            ! Write interpolated eigenvalues into banddos.hdf (/Local/BS/eigenvalues).            
            ! Minimal band-path k-points for banddos.hdf. Mark the two endpoints as special points
            ! (labeled "S"/"E") so /kpts/specialPointLabels is non-empty for the plotting tools.
            kpts_band%nkpt = nfine
            kpts_band%bk   = kpts_fine
            allocate(kpts_band%wtkpt(nfine))
            kpts_band%wtkpt = 1.0 / nfine

            kpts_band%numSpecialPoints = 2
            allocate(kpts_band%specialPointIndices(2))
            allocate(kpts_band%specialPointNames(2))
            allocate(kpts_band%specialPoints(3,2))
            kpts_band%specialPointIndices(1) = 1
            kpts_band%specialPointNames(1)   = 'S'
            kpts_band%specialPoints(:,1)     = kpts_fine(:,1)
            kpts_band%specialPointIndices(2) = nfine
            kpts_band%specialPointNames(2)   = 'E'
            kpts_band%specialPoints(:,2)     = kpts_fine(:,nfine)

#ifdef CPP_HDF
            banddosLocal      = fi%banddos
            banddosLocal%band = .true.
            ! For now use a default weight for projection
            allocate(band_weight(num_wann, nfine, fi%input%jspins))
            band_weight = 1.0
            call openBandDOSFile(banddosFile_id, fi%input, fi%atoms, fi%cell, kpts_band, &
                                fi%sym, banddosLocal, results%ef)
            call writeBandData(banddosFile_id, kpts_band, 'Local', 'Total', &
                            band_weight, eig_interpol)
            call closeBandDOSFile(banddosFile_id)
            deallocate(band_weight)
#endif
        end if 
        deallocate(eig_interpol)

    end subroutine interpolate_bandstructure

end module m_wannier_interpolate
