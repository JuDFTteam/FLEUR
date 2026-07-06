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
    subroutine interpolate_bandstructure(fi, results, kpts_fine)

        use m_available_solvers
        use m_types_solver

        type(t_fleurinput), intent(in) :: fi
        type(t_results),    intent(in) :: results
        real,       intent(in) :: kpts_fine(:,:)

        integer  :: num_bands, num_wann, jspin, ib, ikpt, info, lwork
        real,    allocatable :: eig_interpol(:,:,:)
        real,    allocatable :: rwork(:)
        complex, allocatable :: H_bloch(:,:,:,:)
        complex, allocatable :: H_interpol(:,:,:,:)
        complex, allocatable :: H_interpol_spin(:,:,:,:)
        complex, allocatable :: work(:)
        complex, allocatable :: U_full(:,:,:,:)
        real :: kpath_coord, dk(3) , dk_cart(3)
        integer :: iout , nfine , ioutEv

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

        ! ---------------------------------------------------------------
        ! Step 1: construct diagonal Bloch Hamiltonian H_k from results%eig
        ! H_k(n,m,k,s) = eig(min_band+n-1, k_ibz, s) * delta(n,m)
        ! ---------------------------------------------------------------
        allocate(H_bloch(num_bands, num_bands, fi%kpts%nkptf, fi%input%jspins))
        H_bloch = cmplx(0.0, 0.0)
        do jspin = 1, fi%input%jspins
            do ikpt = 1, fi%kpts%nkptf
                do ib = 1, num_bands
                    H_bloch(ib, ib, ikpt, jspin) = &
                        cmplx(results%eig(fi%wannierlib%min_band + ib - 1, fi%kpts%bkp(ikpt), jspin), 0.0)
                end do
            end do
        end do

        ! ---------------------------------------------------------------
        ! Step 2: Wannier interpolation via wannier_matrix_interpolate
        ! rotate H into Wannier gauge and Fourier interpolate to kpts_fine
        ! ---------------------------------------------------------------
        allocate(H_interpol(num_wann, num_wann, nfine , fi%input%jspins))
        H_interpol = cmplx(0.0, 0.0)
        allocate(U_full(num_bands, num_wann, fi%kpts%nkptf, 1))

        do jspin = 1, fi%input%jspins
            if (num_bands /= num_wann) then
                ! U_full = U_dis @ U_wann  shape: (num_bands x num_wann)
                do ikpt = 1, fi%kpts%nkptf
                    U_full(:,:,ikpt,1) = matmul(results%U_dis(:,:,ikpt,jspin), &
                                                results%U_mat(:,:,ikpt,jspin))
                end do
            else
                U_full(:,:,:,1) = results%U_mat(:,:,:,jspin)
            end if
            call timestart("Wannier Transformation")
            call wannier_matrix_interpolate(fi, H_bloch(:,:,:,jspin:jspin), &
                                            U_full,                          &
                                            fi%kpts, kpts_fine,              &
                                            H_interpol_spin)
            call timestop("Wannier Transformation")
            H_interpol(:,:,:,jspin) = H_interpol_spin(:,:,:,1)
            deallocate(H_interpol_spin)
        end do
        
        deallocate(U_full)

        ! ---------------------------------------------------------------
        ! Step 3: diagonalize interpolated Hamiltonian at each fine k-point
        ! ---------------------------------------------------------------
      
         ! The Wannier-gauge Hamiltonian H_W(k) = sum_R exp(i k.R) H(R) is complex
        ! Hermitian, not real symmetric, even when fi%sym%invs is true: inversion
        ! symmetry only makes the LAPW Hamiltonian real in a special gauge, not the
        ! Wannier gauge produced by w90. Always diagonalize the complex matrix.
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
                if (allocated(zMat_tmp)) then
                    call zMat_tmp%free()
                    deallocate(zMat_tmp)
                end if 
            end do 
        end do !jspins

        !deallocate(H_interpol)

        ! ---------------------------------------------------------------
        ! Write interpolated eigenvalues to wann_bands.dat
        ! Format: kpath_coord  eig_1  eig_2  ...  (one line per k-point)
        ! ---------------------------------------------------------------
        open(newunit=iout, file='wann_bands.dat', status='replace')
        write(iout,'(a)') '# kpath_coord   eig_1  eig_2  ...'

        open(newunit=ioutEv, file='wann_bands_ev.dat', status='replace')
        write(ioutEV,'(a)') '# kpath_coord   eig_1  eig_2  ...'

        do jspin = 1, fi%input%jspins
            if (fi%input%jspins > 1) write(iout,'(a,i1)') '# spin ', jspin
            kpath_coord = 0.0
            do ikpt = 1, nfine 
                if (ikpt > 1) then
                    dk = kpts_fine(:,ikpt) - kpts_fine(:,ikpt-1)
                    dk_cart = matmul(fi%cell%bmat,dk)
                    kpath_coord = kpath_coord + sqrt(dot_product(dk_cart, dk_cart))
                end if
                write(iout,'(f12.6,*(2x,f14.8))') kpath_coord, eig_interpol(:,ikpt,jspin)
                write(ioutEv,'(f12.6,*(2x,f14.8))') kpath_coord, hartree_to_ev_const*eig_interpol(:,ikpt,jspin)
            end do
        end do

        close(iout)
        close(ioutEV)
        deallocate(eig_interpol)

    end subroutine interpolate_bandstructure


    subroutine interpolate_bandstructure_q(fi, results, kpts_fine)

        use m_available_solvers
        use m_types_solver

        type(t_fleurinput), intent(in) :: fi
        type(t_results),    intent(in) :: results
        real,       intent(in) :: kpts_fine(:,:)

        integer  :: num_bands, num_wann, jspin, ib, ikpt, info, lwork
        real,    allocatable :: eig_interpol(:,:,:)
        real,    allocatable :: rwork(:)
        complex, allocatable :: H_bloch(:,:,:,:,:)
        complex, allocatable :: H_interpol(:,:,:,:)
        complex, allocatable :: H_interpol_spin(:,:,:,:,:)
        complex, allocatable :: work(:)
        complex, allocatable :: U_full(:,:,:,:)
        real :: kpath_coord, dk(3) , dk_cart(3)
        integer :: iout , nfine , ioutEv

        class(t_mat),allocatable :: hmat_tmp
        class(t_mat),allocatable :: zMat_tmp
        class(t_solver), allocatable   :: solver, transform 

        type(t_kpts) :: qpts_coarse

        real :: bqpt(3,1)

        if (.not. fi%wannierlib%l_wannierize) return
        if (.not. allocated(results%U_mat)) &
            call juDFT_error('U_mat not allocated; run wannierization first', &
                             calledby='interpolate_bandstructure')

        num_bands = fi%wannierlib%max_band - fi%wannierlib%min_band + 1
        num_wann  = fi%wannierlib%num_wann
        nfine = size(kpts_fine,2)

        bqpt= 0.0 

        qpts_coarse = fi%kpts

        if (num_bands /= num_wann .and. .not. allocated(results%U_dis)) &
            call juDFT_error('U_dis not allocated; disentanglement data missing', &
                             calledby='interpolate_bandstructure')

        ! ---------------------------------------------------------------
        ! Step 1: construct diagonal Bloch Hamiltonian H_k from results%eig
        ! H_k(n,m,k,s) = eig(min_band+n-1, k_ibz, s) * delta(n,m)
        ! ---------------------------------------------------------------
        allocate(H_bloch(num_bands, num_bands, fi%kpts%nkptf, qpts_coarse%nkptf ,fi%input%jspins))
        H_bloch = cmplx(0.0, 0.0)
        do jspin = 1, fi%input%jspins
            do ikpt = 1, fi%kpts%nkptf
                do ib = 1, num_bands
                    H_bloch(ib, ib, ikpt, : ,jspin) = &
                        cmplx(results%eig(fi%wannierlib%min_band + ib - 1, fi%kpts%bkp(ikpt) , jspin), 0.0)
                end do
            end do
        end do

        ! ---------------------------------------------------------------
        ! Step 2: Wannier interpolation via wannier_matrix_interpolate
        ! rotate H into Wannier gauge and Fourier interpolate to kpts_fine
        ! ---------------------------------------------------------------
        allocate(H_interpol(num_wann, num_wann, nfine, fi%input%jspins))
        H_interpol = cmplx(0.0, 0.0)
        allocate(U_full(num_bands, num_wann, fi%kpts%nkptf, 1))

        do jspin = 1, fi%input%jspins
            if (num_bands /= num_wann) then
                ! U_full = U_dis @ U_wann  shape: (num_bands x num_wann)
                do ikpt = 1, fi%kpts%nkptf
                    U_full(:,:,ikpt,1) = matmul(results%U_dis(:,:,ikpt,jspin), &
                                                results%U_mat(:,:,ikpt,jspin))
                end do
            else
                U_full(:,:,:,1) = results%U_mat(:,:,:,jspin)
            end if
            call timestart("Wannier Transformation")
            call wannier_matrixq_interpolate(fi, H_bloch(:,:,:,:,jspin:jspin), &
                                            U_full,                          &
                                            fi%kpts, kpts_fine,              &
                                            H_interpol_spin,qpts_coarse,bqpt)
            call timestop("Wannier Transformation")
            H_interpol(:,:,:,jspin) = H_interpol_spin(:,:,:,1,1)
            deallocate(H_interpol_spin)
        end do
        
        deallocate(U_full)

        ! ---------------------------------------------------------------
        ! Step 3: diagonalize interpolated Hamiltonian at each fine k-point
        ! ---------------------------------------------------------------
      
         ! The Wannier-gauge Hamiltonian H_W(k) = sum_R exp(i k.R) H(R) is complex
        ! Hermitian, not real symmetric, even when fi%sym%invs is true: inversion
        ! symmetry only makes the LAPW Hamiltonian real in a special gauge, not the
        ! Wannier gauge produced by w90. Always diagonalize the complex matrix.
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
                if (allocated(zMat_tmp)) then
                    call zMat_tmp%free()
                    deallocate(zMat_tmp)
                end if 
            end do 
        end do !jspins

        !deallocate(H_interpol)

        ! ---------------------------------------------------------------
        ! Write interpolated eigenvalues to wann_bands.dat
        ! Format: kpath_coord  eig_1  eig_2  ...  (one line per k-point)
        ! ---------------------------------------------------------------
        open(newunit=iout, file='wann_bands_q.dat', status='replace')
        write(iout,'(a)') '# kpath_coord   eig_1  eig_2  ...'

        open(newunit=ioutEv, file='wann_bands_ev_q.dat', status='replace')
        write(ioutEV,'(a)') '# kpath_coord   eig_1  eig_2  ...'

        do jspin = 1, fi%input%jspins
            if (fi%input%jspins > 1) write(iout,'(a,i1)') '# spin ', jspin
            kpath_coord = 0.0
            do ikpt = 1, nfine 
                if (ikpt > 1) then
                    dk = kpts_fine(:,ikpt) - kpts_fine(:,ikpt-1)
                    dk_cart = matmul(fi%cell%bmat,dk)
                    kpath_coord = kpath_coord + sqrt(dot_product(dk_cart, dk_cart))
                end if
                write(iout,'(f12.6,*(2x,f14.8))') kpath_coord, eig_interpol(:,ikpt,jspin)
                write(ioutEv,'(f12.6,*(2x,f14.8))') kpath_coord, hartree_to_ev_const*eig_interpol(:,ikpt,jspin)
            end do
        end do

        close(iout)
        close(ioutEV)
        deallocate(eig_interpol)

    end subroutine interpolate_bandstructure_q

    subroutine test_q_inter(fi, results, kpts_fine)

        use m_available_solvers
        use m_types_solver
        use m_constants

        type(t_fleurinput), intent(in) :: fi
        type(t_results),    intent(in) :: results
        real,       intent(in) :: kpts_fine(:,:)

        integer  :: num_bands, num_wann, jspin, ib, ikpt, info, lwork
        real,    allocatable :: eig_interpol(:,:,:)
        real,    allocatable :: rwork(:)
        complex, allocatable :: H_bloch(:,:,:,:,:)
        complex, allocatable :: H_interpol(:,:,:,:)
        complex, allocatable :: H_interpol_spin(:,:,:,:,:)
        complex, allocatable :: work(:)
        complex, allocatable :: U_full(:,:,:,:)
        real :: kpath_coord, dk(3) , dk_cart(3)
        integer :: iout , nfine , ioutEv, iqpt

        class(t_mat),allocatable :: hmat_tmp
        class(t_mat),allocatable :: zMat_tmp
        class(t_solver), allocatable   :: solver, transform 

        type(t_kpts) :: qpts_coarse
        type(t_fleurinput) :: fi_local

        real :: bqpt(3,1)

        if (.not. fi%wannierlib%l_wannierize) return
        if (.not. allocated(results%U_mat)) &
            call juDFT_error('U_mat not allocated; run wannierization first', &
                             calledby='interpolate_bandstructure')

        num_bands = 9 
        num_wann  = 9 

        fi_local = fi 

        fi_local%wannierlib%num_wann= 9 
        fi_local%wannierlib%num_bands = 9 
        nfine = size(kpts_fine,2)

        bqpt= 1.0/3.0

        qpts_coarse = fi%kpts

        if (num_bands /= num_wann .and. .not. allocated(results%U_dis)) &
            call juDFT_error('U_dis not allocated; disentanglement data missing', &
                             calledby='interpolate_bandstructure')

        ! ---------------------------------------------------------------
        ! Step 1: construct diagonal Bloch Hamiltonian H_k from results%eig
        ! H_k(n,m,k,s) = eig(min_band+n-1, k_ibz, s) * delta(n,m)
        ! ---------------------------------------------------------------
        allocate(H_bloch(num_bands, num_bands, fi%kpts%nkptf, qpts_coarse%nkptf ,fi%input%jspins))
        H_bloch = cmplx(0.0, 0.0)
        do jspin = 1, fi%input%jspins
            do iqpt = 1 , qpts_coarse%nkptf
                do ikpt = 1, fi%kpts%nkptf
                    do ib = 1, num_bands
                        H_bloch(ib, ib, ikpt, iqpt ,jspin) = 2*ib + sin(tpi_const*(fi%kpts%bkf(1,ikpt)+qpts_coarse%bkf(1,iqpt)))*sin(tpi_const*(fi%kpts%bkf(2,ikpt)+qpts_coarse%bkf(2,iqpt))) &
                                                                    + ib* sin(tpi_const*(fi%kpts%bkf(3,ikpt)+qpts_coarse%bkf(3,iqpt))) + cos(tpi_const*(qpts_coarse%bkf(1,iqpt))  )
                    end do
                end do
            end do 
        end do


        call save_npy("Bloch_funktion.npy", H_bloch)
        ! ---------------------------------------------------------------
        ! Step 2: Wannier interpolation via wannier_matrix_interpolate
        ! rotate H into Wannier gauge and Fourier interpolate to kpts_fine
        ! ---------------------------------------------------------------
        allocate(H_interpol(num_wann, num_wann, nfine, fi%input%jspins))
        H_interpol = cmplx(0.0, 0.0)
        allocate(U_full(num_bands, num_wann, fi%kpts%nkptf, 1))
        U_full = 0.0 
        do jspin = 1, fi%input%jspins
            do ib = 1 , num_bands
                U_full(ib,ib,:,:) = 1.0 
            end do 
            
            call timestart("Wannier Transformation")
            call wannier_matrixq_interpolate(fi, H_bloch(:,:,:,:,jspin:jspin), &
                                            U_full,                          &
                                            fi%kpts, kpts_fine,              &
                                            H_interpol_spin,qpts_coarse,bqpt)
            call timestop("Wannier Transformation")
            H_interpol(:,:,:,jspin) = H_interpol_spin(:,:,:,1,1)
            deallocate(H_interpol_spin)
        end do

        call save_npy("U_mat.npy",U_full)
        call save_npy("interpoled.npy",h_interpol)
        call save_npy("k_interpol.py",kpts_fine)
        call save_npy("k_mesh.py",fi%kpts%bkf)


        
        deallocate(U_full)

    end subroutine test_q_inter

end module m_wannier_interpolate
