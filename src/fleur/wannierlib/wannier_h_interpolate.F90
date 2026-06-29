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

    implicit none

contains
    subroutine interpolate_bandstructure(fi, results, kpts_fine)

        type(t_fleurinput), intent(in) :: fi
        type(t_results),    intent(in) :: results
        real,       intent(in) :: kpts_fine(:,:)

        integer  :: num_bands, num_wann, jspin, ib, ikpt, info, lwork
        real,    allocatable :: eig_interpol(:,:,:)
        real,    allocatable :: rwork(:)
        complex, allocatable :: H_bloch(:,:,:,:)
        complex, allocatable :: H_interpol(:,:,:,:)
        complex, allocatable :: H_interpol_spin(:,:,:,:)
        complex, allocatable :: hmat_tmp(:,:)
        complex, allocatable :: work(:)
        complex, allocatable :: U_full(:,:,:,:)
        real :: kpath_coord, dk(3)
        integer :: iout , nfine 


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

            call wannier_matrix_interpolate(fi, H_bloch(:,:,:,jspin:jspin), &
                                            U_full,                          &
                                            fi%kpts, kpts_fine,              &
                                            H_interpol_spin)
            H_interpol(:,:,:,jspin) = H_interpol_spin(:,:,:,1)
            deallocate(H_interpol_spin)
        end do

        deallocate(U_full)

        ! ---------------------------------------------------------------
        ! Step 3: diagonalize interpolated Hamiltonian at each fine k-point
        ! ---------------------------------------------------------------
        lwork = max(1, 2*num_wann - 1)
        allocate(work(lwork))
        allocate(rwork(max(1, 3*num_wann - 2)))
        allocate(hmat_tmp(num_wann, num_wann))
        allocate(eig_interpol(num_wann, nfine, fi%input%jspins))

        do jspin = 1, fi%input%jspins
            do ikpt = 1, nfine 
                hmat_tmp = H_interpol(:,:,ikpt,jspin)
                call cheev('N', 'U', num_wann, hmat_tmp, num_wann, &
                           eig_interpol(:,ikpt,jspin), work, lwork, rwork, info)
                if (info /= 0) &
                    call juDFT_error('cheev failed in interpolate_bandstructure', &
                                     calledby='interpolate_bandstructure')
            end do
        end do

        deallocate(H_interpol, hmat_tmp, work, rwork)

        ! ---------------------------------------------------------------
        ! Write interpolated eigenvalues to wann_bands.dat
        ! Format: kpath_coord  eig_1  eig_2  ...  (one line per k-point)
        ! ---------------------------------------------------------------
        open(newunit=iout, file='wann_bands.dat', status='replace')
        write(iout,'(a)') '# kpath_coord   eig_1  eig_2  ...'

        do jspin = 1, fi%input%jspins
            if (fi%input%jspins > 1) write(iout,'(a,i1)') '# spin ', jspin
            kpath_coord = 0.0
            do ikpt = 1, nfine 
                if (ikpt > 1) then
                    dk = kpts_fine(:,ikpt) - kpts_fine(:,ikpt-1)
                    kpath_coord = kpath_coord + sqrt(dot_product(dk, dk))
                end if
                write(iout,'(f12.6,*(2x,f14.8))') kpath_coord, eig_interpol(:,ikpt,jspin)
            end do
        end do

        close(iout)
        deallocate(eig_interpol)

    end subroutine interpolate_bandstructure

end module m_wannier_interpolate
