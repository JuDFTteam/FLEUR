!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_dfpt_NAC
    
    USE m_juDFT
    USE m_constants
    USE m_types

    implicit none

contains

    subroutine dfpt_NAC(fi,dyn_mat_NAC,qnorm_in)

    type(t_fleurinput), intent(in)     :: fi
    complex, intent(inout)             :: dyn_mat_NAC(:,:)
    real, intent(in),optional          :: qnorm_in(3)          
    real                               :: qnorm_vec(3),dielten(3,3),tempval_denom,tempval_BEC_fixed,tempval_NAC,qnorm_vec_ext(3)
    real, allocatable                  :: borneffcharge(:,:,:)
    integer                            ::  iDir,iDir1,iType,iType_row,iDir_row
    logical                            :: l_present

    allocate(borneffcharge(fi%atoms%ntype,3,3))

    call read_BEC(fi,borneffcharge)
    print*,"test"
    if (present(qnorm_in)) then
        qnorm_vec = qnorm_in
    else
        qnorm_vec(:) = [1.0,1.0,0.0]/norm2([1.0,1.0,0.0])![0.5,0.0,0.5]/norm2([0.5,0.0,0.5])![0.577350269189626,0.577350269189626,0.577350269189626]
    end if
    print*,"qnorm_vec",qnorm_vec
    tempval_denom = 0.0
    qnorm_vec_ext =0.0
    if (.true.) then
        !print*,
        qnorm_vec_ext = matmul(qnorm_vec,fi%cell%bmat)
        qnorm_vec = qnorm_vec_ext
    end if
    print*,"qnorm_vec_ext",qnorm_vec_ext
    print*,"qnorm_vec",qnorm_vec
    call read_dielten(fi,dielten)

    do iDir=1,3
        do iDir1=1,3
            tempval_denom = tempval_denom +qnorm_vec(iDir)*dielten(iDir,iDir1)*qnorm_vec(iDir1)
        end do
    end do

    do itype=1,fi%atoms%ntype
        do iDir=1,3
            do iType_row=1,fi%atoms%ntype
                do iDir_row=1,3
                    tempval_NAC=0.0
                    tempval_NAC = tempval_NAC + sum(borneffcharge(itype,iDir,:)*qnorm_vec(:))*sum(borneffcharge(iType_row,iDir_row,:)*qnorm_vec(:))/tempval_denom
                    dyn_mat_NAC(3 *(itype-1)+iDir,3 *(iType_row-1)+iDir_row) = (fpi_const/fi%cell%omtil)*tempval_NAC*cmplx(1.,0.)
                end do
            end do
        end do
    end do

    end subroutine dfpt_NAC

    subroutine read_BEC(fi,borneffcharge)

        type(t_fleurinput), intent(in)     :: fi
        real, allocatable, intent(inout)   :: borneffcharge(:,:,:)
        integer                            :: iDType,iDir
        character(len=200)                 :: dummy_line

        borneffcharge(:,:,:)=0.0
        open(111, file="born_eff_charge", status='old', action='read', form='formatted')

        do iDType = 1, fi%atoms%ntype
            read(111,'(A)') dummy_line

            do iDir = 1, 3
                read(111,*) borneffcharge(iDType, iDir, 1:3)
            end do
        end do

        close(111)
    
    end subroutine read_BEC

    subroutine read_dielten(fi,dielten)

        type(t_fleurinput), intent(in)     :: fi
        real, intent(inout)   :: dielten(:,:)
        integer                            :: iDir

        dielten(:,:)=0.0
        open(110, file="diel_tensor", status='old', action='read', form='formatted')

        do iDir = 1, 3
            read(110,*) dielten( iDir, 1:3)
            !print*,dielten( iDir, 1:3)
        end do

        close(110)
    
    end subroutine read_dielten

    subroutine get_NAC(fi,qpts,dyn_mat_q,iQ)

        type(t_fleurinput), intent(in)     :: fi
        type(t_kpts), intent(in)           :: qpts
        complex, intent(inout)             :: dyn_mat_q(:,:)
        integer, intent(in)                :: iQ
        
        complex, allocatable               :: dyn_mat_NAC(:,:) 
        real                               :: qnorm_vec(3)                    

        allocate(dyn_mat_NAC(3*fi%atoms%ntype,3*fi%atoms%ntype))

        dyn_mat_NAC(:,:) = cmplx(0.0,0.0)
        print*,"sum(dyn_mat_q) 1",sum(dyn_mat_q(:,:))

        print*,"in the routine"
        print*,"iQ",iQ
        if (iQ==1) then
            print*,"qpts%bk(:,iQ+1)",qpts%bk(:,iQ+1)
            qnorm_vec(:)=qpts%bk(:,iQ+1)/norm2(qpts%bk(:,iQ+1))
            print*,"qnorm_vec(:)",qnorm_vec(:)
        else
            print*,"qpts%bk(:,iQ+1)",qpts%bk(:,iQ-1)
            qnorm_vec(:)=qpts%bk(:,iQ-1)/norm2(qpts%bk(:,iQ-1))
            print*,"qnorm_vec(:)",qnorm_vec(:)
        end if
        call dfpt_NAC(fi,dyn_mat_NAC,qnorm_vec)
        !stop
        !dyn_mat_q(:,:) = cmplx(0.0,0.0)
        dyn_mat_q = dyn_mat_q +dyn_mat_NAC
        print*,"dyn_mat_NAC",dyn_mat_NAC
        print*,"sum(dyn_mat_q)",sum(dyn_mat_q(:,:))

    end subroutine get_NAC

    subroutine get_NAC_ewald(fi,qpts,stars,dyn_mat_NAC,qvec,iQ)

        type(t_fleurinput), intent(in)     :: fi
        type(t_kpts), intent(in)           :: qpts
        type(t_stars) , intent(in)         :: stars
        complex, intent(inout)             :: dyn_mat_NAC(:,:)
        real, intent(in)                   :: qvec(3)
        integer, intent(in)                :: iQ

        complex, allocatable               :: dyn_mat_NAC_gq(:,:)
        complex                            :: grid_phase_sum, tau_phase
        real                               :: gqvec(3), gauss_fact, phas, phas_tau, tau_diff(3)
        integer                            :: istar, ft_lim(2,3), iz, iy, ix
        integer                            :: itype_col, iType_row, iDir_col, iDir_row
        integer                            :: idx_col, idx_row, star_max



        !print*,"fi%sym%nop",fi%sym%nop
        !stop
        star_max = min(stars%ng3,20)
        allocate(dyn_mat_NAC_gq(3*fi%atoms%ntype,3*fi%atoms%ntype))
        dyn_mat_NAC = (0.0,0.0)

        ft_lim(2,:) = qpts%nkpt3(:)/2
        ft_lim(1,:) = ft_lim(2,:) - qpts%nkpt3(:) + 1

        do istar = 1, star_max
            !print*,"istar",istar
            gqvec = stars%kv3(:,istar) + qvec
            if (norm2(gqvec) .lt. 1e-8) cycle

            call dfpt_NAC(fi,dyn_mat_NAC_gq,gqvec)
            call get_gaussian(fi,gqvec,gauss_fact)
            !print*,"gauss_fact", gauss_fact

            grid_phase_sum = (0.0,0.0)
            do iz = ft_lim(1,3), ft_lim(2,3)
                do iy = ft_lim(1,2), ft_lim(2,2)
                    do ix = ft_lim(1,1), ft_lim(2,1)
                        phas = tpi_const * (gqvec(1)*real(ix) + gqvec(2)*real(iy) + gqvec(3)*real(iz))
                        grid_phase_sum = grid_phase_sum + cmplx(cos(phas), sin(phas))
                    end do
                end do
            end do

            do itype_col = 1, fi%atoms%ntype
                do iType_row = 1, fi%atoms%ntype
                    tau_diff = fi%atoms%taual(:,itype_col) - fi%atoms%taual(:,iType_row)
                    phas_tau = tpi_const * dot_product(gqvec, tau_diff)
                    tau_phase = cmplx(cos(phas_tau), sin(phas_tau))

                    do iDir_col = 1, 3
                        idx_col = 3 * (itype_col - 1) + iDir_col
                        do iDir_row = 1, 3
                            idx_row = 3 * (iType_row - 1) + iDir_row
                            dyn_mat_NAC(idx_col, idx_row) = dyn_mat_NAC(idx_col, idx_row) + &
                                dyn_mat_NAC_gq(idx_col, idx_row) * gauss_fact * grid_phase_sum !* tau_phase
                        end do
                    end do
                end do
            end do
        end do
        !stop

        deallocate(dyn_mat_NAC_gq)

    end subroutine get_NAC_ewald

    subroutine get_NAC_ewald_r(fi,qpts,stars,dyn_mat_NAC_r,qvec,iQ)

        type(t_fleurinput), intent(in)     :: fi
        type(t_kpts), intent(in)           :: qpts
        type(t_stars), intent(in)          :: stars
        complex, intent(inout)             :: dyn_mat_NAC_r(:,:,:)
        real, intent(in)                   :: qvec(3)
        integer, intent(in)                :: iQ

        complex, allocatable               :: dyn_mat_NAC_gq(:,:)
        complex                            :: cell_phase, tau_phase
        real                               :: gqvec(3), gauss_fact, phas, phas_tau, tau_diff(3)
        integer                            :: istar, ft_lim(2,3), iz, iy, ix, iGrid
        integer                            :: itype_col, iType_row, iDir_col, iDir_row
        integer                            :: idx_col, idx_row, star_max


        !print*,"fi%sym%nop",fi%sym%nop
        !stop
        star_max = min(stars%ng3,20)
        allocate(dyn_mat_NAC_gq(3*fi%atoms%ntype,3*fi%atoms%ntype))
        dyn_mat_NAC_r = (0.0,0.0)

        ft_lim(2,:) = qpts%nkpt3(:)/2
        ft_lim(1,:) = ft_lim(2,:) - qpts%nkpt3(:) + 1

        iGrid = 0
        do iz = ft_lim(1,3), ft_lim(2,3)
            do iy = ft_lim(1,2), ft_lim(2,2)
                do ix = ft_lim(1,1), ft_lim(2,1)
                    iGrid = iGrid + 1
                    do istar = 1, star_max
                        gqvec = stars%kv3(:,istar) + qvec
                        if (norm2(gqvec) .lt. 1e-8) cycle

                        call dfpt_NAC(fi,dyn_mat_NAC_gq,gqvec)
                        call get_gaussian(fi,gqvec,gauss_fact)
                        !print*,"gauss_fact", gauss_fact

                        phas = tpi_const * (gqvec(1)*real(ix) + gqvec(2)*real(iy) + gqvec(3)*real(iz))
                        cell_phase = cmplx(cos(phas), sin(phas))

                        do itype_col = 1, fi%atoms%ntype
                            do iType_row = 1, fi%atoms%ntype
                                tau_diff = fi%atoms%taual(:,itype_col) - fi%atoms%taual(:,iType_row)
                                phas_tau = tpi_const * dot_product(gqvec, tau_diff)
                                tau_phase = cmplx(cos(phas_tau), sin(phas_tau))

                                do iDir_col = 1, 3
                                    idx_col = 3 * (itype_col - 1) + iDir_col
                                    do iDir_row = 1, 3
                                        idx_row = 3 * (iType_row - 1) + iDir_row
                                        dyn_mat_NAC_r(iGrid, idx_col, idx_row) = dyn_mat_NAC_r(iGrid, idx_col, idx_row) + &
                                            dyn_mat_NAC_gq(idx_col, idx_row) * gauss_fact * cell_phase !* tau_phase
                                    end do
                                end do
                            end do
                        end do
                    end do
                    !print*,"sum(dyn_mat_NAC_r(iGrid,:,:))",sum(dyn_mat_NAC_r(iGrid,:,:))
                end do
            end do
        end do

        deallocate(dyn_mat_NAC_gq)

    end subroutine get_NAC_ewald_r

    subroutine get_gaussian(fi,gqvec,gauss_fact)

        type(t_fleurinput), intent(in)     :: fi
        real, intent(in)                    :: gqvec(3)
        real, intent(inout)                 :: gauss_fact

        integer                         :: iDir, iDir1
        real                            :: tempval_exp, para_ewald, dielten(3,3)

        dielten(:,:) = 0.0
        tempval_exp = 0.0
        para_ewald = 0.00000005
        open(110, file="diel_tensor", status='old', action='read', form='formatted')

        do iDir = 1, 3
            read(110,*) dielten(iDir, 1:3)
        end do
        close(110)

        do iDir = 1, 3
            do iDir1 = 1, 3
                tempval_exp = tempval_exp + gqvec(iDir) * dielten(iDir, iDir1) * gqvec(iDir1)
            end do
        end do

        tempval_exp = -tempval_exp/(4*para_ewald)
        gauss_fact = exp(tempval_exp)

    end subroutine get_gaussian

end module
