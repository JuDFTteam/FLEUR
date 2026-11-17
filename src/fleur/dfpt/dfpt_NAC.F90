!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_dfpt_NAC
    
    USE m_juDFT
    USE m_constants
    USE m_inv3
    USE m_dfpt_dynmat_eig, only: hermitize_scale_dynmat
    USE m_dfpt_data_io, only: read_dfpt_bec, read_dfpt_dielten
    USE m_types_atoms
    USE m_types_fleurinput

    implicit none
    PRIVATE
    PUBLIC :: dfpt_nac, get_nac_ewald
    PUBLIC :: subtract_nac, add_nac

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

    call read_dfpt_bec(fi%atoms,borneffcharge)
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
    call read_dfpt_dielten(dielten)

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

    subroutine get_NAC_ewald(fi,bec,dielten,qvec,dyn_mat_LR,qdir)

        type(t_fleurinput), intent(in)  :: fi
        real,               intent(in)  :: bec(:,:,:)      ! (ntype, displacement, field), cartesian
        real,               intent(in)  :: dielten(3,3)
        real,               intent(in)  :: qvec(3)         ! internal coordinates
        complex,            intent(out) :: dyn_mat_LR(:,:) ! (3nat, 3nat)
        real, optional,     intent(in)  :: qdir(3)         ! internal coordinates

        real, parameter :: expo_max = 14.0   ! terms with K.eps.K/4Lambda above this are dropped

        real    :: eps(3,3), epsinv(3,3), det, lambda, kmax, K(3), KepsK, gq(3), phas
        real    :: blk(3*fi%atoms%nat,3*fi%atoms%nat)
        integer :: nG(3), ix, iy, iz, iAtom, jAtom, ir, jc

        eps = (dielten + transpose(dielten))/2.0
        call inv3(eps,epsinv,det)
        if (eps(1,1) <= 0.0 .or. eps(1,1)*eps(2,2) - eps(1,2)**2 <= 0.0 .or. det <= 0.0) &
            call juDFT_error("Dielectric tensor is not positive definite.", calledby="get_NAC_ewald")

        ! Ewald parameter: the shortest |b_i|^2
        lambda = minval(sum(fi%cell%bmat**2, dim=2))
        ! smallest eigenvalue of eps >= 1/|eps^-1|_F bounds |K| inside the cutoff
        kmax = sqrt(4.0*lambda*expo_max*norm2(epsinv))
        do ix = 1, 3
            nG(ix) = ceiling(kmax*norm2(fi%cell%amat(:,ix))/tpi_const + abs(qvec(ix))) + 1
        end do

        dyn_mat_LR = cmplx(0.0,0.0)
        do iz = -nG(3), nG(3)
            do iy = -nG(2), nG(2)
                do ix = -nG(1), nG(1)
                    gq = qvec + real([ix,iy,iz])
                    K = matmul(gq, fi%cell%bmat)
                    KepsK = dot_product(K, matmul(eps,K))
                    if (KepsK < 1e-10 .or. KepsK/(4.0*lambda) > expo_max) cycle

                    call nac_kblock(fi%atoms, bec, eps, K, blk)
                    blk = blk*exp(-KepsK/(4.0*lambda))
                    do jAtom = 1, fi%atoms%nat
                        jc = 3*(jAtom-1)
                        do iAtom = 1, fi%atoms%nat
                            ir = 3*(iAtom-1)
                            phas = tpi_const*dot_product(gq, fi%atoms%taual(:,jAtom) - fi%atoms%taual(:,iAtom))
                            dyn_mat_LR(ir+1:ir+3,jc+1:jc+3) = dyn_mat_LR(ir+1:ir+3,jc+1:jc+3) &
                                                            + blk(ir+1:ir+3,jc+1:jc+3)*cmplx(cos(phas),sin(phas))
                        end do
                    end do
                end do
            end do
        end do

        if (present(qdir)) then
            K = matmul(qdir, fi%cell%bmat)
            if (norm2(K) > 1e-8) then
                call nac_kblock(fi%atoms, bec, eps, K, blk)
                dyn_mat_LR = dyn_mat_LR + blk
            end if
        end if

        dyn_mat_LR = dyn_mat_LR*fpi_const/fi%cell%omtil

    end subroutine get_NAC_ewald

    subroutine nac_kblock(atoms,bec,eps,K,blk)
        !! (K Z_kappa)_alpha (K Z_kappa')_beta / (K eps K) for one cartesian K, as a 3nat x 3nat block.

        type(t_atoms), intent(in)  :: atoms
        real,          intent(in)  :: bec(:,:,:), eps(3,3), K(3)
        real,          intent(out) :: blk(:,:)

        real    :: kz(3,atoms%nat), v(3*atoms%nat)
        integer :: iAtom

        do iAtom = 1, atoms%nat
            kz(:,iAtom) = matmul(bec(atoms%itype(iAtom),:,:), K)
        end do
        v = reshape(kz, [3*atoms%nat])
        blk = spread(v, 2, size(v))*spread(v, 1, size(v))/dot_product(K, matmul(eps,K))

    end subroutine nac_kblock

    subroutine subtract_NAC(fi,qbk,dyn_mat)
        !! Removes the dipole-dipole part from the stored coarse-mesh dynamical matrices,
        !! so that only the analytic part is interpolated (Gonze-Lee).

        type(t_fleurinput), intent(in)    :: fi
        real,               intent(in)    :: qbk(:,:)       ! (3, nq), internal coordinates
        complex,            intent(inout) :: dyn_mat(:,:,:) ! (3nat, 3nat, nq), mass-scaled

        real, allocatable :: bec(:,:,:)
        real              :: dielten(3,3)
        complex           :: dyn_mat_LR(3*fi%atoms%nat,3*fi%atoms%nat)
        integer           :: iQ

        allocate(bec(fi%atoms%ntype,3,3))
        call read_dfpt_bec(fi%atoms, bec)
        call read_dfpt_dielten(dielten)
        do iQ = 1, size(qbk,2)
            call get_NAC_ewald(fi, bec, dielten, qbk(:,iQ), dyn_mat_LR)
            call hermitize_scale_dynmat(fi%atoms, dyn_mat_LR)
            dyn_mat(:,:,iQ) = dyn_mat(:,:,iQ) - dyn_mat_LR
        end do

    end subroutine subtract_NAC

    subroutine add_NAC(fi,qbk,dyn_mat,l_path)
        !! Adds the dipole-dipole part to interpolated dynamical matrices. With l_path, a Gamma
        !! point gets the LO-TO splitting along the direction of the neighbouring path point.

        type(t_fleurinput), intent(in)    :: fi
        real,               intent(in)    :: qbk(:,:)       ! (3, nq), internal coordinates
        complex,            intent(inout) :: dyn_mat(:,:,:) ! (3nat, 3nat, nq), not mass-scaled
        logical,            intent(in)    :: l_path

        real, allocatable :: bec(:,:,:)
        real              :: dielten(3,3), qdir(3)
        complex           :: dyn_mat_LR(3*fi%atoms%nat,3*fi%atoms%nat)
        integer           :: iQ

        allocate(bec(fi%atoms%ntype,3,3))
        call read_dfpt_bec(fi%atoms, bec)
        call read_dfpt_dielten(dielten)
        do iQ = 1, size(qbk,2)
            if (l_path .and. gamma_direction(qbk, iQ, qdir)) then
                call get_NAC_ewald(fi, bec, dielten, qbk(:,iQ), dyn_mat_LR, qdir)
            else
                call get_NAC_ewald(fi, bec, dielten, qbk(:,iQ), dyn_mat_LR)
            end if
            dyn_mat(:,:,iQ) = dyn_mat(:,:,iQ) + dyn_mat_LR
        end do

    end subroutine add_NAC

    logical function gamma_direction(bk, iQ, qdir)
        !! True if path point iQ is a Gamma point; qdir is then the direction
        !! to the next path point, or to the previous one at the end of the path.

        real,    intent(in)  :: bk(:,:)
        integer, intent(in)  :: iQ
        real,    intent(out) :: qdir(3)

        qdir = 0.0
        gamma_direction = .false.
        if (norm2(bk(:,iQ) - nint(bk(:,iQ))) > 1e-8) return
        if (iQ < size(bk,2)) qdir = bk(:,iQ+1) - bk(:,iQ)
        if (norm2(qdir) < 1e-8 .and. iQ > 1) qdir = bk(:,iQ-1) - bk(:,iQ)
        gamma_direction = norm2(qdir) > 1e-8

    end function gamma_direction

end module
