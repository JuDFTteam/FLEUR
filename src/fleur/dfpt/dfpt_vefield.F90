!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_dfpt_vefield
    
    use m_ylm
    use m_sphbes
    use m_constants
    use m_dfpt_gradient
    use m_types_cell
    use m_inv3
    use m_phasy1
    use m_npy
    use m_types_atoms
    use m_types_dfpt
    use m_types_potden
    use m_types_sphhar
    use m_types_stars
    use m_types_sym
    implicit none
    private
    public :: dfpt_vefield

contains
    
    subroutine dfpt_vefield(dfpt,starsq,atoms,sym,sphhar,cell,dfptvefield,iDir)
        !Currently only spin=1





        type(t_dfpt),     intent(in)               :: dfpt
        TYPE(t_stars),      intent(in)               :: starsq
        type(t_atoms),      intent(in)               :: atoms
        type(t_sym),        intent(in)               :: sym
        type(t_sphhar),     intent(in)               :: sphhar
        type(t_cell),       intent(in)               :: cell
        type(t_potden),     intent(inout)            :: dfptvefield
        integer,             intent(in)              :: iDir

        complex,allocatable                          :: v1efield_mt(:,:,:,:)
        real,allocatable                             :: resultreal(:,:,:,:)
        real,allocatable                             :: resultimag(:,:,:,:)

        complex, allocatable                         :: ylm(:), ylm_m(:)
        complex                                      :: pref
        real                                         :: qlim, qcart(3), qnormvec(3), qvec_ext(3),qvec_int(3),det,inv_bmat(3,3)
        integer                                      :: n, lmax, l, iop, m, ll1, lm, i, imax
        real, allocatable                            :: sbf(:,:)
        !real                              :: sbf(0:atoms%lmaxd)
        real                                         :: sumreal, sumimag
        integer                                         :: int
        complex                                      :: phas
        complex                           :: pylm(( atoms%lmaxd + 1 ) ** 2, atoms%ntype)
        character(len=20)                   :: densave_string


        qlim = dfpt%qlim
        qcart = matmul(starsq%center, cell%bmat)
        pref = ImagUnit/(2*qcart(iDir))

        !interstitial region
        dfptvefield%pw(:,1) = 0.0
        dfptvefield%pw(1,1) = pref
        
        
        !MT-region

        allocate(v1efield_mt(atoms%jmtd,sphhar%nlhd +1,atoms%ntype,1))
        v1efield_mt(:,:,:,:) = 0
        allocate(resultreal(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,1))
        allocate(resultimag(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,1))
        
        qnormvec = 0.0
        qnormvec(iDir) = 1
        qvec_ext = qlim*qnormvec
        call inv3(cell%bmat,inv_bmat(:,:),det)
        qvec_int = matmul(qvec_ext,transpose(inv_bmat))
        !print*,"starsq%center",starsq%center
        call phasy1( atoms, starsq, sym, cell, 1, pylm )
        !call save_npy('pylm.npy', pylm)
        do n =1, atoms%ntype
            lmax = atoms%lmax(n)
            imax = atoms%jri(n)
            allocate(sbf(0:lmax,imax))
            do i=1, imax
                call sphbes(lmax,qlim*atoms%rmsh(i,n),sbf(:,i))
            end do
            do l = 0,lmax
                ll1 = l*(l+1)+1
                do m =-l,l
                    lm = ll1 + m
                    v1efield_mt(1:imax,lm,n,1) = pref*sbf(l,1:imax)*(pylm(lm,n))
                end do
            end do
            deallocate(sbf)
        end do

        !go to lattice harmonics:
        call sh_to_lh(sym, atoms, sphhar, 1, 3, v1efield_mt(:, :, :, :), resultreal(:,:,:,:), resultimag(:,:,:,:))

        ! return final potential:
        !dfptvefield%pw(:,:) = dfptvefield%pw(:,:)
        dfptvefield%mt(:,:,:,1) = resultreal(:,:,:,1) 
        dfptvefield%mtIm(:,:,:,1) = resultimag(:,:,:,1) 

        !do n = 1,atoms%ntype

        !write(densave_string,"(a,i0,a)")"vext1_it1_pw_",iDir,".npy"
        !call save_npy(densave_string, dfptvefield%pw(:,1))
        !write(densave_string,"(a,i0,a)")"vext1_it1_mt_",iDir,".npy"
        !call save_npy(densave_string, dfptvefield%mt(:,:,:,1))
        !write(densave_string,"(a,i0,a)")"vext1_it1_mtIm_",iDir,".npy"
        !call save_npy(densave_string, dfptvefield%mtIm(:,:,:,1))
        !stop
        !print*,"max(dfptvefield%mtIm(:,:,:,1))",dfptvefield%mtIm(:,:,:,1)
    end subroutine dfpt_vefield



end module m_dfpt_vefield
