module m_dfpt_vefield
    
    implicit none

    contains
    
        subroutine dfpt_vefield(juphon,starsq,atoms,sym,sphhar,cell,dfptvefield,dfptvefieldimag,iDir,q_sign)
            use m_types
            use m_ylm
            use m_sphbes
            USE m_constants
            use m_dfpt_gradient
            USE m_types_cell
            use m_inv3
            use m_phasy1



            type(t_juphon),     intent(in)               :: juphon
            TYPE(t_stars),      intent(in)               :: starsq
            type(t_atoms),      intent(in)               :: atoms
            type(t_sym),        intent(in)               :: sym
            type(t_sphhar),     intent(in)               :: sphhar
            type(t_cell),       intent(in)               :: cell
            type(t_potden),     intent(inout)            :: dfptvefield, dfptvefieldimag
            integer,             intent(in)              :: iDir 
            integer,            intent(in)               :: q_sign

            complex,allocatable                          :: v1efield_mt(:,:,:,:)
            real,allocatable                             :: resultreal(:,:,:,:)
            real,allocatable                             :: resultimag(:,:,:,:)

            complex, allocatable                         :: ylm(:), ylm_m(:)
            complex                                      :: pref
            real                                         :: qlim, qnormvec(3), qvec_ext(3),qvec_int(3),det,inv_bmat(3,3) 
            integer                                      :: n, lmax, l, iop, m, ll1, lm, i, imax
            real, allocatable                            :: sbf(:,:)
            !real                              :: sbf(0:atoms%lmaxd)
            real                                         :: sumreal, sumimag
            TYPE(t_stars)                                  :: starsq_m
            integer                                         :: int
            complex                                      :: phas
            complex                           :: pylm(( atoms%lmaxd + 1 ) ** 2, atoms%ntype),pylm_m(( atoms%lmaxd + 1 ) ** 2, atoms%ntype)


            !print*,'q_sign',q_sign
            qlim = juphon%qlim 

            !interstitial region

            dfptvefield%pw(:,1) = 0.0
            dfptvefield%pw(1,1) = q_sign*cmplx(0.0,1/(2*qlim))
            
            
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
            starsq_m = starsq
            starsq_m%center = -starsq%center
            call phasy1( atoms, starsq, sym, cell, 1, pylm )
            !call phasy1( atoms, starsq_m, sym, cell, 1, pylm_m )            !for minus solution
            do n =1, atoms%ntype
                lmax = atoms%lmax(n)
                imax = atoms%jri(n)
                allocate(sbf(0:lmax,imax))
                do i=1, imax
                    call sphbes(lmax,qlim*atoms%rmsh(i,n),sbf(:,i))
                end do
                do l = 0,lmax
                    pref =  q_sign*fpi_const*(ImagUnit**(l+1))/2
                    ll1 = l*(l+1)+1
                    do m =-l,l
                        lm = ll1 + m 
                        v1efield_mt(1:imax,lm,n,1) = ImagUnit/(qlim*2)*sbf(l,1:imax)*(pylm(lm,n))
                    end do
                end do
                deallocate(sbf)
            end do

            !go to lattice harmonics:
            call sh_to_lh(sym, atoms, sphhar, 1, 3, v1efield_mt(:, :, :, :), resultreal(:,:,:,:), resultimag(:,:,:,:))

            ! return final potential:
            !dfptvefield%pw(:,:) = dfptvefield%pw(:,:)
            dfptvefield%mt(:,:,:,1) = resultreal(:,:,:,1) 
            dfptvefieldimag%mt(:,:,:,1) = resultimag(:,:,:,1) 
            !print*,"max(dfptvefieldimag%mt(:,:,:,1))",dfptvefieldimag%mt(:,:,:,1)
        end subroutine dfpt_vefield

end module m_dfpt_vefield
