module m_dfpt_vefield
    
    implicit none

    contains
    
        subroutine dfpt_vefield(juphon,atoms,sym,sphhar,cell,dfptvefield,dfptvefieldimag,iDir)
            use m_types
            use m_ylm
            use m_sphbes
            USE m_constants
            use m_dfpt_gradient
            USE m_types_cell
            use m_inv3


            type(t_juphon),     intent(in)               :: juphon
            type(t_atoms),      intent(in)               :: atoms
            type(t_sym),        intent(in)               :: sym
            type(t_sphhar),     intent(in)               :: sphhar
            type(t_cell),       intent(in)               :: cell
            type(t_potden),     intent(inout)            :: dfptvefield, dfptvefieldimag
            integer,             intent(in)              :: iDir 

            complex,allocatable                          :: v1efield_mt(:,:,:,:)
            real,allocatable                             :: resultreal(:,:,:,:)
            real,allocatable                             :: resultimag(:,:,:,:)

            complex, allocatable                         :: ylm(:)
            complex                                      :: pref
            real                                         :: qlim, qnormvec(3), qnormvecintern(3) 
            integer                                      :: n, lmax, l, iop, m, ll1, lm, i, imax
            real, allocatable                            :: sbf(:,:)
            real                                         :: sumreal, sumimag



            qlim = juphon%qlim 

            !interstitial region

            dfptvefield%pw(:,1) = 0.0
            dfptvefield%pw(1,1) = cmplx(0.0,1/juphon%qlim)
            
            !MT-region

            allocate(v1efield_mt(atoms%jmtd,sphhar%nlhd +1,atoms%ntype,1))
            v1efield_mt(:,:,:,:) = 0
            allocate(resultreal(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,1))
            allocate(resultimag(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,1))
            
            qnormvec = 0.0
            qnormvec(iDir) = 1

            do n =1, atoms%ntype
                lmax = atoms%lmax(n)
                imax = atoms%jri(n)
                allocate(ylm((lmax+1)**2))
                ylm = 0.0
                call ylm4(lmax,qnormvec(:),ylm(:))
                allocate(sbf(0:lmax,imax))
                do i=1, imax
                    call sphbes(lmax,qlim*atoms%rmsh(i,n),sbf(:,i))
                end do
                do l = 0,lmax
                    pref = fpi_const*(ImagUnit**(l+1))
                    ll1 = l*(l+1)+1
                    do m =-l,l
                        lm = ll1 + m 
                        v1efield_mt(1:imax,lm,n,1) = pref*sbf(l,1:imax)*conjg(ylm(lm))*1/qlim
                    end do
                end do
                deallocate(sbf,ylm)
            end do

            !go to lattice harmonics:
            call sh_to_lh(sym, atoms, sphhar, 1, 3, v1efield_mt(:, :, :, :), resultreal(:,:,:,:), resultimag(:,:,:,:))

            ! return final potential:
            dfptvefield%mt(:,:,:,1) = resultreal(:,:,:,1) 
            dfptvefieldimag%mt(:,:,:,1) = resultimag(:,:,:,1) 

        end subroutine dfpt_vefield

end module m_dfpt_vefield
