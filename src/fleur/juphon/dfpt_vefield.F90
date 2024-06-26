module m_dfpt_vefield
    
    implicit none

    contains
    
        subroutine dfpt_vefield(juphon,atoms,sym,sphhar,cell,dfptvefield,dfptvefieldimag)
            use m_types
            use m_ylm
            use m_sphbes
            USE m_constants
            use m_dfpt_gradient
            USE m_types_cell
            use m_inv3
            !use m_SphBessel


            type(t_juphon),     intent(in)               :: juphon
            type(t_atoms),      intent(in)               :: atoms
            type(t_sym),        intent(in)               :: sym
            TYPE(t_sphhar),    INTENT(IN)                :: sphhar
            TYPE(t_cell),      INTENT(IN)                :: cell
            type(t_potden),     intent(inout)            :: dfptvefield, dfptvefieldimag
            !real,               intent(out)              :: dfptvCoulreal(:,0:,:)!(atoms%jmtd,0:sphharl%nlhd,atoms%ntype)   
            !real,               intent(out)              :: dfptvCoulimag(:,0:,:)

            complex,allocatable                          :: v1efield_mt(:,:,:,:)
            real,allocatable                             :: resultreal(:,:,:,:)
            real,allocatable                             :: resultimag(:,:,:,:)

            !type(t_potden),     intent(inout)            :: v1efield_mt
            complex, allocatable                         :: ylm(:)
            complex                                      :: pref
            real                                         :: qlim, qnormvec(3), qnormvecintern(3)
            real                                         :: inv_bmat(3,3), determinant
            integer                                      :: n, lmax, l, iop, m, ll1, lm, i, imax
            !real, allocatable, dimension(:,:)            :: il, kl     
            real, allocatable                            :: sbf(:,:)
            real                                        :: sumreal, sumimag

            !print*, "dfptvefield%mt", shape(dfptvefield%mt)
            !print*, "dfptvefield%pw", shape(dfptvefield%pw)
            !print*, "dfptvefieldimag%mt", shape(dfptvefieldimag%mt)
            !print*, "dfptvefieldimag%pw", shape(dfptvefieldimag%pw)


            !interstitial region
            dfptvefield%pw(:,1) = 0.0
            !print*,'dfptvefieldimag%pw(1,ispin)',dfptvefieldimag%pw(:,1)
            dfptvefield%pw(1,1) = cmplx(0.0,1/juphon%qlim)
            call save_npy("vefield-pw.npy",dfptvefield%pw(:,1))
            !print*,'dfptvefieldimag%pw(1,ispin)',dfptvefieldimag%pw(:,1)
            !print*,'dfptvefieldimag%pw(1,ispin)',dfptvefieldimag%mt(:,1,1,1)
            !print*, 'stop'
            !stop

            !MT-part
            allocate(v1efield_mt(atoms%jmtd,sphhar%nlhd +1,atoms%ntype,1))
            allocate(resultreal(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,1))
            allocate(resultimag(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,1))
            !print*,sphhar%nlhd
            !print*, "v1efield_mt shape", shape(v1efield_mt)
            !stop
            qlim = juphon%qlim 
            print*, 'qlim', qlim
            qnormvec = 0.0
            qnormvec(3) = 1
            !print*, "cell%bmat",cell%bmat
            !print*, "cell%bbmat",cell%bbmat
            !inv_bmat(:,:) = matmul()
            inv_bmat(:,:) =0.0
            call inv3(cell%bmat,inv_bmat(:,:),determinant)
            qnormvecintern(:) = matmul(qnormvec,cell%bmat(:,:))!matmul(cell%amat,qnormvec)!matmul(qnormvec(:),inv_bmat(:,:))
            print*, 'qnormvecintern(:)',qnormvecintern(:)
            !stop

            print*, 'qnormvec',qnormvec(:)
            !print*, 'qnormvectrans', qnorm
            do n =1, atoms%ntype
                lmax = atoms%lmax(n)
                imax = atoms%jri(n)
                !print*,ylm(:)
                allocate(ylm((lmax+1)**2))
                ylm = 0.0
                !print*,ylm(:)
                call ylm4(lmax,qnormvecintern(:),ylm(:))
                !print*, "lower bound"
                !print*,lbound(ylm,1)
                !print*,ylm(:)
                call save_npy("ylm.npy",ylm(:))
                allocate(sbf(0:lmax,imax))
                !print*, "lower bound"
                !print*,lbound(sbf,1)
                do i=1, imax
                    !print*, "radial mesh"
                    !print*, atoms%rmsh(i,n)*qlim
                    call sphbes(lmax,qlim*atoms%rmsh(i,n),sbf(:,i))
                    !print*, sbf(:,i)
                    !print*, atoms%rmsh(i,n)
                end do
                call save_npy("sphbes.npy",sbf(:,:))
                v1efield_mt(:,:,:,:) = 0
                do l = 0,lmax
                    pref = fpi_const*(ImagUnit**(l+1))
                    !print*,pref
                    ll1 = l*(l+1)+1
                    !print*, "l ",l
                    do m =-l,l
                        lm = ll1 + m 
                        v1efield_mt(1:imax,lm,n,1) = pref*sbf(l,1:imax)*conjg(ylm(lm))*1/qlim!pref*sbf(l,1:imax)*conjg(ylm(lm))*1/qlim  
                    end do
                end do
                deallocate(sbf,ylm)
            end do
            call save_npy("vefield-mt.npy",v1efield_mt(:,:,1,1))
            !print*, vr(:,:,:)
            !write(1000,*) vr(:,:,:)
            !call save_npy("v1efield_mt.npy",v1efield_mt(:,:,1,1))

            !go to lattice harmonics: save in seperate array?
            
            call sh_to_lh(sym, atoms, sphhar, 1, 3, v1efield_mt(:, :, :, :), resultreal(:,:,:,:), resultimag(:,:,:,:))

            ! return final potential:
            dfptvefield%mt(:,:,:,1) = resultreal(:,:,:,1) 
            dfptvefieldimag%mt(:,:,:,1) = resultimag(:,:,:,1) 

        end subroutine dfpt_vefield

end module m_dfpt_vefield
