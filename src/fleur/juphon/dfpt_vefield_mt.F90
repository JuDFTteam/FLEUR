    module m_dfpt_vefield_mt
    
    implicit none

    contains
    
        subroutine dfpt_vefield_mt(juphon,atoms,sym,sphhar,dfptvCoulreal,dfptvCoulimag)
            use m_types
            use m_ylm
            use m_sphbes
            USE m_constants
            use m_dfpt_gradient
            !use m_SphBessel

            type(t_juphon),     intent(in)               :: juphon
            type(t_atoms),      intent(in)               :: atoms
            type(t_sym),        intent(in)               :: sym
            TYPE(t_sphhar),    INTENT(IN)                :: sphhar
            real,               intent(out)              :: dfptvCoulreal(:,0:,:)!(atoms%jmtd,0:sphharl%nlhd,atoms%ntype)   
            real,               intent(out)              :: dfptvCoulimag(:,0:,:)

            complex,allocatable                          :: vcoul(:,:,:,:)
            real,allocatable                             :: resultreal(:,:,:,:)
            real,allocatable                             :: resultimag(:,:,:,:)

            !type(t_potden),     intent(inout)            :: vCoul
            complex, allocatable                         :: ylm(:)
            complex                                      :: pref
            real                                         :: qlim, qnormvec(3)
            integer                                      :: n, lmax, l, iop, m, ll1, lm, i, imax
            !real, allocatable, dimension(:,:)            :: il, kl     
            real, allocatable                            :: sbf(:,:)
            real                                        :: sumreal, sumimag

            
            allocate(vcoul(atoms%jmtd,sphhar%nlhd +1,atoms%ntype,1))
            allocate(resultreal(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,1))
            allocate(resultimag(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,1))
            !print*,sphhar%nlhd
            !print*, "vcoul shape", shape(vcoul)
            !stop
            qlim = juphon%qlim 
            print*, 'qlim', qlim
            qnormvec = 0.0
            qnormvec(1) = 1
            !print*,"dimmension"
            !print*,shape(vr(:,0:,:))
            !print*, "sym"
            !print*, sym%nop
            !print*, qnormvec
            !Calculate spherical harmonics:
            !do n = 1, atoms%ntype
                !lmax = atoms%lmax(n)
                !CALL ylm4(atoms%lmaxd, rg, ylm(:,iOp))!kee
                !print*, "test"
                !print*, (lmax+1)**2
                !print*,ylm(:)
                !allocate(ylm((lmax+1)**2))
                !print*,ylm(:)
                !allocate(ylm(((lmax+1)**2,sym%nop)))
                !call ylm4(lmax,qnormvec,ylm(:))
                !print*,shape(ylm(:))
                !print*, ylm(:)
                !print*, lmax
                !do l = 0,lmax
                !    ll1 = l*(l+1)+1
                !    do m =-l,l
                !        lm = ll1 + m 
                        !print*, lm

                    !print *, l
                    !vr(:,0:,n) =            i/qlim
                    !end do
                !end do
               ! deallocate(ylm)
            !end do
            !print*, "between routines"
            !print*, ylm(:)
            !print*,"now bessel functions"
            !print*, shape(il(:,:))
            !do n = 1, atoms%ntype
            !    call sphbes( atoms%lmax(n), qlim, sbf ) !adjust to different atom types, argument of bessel function not finished adjust to radial grid
            !    print*, sbf
                !print*, shape(sbf)
            !end do
            !determine final coefficient:
            !print*, "complete version"
            !allocate(vcoul(atoms%jmtd,1,atoms%ntype))
            !print*, shape(vr)
            !print*, shape(vcoul)
            !print*, "test"
            do n =1, atoms%ntype
                lmax = atoms%lmax(n)
                imax = atoms%jri(n)
                !print*,ylm(:)
                allocate(ylm((lmax+1)**2))
                ylm = 0.0
                !print*,ylm(:)
                call ylm4(lmax,qnormvec,ylm(:))
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
                vcoul(:,:,:,:) = 0
                do l = 0,lmax
                    pref = fpi_const*(ImagUnit**(l+1))
                    !print*,pref
                    ll1 = l*(l+1)+1
                    !print*, "l ",l
                    do m =-l,l
                        lm = ll1 + m 
                        !print*, lm
                        !print*,"shapes"
                        !print*,sbf(l+1,1:imax)
                        !print*, "sbf", shape(sbf(:,1:imax))
                        !print*, "vr", shape(vr(1:imax,:,n))
                        !print*, "l+1",l+1
                        !print*, "lm",lm
                        !print*, "sbf",sbf(l,1:imax)
                        vcoul(1:imax,lm,n,1) = pref*sbf(l,1:imax)*conjg(ylm(lm))*1/qlim!pref*sbf(l,1:imax)*conjg(ylm(lm))*1/qlim  
                        !print*, "vr",vr(:,lm-1,1)
                    end do
                end do
                deallocate(sbf,ylm)
            end do
            !print*, vr(:,:,:)
            !write(1000,*) vr(:,:,:)
            call save_npy("vcoul.npy",vcoul(:,:,1,1))
            !print*,vcoul(1,:9,1,1)
            !print*, lbound(vr,1)
            !print*, lbound(vr,2)
            !print*, lbound(vr,3)
            !print*, vr(1,0,1)
            !print*, vr(:,63,1)
            !print*, vr(:,64,1)
            !print*, shape(vr(:,:,1))
            !print*, vcoul(:,2,:,1)
            


            !go to lattice harmonics: save in seperate array?
            
            call sh_to_lh(sym, atoms, sphhar, 1, 3, vcoul(:, :, :, :), resultreal(:,:,:,:), resultimag(:,:,:,:))
            !print*, "real part:"
            !print*, abs(resultreal(2,:,:,1))
            !print*, "imag part:"
            !print*, abs(resultimag(2,:,:,1))
            !print*, 'shape', shape(resultreal(:,:,:,:))
            !print*, 'vcoul(1,:9,1,1)'
            !print*,vcoul(1,:9,1,1)
            !print*, 'resultreal'
            !print*, resultreal(1,0:8,1,1)
            !print*, 'resultimag'
            !print*, resultimag(1,0:8,1,1)
            !call save_npy("resultreal.npy",resultreal(:,:,1,1))
            !call save_npy("resultimag.npy",resultimag(:,:,1,1))  


            !sumreal = sum(abs(resultreal(:,:,:,1)))
            !sumimag = sum(abs(resultimag(:,:,:,1)))
            !print*,'sumreal', sumreal
            !print*,'sumimag', sumimag

            !allocate(ylm((lmax+1)**2, lapw%nv(igSpin)), stat=ierr)
            
            !q_norm = ()
            !call ylm4_batched(0,0,)
            !print*, 


            ! return final potential:
            dfptvCoulreal(:,:,:) = resultreal(:,:,:,1) 
            dfptvCoulimag(:,:,:) = resultimag(:,:,:,1) 
            
            !print*, "test"
            !print*, "shape", shape(sphhar%clnu(:,:,:))
            !print*, "phase stuff", sphhar%clnu(:,:,:)
            !print*, "Stop"
            !stop
            !print*, juphon%qlim
            !print*,shape(vr)! vCoul(:,0:,:)
            !lmax = 9
            !vr(:,0:,n) =            i/qlim

            
            !vr(:,0:,n) =            i/qlim
        end subroutine dfpt_vefield_mt

end module m_dfpt_vefield_mt
