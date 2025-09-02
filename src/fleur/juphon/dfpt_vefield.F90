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
            use m_npy




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
            character(len=20)                   :: densave_string


            !print*,'q_sign',q_sign
            qlim = juphon%qlim 

            !interstitial region
            dfptvefield%pw(:,1) = 0.0
            dfptvefield%pw(1,1) = cmplx(0.0,1/(2*qlim))!
            
            
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
            !print*,"starsq%center",starsq%center
            call phasy1( atoms, starsq, sym, cell, 1, pylm )
            call phasy1( atoms, starsq_m, sym, cell, 1, pylm_m )            !for minus solution
            !call save_npy('pylm.npy', pylm)
            do n =1, atoms%ntype
                lmax = atoms%lmax(n)
                imax = atoms%jri(n)
                allocate(sbf(0:lmax,imax))
                do i=1, imax
                    call sphbes(lmax,qlim*atoms%rmsh(i,n),sbf(:,i))
                end do
                do l = 0,lmax
                    !pref =  -q_sign*fpi_const*(ImagUnit**(l+1))/2!q_sign*
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

            !do n = 1,atoms%ntype

            !write(densave_string,"(a,i0,a)")"vext1_it1_pw_",iDir,".npy"
            !call save_npy(densave_string, dfptvefield%pw(:,1))
            !write(densave_string,"(a,i0,a)")"vext1_it1_mt_",iDir,".npy"
            !call save_npy(densave_string, dfptvefield%mt(:,:,:,1))
            !write(densave_string,"(a,i0,a)")"vext1_it1_mtIm_",iDir,".npy"
            !call save_npy(densave_string, dfptvefieldimag%mt(:,:,:,1))
            !stop
            !print*,"max(dfptvefieldimag%mt(:,:,:,1))",dfptvefieldimag%mt(:,:,:,1)
        end subroutine dfpt_vefield


        subroutine dfpt_vefield_int(fi,stars,sphhar,fmpi,rho,q_sign)
            
            use m_types
            use m_dfpt_potden_offset
            use m_make_stars
            type(t_fleurinput), intent(in)     :: fi
            TYPE(t_stars),      INTENT(IN)     :: stars
            type(t_sphhar),    intent(in)      :: sphhar
            type(t_potden), intent(inout)         :: rho

            type(t_mpi),        intent(in)     :: fmpi
            integer,      intent(in)           :: q_sign
            type(t_potden)                     :: vExt1, vExt1Im
            type(t_stars)                      :: starsq_vext
            integer                            :: iDir_col
            real                               :: qvec_int(3)
            complex                            ::offset_out

            iDir_col =  1
            qvec_int= q_sign*fi%juPhon%qvec_efield(iDir_col,:)

            CALL starsq_vext%reset_stars()

            CALL make_stars(starsq_vext, fi%sym, fi%atoms, fi%vacuum, sphhar, fi%input, fi%cell, fi%noco, fmpi, qvec_int, 1, iDir_col,fi%juPhon%l_efield)
            starsq_vext%ufft = stars%ufft

            call vExt1%init(starsq_vext, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
            call vExt1Im%init(starsq_vext, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)
            call dfpt_vefield(fi%juPhon,starsq_vext,fi%atoms,fi%sym,sphhar,fi%cell,vExt1,vExt1Im,iDir_col,q_sign)
            

            !print*,"atoms",fi%atoms%ntype
            !print*,"shape(vExt1Im%mt)",shape(vExt1Im%mt)
            !print*,"shape(vExt1Im%mt)",shape(vExt1Im%mt)

            print*,"Vext1_pw",sum(vExt1%pw(:,:))
            print*,"vext1mtIm_1",sum(vExt1Im%mt(:,:,1,1))
            print*,"vext1mt_1",sum(vExt1%mt(:,:,1,1))

            !print*,"vext1mtIm_2",sum(vExt1Im%mt(:,:,2,1))
            !print*,"vext1mt_2",sum(vExt1%mt(:,:,2,1))

            !print*,"vext1mtIm",sum(vExt1Im%mt(:,:,:,1))
            !print*,"vext1mt",sum(vExt1%mt(:,:,:,1))
            !vExt1%mt = 0.0
            !vExt1Im%mt = 0.0
            offset_out = 0.0
            call  dfpt_potden_offset(1,fmpi,starsq_vext,fi%cell,fi%atoms,vExt1,vExt1Im,.FALSE.,.FALSE.,offset_out)
            vExt1Im%mt = 0.0
            !call  dfpt_potden_offset(1,fmpi,stars,fi%cell,fi%atoms,rho,vExt1Im,.FALSE.,.TRUE.,offset_out)
            print*,"offset_out",offset_out

        end subroutine dfpt_vefield_int


        subroutine dfpt_vefield_realspace_MT(fi,stars,sphhar,fmpi,rho,q_sign)

            use m_types
            use m_dfpt_potden_offset
            use m_make_stars
            USE m_ylm
            !use 

            type(t_fleurinput), intent(in)     :: fi
            TYPE(t_stars),      INTENT(IN)     :: stars
            type(t_sphhar),    intent(in)      :: sphhar
            type(t_potden), intent(inout)         :: rho

            type(t_mpi),        intent(in)     :: fmpi
            integer,      intent(in)           :: q_sign
            type(t_potden)                     :: vExt1, vExt1Im
            type(t_stars)                      :: starsq_vext
            integer                            :: iDir_col,n,lmax,imax
            real                               :: qvec_int(3)
            complex                            :: offset_out
            complex, allocatable               :: ylm(:)


            iDir_col =  1
            qvec_int(:) =0.0
            qvec_int= q_sign*fi%juPhon%qvec_efield(iDir_col,:)
            print*,"qvec_int",qvec_int
            print*,"fi%juPhon%qvec_efield(1,:)",fi%juPhon%qvec_efield(1,:)

            CALL starsq_vext%reset_stars()

            CALL make_stars(starsq_vext, fi%sym, fi%atoms, fi%vacuum, sphhar, fi%input, fi%cell, fi%noco, fmpi, qvec_int, 1, iDir_col,fi%juPhon%l_efield)
            starsq_vext%ufft = stars%ufft
            print*,"starsq_vext%center",starsq_vext%center
            call vExt1%init(starsq_vext, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
            call vExt1Im%init(starsq_vext, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)
            call dfpt_vefield(fi%juPhon,starsq_vext,fi%atoms,fi%sym,sphhar,fi%cell,vExt1,vExt1Im,iDir_col,1)
            print*,"wefdsf"
            do n =1, fi%atoms%ntype
                lmax = fi%atoms%lmax(n)
                imax = fi%atoms%jri(n)

                allocate(ylm((lmax+1)**2))

                !CALL ylm4(lmax, , ylm(:))



                deallocate(ylm)


            end do

            


        

        end subroutine dfpt_vefield_realspace_MT



end module m_dfpt_vefield
