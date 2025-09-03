module m_dfpt_born_effcharge
    use m_types
    USE m_make_stars
    use m_convol
    use m_dfpt_dynmat
    use m_dfpt_vefield
    !use m_dfpt_potden_offset
    USE m_intgr, ONLY: intgr3
    USE m_step_function
    use m_inv3
    USE m_dfpt_vgen


    implicit none 

    contains 

        subroutine dfpt_born_eff_charge_element_nef(fi,stars,starsq,sphhar,fmpi,rho,denIn1,denIn1Im,grRho,BEC_element,BEC_contributions_element,iDir_den,q_sign,hybdat,xcpot,nococonv,vTot)


            type(t_fleurinput), intent(in)     :: fi
            type(t_sphhar),    intent(in)      :: sphhar
            TYPE(t_stars),      INTENT(IN)     :: stars, starsq
            type(t_potden), intent(inout)         :: rho
            type(t_potden), intent(inout)         :: denIn1,denIn1Im
            type(t_potden),optional, intent(in)     ::grRho
            TYPE(t_mpi),        intent(in)     :: fmpi
            complex, intent(inout)             :: BEC_element(:,:)
            complex, intent(inout)             :: BEC_contributions_element(:,:,:)
            integer, intent(in)                :: iDir_den,q_sign

            TYPE(t_hybdat),     INTENT(INOUT) :: hybdat
            CLASS(t_xcpot),     INTENT(IN)    :: xcpot
            TYPE(t_nococonv),   INTENT(IN)    :: nococonv
            TYPE(t_potden),     INTENT(IN)    :: vTot

            type(t_potden)                     :: vExt1pho, vExt1Impho
            type(t_stars)                      :: starsq_vextpho
            TYPE(t_fftgrid)                    :: fftgrid_dummy
            complex, allocatable               :: pwwq2(:),pww2(:),denIn1_pw(:),rho_pw(:),theta1full0(:,:,:),theta1full(:,:,:),theta1_pw0(:,:,:),theta1_pw(:,:,:)
            complex                            :: offset_out,tempval_pw,tempval_mt,tempval_SF_IR,tempval_grrho,sigma_loc(2)
            real, allocatable                  :: denIn1_mt(:,:,:),denIn1_mt_Im(:,:,:),rho_mt(:,:,:), grRho_mt(:,:,:)!,mt_r(:,:),mt_Im(:,:)
            integer                            :: iType,iDir,iDtype
            real                               :: qvec_int(3),int_mt_r,int_mt_Im
            type(t_potden)                     :: den1_dummy,den1Im_dummy

            TYPE(t_potden)  :: rho_loc0
            TYPE(t_juphon)                     :: juphon_int

            
            juphon_int = fi%juPhon
            juphon_int%l_efield = .FALSE.
            juphon_int%l_borneffcharge = .FALSE.
            juphon_int%l_phonon = .TRUE.
            !CALL dfpt_vefield_int(fi,stars,sphhar,fmpi,rho,-1)
            !CALL dfpt_vefield_int(fi,stars,sphhar,fmpi,rho,-1)
            !stop

            !write(*,*) 'Born effective charge'
            !print*,"shape(BEC_row(:))",shape(BEC_row)
            !print*,BEC_row(:)
            !print*,"shape(denIn1%mt)",shape(denIn1%mt)
            
            ALLOCATE(pwwq2(starsq%ng3))
            ALLOCATE(denIn1_pw(starsq%ng3))
            ALLOCATE(rho_pw(stars%ng3),pww2(stars%ng3))
            ALLOCATE(denIn1_mt(fi%atoms%jmtd,0:sphhar%nlhd,fi%atoms%ntype),denIn1_mt_Im(fi%atoms%jmtd,0:sphhar%nlhd,fi%atoms%ntype))
            ALLOCATE(rho_mt(fi%atoms%jmtd,0:sphhar%nlhd,fi%atoms%ntype))
            ALLOCATE(theta1full0(0:27*stars%mx1*stars%mx2*stars%mx3-1,fi%atoms%ntype,3))
            ALLOCATE(theta1full(0:27*starsq%mx1*starsq%mx2*starsq%mx3-1,fi%atoms%ntype,3))
            ALLOCATE(theta1_pw0(stars%ng3,fi%atoms%ntype,3))
            ALLOCATE(theta1_pw(starsq%ng3,fi%atoms%ntype,3))
            !ALLOCATE(mt_r(fi%atoms%jmtd,fi%atoms%ntype))
            !ALLOCATE(mt_Im(fi%atoms%jmtd,fi%atoms%ntype))

            !remove spin dependence
            denIn1_pw  = (denIn1%pw(:,1)+denIn1%pw(:,fi%input%jspins))/(3.0-fi%input%jspins)
            denIn1_mt = (denIn1%mt(:,0:,:,1)+denIn1%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)
            denIn1_mt_Im = (denIn1Im%mt(:,0:,:,1)+denIn1Im%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)



            CALL rho_loc0%copyPotDen(denIn1)
            CALL rho_loc0%resetPotDen()
            
            do iDtype = 1,fi%atoms%ntype
                do iDir =1,3
                    CALL starsq_vextpho%reset_stars()
                    print*,"shape(q_vec)",shape(fi%juPhon%qvec_efield)
                    qvec_int= 0.0!fi%juPhon%qvec_efield(1,:)
                    CALL make_stars(starsq_vextpho, fi%sym, fi%atoms, fi%vacuum, sphhar, fi%input, fi%cell, fi%noco, fmpi, qvec_int, iDtype, iDir)
                    starsq_vextpho%ufft = stars%ufft

                    CALL den1_dummy%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_DEN, l_dfpt=.TRUE.)
                    CALL den1Im_dummy%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_DEN, l_dfpt=.FALSE.)

                    CALL vExt1pho%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
                    CALL vExt1Impho%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)
                    
                    sigma_loc = cmplx(0.0,0.0)
                    print*,"right before it"
                    CALL dfpt_vgen(hybdat,fi%field,fi%input,xcpot,fi%atoms,sphhar,stars,fi%vacuum,fi%sym,&
                                juphon_int,fi%cell,fmpi,fi%noco,nococonv,rho_loc0,vTot,&
                                starsq,den1Im_dummy,vExt1pho,.FALSE.,vExt1Impho,den1_dummy,iDtype,iDir,[1,1],sigma_loc,l_vextpho=.TRUE.)
                    
                    print*,"sum(vExt1pho%pw(:,1))", sum(vExt1pho%pw(:,1))
                    !stop
                    tempval_pw = CMPLX(0.0,0.0)
                    ! call make_stars
                    !qvec_int= fi%juPhon%qvec_efield(iQ,:)
                    !print*,"qvec_int",qvec_int
                    !CALL starsq_vext%reset_stars()
                    !CALL make_stars(starsq_vext, fi%sym, fi%atoms, fi%vacuum, sphhar, fi%input, fi%cell, fi%noco, fmpi, qvec_int, 1, iQ,.true.)
                    !starsq_vext%ufft = starsq%ufft !why again
                    !call vExt1%init(starsq_vext, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
                    !call vExt1Im%init(starsq_vext, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)
                    !call dfpt_vefield(fi%juPhon,starsq_vext,fi%atoms,fi%sym,sphhar,fi%cell,vExt1,vExt1Im,iQ,q_sign)
                    !interstitial
                    pwwq2 = CMPLX(0.0,0.0)
                    !print*,"I make it this far"
                    !print*,"sum(vExt1pho%pw(:,1))",sum(vExt1pho%pw(:,1))
                    CALL dfpt_convol_big(1, starsq, stars, vExt1pho%pw(:,1), CMPLX(1.0,0.0)*stars%ufft, pwwq2)
                    print*,"pwwq2",sum(pwwq2)
                    !stop
                    !call save_npy("pwwq2.npy",pwwq2)
                    !stop
                    CALL dfpt_int_pw(starsq, fi%cell, denIn1_pw, pwwq2, tempval_pw)
                    BEC_element(iDtype,iDir) = BEC_element(iDtype,iDir) + tempval_pw 
                    !BEC_contributions_element(1) = BEC_contributions_element(1) + tempval_pw
                    !print*,"I make it this far 2"
                    !stop
                    do iType = 1, fi%atoms%ntype
                        tempval_mt = CMPLX(0.0,0.0) 
                        call dfpt_int_mt(fi%atoms, sphhar, fi%sym, iType, denIn1_mt, denIn1_mt_Im, vExt1pho%mt(:,0:,:,1), vExt1Impho%mt(:,0:,:,1), tempval_mt)
                        BEC_element(iDtype,iDir) =  BEC_element(iDtype,iDir) + tempval_mt
                        print*,"BEC_element(iDtype,iDir) ",BEC_element(iDtype,iDir) 
                        BEC_contributions_element(iDtype,iDir,1) =  BEC_contributions_element(iDtype,iDir,1) + tempval_mt
                        BEC_contributions_element(iDtype,iDir,1+iType) =  BEC_contributions_element(iDtype,iDir,1+iType) + tempval_mt
                    end do
                    !stop
                end do
            end do
            

            write(*,*) 'Determine integrals here'
        
        end subroutine dfpt_born_eff_charge_element_nef


        subroutine dfpt_born_eff_charge_element(fi,stars,starsq,sphhar,fmpi,rho,denIn1,denIn1Im,grRho,BEC_element,BEC_contributions_element,iDir_den,iDType,iQ,q_sign)


            type(t_fleurinput), intent(in)     :: fi
            type(t_sphhar),    intent(in)      :: sphhar
            TYPE(t_stars),      INTENT(IN)     :: stars, starsq
            type(t_potden), intent(inout)         :: rho
            type(t_potden), intent(inout)         :: denIn1,denIn1Im
            type(t_potden),optional, intent(in)     ::grRho
            TYPE(t_mpi),        intent(in)     :: fmpi
            complex, intent(inout)             :: BEC_element
            complex, intent(inout)             :: BEC_contributions_element(:)
            integer, intent(in)                :: iDir_den,iQ,q_sign,iDType

            type(t_potden)                     :: vExt1, vExt1Im
            type(t_stars)                      :: starsq_vext
            TYPE(t_fftgrid)                    :: fftgrid_dummy
            complex, allocatable               :: pwwq2(:),pww2(:),denIn1_pw(:),rho_pw(:),theta1full0(:,:,:),theta1full(:,:,:),theta1_pw0(:,:,:),theta1_pw(:,:,:)
            complex                            :: offset_out,tempval_pw,tempval_mt,tempval_SF_IR,tempval_grrho
            real, allocatable                  :: denIn1_mt(:,:,:),denIn1_mt_Im(:,:,:),rho_mt(:,:,:), grRho_mt(:,:,:)!,mt_r(:,:),mt_Im(:,:)
            integer                            :: iType,iDir
            real                               :: qvec_int(3),int_mt_r,int_mt_Im
            

            !CALL dfpt_vefield_int(fi,stars,sphhar,fmpi,rho,-1)
            !CALL dfpt_vefield_int(fi,stars,sphhar,fmpi,rho,-1)
            !stop

            !write(*,*) 'Born effective charge'
            !print*,"shape(BEC_row(:))",shape(BEC_row)
            !print*,BEC_row(:)
            !print*,"shape(denIn1%mt)",shape(denIn1%mt)
            
            ALLOCATE(pwwq2(starsq%ng3))
            ALLOCATE(denIn1_pw(starsq%ng3))
            ALLOCATE(rho_pw(stars%ng3),pww2(stars%ng3))
            ALLOCATE(denIn1_mt(fi%atoms%jmtd,0:sphhar%nlhd,fi%atoms%ntype),denIn1_mt_Im(fi%atoms%jmtd,0:sphhar%nlhd,fi%atoms%ntype))
            ALLOCATE(rho_mt(fi%atoms%jmtd,0:sphhar%nlhd,fi%atoms%ntype))
            ALLOCATE(theta1full0(0:27*stars%mx1*stars%mx2*stars%mx3-1,fi%atoms%ntype,3))
            ALLOCATE(theta1full(0:27*starsq%mx1*starsq%mx2*starsq%mx3-1,fi%atoms%ntype,3))
            ALLOCATE(theta1_pw0(stars%ng3,fi%atoms%ntype,3))
            ALLOCATE(theta1_pw(starsq%ng3,fi%atoms%ntype,3))
            !ALLOCATE(mt_r(fi%atoms%jmtd,fi%atoms%ntype))
            !ALLOCATE(mt_Im(fi%atoms%jmtd,fi%atoms%ntype))

            !remove spin dependence
            denIn1_pw  = (denIn1%pw(:,1)+denIn1%pw(:,fi%input%jspins))/(3.0-fi%input%jspins)
            denIn1_mt = (denIn1%mt(:,0:,:,1)+denIn1%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)
            denIn1_mt_Im = (denIn1Im%mt(:,0:,:,1)+denIn1Im%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)

            rho_pw = -(rho%pw(:,1)+rho%pw(:,fi%input%jspins))/(3.0-fi%input%jspins) !check if minus is NOT needed
            rho_mt = (rho%mt(:,0:,:,1)+rho%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)
            
            !print*,"shape(denIn1_mt)",shape(denIn1_mt)



            !print*,"shape(denIn1)",shape(denIn1%mt)
            !print*,"shape(denIn1Im%mt)",shape(denIn1Im%mt)
            
            !denIn1Im%mt = 0.0
            !call  dfpt_potden_offset(1,fmpi,starsq,fi%cell,fi%atoms,rho,denIn1Im,.FALSE.,.TRUE.,offset_out)
            !print*,"offset_out",offset_out
            !print*, "irank", fmpi%irank
            !print*,"juphon%qvec_efield",fi%juphon%qvec_efield
            !stop
            int_mt_r = 0.0
            int_mt_Im = 0.0

            tempval_pw = CMPLX(0.0,0.0)
            ! call make_stars
            qvec_int= fi%juPhon%qvec_efield(iQ,:)
            !print*,"qvec_int",qvec_int
            CALL starsq_vext%reset_stars()
            CALL make_stars(starsq_vext, fi%sym, fi%atoms, fi%vacuum, sphhar, fi%input, fi%cell, fi%noco, fmpi, qvec_int, 1, iQ,.true.)
            starsq_vext%ufft = starsq%ufft !why again
            call vExt1%init(starsq_vext, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
            call vExt1Im%init(starsq_vext, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)
            call dfpt_vefield(fi%juPhon,starsq_vext,fi%atoms,fi%sym,sphhar,fi%cell,vExt1,vExt1Im,iQ,q_sign)
            !interstitial
            pwwq2 = CMPLX(0.0,0.0)
            print*,"I make it this far"
            CALL dfpt_convol_big(1, starsq, stars, vExt1%pw(:,1), CMPLX(1.0,0.0)*stars%ufft, pwwq2)
            call save_npy("pwwq2.npy",pwwq2)
            CALL dfpt_int_pw(starsq, fi%cell, denIn1_pw, pwwq2, tempval_pw)
            !print*,"I make it this far 1/2"
            !print*,"BEC_element",BEC_element
            !print*,"tempval_pw",tempval_pw
            BEC_element = BEC_element + tempval_pw 
            BEC_contributions_element(1) = BEC_contributions_element(1) + tempval_pw
            !print*,"I make it this far 2"
            do iType = 1, fi%atoms%ntype
                tempval_mt = CMPLX(0.0,0.0) 
                call dfpt_int_mt(fi%atoms, sphhar, fi%sym, iType, denIn1_mt, denIn1_mt_Im, vExt1%mt(:,0:,:,1), vExt1Im%mt(:,0:,:,1), tempval_mt)
                BEC_element =  BEC_element + tempval_mt
                BEC_contributions_element(2) =  BEC_contributions_element(2) + tempval_mt
                BEC_contributions_element(4+iType) =  BEC_contributions_element(4+iType) + tempval_mt
            end do

            !surface integral MT
            if (iQ .eq. iDir_den) then
                print*,"shape(rho%mt)",shape(rho%mt)
                call intgr3(rho%mt(:,0,iDType,1),fi%atoms%rmsh(1,iDType),fi%atoms%dx(iDType),fi%atoms%jri(iDType),int_mt_r)
                denIn1Im%mt = 0.0
                !call dfpt_potden_offset(1,fmpi,starsq,fi%cell,fi%atoms,rho,denIn1Im,.FALSE.,.TRUE.,offset_out)
                !call  dfpt_potden_offset(1,fmpi,starsq,fi%cell,fi%atoms,rho,denIn1Im,.FALSE.,.TRUE.,offset_out)
                !print*,"offset_out",offset_out
                !CALL intgr3(denIn1_mt_Im(:,0,iDType),fi%atoms%rmsh(1,iDType),fi%atoms%dx(iDType),fi%atoms%jri(iDType),int_mt_Im)
                BEC_element = BEC_element-int_mt_r!+int_mt_Im)
                BEC_contributions_element(3) =  BEC_contributions_element(3)  -int_mt_r!+int_mt_Im)

            end if 
            !surface integral IR
            pww2 = CMPLX(0.0,0.0)
            theta1full0 = CMPLX(0.0,0.0)
            !CALL stepf_analytical(fi%sym, stars, fi%atoms, fi%input, fi%cell, fmpi, fftgrid_dummy, [0.0,0.0,0.0], iDType, iDir_den, 1, theta1full0(0:,:,:))
            CALL stepf_analytical(fi%sym, stars, fi%atoms, fi%input, fi%cell, fmpi, fftgrid_dummy,fi%juPhon%qvec_efield(iQ,:), iDType, iDir_den, 1, theta1full(0:,:,:))   
            !fftgrid_dummy%grid = theta1full0(0:,1,1)
            !CALL fftgrid_dummy%takeFieldFromGrid(stars, theta1_pw0(:))
            !theta1_pw0(:) = theta1_pw0(:) * 3 * stars%mx1 * 3 * stars%mx2 * 3 * stars%mx3
            !CALL fftgrid_dummy%perform_fft(forward=.false.)
            !theta1full0(0:,1,1) = fftgrid_dummy%grid
            
            DO iType = 1, fi%atoms%ntype
                DO iDir = 1, 3
                    fftgrid_dummy%grid = theta1full(0:, iType, iDir)
                    CALL fftgrid_dummy%takeFieldFromGrid(starsq, theta1_pw(:, iType, iDir))
                    theta1_pw(:, iType, iDir) = theta1_pw(:, iType, iDir) * 3 * starsq%mx1 * 3 * starsq%mx2 * 3 * starsq%mx3
                    CALL fftgrid_dummy%perform_fft(forward=.false.)
                    theta1full(0:, iType, iDir) = fftgrid_dummy%grid
                END DO
             END DO
            
            tempval_SF_IR = CMPLX(0.0,0.0)

            CALL dfpt_convol_big(1, stars, stars, rho_pw, theta1full(0:,iDtype,iDir_den), pww2)
            CALL dfpt_int_pw(stars, fi%cell, pww2, vExt1%pw(:,1), tempval_SF_IR)
            BEC_element= BEC_element+ tempval_SF_IR
            BEC_contributions_element(4) =  BEC_contributions_element(4) + tempval_SF_IR

            !calculate gradient integral:
            print*,"doing gradient shit"
            tempval_grrho = 0.0 
            grRho_mt = (grRho%mt(:,0:,:,1)+grRho%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)
            !print*,"grho",grRho_mt
            !print*,"sum(grho)",sum(grRho_mt)
            !print*,"iDType",iDType
            call dfpt_int_mt(fi%atoms, sphhar, fi%sym, iDType, grRho_mt, grRho_mt*0, vExt1%mt(:,0:,:,1), vExt1Im%mt(:,0:,:,1), tempval_grrho)
            !print*,"tempval_grrho",tempval_grrho
            !print*,"BEC_contributions_element(4+fi%atoms%ntype)",BEC_contributions_element(5+fi%atoms%ntype)
            BEC_contributions_element(5+fi%atoms%ntype) = BEC_contributions_element(5+fi%atoms%ntype) + tempval_grrho
            !print*,"BEC_contributions_element(4+fi%atoms%ntype)",BEC_contributions_element(5+fi%atoms%ntype)
            !stop

            write(*,*) 'Determine integrals here'
        
        end subroutine dfpt_born_eff_charge_element

        subroutine dfpt_born_eff_charge_final(fi,born_eff_charge,born_eff_charge_contributions)

            type(t_fleurinput), intent(in)    :: fi
            complex, intent(inout)   :: born_eff_charge(:,:,:)
            complex, intent(inout)   :: born_eff_charge_contributions(:,:,:,:)
            integer                  :: iDType,iDir, j 
            complex                  :: dielten_iden(3,3) 
            character(len=20)         :: atom_string


            born_eff_charge(:,:,:) = -born_eff_charge(:,:,:)
            
            atom_string = 'atom No:'

            open( 112, file="born_eff_charge", status='replace', action='write', form='formatted')
            write(*,*) '-------------------------' 
            write(*,*) "Born Effective Charge" 
            do iDType = 1, fi%atoms%ntype
                write(*,'(A,I4)') atom_string, iDType
                write(112,'(A,I4)') atom_string, iDType
                do iDir = 1,3
                do j = 1,2
                    write(112,'(2es16.8)', ADVANCE='NO') born_eff_charge(iDType,iDir,j) 
                    write(112, '(A)', ADVANCE='NO') ' ' 
                    write(*,'(2es16.8)', ADVANCE='NO') born_eff_charge(iDType,iDir,j)
                    write(*, '(A)', ADVANCE='NO') ' ' 
                end do
                write(112,'(2es16.8)')born_eff_charge(iDType,iDir,3)
                write(*,'(2es16.8)')born_eff_charge(iDType,iDir,3)
                end do
            end do
            close(112)
            write(*,*) '-------------------------' 

            !born_eff_charge_IR(:,:,:) = -born_eff_charge_IR(:,:,:)
            
            atom_string = 'atom No:'

            !atom_string = 'atom No:'
            IF (fi%juPhon%l_efield) THEN
                call write_born_effective_charge(112, "born_eff_charge", born_eff_charge, atom_string, fi)
                call write_born_effective_charge(112, "born_eff_charge_MT", -born_eff_charge_contributions(:,:,:,1), atom_string, fi)
                call write_born_effective_charge(112, "born_eff_charge_MT_1", -born_eff_charge_contributions(:,:,:,2), atom_string, fi)
                call write_born_effective_charge(112, "born_eff_charge_MT_2", -born_eff_charge_contributions(:,:,:,3), atom_string, fi)
            ELSE
                call write_born_effective_charge(112, "born_eff_charge", born_eff_charge, atom_string, fi)
                call write_born_effective_charge(113, "born_eff_charge_IR", -born_eff_charge_contributions(:,:,:,1), atom_string, fi)
                call write_born_effective_charge(114, "born_eff_charge_MT", -born_eff_charge_contributions(:,:,:,2), atom_string, fi)
                call write_born_effective_charge(115, "born_eff_charge_SF-MT", -born_eff_charge_contributions(:,:,:,3), atom_string, fi)
                call write_born_effective_charge(116, "born_eff_charge_SF-IR", -born_eff_charge_contributions(:,:,:,4), atom_string, fi)
                call write_born_effective_charge(117, "born_eff_charge_MT_1", -born_eff_charge_contributions(:,:,:,5), atom_string, fi)
                call write_born_effective_charge(118, "born_eff_charge_MT_2", -born_eff_charge_contributions(:,:,:,6), atom_string, fi)
                call write_born_effective_charge(119, "born_eff_charge_grRho", -born_eff_charge_contributions(:,:,:,7), atom_string, fi)
            END IF
        
        end subroutine dfpt_born_eff_charge_final

        subroutine write_born_effective_charge(file_id, filename, BEC_contribution, atom_string, fi)

            integer, intent(in) :: file_id
            character(len=*), intent(in) :: filename
            complex, intent(in) :: BEC_contribution(:,:,:)
            character(len=*), intent(in) :: atom_string
            type(t_fleurinput), intent(in) :: fi
        
            integer :: iDType, iDir, j
        
            open(file_id, file=filename, status='replace', action='write', form='formatted')
            write(*,*) '-------------------------'
            write(*,*) filename
            
            do iDType = 1, fi%atoms%ntype
                write(*,'(A,I4)') atom_string, iDType
                write(file_id,'(A,I4)') atom_string, iDType
        
                do iDir = 1, 3
                    do j = 1, 2
                        write(file_id, '(2es16.8)', ADVANCE='NO') BEC_contribution(iDType, iDir, j)
                        write(file_id, '(A)', ADVANCE='NO') ' '
                        write(*, '(2es16.8)', ADVANCE='NO') BEC_contribution(iDType, iDir, j)
                        write(*, '(A)', ADVANCE='NO') ' '
                    end do
        
                    write(file_id, '(2es16.8)') BEC_contribution(iDType, iDir, 3)
                    write(*, '(2es16.8)') BEC_contribution(iDType, iDir, 3)
                end do
            end do
        
            close(file_id)
            write(*,*) '-------------------------'
        end subroutine write_born_effective_charge
        
end module 