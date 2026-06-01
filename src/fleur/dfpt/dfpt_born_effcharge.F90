!--------------------------------------------------------------------------------
! Copyright (c) 2024 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_dfpt_born_effcharge
    use m_types
    USE m_make_stars
    use m_convol
    use m_dfpt_dynmat
    use m_dfpt_vefield
    USE m_intgr, ONLY: intgr3
    USE m_step_function
    use m_inv3
    USE m_dfpt_vgen
    USE m_constants
    USE m_checkdopall
    USE m_types_fftGrid
    




    implicit none 

contains 

    subroutine dfpt_born_eff_charge_element_nef(sternheimerJob,fi,stars,starsq,sphhar,fmpi,rho,denIn1,denIn1Im,grRho,BEC_element,BEC_contributions_element,iDir_den,q_sign,hybdat,xcpot,nococonv,vTot)


        type(t_sternheimerJob), intent(in) :: sternheimerJob 
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
        complex                            :: offset_out,tempval_pw,tempval_mt,tempval_SF_IR,tempval_grrho
        real, allocatable                  :: denIn1_mt(:,:,:),denIn1_mt_Im(:,:,:),rho_mt(:,:,:), grRho_mt(:,:,:)!,mt_r(:,:),mt_Im(:,:)
        integer                            :: iType,iDir,iDtype
        real                               :: qvec_int(3),int_mt_r,int_mt_Im
        type(t_potden)                     :: den1_dummy,den1Im_dummy

        TYPE(t_potden)  :: rho_loc0

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
                print*,"shape(q_vec)",shape(fi%dfpt%qvec_efield)
                qvec_int= 0.0!fi%dfpt%qvec_efield(1,:)
                CALL make_stars(starsq_vextpho, fi%sym, fi%atoms, fi%vacuum, sphhar, fi%input, fi%cell, fi%noco, fmpi, qvec_int, iDtype, iDir)
                starsq_vextpho%ufft = stars%ufft

                CALL den1_dummy%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_DEN, l_dfpt=.TRUE.)
                CALL den1Im_dummy%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_DEN, l_dfpt=.FALSE.)

                CALL vExt1pho%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
                CALL vExt1Impho%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)
                print*,"right before it"
                CALL dfpt_vgen(sternheimerJob,hybdat,fi%field,fi%input,xcpot,fi%atoms,sphhar,stars,fi%vacuum,fi%sym,&
                            fi%dfpt,fi%cell,fmpi,fi%noco,nococonv,rho_loc0,vTot,&
                            starsq,den1Im_dummy,vExt1pho,.FALSE.,vExt1Impho,den1_dummy,iDtype,iDir,[1,1],l_vextpho=.TRUE.)
                
                print*,"sum(vExt1pho%pw(:,1))", sum(vExt1pho%pw(:,1))
                !stop
                tempval_pw = CMPLX(0.0,0.0)
                ! call make_stars
                !qvec_int= fi%dfpt%qvec_efield(iQ,:)
                !print*,"qvec_int",qvec_int
                !CALL starsq_vext%reset_stars()
                !CALL make_stars(starsq_vext, fi%sym, fi%atoms, fi%vacuum, sphhar, fi%input, fi%cell, fi%noco, fmpi, qvec_int, 1, iQ,.true.)
                !starsq_vext%ufft = starsq%ufft !why again
                !call vExt1%init(starsq_vext, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
                !call vExt1Im%init(starsq_vext, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)
                !call dfpt_vefield(fi%dfpt,starsq_vext,fi%atoms,fi%sym,sphhar,fi%cell,vExt1,vExt1Im,iQ,q_sign)
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
        type(t_potden), intent(in)         :: rho
        ! type(t_potden), intent(inout)         :: rho_core
        type(t_potden), intent(in)         :: denIn1,denIn1Im
        type(t_potden),optional, intent(in)     ::grRho
        TYPE(t_mpi),        intent(in)     :: fmpi
        !complex, intent(in)                 :: grrho_val(:,:,:,:)
        complex, intent(inout)             :: BEC_element
        complex, intent(inout)             :: BEC_contributions_element(:)
        integer, intent(in)                :: iDir_den,iQ,q_sign,iDType

        type(t_potden)                     :: vExt1, vExt1Im
        type(t_stars)                      :: starsq_vext
        TYPE(t_fftgrid)                    :: fftgrid_dummy
        complex, allocatable               :: pwwq2(:),denIn1_pw(:),rho_pw(:),theta1full(:,:,:),theta1_pw0(:,:,:),theta1_pw(:,:,:)
        complex                            :: offset_out,tempval_pw,tempval_mt,tempval_SF_IR,tempval_grrho,tempval_SF
        real, allocatable                  :: denIn1_mt(:,:,:),denIn1_mt_Im(:,:,:),rho_mt(:,:,:), grRho_mt(:,:,:)!,mt_r(:,:),mt_Im(:,:)
        integer                            :: iType,iDir
        real                               :: qvec_int(3),int_mt_r,int_mt_Im
        

        
        ALLOCATE(pwwq2(starsq%ng3))
        ALLOCATE(denIn1_pw(starsq%ng3))
        ALLOCATE(rho_pw(stars%ng3))
        ALLOCATE(denIn1_mt(fi%atoms%jmtd,0:sphhar%nlhd,fi%atoms%ntype),denIn1_mt_Im(fi%atoms%jmtd,0:sphhar%nlhd,fi%atoms%ntype))
        ALLOCATE(rho_mt(fi%atoms%jmtd,0:sphhar%nlhd,fi%atoms%ntype))
        ALLOCATE(theta1full(0:27*starsq%mx1*starsq%mx2*starsq%mx3-1,fi%atoms%ntype,3))
        !ALLOCATE(theta1_pw0(stars%ng3,fi%atoms%ntype,3))
        ALLOCATE(theta1_pw(starsq%ng3,fi%atoms%ntype,3))


        !remove spin dependence
        denIn1_pw  = (denIn1%pw(:,1)+denIn1%pw(:,fi%input%jspins))/(3.0-fi%input%jspins)
        denIn1_mt = (denIn1%mt(:,0:,:,1)+denIn1%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)
        denIn1_mt_Im = (denIn1Im%mt(:,0:,:,1)+denIn1Im%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)


        IF (.FALSE.) denIn1_mt(:,0:,iDType) = &
                            denIn1_mt(:,0:,iDType) - &
                            (grRho%mt(:,0:,iDType,1)+grRho%mt(:,0:,iDType,fi%input%jspins))/(3.0-fi%input%jspins)
        ! This is some debugging stuff 
        !IF (.FALSE.) denIn1%mt(:,0:,iDType,1) = &
        !                    denIn1%mt(:,0:,iDType,1) - &
        !                    (grRho%mt(:,0:,iDType,1)+grRho%mt(:,0:,iDType,fi%input%jspins))/(3.0-fi%input%jspins)
        rho_pw = (rho%pw(:,1)+rho%pw(:,fi%input%jspins))/(3.0-fi%input%jspins)
        rho_mt = (rho%mt(:,0:,:,1)+rho%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)
        
        int_mt_r = 0.0
        int_mt_Im = 0.0

        tempval_pw = CMPLX(0.0,0.0)

        call vExt1%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
        call vExt1Im%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)
        call dfpt_vefield(fi%dfpt,starsq,fi%atoms,fi%sym,sphhar,fi%cell,vExt1,vExt1Im,iQ,q_sign)

        !interstitial
        pwwq2 = CMPLX(0.0,0.0)
        CALL dfpt_convol_big(1, starsq, stars, vExt1%pw(:,1), CMPLX(1.0,0.0)*stars%ufft, pwwq2)
        CALL dfpt_int_pw(starsq, fi%cell, denIn1_pw, pwwq2, tempval_pw)
        !call dfpt_checkdopall(fi, sphhar, starsq,vExt1,vExt1Im,denIn1,denIn1Im) 
        !print*,"convol(theta*vext1)*n1",tempval_pw
        BEC_element = BEC_element + 2*tempval_pw 
        BEC_contributions_element(1) = BEC_contributions_element(1) + tempval_pw
        BEC_contributions_element(8) =  BEC_contributions_element(8) +2*tempval_pw 

        !IR integral differently:
        !use that Fourier expression of Step function
        !tempval_pw = CMPLX(0.0,0.0)
        !pwwq2 = CMPLX(0.0,0.0)
        !pwwq2 =vExt1%pw(1,1)*stars%ustep
        !CALL dfpt_int_pw(starsq, fi%cell, denIn1_pw, pwwq2, tempval_pw)
        !print*,"using u step",tempval_pw
        !change order of convolution
        !tempval_pw = CMPLX(0.0,0.0)
        !pwwq2 = CMPLX(0.0,0.0)
        !CALL dfpt_convol_big(1, starsq, stars,denIn1_pw, CMPLX(1.0,0.0)*stars%ufft, pwwq2)
        !CALL dfpt_int_pw(starsq, fi%cell,  vExt1%pw(:,1), pwwq2, tempval_pw)
        !print*,"convol(theta*n1)*vext1",tempval_pw

        !tempval_pw = CMPLX(0.0,0.0)
        !tempval_pw = tempval_pw + fi%cell%omtil*conjg(vExt1%pw(1,1))*pwwq2(1)
        !print*,"use that only G=0 zero contributes:", tempval_pw


        !stop
        !MT integral:   
        do iType = 1, fi%atoms%ntype
            tempval_mt = CMPLX(0.0,0.0) 
            call dfpt_int_mt(fi%atoms, sphhar, fi%sym, iType, denIn1_mt, denIn1_mt_Im, vExt1%mt(:,0:,:,1), vExt1Im%mt(:,0:,:,1), tempval_mt)
            BEC_element =  BEC_element + 2*tempval_mt
            BEC_contributions_element(2) =  BEC_contributions_element(2) + tempval_mt
            BEC_contributions_element(8+iType) =  BEC_contributions_element(8+iType) + tempval_mt
            !print*,"fi%atoms%ntype",fi%atoms%ntype
            BEC_contributions_element(8) =  BEC_contributions_element(8) + 2*tempval_mt
        end do
        
        !surface integral MT
        if (iQ .eq. iDir_den) then
            int_mt_r =0.0
            call intgr3(rho%mt(:,0,iDType,1),fi%atoms%rmsh(1,iDType),fi%atoms%dx(iDType),fi%atoms%jri(iDType),int_mt_r)
            !call intgr3(rho%mt(:,0,iDType,1)-rho_core%mt(:,0,iDType,1),fi%atoms%rmsh(1,iDType),fi%atoms%dx(iDType),fi%atoms%jri(iDType),int_mt_r)

            !print*,"offset_out",offset_out
            !CALL intgr3(denIn1_mt_Im(:,0,iDType),fi%atoms%rmsh(1,iDType),fi%atoms%dx(iDType),fi%atoms%jri(iDType),int_mt_Im)
            BEC_element = BEC_element+fi%atoms%zatom(iDType)-int_mt_r*sfp_const!+int_mt_Im)
            !BEC_element = BEC_element+fi%atoms%econf(iDType)%valence_electrons-int_mt_r*sfp_const!+int_mt_Im)
            BEC_contributions_element(3) =  BEC_contributions_element(3)  -int_mt_r*sfp_const!+int_mt_Im)
            BEC_contributions_element(7) =  BEC_contributions_element(7)  -int_mt_r*sfp_const
            BEC_contributions_element(6) =  BEC_contributions_element(6)  +fi%atoms%zatom(iDType)
            BEC_contributions_element(8) = BEC_contributions_element(8)+fi%atoms%zatom(iDType)
            !BEC_contributions_element(8) =  BEC_contributions_element(8)  +fi%atoms%econf(iDType)%valence_electrons


        end if 
        !surface integral IR
        pwwq2 = CMPLX(0.0,0.0)
        theta1full = CMPLX(0.0,0.0)
        !CALL stepf_analytical(fi%sym, stars, fi%atoms, fi%input, fi%cell, fmpi, fftgrid_dummy, [0.0,0.0,0.0], iDType, iDir_den, 1, theta1full0(0:,:,:))
        CALL stepf_analytical(fi%sym, starsq, fi%atoms, fi%input, fi%cell, fmpi, fftgrid_dummy,q_sign*fi%dfpt%qvec_efield(:,iQ), iDType, iDir_den, 1, theta1full)   
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
        !print*,"sum(rho_pw)",sum(rho_pw)
        !print*,"sum(theta1full(0:,iDType,iDir_den))",sum(theta1full(0:,iDType,iDir_den))
        CALL dfpt_convol_big(2, stars, starsq, rho_pw, theta1full(0:,iDType,iDir_den), pwwq2)
        !print*,"sum(pwwq2)",sum(pwwq2)
        CALL dfpt_int_pw(starsq, fi%cell,pwwq2, vExt1%pw(:,1),   tempval_SF_IR)
        !print*,"vExt1%pw(:,1)",sum(vExt1%pw(:,1))
        !print*,"tempval_SF_IR",tempval_SF_IR
        !stop
        BEC_element = BEC_element+ 2*tempval_SF_IR
        BEC_contributions_element(4) =  BEC_contributions_element(4) + tempval_SF_IR
        BEC_contributions_element(7) =  BEC_contributions_element(7) +2*tempval_SF_IR
        
        write(9989,FMT=8000) "IR Theta1 rho V1ext new              ", tempval_SF_IR
        tempval_SF = cmplx(0.0,0.0)
        CALL dfpt_int_mt_sf(fi%atoms, sphhar, fi%sym, iDir_den, iDType, rho_mt, vExt1%mt(:,0:,:,1), vExt1Im%mt(:,0:,:,1), tempval_SF)
        !dyn_row_HF(col_index) = dyn_row_HF(col_index) + tempval
        write(9989,FMT=8000) "SF rho Vext1 efield                 ", tempval_SF
        !write(9989,'(A,1X,G0)') "SF rho Vext1 efield", tempval_SF

        !calculate gradient integral:
        tempval_grrho = 0.0 
        !print*,"grrho_val imag", sum(grrho_val)
        !stop
        grRho_mt =(grRho%mt(:,0:,:,1)+grRho%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)
        !grRho_mt =(grrho_val(:,0:,:,1)+grrho_val(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)
        call dfpt_int_mt(fi%atoms, sphhar, fi%sym, iDType, grRho_mt, grRho_mt*0, vExt1%mt(:,0:,:,1), vExt1Im%mt(:,0:,:,1), tempval_grrho)
        !print*,"tempval_grrho",tempval_grrho
        !stop
        BEC_contributions_element(5) = BEC_contributions_element(5) + tempval_grrho
        BEC_contributions_element(7) =  BEC_contributions_element(7)  +2*tempval_grrho
        BEC_contributions_element(8) = BEC_contributions_element(8)-2*tempval_grrho
        
        8000  FORMAT (a,2E19.8E2)
    end subroutine dfpt_born_eff_charge_element



    subroutine dfpt_born_eff_charge_final(fi,born_eff_charge,born_eff_charge_contributions)

         USE m_xmlOutput

        type(t_fleurinput), intent(in)    :: fi
        complex, intent(inout)   :: born_eff_charge(:,:,:)
        complex, intent(inout)   :: born_eff_charge_contributions(:,:,:,:)
        integer                  :: iDType,iDir, j 
        complex                  :: dielten_iden(3,3) 
        character(len=20)        :: atom_string
        character(len=20)        :: filename
        integer                  ::i, file_int
        CHARACTER(LEN=20)        :: attributes(2)


        !born_eff_charge(:,:,:) = -born_eff_charge(:,:,:)
        
        atom_string = 'atom No:'

        open( 111, file="born_eff_charge", status='replace', action='write', form='formatted')
        write(*,*) '-------------------------' 
        write(*,*) "Born Effective Charge" 
        do iDType = 1, fi%atoms%ntype
            !write(*, '(A,I4,1X,A)', ADVANCE='NO') atom_string, iDType,fi%atoms%speciesName(iDType)
            write(*,'(A,I4,1X,A)') atom_string, iDType,fi%atoms%speciesName(iDType)
            write(111,'(A,I4,1X,A)') atom_string, iDType,fi%atoms%speciesName(iDType)
            do iDir = 1,3
                do j = 1,2
                    write(111,'(es16.8)', ADVANCE='NO') real(born_eff_charge(iDType,iDir,j)) 
                    write(111, '(A)', ADVANCE='NO') ' ' 
                    write(*,'(es16.8)', ADVANCE='NO') real(born_eff_charge(iDType,iDir,j))
                    write(*, '(A)', ADVANCE='NO') ' ' 
                end do
            write(111,'(es16.8)')real(born_eff_charge(iDType,iDir,3))
            write(*,'(es16.8)')real(born_eff_charge(iDType,iDir,3))
            end do
        end do
        close(111)
        write(*,*) '-------------------------' 

        !save in out.xml
        CALL openXMLElementNoAttributes('Phonons')
        CALL openXMLElementNoAttributes('efield')
        do iDType = 1, fi%atoms%ntype 
            do iDir = 1,3
                attributes = ''
                WRITE(attributes(1),'(i0)') iDType
                WRITE(attributes(2),'(i0)') iDir
                CALL writeXMLElementPoly('borneffcharge',(/ 'iDtype' , 'iDir  '/), attributes,real(born_eff_charge(iDType,iDir,:)))
            end do
        end do
        CALL closeXMLElement('efield')
        CALL closeXMLElement('Phonons')
        
        atom_string = 'atom No:'
        ! Turn off debugging stuff
        IF (.FALSE.) THEN
            IF (fi%dfpt%l_efield) THEN
                call write_born_effective_charge(112, "born_eff_charge", born_eff_charge, atom_string, fi)
                call write_born_effective_charge(112, "born_eff_charge_MT", -born_eff_charge_contributions(:,:,:,1), atom_string, fi)
                call write_born_effective_charge(112, "born_eff_charge_MT_1", -born_eff_charge_contributions(:,:,:,2), atom_string, fi)
                call write_born_effective_charge(112, "born_eff_charge_MT_2", -born_eff_charge_contributions(:,:,:,3), atom_string, fi)
            ELSE
                call write_born_effective_charge(112, "born_eff_charge_cmplx", born_eff_charge, atom_string, fi)
                call write_born_effective_charge(113, "born_eff_charge_IR", born_eff_charge_contributions(:,:,:,1), atom_string, fi)
                call write_born_effective_charge(114, "born_eff_charge_MT", born_eff_charge_contributions(:,:,:,2), atom_string, fi)
                call write_born_effective_charge(115, "born_eff_charge_SF-MT", born_eff_charge_contributions(:,:,:,3), atom_string, fi)
                call write_born_effective_charge(116, "born_eff_charge_SF-IR", born_eff_charge_contributions(:,:,:,4), atom_string, fi)
                call write_born_effective_charge(117, "born_eff_charge_grRho", born_eff_charge_contributions(:,:,:,5), atom_string, fi)
                call write_born_effective_charge(118, "born_eff_charge_Zbare", born_eff_charge_contributions(:,:,:,6), atom_string, fi)
                call write_born_effective_charge(119, "born_eff_charge_SF", born_eff_charge_contributions(:,:,:,7), atom_string, fi)
                call write_born_effective_charge(120, "born_eff_charge_pulay_only", born_eff_charge_contributions(:,:,:,8), atom_string, fi)
                DO i =1,fi%atoms%nat
                    file_int = 120+i
                    write(filename, '("born_eff_charge_MT_",i0)') i
                    !top
                    call write_born_effective_charge(120, filename, born_eff_charge_contributions(:,:,:,8+i), atom_string, fi)
                    !call write_born_effective_charge(121, "born_eff_charge_MT_2", born_eff_charge_contributions(:,:,:,9), atom_string, fi)
                END DO
            END IF
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
        !write(*,*) '-------------------------'
        !write(*,*) filename
        
        do iDType = 1, fi%atoms%ntype
            !write(*,'(A,I4)') atom_string, iDType
            write(file_id,'(A,I4,1X,A)') atom_string, iDType,fi%atoms%speciesName(iDType)
    
            do iDir = 1, 3
                do j = 1, 2
                    write(file_id, '(2es16.8)', ADVANCE='NO') BEC_contribution(iDType, iDir, j)
                    write(file_id, '(A)', ADVANCE='NO') ' '
                    !write(*, '(2es16.8)', ADVANCE='NO') BEC_contribution(iDType, iDir, j)
                    !write(*, '(A)', ADVANCE='NO') ' '
                end do
    
                write(file_id, '(2es16.8)') BEC_contribution(iDType, iDir, 3)
                !write(*, '(2es16.8)') BEC_contribution(iDType, iDir, 3)
            end do
        end do
    
        close(file_id)
        !write(*,*) '-------------------------'
    end subroutine write_born_effective_charge

        subroutine dfpt_checkdopall(fi, sphhar, starsq,vExt1,vExt1Im,denIn1,denIn1Im)
        
        type(t_fleurinput), intent(in)     :: fi
        type(t_sphhar),    intent(in)      :: sphhar
        TYPE(t_stars),      INTENT(IN)     ::  starsq
        type(t_potden), intent(in)         :: vExt1, vExt1Im
        type(t_potden), intent(in)         :: denIn1,denIn1Im

        type(t_potden)                     :: product,productIm
        TYPE(t_fftgrid) :: fftgrid_pot, fftgrid_den, fftgrid_product

        WRITE (oUnit,*) "potential perturbation checkdopall"
        CALL checkDOPALL(fi%input, sphhar, starsq,fi%atoms, fi%sym, fi%vacuum, fi%cell,vExt1,1, vExt1Im,'vext1')
        
        WRITE (oUnit,*) "density response checkdopall"
        CALL checkDOPALL(fi%input, sphhar, starsq,fi%atoms, fi%sym, fi%vacuum, fi%cell,denIn1,1, denIn1Im,'den1')

        product =vExt1
        productIm =vExt1
        CALL product%resetPotDen()
        CALL productIm%resetPotDen()
        print*,"vext1%pw(:,1)",sum(vext1%pw(:,1))
        print*,"denIn1%pw(:,1)",sum(denIn1%pw(:,1))
        print*,"product%pw(:,1)",sum(product%pw(:,1))
        !Do the PW stuff
        CALL fftgrid_pot%init((/3*starsq%mx1,3*starsq%mx2,3*starsq%mx3/))
        CALL fftgrid_den%init((/3*starsq%mx1,3*starsq%mx2,3*starsq%mx3/))
        !CALL fftgrid_product%init((/3*starsq%mx1,3*starsq%mx2,3*starsq%mx3/))

        CALL fftgrid_pot%putFieldOnGrid(starsq,vExt1%pw(:,1))
        !stop
        CALL fftgrid_pot%perform_fft(forward=.false.)

        CALL fftgrid_den%putFieldOnGrid(starsq,denIn1%pw(:,1))
        !stop
        CALL fftgrid_den%perform_fft(forward=.false.)

        print*,"fftgrid_pot%grid",(fftgrid_pot%grid)
        print*,"fftgrid_pot%grid",sum(fftgrid_pot%grid)
        print*,"fftgrid_den%grid",sum(fftgrid_den%grid)
        !stop
        fftgrid_pot%grid =(fftgrid_pot%grid)*(fftgrid_den%grid)
        print*,"fftgrid_pot%grid",sum(fftgrid_pot%grid)
        CALL fftgrid_pot%perform_fft(forward=.true.)
        CALL fftgrid_pot%takeFieldFromGrid(starsq,product%pw(:,1))
        product%pw(:,1) = product%pw(:,1)*starsq%nstr
        print*,",product%pw(:,1)",sum(product%pw(:,1))
        WRITE (oUnit,*) "trying new stuff"
        !Do MT shit
        product%mt= 0.0
        productIm%mt= 0.0
        !CALL product_MT(fi,sphhar,product%mt,productIm%mt)
        !stop
        CALL checkDOPALL(fi%input, sphhar, starsq,fi%atoms, fi%sym, fi%vacuum, fi%cell,product,1, product,'v1n1')

        stop
    
    end subroutine dfpt_checkdopall


        
end module 