module m_dfpt_born_effcharge
    use m_types
    USE m_make_stars
    use m_convol
    use m_dfpt_dynmat
    use m_dfpt_vefield
    use m_dfpt_potden_offset
    USE m_intgr, ONLY: intgr3

    use m_inv3


    implicit none 

    contains 
        subroutine dfpt_born_eff_charge_element(fi,stars,starsq,sphhar,fmpi,rho,denIn1,denIn1Im,BEC_element,iDir_den,iDType,iQ,q_sign)


            type(t_fleurinput), intent(in)     :: fi
            type(t_sphhar),    intent(in)      :: sphhar
            TYPE(t_stars),      INTENT(IN)     :: stars, starsq
            type(t_potden), intent(inout)         :: rho
            type(t_potden), intent(inout)         :: denIn1,denIn1Im
            TYPE(t_mpi),        intent(in)     :: fmpi
            complex, intent(inout)             :: BEC_element
            integer, intent(in)                :: iDir_den,iQ,q_sign

            type(t_potden)                     :: vExt1, vExt1Im
            type(t_stars)                      :: starsq_vext
            complex, allocatable               :: pwwq2(:),tempval_pw,tempval_mt, denIn1_pw(:)
            real, allocatable                  :: denIn1_mt(:,:,:),denIn1_mt_Im(:,:,:)!,mt_r(:,:),mt_Im(:,:)
            integer                            :: iType,iDType
            real                               :: qvec_int(3),int_mt_r,int_mt_Im
            complex                            :: offset_out



            !write(*,*) 'Born effective charge'
            !print*,"shape(BEC_row(:))",shape(BEC_row)
            !print*,BEC_row(:)
            !print*,"shape(denIn1%mt)",shape(denIn1%mt)
            
            ALLOCATE(pwwq2(starsq%ng3))
            ALLOCATE(denIn1_pw(starsq%ng3))
            ALLOCATE(denIn1_mt(fi%atoms%jmtd,0:sphhar%nlhd,fi%atoms%ntype),denIn1_mt_Im(fi%atoms%jmtd,0:sphhar%nlhd,fi%atoms%ntype))
            !ALLOCATE(mt_r(fi%atoms%jmtd,fi%atoms%ntype))
            !ALLOCATE(mt_Im(fi%atoms%jmtd,fi%atoms%ntype))

            !remove spin dependence
            denIn1_pw  = (denIn1%pw(:,1)+denIn1%pw(:,fi%input%jspins))/(3.0-fi%input%jspins)
            denIn1_mt = (denIn1%mt(:,0:,:,1)+denIn1%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)
            denIn1_mt_Im = (denIn1Im%mt(:,0:,:,1)+denIn1Im%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)

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
            
            CALL dfpt_int_pw(starsq, fi%cell, denIn1_pw, pwwq2, tempval_pw)
            !print*,"I make it this far 1/2"
            !print*,"BEC_element",BEC_element
            !print*,"tempval_pw",tempval_pw
            BEC_element = BEC_element + tempval_pw 
            !print*,"I make it this far 2"
            do iType = 1, fi%atoms%ntype
                tempval_mt = CMPLX(0.0,0.0) 
                call dfpt_int_mt(fi%atoms, sphhar, fi%sym, iType, denIn1_mt, denIn1_mt_Im, vExt1%mt(:,0:,:,1), vExt1Im%mt(:,0:,:,1), tempval_mt)
                BEC_element =  BEC_element + tempval_mt
            end do
            !print*,"I make it this far 3"
            !surface integral
            if (iQ .eq. iDir_den) then
                call intgr3(denIn1_mt(:,0,iDType),fi%atoms%rmsh(1,iDType),fi%atoms%dx(iDType),fi%atoms%jri(iDType),int_mt_r)
                CALL intgr3(denIn1_mt_Im(:,0,iDType),fi%atoms%rmsh(1,iDType),fi%atoms%dx(iDType),fi%atoms%jri(iDType),int_mt_Im)
                BEC_element = BEC_element-(int_mt_r+int_mt_Im)
            end if 




            !print*,"BEC_row", BEC_row
            write(*,*) 'Determine integrals here'
            




        end subroutine dfpt_born_eff_charge_element

        subroutine dfpt_born_eff_charge_final(fi,born_eff_charge)

            type(t_fleurinput), intent(in)    :: fi
            complex, intent(inout)   :: born_eff_charge(:,:,:)
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
        end subroutine dfpt_born_eff_charge_final
        
end module 