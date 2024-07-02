module m_dfpt_dielecten
    use m_types
    use m_dfpt_vefield
    use m_convol
    use m_dfpt_dynmat

    implicit none

    contains

        subroutine dfpt_dielecten_row(fi,stars,starsq,sphhar,denIn1,denIn1Im,dieltensor_row)

            type(t_fleurinput), intent(in)     :: fi
            type(t_sphhar),    intent(in)      :: sphhar
            TYPE(t_stars),      INTENT(IN)     :: stars, starsq
            type(t_potden), intent(in)         :: denIn1,denIn1Im
            complex, intent(inout)             :: dieltensor_row(:)
            
            type(t_potden)                     :: vExt1, vExt1Im

            
            complex, allocatable               :: pwwq2(:),tempval_pw,tempval_mt, denIn1_pw(:), dieltensor_HF(:)
            real, allocatable                  :: denIn1_mt(:,:,:),denIn1_mt_Im(:,:,:)
            integer                            :: iDtype_col,iDir_col,col_index,iType


            !CALL vExt1%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
            !CALL vExt1Im%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)
            ALLOCATE(pwwq2(starsq%ng3))
            ALLOCATE(denIn1_pw(starsq%ng3))
            ALLOCATE(denIn1_mt(fi%atoms%jmtd,0:sphhar%nlhd,fi%atoms%ntype),denIn1_mt_Im(fi%atoms%jmtd,0:sphhar%nlhd,fi%atoms%ntype))
            ALLOCATE(dieltensor_HF(SIZE(dieltensor_row)))
            print*, "Im in dfpt_dielecten"
            print*, "dieltensor_row", dieltensor_HF
            denIn1_pw  = (denIn1%pw(:,1)+denIn1%pw(:,fi%input%jspins))/(3.0-fi%input%jspins)
            denIn1_mt = (denIn1%mt(:,0:,:,1)+denIn1%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)
            denIn1_mt_Im = (denIn1Im%mt(:,0:,:,1)+denIn1Im%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)
            do iDtype_col = 1, fi%atoms%ntype
                print*,"iDtype_col",iDtype_col 
                do iDir_col = 1, 3
                    print*,"iDir_col",iDir_col
                    
                    ! \rho(1)V_{ext}(1) integral (HF)
                    !interstitial
                    tempval_pw = CMPLX(0.0,0.0)
                    col_index = 3 * (iDtype_col - 1) + iDir_col
                    call vExt1%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
                    call vExt1Im%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)
                    call dfpt_vefield(fi%juPhon,fi%atoms,fi%sym,sphhar,fi%cell,vExt1,vExt1Im,iDir_col)
                    !print*,'shape Vext1', shape(vExt1%pw)
                    !print*,'shape Vext1im', shape(vExt1Im%pw)
                    ! IR integral:
                    pwwq2 = CMPLX(0.0,0.0)
                    CALL dfpt_convol_big(1, starsq, stars, vExt1%pw(:,1), CMPLX(1.0,0.0)*stars%ufft, pwwq2)
                    CALL dfpt_int_pw(starsq, fi%cell, denIn1_pw, pwwq2, tempval_pw)
                    print*, 'tempval_mt',tempval_pw
                    dieltensor_HF(col_index) = dieltensor_HF(col_index) + tempval_pw
                    print*, "dieltensor_row",dieltensor_row(:)


                    !Muffin-tin 
                    do iType = 1, fi%atoms%ntype
                        tempval_mt = CMPLX(0.0,0.0)               
                        call dfpt_int_mt(fi%atoms, sphhar, fi%sym, iType, denIn1_mt, denIn1_mt_Im, vExt1%mt(:,0:,:,1), vExt1Im%mt(:,0:,:,1), tempval_mt)
                        print*, 'tempval_mt',tempval_mt
                        dieltensor_HF(col_index) = dieltensor_HF(col_index) + tempval_mt
                    end do





                    !\sum f(1)\epsilon(1)

                    !interstitial

                    !Muffin-tin
                



                end do
            end do

            !stop
        end subroutine dfpt_dielecten_row
end module 