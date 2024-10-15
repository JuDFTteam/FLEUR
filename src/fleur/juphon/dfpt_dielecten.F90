module m_dfpt_dielecten
    use m_types
    use m_dfpt_vefield
    use m_convol
    use m_dfpt_dynmat
    use m_npy

    implicit none

    contains

        subroutine dfpt_dielecten_HF_int(fi,stars,starsq,sphhar,fmpi,denIn1,denIn1Im,results,results1,dieltensor_row)


            
            type(t_fleurinput), intent(in)     :: fi
            type(t_sphhar),    intent(in)      :: sphhar
            TYPE(t_stars),      INTENT(IN)     :: stars, starsq
            type(t_potden), intent(in)         :: denIn1,denIn1Im
            type(t_results), intent(in)        :: results, results1
            TYPE(t_mpi),        intent(in)     :: fmpi
            complex, intent(inout)             :: dieltensor_row(:)



            type(t_potden)                     :: vExt1, vExt1Im
            complex, allocatable               :: pwwq2(:),tempval_pw,tempval_mt, denIn1_pw(:), dieltensor_HF(:),dieltensor_occu(:),dieltensor_(:)
            real, allocatable                  :: denIn1_mt(:,:,:),denIn1_mt_Im(:,:,:) 
            integer                            :: iDir_col,iType
            
            ALLOCATE(pwwq2(starsq%ng3))
            ALLOCATE(denIn1_pw(starsq%ng3))
            ALLOCATE(denIn1_mt(fi%atoms%jmtd,0:sphhar%nlhd,fi%atoms%ntype),denIn1_mt_Im(fi%atoms%jmtd,0:sphhar%nlhd,fi%atoms%ntype))
            ALLOCATE(dieltensor_HF(SIZE(dieltensor_row)),dieltensor_occu(SIZE(dieltensor_row)))
            dieltensor_HF(:) = CMPLX(0.0,0.0)
            !print*, "Im in dfpt_dielecten"

            !remove spin dependence
            denIn1_pw  = (denIn1%pw(:,1)+denIn1%pw(:,fi%input%jspins))/(3.0-fi%input%jspins)
            denIn1_mt = (denIn1%mt(:,0:,:,1)+denIn1%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)
            denIn1_mt_Im = (denIn1Im%mt(:,0:,:,1)+denIn1Im%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)
            
            ! \rho(1)V_{ext}(1) integral (HF)
            
            do iDir_col = 1, 3   
                
                !interstitial
                tempval_pw = CMPLX(0.0,0.0)
                call vExt1%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
                call vExt1Im%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)
                call dfpt_vefield(fi%juPhon,fi%atoms,fi%sym,sphhar,fi%cell,vExt1,vExt1Im,iDir_col)

                ! IR integral:
                pwwq2 = CMPLX(0.0,0.0)
                CALL dfpt_convol_big(1, starsq, stars, vExt1%pw(:,1), CMPLX(1.0,0.0)*stars%ufft, pwwq2)
                CALL dfpt_int_pw(starsq, fi%cell, denIn1_pw, pwwq2, tempval_pw)
                dieltensor_HF(iDir_col) = dieltensor_HF(iDir_col) + tempval_pw


                !Muffin-tin 
                do iType = 1, fi%atoms%ntype
                    tempval_mt = CMPLX(0.0,0.0)               
                    call dfpt_int_mt(fi%atoms, sphhar, fi%sym, iType, denIn1_mt, denIn1_mt_Im, vExt1%mt(:,0:,:,1), vExt1Im%mt(:,0:,:,1), tempval_mt)
                    dieltensor_HF(iDir_col) = dieltensor_HF(iDir_col) + tempval_mt
                end do
            end do

            dieltensor_row(:)= dieltensor_HF(:)
            
        end subroutine dfpt_dielecten_HF_int

        subroutine dfpt_dielecten_final(fi, dielecten)

            type(t_fleurinput), intent(in)    :: fi
            complex, intent(inout)   :: dielecten(:,:)
            integer                  :: iDir, j 
            complex                  :: dielten_iden(3,3) 

            dielten_iden(:,:) =0
            DO j = 1,3
                dielten_iden(j,j) = CMPLX(1,0)
             END DO
            dielecten(:,:) = dielten_iden(:,:) - (fpi_const/fi%cell%omtil)*dielecten(:,:)
            open( 110, file="diel_tensor", status='replace', action='write', form='formatted')
            write(*,*) '-------------------------' 
            write(*,*) "Dielectric tensor" 
            do iDir = 1,3
               do j = 1,2
                  write(110,'(2es16.8)', ADVANCE='NO') dielecten(iDir,j) 
                  write(110, '(A)', ADVANCE='NO') ' ' 
                  write(*,'(2es16.8)', ADVANCE='NO') dielecten(iDir,j)
                  write(*, '(A)', ADVANCE='NO') ' ' 
               end do
               write(110,'(2es16.8)')dielecten(iDir,3)
               write(*,'(2es16.8)')dielecten(iDir,3)
            end do
            close(110)
            write(*,*) '-------------------------' 
        end subroutine dfpt_dielecten_final

end module 