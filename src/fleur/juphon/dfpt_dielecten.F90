module m_dfpt_dielecten
    use m_types
    use m_dfpt_vefield
    use m_convol
    use m_dfpt_dynmat

    implicit none

    contains

        subroutine dfpt_dielecten_row_HF(fi,stars,starsq,sphhar,fmpi,denIn1,denIn1Im,results,results1,no_row,dieltensor_row)


            
            type(t_fleurinput), intent(in)     :: fi
            type(t_sphhar),    intent(in)      :: sphhar
            TYPE(t_stars),      INTENT(IN)     :: stars, starsq
            type(t_potden), intent(in)         :: denIn1,denIn1Im
            type(t_results), intent(in)        :: results, results1
            TYPE(t_mpi),        intent(in)     :: fmpi
            complex, intent(inout)             :: dieltensor_row(:)
            integer, intent(in)                :: no_row
            
            type(t_potden)                     :: vExt1, vExt1Im

            
            complex, allocatable               :: pwwq2(:),tempval_pw,tempval_mt, denIn1_pw(:), dieltensor_HF(:),dieltensor_occu(:),dieltensor_(:)
            real, allocatable                  :: denIn1_mt(:,:,:),denIn1_mt_Im(:,:,:), we1(:),eig1(:), we1_data(:,:,:,:),eig1_data(:,:,:,:)
            integer                            :: iDtype_col,iDir_col,col_index,iType,jsp,nk_i, nk,len_kpoints


            !CALL vExt1%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
            !CALL vExt1Im%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)
            ALLOCATE(pwwq2(starsq%ng3))
            ALLOCATE(denIn1_pw(starsq%ng3))
            ALLOCATE(denIn1_mt(fi%atoms%jmtd,0:sphhar%nlhd,fi%atoms%ntype),denIn1_mt_Im(fi%atoms%jmtd,0:sphhar%nlhd,fi%atoms%ntype))
            ALLOCATE(dieltensor_HF(SIZE(dieltensor_row)),dieltensor_occu(SIZE(dieltensor_row)))
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
                end do
            end do 
            dieltensor_row(:)= dieltensor_HF(:)



            !\sum f(1)\epsilon(1)

            !interstitial

            !Muffin-tin
            print*,"shape(results1%w_iks(:,:,:))",shape(results%w_iks(:,:,:))
            print*,"shape(results1%eig(:,:,:))",shape(results1%eig(:,:,:))
            print*,"fi%input%neig",fi%input%neig
            !print*, "fmpi%k_list",shape(fmpi%)
            !len_kpoints = shape(fmpi%k_list,1)
            allocate(we1_data(fi%input%neig,size(fmpi%k_list), MERGE(1,fi%input%jspins,fi%noco%l_noco),3*fi%atoms%ntype))
            allocate(eig1_data(fi%input%neig,size(fmpi%k_list), MERGE(1,fi%input%jspins,fi%noco%l_noco),3*fi%atoms%ntype))
            !print*,"shape(we1_data)",shape(we1_data)
            !allocate(we1_data(,fmpi%k_list, MERGE(1,fi%input%jspins,fi%noco%l_noco),fi%atoms%ntype))
            we1_data(:,:,:,no_row) = results1%w_iks
            eig1_data(:,:,:,no_row) = results1%eig
            
            !print*,shape(we1_data)
            !stop
            !print*,(results1%eig(:,:,:))
            !print*,kind(results%w_iks(:,:,:))
            DO jsp = 1, MERGE(1,fi%input%jspins,fi%noco%l_noco)
                DO nk_i = 1,size(fmpi%k_list)
                    nk = fmpi%k_list(nk_i)
                    print*, "nk", nk
                    print*,"Doing the shit"
                    print*,size(fmpi%k_list)
                    !print*, shape(results%w_iks)
                    !print*, shape(results%eig)
                    !we  = results%w_iks(:,nk,jsp)
                    we1 = results1%w_iks(:,nk,jsp)
                    print*, "shape(we1)",shape(we1) 
                    !eig = results%eig(:,nk,jsp)
                    eig1 = results1%eig(:,nk,jsp)
                    print*,eig1
                    print*,"shape(eig1)",shape(eig1)
                    print*,kind(eig1)
                END DO
            END DO
            call save_npy("we1_data.npy",we1_data(:,:,:,:))
            !print*, "no_row",no_row
            !if (no_row == 3*fi%atoms%ntype) then
            !    print*,"end of rows"
            !    do 
            !        do jsp = 1, MERGE(1,fi%input%jspins,fi%noco%l_noco)
            !            do nk_i = 1,size(fmpi%k_list)
            !                nk = fmpi%k_list(nk_i)
                            


            !            end do
            !        end do
            !    end do 
            !end if 
            !print*,"we1_data(:,1,1,:)",we1_data(30,1100,1,:)

            
        end subroutine dfpt_dielecten_row_HF

        subroutine dfpt_dielecten_occ1(fi,fmpi,results1,we1_data,eig1_data,diel_tensor_occ1,no_row)

            type(t_fleurinput), intent(in)     :: fi
            type(t_results), intent(in)        :: results1
            type(t_mpi),        intent(in)     :: fmpi
            real, intent(inout)                ::  we1_data(:,:,:,:),eig1_data(:,:,:,:)
            real, intent(inout)                :: diel_tensor_occ1(:,:)
            integer, intent(in)                :: no_row
            integer                           :: ten_row, ten_col, jsp, nk_i,nk, nband
            real                             :: temp_val


            diel_tensor_occ1(:,:) = 0.0

            we1_data(:,:,:,no_row) = results1%w_iks
            eig1_data(:,:,:,no_row) = results1%eig
            if (no_row == 3*fi%atoms%ntype) then
                print*,"end of rows"
                do ten_row = 1,3*fi%atoms%ntype
                    do ten_col = 1,3*fi%atoms%ntype
                        temp_val = 0.0
                        do jsp = 1, merge(1,fi%input%jspins,fi%noco%l_noco)
                            do nk_i = 1, size(fmpi%k_list)
                                nk = fmpi%k_list(nk_i)
                                do nband = 1,fi%input%neig
                                    temp_val = temp_val + we1_data(nband,nk,jsp,ten_row)*eig1_data(nband,nk,jsp,ten_col)
                                    print*,'temp_val',temp_val
                                end do
                            end do
                        end do
                        diel_tensor_occ1(ten_row,ten_col) = temp_val
                    end do
                end do
                print*,"diel_tensor_occ1",diel_tensor_occ1(:,:)
                print*,"test"
                stop
            end if 
        end subroutine dfpt_dielecten_occ1

end module 