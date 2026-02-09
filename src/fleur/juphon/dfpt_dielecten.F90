!--------------------------------------------------------------------------------
! Copyright (c) 2024 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_dfpt_dielecten
    use m_types
    use m_dfpt_vefield
    use m_convol
    use m_dfpt_dynmat
    use m_npy
    USE m_make_stars
    use m_inv3


    implicit none

contains

    subroutine dfpt_dielecten_HF_int(fi,stars,starsq,sphhar,fmpi,denIn1,denIn1Im,results,results1,dieltensor_row,rho,iDir_den,q_sign)


        
        type(t_fleurinput), intent(in)     :: fi
        type(t_sphhar),    intent(in)      :: sphhar
        TYPE(t_stars),      INTENT(IN)     :: stars, starsq
        type(t_potden), intent(inout)         :: denIn1,denIn1Im
        type(t_results), intent(in)        :: results, results1
        TYPE(t_mpi),        intent(in)     :: fmpi
        complex, intent(inout)             :: dieltensor_row(:)
        integer, intent(in)                :: q_sign
        type(t_potden), intent(inout)         :: rho



        type(t_potden)                     :: vExt1, vExt1Im
        type(t_stars)                      :: starsq_vext
        complex, allocatable               :: pwwq2(:),tempval_pw,tempval_mt, denIn1_pw(:), dieltensor_HF(:)
        real, allocatable                  :: denIn1_mt(:,:,:),denIn1_mt_Im(:,:,:) 
        integer                            :: iDir_col,iType
        integer                            :: iDir_den
        complex                            :: diel_tensor_int_IR(3),diel_tensor_int_MT(3)
        complex, allocatable                :: diel_tensor_int_MT_atom(:,:)
        character(len=20)                  :: status_string,densave_string,vextsave_string
        character(len=20)                   :: filename_string
        real                                         ::qvec_int(3)
        complex                            ::offset_out
                                
        
        ALLOCATE(pwwq2(starsq%ng3))
        ALLOCATE(denIn1_pw(starsq%ng3))
        ALLOCATE(denIn1_mt(fi%atoms%jmtd,0:sphhar%nlhd,fi%atoms%ntype),denIn1_mt_Im(fi%atoms%jmtd,0:sphhar%nlhd,fi%atoms%ntype))
        ALLOCATE(dieltensor_HF(SIZE(dieltensor_row)))
        allocate(diel_tensor_int_MT_atom(3,fi%atoms%ntype))
        dieltensor_HF(:) = CMPLX(0.0,0.0)

        !remove spin dependence
        denIn1_pw  = (denIn1%pw(:,1)+denIn1%pw(:,fi%input%jspins))/(3.0-fi%input%jspins)
        denIn1_mt = (denIn1%mt(:,0:,:,1)+denIn1%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)
        denIn1_mt_Im = (denIn1Im%mt(:,0:,:,1)+denIn1Im%mt(:,0:,:,fi%input%jspins))/(3.0-fi%input%jspins)

        ! \rho(1)V_{ext}(1) integral (HF)
        diel_tensor_int_IR = CMPLX(0.0,0.0)
        diel_tensor_int_MT = CMPLX(0.0,0.0)
        diel_tensor_int_MT_atom = cmplx(0.0,0.0)



        do iDir_col = 1, 3   
            ! call make_stars 
            tempval_pw = CMPLX(0.0,0.0)
            qvec_int=fi%juPhon%qvec_efield(:,iDir_col)
            CALL starsq_vext%reset_stars()
            CALL make_stars(starsq_vext, fi%sym, fi%atoms, fi%vacuum, sphhar, fi%input, fi%cell, fi%noco, fmpi, qvec_int, 1, iDir_col,fi%juPhon%l_efield)
            starsq_vext%ufft = starsq%ufft

            call vExt1%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
            call vExt1Im%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)
            call dfpt_vefield(fi%juPhon,starsq,fi%atoms,fi%sym,sphhar,fi%cell,vExt1,vExt1Im,iDir_col,1)


            if (iDir_col .eq. iDir_den) then
                !interstitial
                pwwq2 = CMPLX(0.0,0.0)
                
                CALL dfpt_convol_big(1, starsq , stars, vExt1%pw(:,1), CMPLX(1.0,0.0)*stars%ufft, pwwq2)!starsq
                call save_npy("pwwq2.npy",pwwq2)
                !CALL dfpt_int_pw(starsq, fi%cell, vExt1%pw(:,1), pwwq2, tempval_pw)!denIn1_pw
                CALL dfpt_int_pw(starsq, fi%cell, denIn1_pw, pwwq2, tempval_pw)!denIn1_pw
                
                dieltensor_HF(iDir_col) = dieltensor_HF(iDir_col) + tempval_pw 
                diel_tensor_int_IR(iDir_col) =  diel_tensor_int_IR(iDir_col) + tempval_pw


                !Muffin-tin 
                do iType = 1, fi%atoms%nat
                    tempval_mt = CMPLX(0.0,0.0) 
        
                    call dfpt_int_mt(fi%atoms, sphhar, fi%sym, iType, denIn1_mt, denIn1_mt_Im, vExt1%mt(:,0:,:,1), vExt1Im%mt(:,0:,:,1), tempval_mt)!denIn1_mt
                    dieltensor_HF(iDir_col) = dieltensor_HF(iDir_col) + tempval_mt
                    diel_tensor_int_MT(iDir_col) =  diel_tensor_int_MT(iDir_col) + tempval_mt
                    diel_tensor_int_MT_atom(iDir_col,iType) = diel_tensor_int_MT_atom(iDir_col,iType) + tempval_mt
                    
                end do
            end if 
        end do
        IF (fmpi%irank==0) THEN
            status_string = 'old'
            if (iDir_den == 1) status_string = 'replace'
            open(113, file="diel_tensor_int", status=status_string, position='append', action='write', form='formatted')
            write(113,'(6es16.8)') dieltensor_HF(:)
            close(113)
            open(111, file="diel_tensor_int_IR", status=status_string, position='append',action='write', form='formatted')
            write(111,'(6es16.8)') diel_tensor_int_IR(:)
            close(111)
            open(112, file="diel_tensor_int_MT", status=status_string, position='append', action='write', form='formatted')
            write(112,'(6es16.8)') diel_tensor_int_MT(:)
            close(112)
            if (fi%atoms%ntype >1) then
                do iType = 1,fi%atoms%ntype
                    write(filename_string, '(A,I0)') "diel_tensor_int_MT_", iType
                    open(112+itype, file=filename_string, status=status_string, position='append', action='write', form='formatted')
                    write(112+itype,'(6es16.8)') diel_tensor_int_MT_atom(:,iType)
                    close(112+itype)
                end do
            end if 
            END IF 


        dieltensor_row(:)= dieltensor_HF(:)
        
    end subroutine dfpt_dielecten_HF_int


    subroutine dfpt_dielecten_final_new(fi, dielecten)

        !USE m_inv3
        USE m_xmlOutput
        type(t_fleurinput), intent(in)    :: fi
        complex, intent(inout)   :: dielecten(:,:)
        integer                  :: iDir, j 
        complex                  :: dielten_iden(3,3) 
        complex                  :: dielecten_inv(3,3)
        real                     :: det,dielecten_r(3,3)
        CHARACTER(LEN=20)        :: attributes(1)


        dielecten_inv(:,:)= 0
        DO j = 1,3
            dielecten_inv(j,j) = CMPLX(1,0) + (4*fpi_const/fi%cell%omtil)*real(dielecten(j,j))
            END DO
        dielecten(:,:) = cmplx(0.0,0.0)
        call inv3(real(dielecten_inv),dielecten_r,det)
        dielecten = dielecten_r
        open( 110, file="diel_tensor", status='replace', action='write', form='formatted')
        write(*,*) '-------------------------' 
        write(*,*) "High Frequency Dielectric tensor" 
        do iDir = 1,3
            do j = 1,2
                write(110,'(2es16.8)', ADVANCE='NO') real(dielecten(iDir,j)) 
                write(110, '(A)', ADVANCE='NO') ' ' 
                write(*,'(2es16.8)', ADVANCE='NO') real(dielecten(iDir,j))
                write(*, '(A)', ADVANCE='NO') ' ' 
            end do
            write(110,'(2es16.8)')real(dielecten(iDir,3))
            write(*,'(2es16.8)')real(dielecten(iDir,3))
        end do
        close(110)
        write(*,*) '-------------------------' 

        !save in out.xml
        CALL openXMLElementNoAttributes('Phonons')
        CALL openXMLElementNoAttributes('efield')
        attributes = ''
        WRITE(attributes(1),'(f15.8)') fi%juPhon%qlim
        CALL writeXMLElementPoly('dielconst',(/'qlim'/), attributes,real(dielecten(:,1)))
        CALL closeXMLElement('efield')
        CALL closeXMLElement('Phonons')

    end subroutine dfpt_dielecten_final_new


    subroutine dfpt_dielecten_final_old(fi, dielecten)

        !USE m_inv3
        USE m_xmlOutput
        type(t_fleurinput), intent(in)    :: fi
        complex, intent(inout)   :: dielecten(:,:)
        integer                  :: iDir, j 
        complex                  :: dielten_iden(3,3) 
        complex                  :: dielecten_inv(3,3)
        real                     :: det,dielecten_r(3,3)
        CHARACTER(LEN=20)        :: attributes(1)

        dielecten_inv(:,:)= 0
        DO j = 1,3
            dielecten_inv(j,j) = CMPLX(1,0) - (4*fpi_const/fi%cell%omtil)*real(dielecten(j,j))
            END DO
        dielecten = dielecten_inv
        open( 110, file="diel_tensor_old", status='replace', action='write', form='formatted')
        write(*,*) '-------------------------' 
        write(*,*) "High Frequency Dielectric tensor (old formula)" 
        do iDir = 1,3
            do j = 1,2
                write(110,'(2es16.8)', ADVANCE='NO') real(dielecten(iDir,j)) 
                write(110, '(A)', ADVANCE='NO') ' ' 
                write(*,'(2es16.8)', ADVANCE='NO') real(dielecten(iDir,j))
                write(*, '(A)', ADVANCE='NO') ' ' 
            end do
            write(110,'(2es16.8)')real(dielecten(iDir,3))
            write(*,'(2es16.8)')real(dielecten(iDir,3))
        end do
        close(110)
        write(*,*) '-------------------------' 

    end subroutine dfpt_dielecten_final_old

end module 