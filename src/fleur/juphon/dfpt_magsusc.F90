!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_dfpt_magsusc
    use m_types
    use m_dfpt_vbfield
    use m_convol
    use m_dfpt_dynmat
    use m_npy
    USE m_make_stars
    use m_inv3
    USE m_checkdopall




    implicit none

contains

    subroutine dfpt_magnetic_susc(fi,stars,starsq,sphhar,fmpi,denIn1,denIn1Im,magnetic_susc)


        
        type(t_fleurinput), intent(in)     :: fi
        type(t_sphhar),    intent(in)      :: sphhar
        type(t_stars),     intent(in)      :: stars, starsq
        type(t_potden), intent(in)         :: denIn1,denIn1Im
        type(t_mpi),        intent(in)     :: fmpi
        complex, intent(inout)             :: magnetic_susc



        type(t_potden)                     :: vExt1, vExt1Im
        type(t_stars)                      :: starsq_vext
        complex, allocatable               :: pwwq2(:,:),tempval_pw,tempval_mt, denIn1_pw(:), dieltensor_HF(:)
        real, allocatable                  :: denIn1_mt(:,:,:),denIn1_mt_Im(:,:,:) 
        integer                            :: iDir_col,iType,iSpin
        complex                            :: diel_tensor_int_IR(3),diel_tensor_int_MT(3)
        complex, allocatable               :: diel_tensor_int_MT_atom(:,:)
        character(len=20)                  :: status_string,densave_string,vextsave_string
        character(len=20)                  :: filename_string
        real                               ::qvec_int(3)
        complex                            ::offset_out,magnetic_susc_final
 
        
        ALLOCATE(pwwq2(starsq%ng3,2))
        magnetic_susc = CMPLX(0.0,0.0)


        ! \rho(1)V_{ext}(1) integral (HF)
        diel_tensor_int_IR = CMPLX(0.0,0.0)
        diel_tensor_int_MT = CMPLX(0.0,0.0)
        diel_tensor_int_MT_atom = cmplx(0.0,0.0)


        do iSpin=1,2 
            do iDir_col = 1, 1   
                tempval_pw = CMPLX(0.0,0.0)
                call vExt1%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.TRUE.)
                call vExt1Im%init(starsq, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT, l_dfpt=.FALSE.)
                call dfpt_vbfield(fi%input,starsq,fi%noco,fi%atoms,vExt1,vExt1Im)
                call checkDOPAll(fi%input, sphhar, starsq,fi%atoms, fi%sym, fi%vacuum, fi%cell,vExt1,iSpin,vExt1Im) 
                !interstitial
                pwwq2 = CMPLX(0.0,0.0)

                CALL dfpt_convol_big(1, starsq , stars, vExt1%pw(:,iSpin), CMPLX(1.0,0.0)*stars%ufft, pwwq2(:,iSpin))
                CALL dfpt_int_pw(starsq, fi%cell, denIn1%pw(:,iSpin), pwwq2(:,iSpin), tempval_pw) 
                magnetic_susc =  magnetic_susc+ tempval_pw


                !Muffin-tin 
                do iType = 1, fi%atoms%nat
                    tempval_mt = CMPLX(0.0,0.0) 
                    call dfpt_int_mt(fi%atoms, sphhar, fi%sym, iType, denIn1%mt(:,:,:,iSpin), denIn1Im%mt(:,:,:,iSpin), vExt1%mt(:,0:,:,iSpin), -vExt1Im%mt(:,0:,:,iSpin), tempval_mt)!denIn1_mt ! same minus sign for potential due to mixed use of energy and spin sign
                    !print*,"tempval_mt",tempval_mt
                    magnetic_susc =  magnetic_susc+ tempval_mt
                    !print*,"magnetic_susc",magnetic_susc
                end do
            end do
        end do

        magnetic_susc = -1*magnetic_susc !technically there should be a factor 2 here
       ! IF (fmpi%irank==0) write(9989,*) "magnetic response", magnetic_susc
        IF (fmpi%irank==0) write(oUnit,*) "magnetic response (real): ", real(magnetic_susc), ", (imag): ", aimag(magnetic_susc)
        magnetic_susc_final =6.6918e-4*magnetic_susc/fi%cell%omtil
        IF (fmpi%irank==0) print*,"magnetic susceptibility in SI per vol:",real(magnetic_susc_final)
        IF (fmpi%irank==0) write(oUnit,*) "magnetic susceptibility in SI per vol: ", real(magnetic_susc_final), ", (imag): ", aimag(magnetic_susc_final)
        magnetic_susc_final =4.74891e-6*magnetic_susc/fi%atoms%nat
        IF (fmpi%irank==0) write(oUnit,*) "magnetic susceptibility in cgs per mol: ", real(magnetic_susc_final), ", (imag): ", aimag(magnetic_susc_final)
        IF (fmpi%irank==0) print*,"magnetic susceptibility in cgs per mol:",real(magnetic_susc_final)

    end subroutine dfpt_magnetic_susc


end module m_dfpt_magsusc