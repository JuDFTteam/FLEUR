!--------------------------------------------------------------------------------
! Copyright (c) 2024 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_dfpt_NAC
    
    USE m_juDFT
    USE m_constants
    USE m_types

    implicit none

contains

    subroutine dfpt_NAC(fi,iDtype,iDir_den,dyn_mat_NAC_row)

    type(t_fleurinput), intent(in)     :: fi
    integer, intent(in)                :: iDtype,iDir_den 
    complex, intent(inout)             :: dyn_mat_NAC_row(:,:)
    real                    :: borneffcharge_dumb(2,3,3),dielten(3,3),dielten_dumb(3,3),tempval_denom,tempval_BEC_fixed,tempval_NAC
    real, allocatable       :: borneffcharge(:,:,:)
    integer                 :: qnorm_vec(3), iDir,iDir1,iType,iType_row,iDir_row

    borneffcharge_dumb(:,:,:) =0.0
    dielten_dumb(:,:) =0.0
    !get data for SiC:
    do iDir=1,3
        dielten_dumb(iDir,iDir) = 7.053
        borneffcharge_dumb(1,iDir,iDir)=2.738
        borneffcharge_dumb(2,iDir,iDir)=-2.738
    end do
    allocate(borneffcharge(fi%atoms%ntype,3,3))
    !print*,"dielten",dielten
    !print*,"borneffcharge1",borneffcharge(1,:,:)
    !print*,"borneffcharge2",borneffcharge(2,:,:)
    call read_BEC(fi,borneffcharge)
    !borneffcharge= borneffcharge_dumb
    qnorm_vec(:) = [1,0,0]
    tempval_denom = 0.0
    !print*,"diff", 7.053*borneffcharge-borneffcharge_dumb
    call read_dielten(fi,dielten)
    !call correct_BEC(dielten,borneffcharge)

    !print*,"borneffcharge1",borneffcharge(1,:,:)
    !print*,"borneffcharge2",borneffcharge(2,:,:)
    !print*,"dielten",dielten(:,:)
    !stop

    !dielten_dumb = dielten
    !print*,"diff", dielten-dielten_dumb


    !stop
    do iDir=1,3
        do iDir1=1,3
            tempval_denom = tempval_denom +qnorm_vec(iDir)*dielten(iDir,iDir1)*qnorm_vec(iDir1)
            !Isn't this actually completely redundant?
        end do
    end do
    !print*,"tempval_denom",tempval_denom
    !tempval_BEC_fixed = sum(borneffcharge(iDtype,iDir_den,:)*qnorm_vec(:))
    do itype=1,fi%atoms%ntype
        do iDir=1,3
            do iType_row=1,fi%atoms%ntype
                do iDir_row=1,3
                    tempval_NAC=0.0
                    tempval_NAC = tempval_NAC + sum(borneffcharge(itype,iDir,:)*qnorm_vec(:))*sum(borneffcharge(iType_row,iDir_row,:)*qnorm_vec(:))/tempval_denom
                    dyn_mat_NAC_row(3 *(itype-1)+iDir,3 *(iType_row-1)+iDir_row) = (fpi_const/fi%cell%omtil)*tempval_NAC*cmplx(1.,0.)
                end do
            end do
        end do
    end do
    !dyn_mat_NAC_row(:) =0.0

    !stop


    end subroutine dfpt_NAC

    subroutine read_BEC(fi,borneffcharge)

        type(t_fleurinput), intent(in)     :: fi
        real, allocatable, intent(inout)   :: borneffcharge(:,:,:)
        integer                            :: iDType,iDir
        character(len=200)                 :: dummy_line

        borneffcharge(:,:,:)=0.0
        open(111, file="born_eff_charge", status='old', action='read', form='formatted')

        do iDType = 1, fi%atoms%ntype
            ! Skip header line completely
            read(111,'(A)') dummy_line
            ! Read 3x3 block
            do iDir = 1, 3
                read(111,*) borneffcharge(iDType, iDir, 1:3)
                !print*,borneffcharge(iDType, iDir, 1:3)
            end do
        end do

        close(111)
        !print*,"borneffcharge",borneffcharge(1,:,:)
    
    end subroutine read_BEC

    subroutine read_dielten(fi,dielten)

        type(t_fleurinput), intent(in)     :: fi
        real, intent(inout)   :: dielten(:,:)
        integer                            :: iDir

        dielten(:,:)=0.0
        open(110, file="diel_tensor", status='old', action='read', form='formatted')

        do iDir = 1, 3
            read(110,*) dielten( iDir, 1:3)
            print*,dielten( iDir, 1:3)
        end do

        close(110)
        !print*,"dielten",dielten(:,:)
    
    end subroutine read_dielten

    subroutine correct_BEC(dielten,borneffcharge)

        !type(t_fleurinput), intent(in)     :: fi
        real, intent(in)   :: dielten(:,:)
        real, intent(inout)  :: borneffcharge(:,:,:)
        integer                            :: iDir

        !dielten(:,:)=0.0
        !open(110, file="diel_tensor", status='old', action='read', form='formatted')

        do iDir = 1, 3
            borneffcharge(:,:,iDir) = borneffcharge(:,:,iDir)*dielten(iDir,iDir) !last factor technically not correct --> discuss this
            !read(110,*) dielten( iDir, 1:3)
            !print*,dielten( iDir, 1:3)
        end do

        !close(110)
        !print*,"dielten",dielten(:,:)
    
    end subroutine correct_BEC

end module 