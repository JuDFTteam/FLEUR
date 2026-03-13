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

    subroutine dfpt_NAC(fi,dyn_mat_NAC_row)

    type(t_fleurinput), intent(in)     :: fi
    complex, intent(inout)             :: dyn_mat_NAC_row(:,:)
    real                    :: dielten(3,3),tempval_denom,tempval_BEC_fixed,tempval_NAC
    real, allocatable       :: borneffcharge(:,:,:)
    integer                 :: qnorm_vec(3), iDir,iDir1,iType,iType_row,iDir_row

    allocate(borneffcharge(fi%atoms%ntype,3,3))

    call read_BEC(fi,borneffcharge)
        
    qnorm_vec(:) = [1,0,0]
    tempval_denom = 0.0

    call read_dielten(fi,dielten)

    do iDir=1,3
        do iDir1=1,3
            tempval_denom = tempval_denom +qnorm_vec(iDir)*dielten(iDir,iDir1)*qnorm_vec(iDir1)
        end do
    end do

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



    end subroutine dfpt_NAC

    subroutine read_BEC(fi,borneffcharge)

        type(t_fleurinput), intent(in)     :: fi
        real, allocatable, intent(inout)   :: borneffcharge(:,:,:)
        integer                            :: iDType,iDir
        character(len=200)                 :: dummy_line

        borneffcharge(:,:,:)=0.0
        open(111, file="born_eff_charge", status='old', action='read', form='formatted')

        do iDType = 1, fi%atoms%ntype
            read(111,'(A)') dummy_line

            do iDir = 1, 3
                read(111,*) borneffcharge(iDType, iDir, 1:3)
            end do
        end do

        close(111)
    
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
    
    end subroutine read_dielten


end module 