!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_io_hybrid
  use m_io_matrix
  use m_judft
  use m_types
  implicit none
  !private
  integer,save :: id_olap,id_z,id_v_x,id_coulomb,id_coulomb_spm
  !public:: open_hybrid_io,read_cmt,write_cmt
contains

  SUBROUTINE open_hybrid_io1(l_real)
    implicit none

    LOGICAL,INTENT(IN)          :: l_real
    LOGICAL :: opened=.false.

    if (opened) return
    opened=.true.

    !print *,"Open olap.mat"
    id_olap=OPEN_MATRIX(l_real,lapw_dim_nbasfcn,1,1,"olap.mat")
    !print *,"Open z.mat"
    id_z=OPEN_MATRIX(l_real,lapw_dim_nbasfcn,1,1,"z.mat")
  END SUBROUTINE open_hybrid_io1


  SUBROUTINE open_hybrid_io1b(l_real)
    implicit none

    LOGICAL,INTENT(IN)          :: l_real
    LOGICAL :: opened=.false.

    if (opened) return
    opened=.true.

    !print *,"Open v_x.mat"
    id_v_x=OPEN_MATRIX(l_real,lapw_dim_nbasfcn,1,1,"v_x.mat")
  END SUBROUTINE open_hybrid_io1b


  SUBROUTINE open_hybrid_io2(mpdata,hybrid,input,atoms,l_real)
    IMPLICIT NONE
    type(t_mpdata), intent(in) :: mpdata
    TYPE(t_hybrid),INTENT(IN)   :: hybrid
    TYPE(t_input),INTENT(IN):: input
    TYPE(t_atoms),INTENT(IN)    :: atoms
    LOGICAL,INTENT(IN)          :: l_real
    INTEGER:: irecl_coulomb
    LOGICAL :: opened=.FALSE.



    if (opened) return
    opened=.true.
    OPEN(unit=777,file='cmt',form='unformatted',access='direct',&
         &     recl=input%neig*hybrid%maxlmindx*atoms%nat*16)

#ifdef CPP_NOSPMVEC
    irecl_coulomb = hybrid%maxbasm1 * (hybrid%maxbasm1+1) * 8 / 2
    if (.not.l_real) irecl_coulomb =irecl_coulomb *2
    OPEN(unit=778,file='coulomb',form='unformatted',access='direct', recl=irecl_coulomb)
    id_coulomb=778
#else
    ! if the sparse matrix technique is used, several entries of the
    ! matrix vanish so that the size of each entry is smaller
    irecl_coulomb = ( atoms%ntype*(maxval(hybrid%lcutm1)+1)*(maxval(mpdata%num_radbasfn)-1)**2&
         +   atoms%nat *(maxval(hybrid%lcutm1)+2)*(2*maxval(hybrid%lcutm1)+1)*(maxval(mpdata%num_radbasfn)-1)&
         +   (maxval(mpdata%num_radbasfn)-1)*atoms%nat**2&
         +   ((maxval(hybrid%lcutm1)+1)**2*atoms%nat+maxval(mpdata%n_g))&
         *((maxval(hybrid%lcutm1)+1)**2*atoms%nat+maxval(mpdata%n_g)+1)/2 )*8
    if (.not.l_real) irecl_coulomb =irecl_coulomb *2
    OPEN(unit=778,file='coulomb1',form='unformatted',access='direct', recl=irecl_coulomb)
    id_coulomb_spm=778
#endif
  END SUBROUTINE open_hybrid_io2

  subroutine write_cmt(cmt,nk)
    implicit none
    complex,INTENT(IN):: cmt(:,:,:)
    integer,INTENT(IN):: nk

    write(777,rec=nk) cmt
  end subroutine write_cmt

  subroutine read_cmt(cmt,nk)
    implicit none
    complex,INTENT(OUT):: cmt(:,:,:)
    integer,INTENT(IN):: nk

    read(777,rec=nk) cmt
  end subroutine read_cmt

  subroutine write_coulomb(nk,l_real,coulomb)
    implicit none
    complex,intent(in) :: coulomb(:)
    integer,intent(in) :: nk
    logical,intent(in) :: l_real

    if (l_real) THEN
       write(id_coulomb,rec=nk) real(coulomb)
    else
       write(id_coulomb,rec=nk) coulomb
    end if
  end subroutine write_coulomb

   subroutine write_coulomb_spm_r(nk,coulomb_mt1,coulomb_mt2,coulomb_mt3,coulomb_mtir)
    implicit none
    real,intent(in)    :: coulomb_mt1(:,:,:,:)
    real,intent(in) :: coulomb_mt2(:,:,:,:), coulomb_mt3(:,:,:)
    real,intent(in) :: coulomb_mtir(:)
    integer,intent(in) :: nk

    !print *, "write coulomb",nk,size(coulomb_mt1),size(coulomb_mt2),size(coulomb_mt3),size(coulomb_mtir)
    write(id_coulomb_spm,rec=nk) coulomb_mt1,coulomb_mt2,coulomb_mt3,coulomb_mtir
  end subroutine write_coulomb_spm_r

   subroutine write_coulomb_spm_c(nk,coulomb_mt1,coulomb_mt2,coulomb_mt3,coulomb_mtir)
    implicit none
    real,intent(in)    :: coulomb_mt1(:,:,:,:)
    complex,intent(in) :: coulomb_mt2(:,:,:,:), coulomb_mt3(:,:,:)
    complex,intent(in) :: coulomb_mtir(:)
    integer,intent(in) :: nk

    write(id_coulomb_spm,rec=nk) coulomb_mt1,coulomb_mt2,coulomb_mt3,coulomb_mtir
  end subroutine write_coulomb_spm_c

     subroutine read_coulomb_spm_r(nk,coulomb_mt1,coulomb_mt2,coulomb_mt3,coulomb_mtir)
    implicit none
    real,intent(out)    :: coulomb_mt1(:,:,:,:)
    real,intent(out) :: coulomb_mt2(:,:,:,:), coulomb_mt3(:,:,:)
    real,intent(out) :: coulomb_mtir(:)
    integer,intent(in) :: nk

    !print *, "read coulomb",nk,size(coulomb_mt1),size(coulomb_mt2),size(coulomb_mt3),size(coulomb_mtir)
    read(id_coulomb_spm,rec=nk) coulomb_mt1,coulomb_mt2,coulomb_mt3,coulomb_mtir
  end subroutine read_coulomb_spm_r

   subroutine read_coulomb_spm_c(nk,coulomb_mt1,coulomb_mt2,coulomb_mt3,coulomb_mtir)
    implicit none
    real,intent(out)    :: coulomb_mt1(:,:,:,:)
    complex,intent(out) :: coulomb_mt2(:,:,:,:), coulomb_mt3(:,:,:)
    complex,intent(out) :: coulomb_mtir(:)
    integer,intent(in) :: nk
    read(id_coulomb_spm,rec=nk) coulomb_mt1,coulomb_mt2,coulomb_mt3,coulomb_mtir
  end subroutine read_coulomb_spm_c

  subroutine read_coulomb_r(nk,coulomb)
    implicit none
    real   ,intent(out) :: coulomb(:)
    integer,intent(in) :: nk

    read(id_coulomb,rec=nk) coulomb
  end subroutine read_coulomb_r

  subroutine read_coulomb_c(nk,coulomb)
    implicit none
    complex,intent(out) :: coulomb(:)
    integer,intent(in) :: nk

    read(id_coulomb,rec=nk) coulomb
  end subroutine read_coulomb_c



  subroutine read_olap(mat,rec)
    implicit none
    TYPE(t_mat),INTENT(INOUT):: mat
    INTEGER,INTENT(IN)           :: rec

    CALL read_matrix(mat,rec,id_olap)
  END subroutine read_olap

    subroutine write_olap(mat,rec)
    implicit none
    TYPE(t_mat),INTENT(IN)   :: mat
    INTEGER,INTENT(IN)           :: rec

    CALL write_matrix(mat,rec,id_olap)
  END subroutine write_olap

  subroutine read_z(mat,rec)
    implicit none
    TYPE(t_mat),INTENT(INOUT):: mat
    INTEGER,INTENT(IN)           :: rec
    !print *,"read z:",rec

    CALL read_matrix(mat,rec,id_z)
  END subroutine read_z

    subroutine write_z(mat,rec)
    implicit none
    TYPE(t_mat),INTENT(IN)   :: mat
    INTEGER,INTENT(IN)           :: rec
     !print *,"write z:",rec
   CALL write_matrix(mat,rec,id_z)
  END subroutine write_z

  subroutine read_v_x(mat,rec)
    implicit none
    TYPE(t_mat),INTENT(INOUT):: mat
    INTEGER,INTENT(IN)           :: rec

    CALL read_matrix(mat,rec,id_v_x)
  END subroutine read_v_x

  subroutine write_v_x(mat,rec)
    implicit none
    TYPE(t_mat),INTENT(IN)   :: mat
    INTEGER,INTENT(IN)           :: rec

    CALL write_matrix(mat,rec,id_v_x)
  END subroutine write_v_x




end module m_io_hybrid
