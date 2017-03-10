!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_cusolver
  use m_juDFT
  INTEGER,PARAMETER   :: NGPU_CONST=1
  INTEGER,PARAMETER   :: CUSOLVER_EIG_TYPE_1=1
  CHARACTER,PARAMETER :: CUSOLVER_EIG_MODE_VECTOR='V'
  !**********************************************************
  !     Solve the generalized eigenvalue problem
  !     using the MAGMA library for multiple GPUs
  !**********************************************************
CONTAINS
  SUBROUTINE cusolver_diag(nsize,eig,ne,a_r,b_r,z_r,a_c,b_c,z_c)
    use m_packed_to_full
#ifdef CPP_MAGMA
    !use magma
#endif    
#include"cpp_double.h"
    IMPLICIT NONE

    ! ... Arguments ...

    INTEGER, INTENT (IN) :: nsize
  
    REAL,    INTENT(OUT) :: eig(:)
    INTEGER, INTENT(INOUT) :: ne

    REAL, OPTIONAL,ALLOCATABLE, INTENT (INOUT) :: a_r(:),b_r(:)
    REAL, OPTIONAL,ALLOCATABLE, INTENT (INOUT) :: z_r(:,:)
    COMPLEX, OPTIONAL,ALLOCATABLE, INTENT (INOUT) :: a_c(:),b_c(:)
    COMPLEX, OPTIONAL,ALLOCATABLE, INTENT (INOUT) :: z_c(:,:)

#ifdef CPP_MAGMA

    ! ... Local Variables ..
    INTEGER iind,ind1,ind2,info,lwork,liwork,lrwork,err,i,mout(1)
    REAL eigTemp(nsize)
    LOGICAL:: initialized=.false.

    REAL,    ALLOCATABLE :: rwork(:)
    INTEGER, ALLOCATABLE :: iwork(:)

    REAL, ALLOCATABLE :: largea_r(:,:),largeb_r(:,:)
    COMPLEX, ALLOCATABLE :: largea_c(:,:),largeb_c(:,:)
    COMPLEX,ALLOCATABLE :: work(:)

    INTEGER :: ifail(nsize)
    
    LOGICAL :: l_real
    l_real=present(a_r)

    !**********************************
    !expand from packed to full storage
    !**********************************
    !hamiltonian
    if (l_real) THEN
       call packed_to_full(nsize,a_r,largea_r)
       call packed_to_full(nsize,b_r,largeb_r)
       !deallocate(a_r,b_r)
    ELSE
       call packed_to_full(nsize,a_c,largea_c)
       call packed_to_full(nsize,b_c,largeb_c)
       !deallocate(a_c,b_c)
    Endif

    if (l_real) call juDFT_error("REAL diagonalization not implemented in magma.F90")

    err=cusolverDnCreate(handle)
    !$acc data copyin(largea_c,largeb_c)
    !$acc host_data use_device(largea_c,largeb_c,eigtemp,work)
    err=cusolverDnZhegvd( handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, cublas_Fill_Mode_upper, nsize, largea_c,nsize,largeb_c,nsize, eigtemp, work, lwork, devInfo);
    !$acc end host_data
    !$acc end data
    err=cusolverDnDestroy(handle)

    DO i = 1, ne
       eig(i) = eigTemp(i)
       z_c(:nsize,i)=largea_c(:nsize,i)
    END DO
    !call judft_error("Eigenvectors are not calculated in MAGMA")
#endif
  END SUBROUTINE magma_diag
END MODULE m_magma

module priv_cusolver_interface
contains
  integer(c_int) function cusolverDnCreate(cusolver_Hndl) bind(C,name="cusolverDnCreate")
     
    use iso_c_binding
    implicit none
     
    type(c_ptr)::cusolver_Hndl
     
  end function cusolverDnCreate

  integer(c_int) function cusolverDnDestroy(cusolver_Hndl) bind(C,name="cusolverDnDestroy")
     
    use iso_c_binding
    implicit none
     
    type(c_ptr),value::cusolver_Hndl
     
  end function cusolverDnDestroy

  integer(c_int) function cusolverDnZhegvd( handle, CUSOLVER_EIG_TYPE, CUSOLVER_EIG_MODE, cublas_Fill_Mode, nsize, largea_c,nsize,largeb_c,nsize, eigtemp, work, lwork, devInfo) bind(C,name="cusolverDnZhegvd") 
      
    use iso_c_binding
    implicit none
    
    type(c_ptr),value      :: handle
    integer(c_int),value   :: cusolver_eig_type
    character(c_char),value:: cusolver_eig_mode
    integer(c_int),value   :: cublas_fill_mode
    integer(c_int),value   ::n
    type(c_ptr),value::d_A
    integer(c_int),value::lda
    type(c_ptr),value::Lwork
  end function cusolverDnZhegvd

end module priv_cusolver_interface
