!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_magma
  use m_juDFT
  INTEGER,PARAMETER :: NGPU_CONST=1
  !**********************************************************
  !     Solve the generalized eigenvalue problem
  !     using the MAGMA library for multiple GPUs
  !**********************************************************
CONTAINS
  SUBROUTINE magma_diag(nsize,eig,ne,a_r,b_r,z_r,a_c,b_c,z_c)
#ifdef CPP_MAGMA
    use magma
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
    
    LOGICAL :: l_real
    l_real=present(a_r)

    print *,"MAGMA start"
    IF (.NOT.initialized) THEN
       initialized=.true.
       call magmaf_init()
       print *,"MAGMA init"
    ENDIF

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

    !Query the workspace size 
    allocate(work(1),rwork(1),iwork(1))
    print *,"Magma workspace query"
    call flush()
    call magmaf_zhegvdx(1,'v','i','l',nsize,largea_c,nsize,largeb_c,nsize,&
         0.0,0.0,1,ne,mout,eigTemp,work,-1,rwork,-1,iwork,-1,err)
    lwork=work(1)
    lrwork=rwork(1)
    liwork=iwork(1)
    print*,"MAGMA:",lwork,lrwork,liwork
    deallocate(work,rwork,iwork)
    allocate(work(lwork),rwork(lrwork),iwork(liwork))
    if (err/=0) call juDFT_error("Failed to allocate workspaces",calledby="magma.F90")
    !Now the diagonalization
    print *,"Magma diagonalization"
    print *,nsize,shape(largea_c),shape(eigTemp),ne
    call magmaf_zhegvdx(1,'v','i','l',nsize,largea_c,nsize,largeb_c,nsize,&
         0.0,0.0,1,ne,mout,eigTemp,work,lwork,rwork,lrwork,iwork,liwork,err)
    print*,"MAGMA info:",err
    if (err/=0) call juDFT_error("Magma failed to diagonalize Hamiltonian")
    print *,"MAGMA mout:",mout

    DO i = 1, ne
       eig(i) = eigTemp(i)
       z_c(:nsize,i)=largea_c(:nsize,i)
    END DO
    !call judft_error("Eigenvectors are not calculated in MAGMA")
#endif
  END SUBROUTINE magma_diag
END MODULE m_magma

