!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_geneigprobl
  USE m_juDFT
  !**********************************************************
  !     Solve the generalized eigenvalue problem
  !     Frank Freimuth, November 2006
  !**********************************************************
CONTAINS
  SUBROUTINE geneigprobl(nbasfcn, nsize,neigd,l_J,eig,ne,a_r,b_r,z_r,a_c,b_c,z_c)
#include"cpp_double.h"
    USE m_packed_to_full
    IMPLICIT NONE

    ! ... Arguments ...

    INTEGER, INTENT (IN) :: nbasfcn
    INTEGER, INTENT (IN) :: neigd
    INTEGER, INTENT (IN) :: nsize
    LOGICAL, INTENT (IN) :: l_J

    REAL,    INTENT(OUT) :: eig(:)
    INTEGER, INTENT(OUT) :: ne

    REAL,OPTIONAL, ALLOCATABLE, INTENT (INOUT) :: a_r(:),b_r(:)
    REAL,OPTIONAL, ALLOCATABLE, INTENT (INOUT) :: z_r(:,:)
    COMPLEX,OPTIONAL, ALLOCATABLE, INTENT (INOUT) :: a_c(:),b_c(:)
    COMPLEX,OPTIONAL, ALLOCATABLE, INTENT (INOUT) :: z_c(:,:)


    ! ... Local Variables ..

    INTEGER iind,ind1,ind2,info,lwork,liwork,lrwork,err,i
    INTEGER sizez,iu 
    REAL :: lb,ub
    ! 'sizez' is needed, as some compilers sometimes produce errors,
    ! if the size command is used directly as a lapack argument.  

    REAL toler, eigTemp(nsize)

    REAL,    ALLOCATABLE :: work(:)
    INTEGER, ALLOCATABLE :: iwork(:),isuppz(:)
    REAL, ALLOCATABLE :: largea_r(:,:),largeb_r(:,:)
    COMPLEX, ALLOCATABLE :: largea_c(:,:),largeb_c(:,:)
    COMPLEX,ALLOCATABLE :: cwork(:)

    LOGICAL :: l_real

    l_real=PRESENT(a_r)


    !**********************************
    !expand from packed to full storage: full storage lapack-routines
    !are faster than the packed lapack-routines.
    !**********************************
    !hamiltonian
    IF (l_real) THEN
       call packed_to_full(nsize,a_r,largea_r)
       DEALLOCATE (a_r)
       call packed_to_full(nsize,b_r,largeb_r)
       DEALLOCATE (b_r)
    ELSE
       call packed_to_full(nsize,a_c,largea_c)
       DEALLOCATE (a_c)
       call packed_to_full(nsize,b_c,largeb_c)
       DEALLOCATE (b_c)
    ENDIF



    IF (l_real) then
    CALL CPP_LAPACK_spotrf('U',nsize,largeb_r,nsize,info)
    IF (info /= 0)  CALL juDFT_error("error in spotrf",calledby ="geneigprobl")

    CALL CPP_LAPACK_ssygst(1,'U',nsize,largea_r,nsize,largeb_r,nsize,info)
    IF (info /= 0)  CALL juDFT_error("error in ssygst",calledby ="geneigprobl")

    toler = 2.0*TINY(toler)
    liwork = 10*nsize
    ALLOCATE ( iwork(liwork), stat=err )
    IF (err/=0)  CALL juDFT_error("error allocating iwork",calledby ="geneigprobl")

    lwork = 26*nsize
    ALLOCATE ( work(lwork), stat=err )
    IF (err/=0)  CALL juDFT_error(" error allocating work",calledby ="geneigprobl")
    ALLOCATE ( isuppz(2*nsize), stat=err )
    IF (err /= 0)  CALL juDFT_error("error allocating isuppz",calledby ="geneigprobl")
    IF (ALLOCATED(z_r)) THEN
       IF (.NOT.(SIZE(z_r,1)==nbasfcn.AND.SIZE(z_r,2)==neigd)) DEALLOCATE(z_r)
    ENDIF
    IF (.NOT.ALLOCATED(z_r)) THEN
       ALLOCATE ( z_r(nbasfcn,neigd), stat=err )
       IF (err/=0) THEN
          WRITE(*,*) nbasfcn,neigd,err
          CALL juDFT_error("error allocating z",calledby ="geneigprobl")
       ENDIF
    ENDIF
    sizez= SIZE(z_r,1) 
    iu   = MIN(nsize,neigd)
    IF (l_J) THEN
       CALL CPP_LAPACK_ssyevr('N','I','U',nsize,largea_r, nsize,lb,ub,1,iu,toler,ne,eigTemp,z_r,&
            sizez,isuppz,work,lwork,iwork,liwork,info)
    ELSE
       CALL CPP_LAPACK_ssyevr('V','I','U',nsize,largea_r,nsize,lb,ub,1,iu,toler,ne,eigTemp,z_r,&
            sizez,isuppz,work,lwork,iwork,liwork,info)
    ENDIF
    IF (info /= 0)  CALL juDFT_error("error in ssyevr",calledby ="geneigprobl")
    DEALLOCATE (isuppz,work,iwork)

    CALL CPP_LAPACK_strtrs('U','N','N',nsize,ne,largeb_r, nsize,z_r,sizez,info)
    IF (info /= 0)  CALL juDFT_error("error in strtrs",calledby ="geneigprobl")
 ELSE

    CALL CPP_LAPACK_cpotrf('U',nsize,largeb_c,nsize,info)
    IF (info /= 0)  CALL juDFT_error("error in cpotrf",calledby ="geneigprobl")

    CALL CPP_LAPACK_chegst(1,'U',nsize,largea_c,nsize,largeb_c,nsize,info)
    IF (info /= 0)  CALL juDFT_error(" error in chegst",calledby ="geneigprobl")

    toler = 2.0*TINY(toler)
    liwork = 50*nsize
    ALLOCATE ( iwork(liwork), stat=err )
    IF (err/=0)  CALL juDFT_error("error allocating iwork",calledby ="geneigprobl")

    lwork = 20*nsize
    ALLOCATE( cwork(lwork), stat=err )
    IF (err/=0)  CALL juDFT_error("error allocating cwork",calledby ="geneigprobl")
    ALLOCATE( isuppz(10*nsize), stat=err )
    IF (err/=0)  CALL juDFT_error("error allocating isuppz",calledby ="geneigprobl")

    lrwork = 84*nsize
    ALLOCATE (work(lrwork), stat=err )
    IF (err/=0)  CALL juDFT_error(" error allocating work",calledby ="geneigprobl")
    IF (ALLOCATED(z_c)) THEN
       IF (.NOT.(SIZE(z_c,1)==nbasfcn.AND.SIZE(z_c,2)==neigd)) DEALLOCATE(z_c)
    ENDIF
    IF (.NOT.ALLOCATED(z_c)) THEN
       ALLOCATE ( z_c(nbasfcn,neigd), stat=err )
       IF (err/=0) THEN
          WRITE(*,*) nbasfcn,neigd,err
          CALL juDFT_error("error allocating z",calledby ="geneigprobl")
       ENDIF
    ENDIF
    sizez= SIZE(z_c,1) 
    iu   = MIN(nsize,neigd)
    IF (l_J) THEN
       CALL CPP_LAPACK_cheevr('N','I','U',nsize,largea_c, nsize,lb,ub,1,iu,toler,ne,eigTemp,z_c,&
            sizez,isuppz,cwork,lwork,work,lrwork,iwork,liwork,info)
    ELSE
       CALL CPP_LAPACK_cheevr('V','I','U',nsize,largea_c, nsize,lb,ub,1,iu,toler,ne,eigTemp,z_c,&
            sizez,isuppz,cwork,lwork,work,lrwork,iwork,liwork,info)
    ENDIF
    IF (info /= 0)  CALL juDFT_error("error in cheevr",calledby ="geneigprobl")
    DEALLOCATE ( isuppz )
    DEALLOCATE ( work   )
    DEALLOCATE ( iwork  )
    DEALLOCATE ( cwork  )

    CALL CPP_LAPACK_ctrtrs('U','N','N',nsize,ne,largeb_c, nsize,z_c,sizez,info)
    IF (info /= 0)  CALL juDFT_error("error in ctrtrs",calledby ="geneigprobl")
 ENDIF

 DO i = 1, neigd
    eig(i) = eigTemp(i)
 END DO

END SUBROUTINE geneigprobl
END MODULE m_geneigprobl

