!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_geneigprobl
  use m_juDFT
  !**********************************************************
  !     Solve the generalized eigenvalue problem
  !     Frank Freimuth, November 2006
  !**********************************************************
CONTAINS
  SUBROUTINE geneigprobl(nbasfcn, nsize,neigd,l_J,a,b,z,eig,ne)
#include"cpp_double.h"
    IMPLICIT NONE

    ! ... Arguments ...

    INTEGER, INTENT (IN) :: nbasfcn
    INTEGER, INTENT (IN) :: neigd
    INTEGER, INTENT (IN) :: nsize
    LOGICAL, INTENT (IN) :: l_J

    REAL,    INTENT(OUT) :: eig(:)
    INTEGER, INTENT(OUT) :: ne
#ifdef CPP_F90

#ifdef CPP_INVERSION
    REAL,  INTENT (INOUT) :: a(:),b(:)
    REAL,  INTENT (INOUT) :: z(:,:)
#else
    COMPLEX, INTENT (INOUT)::a(:),b(:)
    COMPLEX, INTENT (INOUT) :: z(:,:)
#endif

#else

#ifdef CPP_INVERSION
    REAL, ALLOCATABLE, INTENT (INOUT) :: a(:),b(:)
    REAL, ALLOCATABLE, INTENT (INOUT) :: z(:,:)
#else
    COMPLEX, ALLOCATABLE, INTENT (INOUT) :: a(:),b(:)
    COMPLEX, ALLOCATABLE, INTENT (INOUT) :: z(:,:)
#endif

#endif

    ! ... Local Variables ..

    INTEGER iind,ind1,ind2,info,lwork,liwork,lrwork,err,i
    INTEGER sizez,iu 
    REAL :: lb,ub
    ! 'sizez' is needed, as some compilers sometimes produce errors,
    ! if the size command is used directly as a lapack argument.  

    REAL toler, eigTemp(nsize)

    REAL,    ALLOCATABLE :: work(:)
    INTEGER, ALLOCATABLE :: iwork(:),isuppz(:)
#ifdef CPP_INVERSION
    REAL, ALLOCATABLE :: largea(:,:),largeb(:,:)
#else
    COMPLEX, ALLOCATABLE :: largea(:,:),largeb(:,:)
    COMPLEX,ALLOCATABLE :: cwork(:)
#endif

    !**********************************
    !expand from packed to full storage: full storage lapack-routines
    !are faster than the packed lapack-routines.
    !**********************************
    !hamiltonian
    ALLOCATE ( largea(nsize,nsize), stat=err )
    IF (err/=0)  CALL juDFT_error("error allocating largea",calledby&
         &     ="geneigprobl")
    largea=0.0
    iind = 0
    DO ind1 = 1, nsize
       DO ind2 = 1, ind1
          iind = iind+1
          largea(ind2,ind1) = a(iind)
       ENDDO
    ENDDO
    !save some storage by deallocation of unused array
#ifndef CPP_F90
    DEALLOCATE (a)
#endif
    !metric
    ALLOCATE ( largeb(nsize,nsize), stat=err )
    IF (err/=0)  CALL juDFT_error("error allocating largeb",calledby ="geneigprobl")
    iind=0
    largeb=0.0
    DO ind1 = 1, nsize
       DO ind2 = 1, ind1
          iind = iind+1
          largeb(ind2,ind1) = b(iind)
       ENDDO
    ENDDO
    !save some storage by deallocation of unused array
#ifndef CPP_F90
    DEALLOCATE (b)
#endif



#ifdef CPP_INVERSION
    CALL CPP_LAPACK_spotrf('U',nsize,largeb,nsize,info)
    IF (info /= 0)  CALL juDFT_error("error in spotrf",calledby ="geneigprobl")

    CALL CPP_LAPACK_ssygst(1,'U',nsize,largea,nsize,largeb,nsize,info)
    IF (info /= 0)  CALL juDFT_error("error in ssygst",calledby ="geneigprobl")

    toler = 2.0*tiny(toler)
    liwork = 10*nsize
    ALLOCATE ( iwork(liwork), stat=err )
    IF (err/=0)  CALL juDFT_error("error allocating iwork",calledby ="geneigprobl")

    lwork = 26*nsize
    ALLOCATE ( work(lwork), stat=err )
    IF (err/=0)  CALL juDFT_error(" error allocating work",calledby ="geneigprobl")
    ALLOCATE ( isuppz(2*nsize), stat=err )
    IF (err /= 0)  CALL juDFT_error("error allocating isuppz",calledby ="geneigprobl")
#ifndef CPP_F90
    IF (allocated(z)) THEN
       IF (.not.(size(z,1)==nbasfcn.and.size(z,2)==neigd)) deallocate(z)
    ENDIF
    IF (.not.allocated(z)) THEN
       ALLOCATE ( z(nbasfcn,neigd), stat=err )
       IF (err/=0) THEN
          write(*,*) nbasfcn,neigd,err
          CALL juDFT_error("error allocating z",calledby ="geneigprobl")
       ENDIF
    ENDIF
#endif
    sizez= size(z,1) 
    iu   = min(nsize,neigd)
#ifndef CPP_F90
    IF (l_J) THEN
       CALL CPP_LAPACK_ssyevr('N','I','U',nsize,largea, nsize,lb,ub,1,iu,toler,ne,eigTemp,z,&
            sizez,isuppz,work,lwork,iwork,liwork,info)
    ELSE
       CALL CPP_LAPACK_ssyevr('V','I','U',nsize,largea,nsize,lb,ub,1,iu,toler,ne,eigTemp,z,&
            sizez,isuppz,work,lwork,iwork,liwork,info)
    ENDIF
#else
    eig = 0.0
    eigTemp = 0.0
#endif
    IF (info /= 0)  CALL juDFT_error("error in ssyevr",calledby ="geneigprobl")
    DEALLOCATE (isuppz,work,iwork)

    CALL CPP_LAPACK_strtrs('U','N','N',nsize,ne,largeb, nsize,z,sizez,info)
    IF (info /= 0)  CALL juDFT_error("error in strtrs",calledby&
         &     ="geneigprobl")
#else

    CALL CPP_LAPACK_cpotrf('U',nsize,largeb,nsize,info)
    IF (info /= 0)  CALL juDFT_error("error in cpotrf",calledby ="geneigprobl")

    CALL CPP_LAPACK_chegst(1,'U',nsize,largea,nsize,largeb,nsize,info)
    IF (info /= 0)  CALL juDFT_error(" error in chegst",calledby ="geneigprobl")

    toler = 2.0*tiny(toler)
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
#ifndef CPP_F90
    IF (allocated(z)) THEN
       IF (.not.(size(z,1)==nbasfcn.and.size(z,2)==neigd)) deallocate(z)
    ENDIF
    IF (.not.allocated(z)) THEN
       ALLOCATE ( z(nbasfcn,neigd), stat=err )
       IF (err/=0) THEN
          write(*,*) nbasfcn,neigd,err
          CALL juDFT_error("error allocating z",calledby ="geneigprobl")
       ENDIF
    ENDIF
#endif
    sizez= size(z,1) 
    iu   = min(nsize,neigd)
#ifndef CPP_F90
    IF (l_J) THEN
       CALL CPP_LAPACK_cheevr('N','I','U',nsize,largea, nsize,lb,ub,1,iu,toler,ne,eigTemp,z,&
            sizez,isuppz,cwork,lwork,work,lrwork,iwork,liwork,info)
    ELSE
#if (1==1)
       CALL CPP_LAPACK_cheevr('V','I','U',nsize,largea, nsize,lb,ub,1,iu,toler,ne,eigTemp,z,&
            sizez,isuppz,cwork,lwork,work,lrwork,iwork,liwork,info)
#else

       CALL CPP_LAPACK_cheevx('V','I','U',nsize,largea, nsize,lb,ub,1,iu,toler,ne,eigTemp,z,&
            sizez,cwork,lwork,work,iwork,isuppz,info)
#endif
    ENDIF
#else
    eig = 0.0
    eigTemp = 0.0
#endif
    IF (info /= 0)  CALL juDFT_error("error in cheevr",calledby ="geneigprobl")
    DEALLOCATE ( isuppz )
    deallocate ( work   )
    deallocate ( iwork  )
    deallocate ( cwork  )

    CALL CPP_LAPACK_ctrtrs('U','N','N',nsize,ne,largeb, nsize,z,sizez,info)
    IF (info /= 0)  CALL juDFT_error("error in ctrtrs",calledby ="geneigprobl")
#endif
    DEALLOCATE ( largea,largeb )

    DO i = 1, neigd
       eig(i) = eigTemp(i)
    END DO

  END SUBROUTINE geneigprobl
END MODULE m_geneigprobl

