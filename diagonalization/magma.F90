MODULE m_magma
  use m_juDFT
  INTEGER,PARAMETER :: NGPU_CONST=1
  !**********************************************************
  !     Solve the generalized eigenvalue problem
  !     using the MAGMA library for multiple GPUs
  !**********************************************************
CONTAINS
  SUBROUTINE magma_diag(nsize,a,b,z,eig,ne)
#ifdef CPP_MAGMA
    use magma
#endif    
#include"cpp_double.h"
    IMPLICIT NONE

    ! ... Arguments ...

    INTEGER, INTENT (IN) :: nsize
  
    REAL,    INTENT(OUT) :: eig(:)
    INTEGER, INTENT(INOUT) :: ne
#ifdef CPP_INVERSION
    REAL, ALLOCATABLE, INTENT (INOUT) :: a(:),b(:)
    REAL, ALLOCATABLE, INTENT (INOUT) :: z(:,:)
#else
    COMPLEX, ALLOCATABLE, INTENT (INOUT) :: a(:),b(:)
    COMPLEX, ALLOCATABLE, INTENT (INOUT) :: z(:,:)
#endif

#ifdef CPP_MAGMA

    ! ... Local Variables ..
    INTEGER iind,ind1,ind2,info,lwork,liwork,lrwork,err,i,mout(1)
    REAL eigTemp(nsize)
    LOGICAL:: initialized=.false.

    REAL,    ALLOCATABLE :: rwork(:)
    INTEGER, ALLOCATABLE :: iwork(:)
#ifdef CPP_INVERSION
    REAL, ALLOCATABLE :: largea(:,:),largeb(:,:)
#else
    COMPLEX, ALLOCATABLE :: largea(:,:),largeb(:,:)
    COMPLEX,ALLOCATABLE :: work(:)
#endif

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
    ALLOCATE ( largea(nsize,nsize), stat=err )
    IF (err/=0)  CALL juDFT_error("error allocating largea",calledby="geneigprobl")
    iind = 0
    DO ind1 = 1, nsize
       DO ind2 = 1, ind1
          iind = iind+1
          largea(ind2,ind1) = a(iind)
       ENDDO
    ENDDO
    !save some storage by deallocation of unused array
    !DEALLOCATE (a)
    !metric
    ALLOCATE ( largeb(nsize,nsize), stat=err )
    IF (err/=0)  CALL juDFT_error("error allocating largeb",calledby ="geneigprobl")
    iind=0
    DO ind1 = 1, nsize
       DO ind2 = 1, ind1
          iind = iind+1
          largeb(ind2,ind1) = b(iind)
       ENDDO
    ENDDO
    !save some storage by deallocation of unused array
    !DEALLOCATE (b)

#ifdef CPP_INVERSION
    call juDFT_error("REAL diagonalization not implemented in magma.F90")
#else
    !Query the workspace size 
    allocate(work(1),rwork(1),iwork(1))
    call magmaf_zhegvdx_2stage_m(NGPU_CONST,1,MagmaVec,MagmaRangeI,MagmaLower,nsize,largea,nsize,largeb,nsize,&
         0.0,0.0,1,ne,mout,eigTemp,work,-1,rwork,-1,iwork,-1,err)
    lwork=work(1)
    lrwork=rwork(1)
    liwork=iwork(1)
    print*,"MAGMA:",lwork,lrwork,liwork
    deallocate(work,rwork,iwork)
    allocate(work(lwork),rwork(lrwork),iwork(liwork))
    if (err/=0) call juDFT_error("Failed to allocate workspaces",calledby="magma.F90")
    !Now the diagonalization
    call magmaf_zhegvdx_2stage_m(NGPU_CONST,1,MagmaVec,MagmaRangeI,MagmaLower,nsize,largea,nsize,largeb,nsize,&
         0.0,0.0,1,ne,mout,eigTemp,work,lwork,rwork,lrwork,iwork,liwork,err)
    print*,"MAGMA info:",err
    if (err/=0) call juDFT_error("Magma failed to diagonalize Hamiltonian")
    print *,"MAGMA mout:",mout
#endif

    DO i = 1, ne
       eig(i) = eigTemp(i)
    END DO
#endif
  END SUBROUTINE magma_diag
END MODULE m_magma

