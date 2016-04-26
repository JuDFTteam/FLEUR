MODULE m_franza
CONTAINS
  SUBROUTINE franza(nbasfcn,neigd, nsize,&
       l_zref,l_J,matind,nred, a,b,gw, z,eig,ne)
    !***********************************************************************
    !
    !     solves the secular equation a*z=eig*b*z
    !
    !     nbase,nval  array dimensions as declared by calling routin
    !                         unchanged on exit
    !     nsize   actual dimension of a, b
    !     a       on entry:   hamiltonian, lower triangle stored row-wise
    !                              = upper triangle stored column-wise 
    !             on exit :   destroyed
    !     b       on entry:   overlap matrix, lower tr. stored row-wise
    !             on exit :   destroyed(overwritten by cholesky factor)
    !     eig     on entry:   ----
    !             on exit :   eigenvalues in ascending rder
    !     z       on exit :   corresponding eigenvectors
    !     ne      on exit :   number of eigenvalues/vectors actually found
    !     nblw    on exit :   number of eigenvalues below e(1)
    !     work,iwork          temporary storage of dimensions
    !                            nrtmp .ge. 8*nv
    !                            nitmp .ge. 5*nv
    !
    !  This subroutine was written by Gustav Bihlmayer 
    !                                 IFF, July 1996
    !
    !***********************************************************************
    ! lapack routines inserted for the cray
    ! please note, that if too many ev's are in your interval, you'll get 
    ! a memory fault. Then change dspevx('V','V' ... to dspevx('V','A'
    ! to see what's going on .... 
    !***********************************************************************
    USE m_juDFT
#include"cpp_double.h"
    IMPLICIT NONE
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN)    :: nbasfcn,neigd
    INTEGER, INTENT (INOUT) :: nsize
    INTEGER, INTENT (OUT)   :: ne
    INTEGER, INTENT (IN)    :: nred,gw
    LOGICAL, INTENT (IN)    :: l_zref,l_J
    !     ..
    !     .. Array Arguments ..
    INTEGER, INTENT (IN)  :: matind(nbasfcn,2)
    REAL,    INTENT (OUT) :: eig(neigd)
#ifdef CPP_INVERSION
    REAL,    INTENT (INOUT):: a(:),b(:)!(matsize)
    REAL,    INTENT (OUT) :: z(nbasfcn,neigd)
#else
    COMPLEX, INTENT (INOUT):: a(:),b(:)
    COMPLEX, INTENT (OUT) :: z(nbasfcn,neigd)
#endif
    !     ..
    !     .. Local Scalars ..
    REAL toler,sq2i
    INTEGER j,i1,i2,j1,j2,k,nsym,jsym,matsz
    INTEGER info,i,iu,ne_a,nsize_a,n
    LOGICAL addstate
    !     ..
    !     .. Local Arrays
#ifdef CPP_INVERSION
    REAL,    ALLOCATABLE :: aa(:), bb(:)
    REAL                 :: zz(nbasfcn,neigd+1)
#else
    COMPLEX, ALLOCATABLE :: aa(:), bb(:), cwork(:)
    COMPLEX              :: zz(nbasfcn,neigd+1)
#endif
    REAL,    ALLOCATABLE ::  work(:), etemp(:)
    INTEGER, ALLOCATABLE ::  iwork(:),ifail(:)
    LOGICAL sort(nbasfcn)
    REAL:: lb,ub

    !     ..
    !     ..
    !     .. External Subroutines ..
#ifdef CPP_INVERSION
    EXTERNAL CPP_LAPACK_spptrf,CPP_LAPACK_sspgst,CPP_LAPACK_sspevx, CPP_LAPACK_stptrs,CPP_BLAS_scopy
#else
    EXTERNAL CPP_LAPACK_cpptrf,CPP_LAPACK_chpgst,CPP_LAPACK_chpevx, CPP_LAPACK_ctptrs,CPP_BLAS_ccopy
#endif
    !     ..
    nsym=1
    IF (l_zref) THEN
       !
       ALLOCATE (aa(nred*(nred+1)/2),bb(nred*(nred+1)/2))
       !
       ! separate H and S matrix in symmetric (aa,bb) and antisymmetric (a,b) part
       !
       nsym=2
       i2=0
       j2=0
       DO i=1,nred
          DO j=1,i 
             i1=(matind(i,1)-1)*matind(i,1)/2+matind(j,1)
             j1=(matind(i,1)-1)*matind(i,1)/2+matind(j,2)
             i2=i2+1
             aa(i2)=a(i1)+a(j1)
             bb(i2)=b(i1)+b(j1)
             IF ((matind(i,1).NE.matind(i,2)).AND.(matind(j,1).NE.matind(j,2))) THEN
                j2=j2+1
                a(j2)=a(i1)-a(j1)
                b(j2)=b(i1)-b(j1)
             ENDIF
          ENDDO
       ENDDO
       nsize = nsize - nred
    ENDIF
    DO jsym=1,nsym

       IF (jsym.EQ.2) THEN
          !
          ! second time calculate symmetric part and store the antisym. EV's and EW's
          !
          matsz=(nred+1)*nred/2
#ifdef CPP_INVERSION
          CALL CPP_BLAS_scopy(matsz,aa,1,a,1)
          CALL CPP_BLAS_scopy(matsz,bb,1,b,1)
#else
          CALL CPP_BLAS_ccopy(matsz,aa,1,a,1)
          CALL CPP_BLAS_ccopy(matsz,bb,1,b,1)
#endif

          ne_a=ne
          k=1
          nsize_a=nsize
          DO i=1,ne
#ifdef CPP_INVERSION
             aa(i)=eig(i)
#else
             aa(i)=CMPLX(eig(i),0.0)
#endif
             DO j=1,nsize
                bb(k)=z(j,i)
                k=k+1
             ENDDO
          ENDDO
          nsize=nred

       ENDIF
       !-gu
       ! --> start with Cholesky factorization of b ( so that b = l * l^t)
       ! --> b is overwritten by l
       !
#ifdef CPP_INVERSION
       CALL CPP_LAPACK_spptrf('U',nsize,b,info)
#else
       CALL CPP_LAPACK_cpptrf('U',nsize,b,info)
#endif
       !
       IF (info.NE.0) THEN
#ifdef CPP_INVERSION
          WRITE (*,*) 'Error in spptrf: info =',info
#else
          WRITE (*,*) 'Error in cpptrf: info =',info
#endif
          CALL juDFT_error("Diagonalization failed",calledby="franza")
       ENDIF
       !
       ! --> now reduce a * z = eig * b * z to the standard form a' * z' = eig * z' 
       ! --> where a' = (l)^-1 * a * (l^t)^-1 and z' = l^t * z
       !
#ifdef CPP_INVERSION
       CALL CPP_LAPACK_sspgst(1,'U',nsize,a,b,info)
#else
       CALL CPP_LAPACK_chpgst(1,'U',nsize,a,b,info)
#endif
       !
       IF (info.NE.0) THEN
#ifdef CPP_INVERSION
          WRITE (6,*) 'Error in sspgst: info =',info
#else
          WRITE (6,*) 'Error in chpgst: info =',info
#endif
          CALL juDFT_error("Diagonalization failed",calledby="franza")
       ENDIF
       !
       ! --> solve a' * z' = eig * z' for eigenvalues eig between lb und ub
       !
       iu = MIN(neigd,nsize)
       IF (l_zref) iu = (neigd-1)/2
       !      toler=2.0*CPP_LAPACK_slamch('S')
       toler=2.0*TINY(sq2i)
       ALLOCATE ( work(8*nbasfcn),iwork(5*nbasfcn),ifail(nbasfcn) )
       ALLOCATE ( etemp(nbasfcn) )
       addstate = gw.NE.0.AND.iu.LT.nsize ! add one state, 
       IF(addstate) iu = iu + 1           ! see below (CF)
#ifdef CPP_INVERSION
       IF(l_J)THEN
          CALL CPP_LAPACK_sspevx('N','I','U',nsize,a,lb,ub,1,iu,toler,ne, etemp,zz,nbasfcn,work,iwork,ifail,info)
       ELSE
          CALL CPP_LAPACK_sspevx('V','I','U',nsize,a,lb,ub,1,iu,toler,ne, etemp,zz,nbasfcn,work,iwork,ifail,info)
       ENDIF
#else
       ALLOCATE ( cwork(2*nbasfcn) )
       IF(l_J)THEN
          CALL CPP_LAPACK_chpevx('N','I','U',nsize,a,lb,ub,1,iu,toler,ne, etemp,zz,nbasfcn,cwork,work,iwork,ifail,info)
       ELSE
          CALL CPP_LAPACK_chpevx('V','I','U',nsize,a,lb,ub,1,iu,toler,ne, etemp,zz,nbasfcn,cwork,work,iwork,ifail,info)
       ENDIF
       DEALLOCATE ( cwork )
#endif
       IF(addstate) THEN ! cut topmost subspace of degenerate states to avoid symmetry breaking (CF)
          iu = ne
          ne = ne - 1
          DO WHILE(ABS(etemp(ne)-etemp(iu)).LT.1d-8)
             ne = ne - 1
          ENDDO
       ENDIF
       z = zz(:,:ne) 
       !
       IF (ne.GT.neigd) THEN
          WRITE(6,*) 'ne=',ne,' > neigd'
          CALL juDFT_error("ne.GT.neigd",calledby="franza")
       ENDIF
       eig(:ne) = etemp(:ne)
       !
       IF(l_J) THEN
          DEALLOCATE ( work,iwork,ifail,etemp )
       ELSE
          IF (info.NE.0) THEN
#ifdef CPP_INVERSION
             WRITE (6,*) 'Error in sspevx: info =',info
#else
             WRITE (6,*) 'Error in chpevx: info =',info
#endif
             WRITE (6,*) 'The following eigenvectors did not converge:'
             WRITE (6,'(30i5)') (ifail(i),i=1,ne)
             CALL juDFT_error("Diagonalization failed",calledby="franza")
          ENDIF
          DEALLOCATE ( work,iwork,ifail,etemp )
          !
          ! --> recover the generalized eigenvectors z by solving z' = l^t * z
          !
#ifdef CPP_INVERSION
          CALL CPP_LAPACK_stptrs('U','N','N',nsize,ne,b,z,nbasfcn,info)
#else
          CALL CPP_LAPACK_ctptrs('U','N','N',nsize,ne,b,z,nbasfcn,info)
#endif
          !
          IF (info.NE.0) THEN
#ifdef CPP_INVERSION
             WRITE (6,*) 'Error in stptrs: info =',info
#else
             WRITE (6,*) 'Error in ctptrs: info =',info
#endif
             CALL juDFT_error("Diagonalization failed",calledby="franza")
          ENDIF
       ENDIF !l_J
       !+gu
    ENDDO
    !
    ! now collect symmetric and antisym. EW's and EV's and sort
    !     
    IF (l_zref) THEN

       k=1
       DO i=1,ne
          DO j=1,nsize
             b(k)=z(j,i)
             k=k+1
          ENDDO
       ENDDO
       !
       ! prepare sort-array: even=.true., odd=.false.
       !
       i=1
       j=1
       eig(ne+1)    = 99.9e9
#ifdef CPP_INVERSION
       aa(ne_a+1) = 99.9e9
#else
       aa(ne_a+1) = CMPLX(99.9e9,0.0)
#endif
       DO k=1,ne+ne_a
          IF (eig(i).LT.REAL(aa(j))) THEN
             sort(k)=.TRUE.
             i=i+1
          ELSE
             sort(k)=.FALSE.
             j=j+1
          ENDIF
       ENDDO
       !
       ! sort EW's and EV's
       !
       i=ne
       j=ne_a
       nsize = nsize + nsize_a
       ne = ne + ne_a
       sq2i=1.0/SQRT(2.0)
       DO k=ne,1,-1
          DO n=1,nsize
#ifdef CPP_INVERSION
             z(n,k)=0.0
#else
             z(n,k)=CMPLX(0.0,0.0)
#endif
          ENDDO
          IF (sort(k)) THEN
             eig(k)=eig(i)
             i1=nred * (i-1)
             DO n=1,nred
                i1=i1+1
#ifdef CPP_INVERSION
                z(matind(n,1),k) = z(matind(n,1),k)+b(i1)*sq2i
                z(matind(n,2),k) = z(matind(n,2),k)+b(i1)*sq2i
#else
                z(matind(n,1),k) = z(matind(n,1),k)+b(i1)*CMPLX(sq2i,0.0)
                z(matind(n,2),k) = z(matind(n,2),k)+b(i1)*CMPLX(sq2i,0.0)
#endif
             ENDDO
             i=i-1
          ELSE
#ifdef CPP_INVERSION
             eig(k)=aa(j)
#else
             eig(k)=REAL(aa(j))
#endif
             j1=nsize_a * (j-1)
             DO n=1,nred
                IF (matind(n,1).NE.matind(n,2)) THEN
                   j1=j1+1
#ifdef CPP_INVERSION
                   z(matind(n,1),k) = z(matind(n,1),k)+bb(j1)*sq2i
                   z(matind(n,2),k) = z(matind(n,2),k)-bb(j1)*sq2i
#else
                   z(matind(n,1),k) = z(matind(n,1),k)+bb(j1) *CMPLX(sq2i,0.0)
                   z(matind(n,2),k) = z(matind(n,2),k)-bb(j1) *CMPLX(sq2i,0.0)
#endif
                ELSE
#ifdef CPP_INVERSION
                   z(matind(n,1),k) = 0.0
#else
                   z(matind(n,1),k) = CMPLX(0.0,0.0)
#endif
                ENDIF
             ENDDO
             j=j-1
          ENDIF
       ENDDO
       DEALLOCATE ( aa,bb )

    ENDIF
    !-gu
  END SUBROUTINE franza
END 
