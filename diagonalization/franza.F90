!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_franza
CONTAINS
  SUBROUTINE franza(nbasfcn,neigd, nsize,&
       l_zref,l_J,matind,nred,gw,eig,ne,hamOvlp,zMat)
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
    USE m_types
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
    TYPE(t_hamOvlp),INTENT(INOUT) :: hamOvlp
    TYPE(t_zMat),INTENT(INOUT)    :: zMat

    !     ..
    !     .. Local Scalars ..
    REAL toler,sq2i
    INTEGER j,i1,i2,j1,j2,k,nsym,jsym,matsz
    INTEGER info,i,iu,ne_a,nsize_a,n
    LOGICAL addstate
    !     ..
    !     .. Local Arrays
    REAL,    ALLOCATABLE :: aa_r(:), bb_r(:)
    REAL,    ALLOCATABLE :: zz_r(:,:)
    COMPLEX, ALLOCATABLE :: aa_c(:), bb_c(:), cwork(:)
    COMPLEX, ALLOCATABLE :: zz_c(:,:)

    REAL,    ALLOCATABLE ::  work(:), etemp(:)
    INTEGER, ALLOCATABLE ::  iwork(:),ifail(:)
    LOGICAL sort(nbasfcn)
    REAL:: lb,ub

    !to select real/complex data
    LOGICAL:: l_real

    l_real=zMat%l_real
!    IF (l_real.AND.PRESENT(a_c)) CALL juDFT_error("BUG in franza, call either with real OR complex data")
    IF (l_real) THEN
       ALLOCATE(zz_r(nbasfcn,neigd+1))
    ELSE
       ALLOCATE(zz_c(nbasfcn,neigd+1))
    ENDIF
    nsym=1
    IF (l_zref) THEN
       !
       ! separate H and S matrix in symmetric (aa,bb) and antisymmetric (a,b) part
       !
       IF (l_real) THEN
          ALLOCATE (aa_r(nred*(nred+1)/2),bb_r(nred*(nred+1)/2))
       ELSE
          ALLOCATE (aa_c(nred*(nred+1)/2),bb_c(nred*(nred+1)/2))
       END IF

       nsym=2
       i2=0
       j2=0
       DO i=1,nred
          DO j=1,i 
             i1=(matind(i,1)-1)*matind(i,1)/2+matind(j,1)
             j1=(matind(i,1)-1)*matind(i,1)/2+matind(j,2)
             i2=i2+1
             IF (l_real) THEN 
                aa_r(i2)=hamOvlp%a_r(i1)+hamOvlp%a_r(j1)
                bb_r(i2)=hamOvlp%b_r(i1)+hamOvlp%b_r(j1)
             ELSE
                aa_c(i2)=hamOvlp%a_c(i1)+hamOvlp%a_c(j1)
                bb_c(i2)=hamOvlp%b_c(i1)+hamOvlp%b_c(j1)
             END IF
             IF ((matind(i,1).NE.matind(i,2)).AND.(matind(j,1).NE.matind(j,2))) THEN
                j2=j2+1
                IF (l_real) THEN 
                   aa_r(j2)=hamOvlp%a_r(i1)-hamOvlp%a_r(j1)
                   bb_r(j2)=hamOvlp%b_r(i1)-hamOvlp%b_r(j1)
                ELSE
                   aa_c(j2)=hamOvlp%a_c(i1)-hamOvlp%a_c(j1)
                   bb_c(j2)=hamOvlp%b_c(i1)-hamOvlp%b_c(j1)
                END IF
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
          IF (l_real) THEN
             CALL CPP_BLAS_scopy(matsz,aa_r,1,hamOvlp%a_r,1)
             CALL CPP_BLAS_scopy(matsz,bb_r,1,hamOvlp%b_r,1)
          ELSE
             CALL CPP_BLAS_ccopy(matsz,aa_c,1,hamOvlp%a_c,1)
             CALL CPP_BLAS_ccopy(matsz,bb_c,1,hamOvlp%b_c,1)
          ENDIF

          ne_a=ne
          k=1
          nsize_a=nsize
          DO i=1,ne
             IF (l_real) THEN
                aa_r(i)=eig(i)
             ELSE
                aa_c(i)=CMPLX(eig(i),0.0)
             ENDIF
             IF (l_real) THEN
                DO j=1,nsize
                   bb_r(k)=zMat%z_r(j,i)
                   k=k+1
                ENDDO
             ELSE
                DO j=1,nsize
                   bb_c(k)=zMat%z_c(j,i)
                   k=k+1
                ENDDO
             ENDIF
          ENDDO
          nsize=nred

       ENDIF
       !-gu
       ! --> start with Cholesky factorization of b ( so that b = l * l^t)
       ! --> b is overwritten by l
       !
       IF (l_real) THEN
          CALL CPP_LAPACK_spptrf('U',nsize,hamOvlp%b_r,info)
       ELSE
          CALL CPP_LAPACK_cpptrf('U',nsize,hamOvlp%b_c,info)
       ENDIF
       !
       IF (info.NE.0) THEN
          WRITE (*,*) 'Error in cpptrf/spptrf: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="franza")
       ENDIF
       !
       ! --> now reduce a * z = eig * b * z to the standard form a' * z' = eig * z' 
       ! --> where a' = (l)^-1 * a * (l^t)^-1 and z' = l^t * z
       !
       IF (l_real) THEN
          CALL CPP_LAPACK_sspgst(1,'U',nsize,hamOvlp%a_r,hamOvlp%b_r,info)
       ELSE
          CALL CPP_LAPACK_chpgst(1,'U',nsize,hamOvlp%a_c,hamOvlp%b_c,info)
       ENDIF
       !
       IF (info.NE.0) THEN
          WRITE (6,*) 'Error in chpgst/sspgst: info =',info
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

       IF (l_real) THEN
          zz_r=0.0
          IF(l_J)THEN
             CALL CPP_LAPACK_sspevx('N','I','U',nsize,hamOvlp%a_r,lb,ub,1,iu,toler,ne, etemp,zz_r,nbasfcn,work,iwork,ifail,info)
          ELSE
             CALL CPP_LAPACK_sspevx('V','I','U',nsize,hamOvlp%a_r,lb,ub,1,iu,toler,ne, etemp,zz_r,nbasfcn,work,iwork,ifail,info)
          ENDIF
       ELSE
          zz_c=0.0
          ALLOCATE ( cwork(2*nbasfcn) )
          IF(l_J)THEN
             CALL CPP_LAPACK_chpevx('N','I','U',nsize,hamOvlp%a_c,lb,ub,1,iu,toler,ne, etemp,zz_c,nbasfcn,cwork,work,iwork,ifail,info)
          ELSE
             CALL CPP_LAPACK_chpevx('V','I','U',nsize,hamOvlp%a_c,lb,ub,1,iu,toler,ne, etemp,zz_c,nbasfcn,cwork,work,iwork,ifail,info)
          ENDIF
          DEALLOCATE ( cwork )
       ENDIF
       IF(addstate) THEN ! cut topmost subspace of degenerate states to avoid symmetry breaking (CF)
          iu = ne
          ne = ne - 1
          DO WHILE(ABS(etemp(ne)-etemp(iu)).LT.1d-8)
             ne = ne - 1
          ENDDO
       ENDIF
       IF (l_real) THEN
          zMat%z_r = zz_r(:,:ne)
       ELSE
          zMat%z_c = zz_c(:,:ne)
       END IF
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
             WRITE (6,*) 'Error in chpexvx/sspevx: info =',info
             WRITE (6,*) 'The following eigenvectors did not converge:'
             WRITE (6,'(30i5)') (ifail(i),i=1,ne)
             CALL juDFT_error("Diagonalization failed",calledby="franza")
          ENDIF
          DEALLOCATE ( work,iwork,ifail,etemp )
          !
          ! --> recover the generalized eigenvectors z by solving z' = l^t * z
          !
          IF (l_real) THEN
             CALL CPP_LAPACK_stptrs('U','N','N',nsize,ne,hamOvlp%b_r,zMat%z_r,nbasfcn,info)
          ELSE
             CALL CPP_LAPACK_ctptrs('U','N','N',nsize,ne,hamOvlp%b_c,zMat%z_c,nbasfcn,info)
          ENDIF
          !
          IF (info.NE.0) THEN
             WRITE (6,*) 'Error in c/stptrs: info =',info
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
          IF (l_real) THEN
             DO j=1,nsize
                hamOvlp%b_r(k)=zMat%z_r(j,i)
                k=k+1
             ENDDO
          ELSE
             DO j=1,nsize
                hamOvlp%b_c(k)=zMat%z_c(j,i)
                k=k+1
             ENDDO
          END IF
       ENDDO
       !
       ! prepare sort-array: even=.true., odd=.false.
       !
       i=1
       j=1
       eig(ne+1)    = 99.9e9
       IF (l_real) THEN
          aa_r(ne_a+1) = 99.9e9
       ELSE
          aa_c(ne_a+1) = CMPLX(99.9e9,0.0)
       ENDIF
       DO k=1,ne+ne_a
          IF (l_real) THEN
             IF (eig(i).LT.REAL(aa_r(j))) THEN
                sort(k)=.TRUE.
                i=i+1
             ELSE
                sort(k)=.FALSE.
                j=j+1
             ENDIF
          ELSE
             IF (eig(i).LT.REAL(aa_c(j))) THEN
                sort(k)=.TRUE.
                i=i+1
             ELSE
                sort(k)=.FALSE.
                j=j+1
             ENDIF
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
             IF (l_real) THEN
                zMat%z_r(n,k)=0.0
             ELSE
                zMat%z_c(n,k)=CMPLX(0.0,0.0)
             ENDIF
          ENDDO
          IF (sort(k)) THEN
             eig(k)=eig(i)
             i1=nred * (i-1)
             DO n=1,nred
                i1=i1+1
                IF (l_real) THEN
                   zMat%z_r(matind(n,1),k) = zMat%z_r(matind(n,1),k)+hamOvlp%b_r(i1)*sq2i
                   zMat%z_r(matind(n,2),k) = zMat%z_r(matind(n,2),k)+hamOvlp%b_r(i1)*sq2i
                ELSE
                   zMat%z_c(matind(n,1),k) = zMat%z_c(matind(n,1),k)+hamOvlp%b_c(i1)*CMPLX(sq2i,0.0)
                   zMat%z_c(matind(n,2),k) = zMat%z_c(matind(n,2),k)+hamOvlp%b_c(i1)*CMPLX(sq2i,0.0)
                ENDIF
             ENDDO
             i=i-1
          ELSE
             IF (l_real) THEN
                eig(k)=aa_r(j)
             ELSE
                eig(k)=REAL(aa_c(j))
             ENDIF
             j1=nsize_a * (j-1)
             DO n=1,nred
                IF (matind(n,1).NE.matind(n,2)) THEN
                   j1=j1+1
                   IF (l_real) THEN
                      zMat%z_r(matind(n,1),k) = zMat%z_r(matind(n,1),k)+bb_r(j1)*sq2i
                      zMat%z_r(matind(n,2),k) = zMat%z_r(matind(n,2),k)-bb_r(j1)*sq2i
                   ELSE
                      zMat%z_c(matind(n,1),k) = zMat%z_c(matind(n,1),k)+bb_c(j1) *CMPLX(sq2i,0.0)
                      zMat%z_c(matind(n,2),k) = zMat%z_c(matind(n,2),k)-bb_c(j1) *CMPLX(sq2i,0.0)
                   ENDIF
                ELSE
                   IF (l_real) THEN
                      zMat%z_r(matind(n,1),k) = 0.0
                   ELSE
                      zMat%z_c(matind(n,1),k) = CMPLX(0.0,0.0)
                   ENDIF
                ENDIF
             ENDDO
             j=j-1
          ENDIF
       ENDDO

    ENDIF
    !-gu
  END SUBROUTINE franza
END 
