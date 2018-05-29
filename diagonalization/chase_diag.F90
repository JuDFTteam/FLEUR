!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!
! @authors: Miriam Hinzen, Gregor Michalicek
!--------------------------------------------------------------------------------
MODULE m_chase_diag
#ifdef CPP_CHASE

IMPLICIT NONE

  interface
    subroutine chase_c( h, n, v, ritzv, nev, nex, deg, tol, mode, opt ) bind( c, name = 'zchase_' )
      use, intrinsic :: iso_c_binding
      complex(c_double_complex)     :: h(n,*), v(n,*)
      integer(c_int)                :: n, deg, nev, nex
      real(c_double)                :: ritzv(*), tol
      character(len=1,kind=c_char)  :: mode, opt
    end subroutine chase_c
  end interface

  interface 
    subroutine chase_r( h, n, v, ritzv, nev, nex, deg, tol, mode, opt ) bind( c, name = 'dchase_' )
      use, intrinsic :: iso_c_binding
      real(c_double_complex)        :: h(n,*), v(n,*)
      integer(c_int)                :: n, deg, nev, nex
      real(c_double)                :: ritzv(*), tol
      character(len=1,kind=c_char)  :: mode, opt
    end subroutine chase_r
  end interface

  CONTAINS

  SUBROUTINE chase_diag(hmat,smat,ne,eig,zmat)

    USE m_types
    USE m_judft
    USE iso_c_binding

    !Simple driver to solve Generalized Eigenvalue Problem using the ChASE library
    IMPLICIT NONE
    TYPE(t_mat),               INTENT(INOUT) :: hmat,smat
    INTEGER,                   INTENT(INOUT) :: ne
    CLASS(t_mat), ALLOCATABLE, INTENT(OUT)   :: zmat
    REAL,                      INTENT(OUT)   :: eig(:)

    INTEGER            :: i, j, nev, nex
    INTEGER            :: info

    REAL(c_double),            ALLOCATABLE  :: eigenvalues(:)

    REAL(c_double),            ALLOCATABLE  :: eigenvectors_r(:,:)
    COMPLEX(c_double_complex), ALLOCATABLE  :: eigenvectors_c(:,:)

    ALLOCATE(t_mat::zmat)
    CALL zmat%alloc(hmat%l_real,hmat%matsize1,ne)

    IF (hmat%l_real) THEN

       ! --> start with Cholesky factorization of b ( so that b = l * l^t)
       ! --> b is overwritten by l
       CALL dpotrf('U',smat%matsize1,smat%data_r,SIZE(smat%data_r,1),info)
       IF (info.NE.0) THEN
          WRITE (*,*) 'Error in dpotrf: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF

       ! --> now reduce a * z = eig * b * z to the standard form a' * z' = eig * z' 
       ! --> where a' = (l)^-1 * a * (l^t)^-1 and z' = l^t * z
       CALL dsygst(1,'U',smat%matsize1,hmat%data_r,SIZE(hmat%data_r,1),smat%data_r,SIZE(smat%data_r,1),info)
       IF (info.NE.0) THEN
          WRITE (6,*) 'Error in dsygst: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF

       ! --> solve a' * z' = eig * z' for eigenvalues eig between lb und ub

       nev = min(ne,hmat%matsize1)

       nex = min(max(nev/4, 45), hmat%matsize1-nev) !dimensioning for workspace

       ALLOCATE(eigenvectors_r(smat%matsize1,nev+nex))
       ALLOCATE(eigenvalues(nev+nex))
       eigenvectors_r = 0.0
       eigenvalues = 0.0

       do j = 1, hmat%matsize1
          do i = 1, j
             hmat%data_r(j,i) = hmat%data_r(i,j)
          end do
       end do

!       if(first_entry_franza) then
          call chase_r(hmat%data_r, hmat%matsize1, eigenvectors_r, eigenvalues, nev, nex, 25, 1e-10, 'R', 'S' )
!       else
!          call chase_r(hmat%data_r, hmat%matsize1, eigenvectors_r, eigenvalues, nev, nex, 25, 1e-10, 'A', 'S' )
!       end if

       ne = nev

       ! --> recover the generalized eigenvectors z by solving z' = l^t * z
       CALL dtrtrs('U','N','N',hmat%matsize1,nev,smat%data_r,eigenvectors_r,zmat%matsize1,info)
       IF (info.NE.0) THEN
          WRITE (6,*) 'Error in dtrtrs: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF

       DO i = 1, ne
          DO j = 1, hmat%matsize1
             zmat%data_r(j,i) = eigenvectors_r(j,i)
          END DO
          eig(i) = eigenvalues(i)
       END DO


       !TODO:  Store eigenvectors array to reuse it in next iteration



       DEALLOCATE(eigenvalues)
       DEALLOCATE(eigenvectors_r)


    ELSE

       ! --> start with Cholesky factorization of b ( so that b = l * l^t)
       ! --> b is overwritten by l
       CALL zpotrf('U',smat%matsize1,smat%data_c,SIZE(smat%data_c,1),info)
       IF (info.NE.0) THEN
          WRITE (*,*) 'Error in zpotrf: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF

       ! --> now reduce a * z = eig * b * z to the standard form a' * z' = eig * z' 
       ! --> where a' = (l)^-1 * a * (l^t)^-1 and z' = l^t * z
       CALL zhegst(1,'U',smat%matsize1,hmat%data_c,SIZE(hmat%data_c,1),smat%data_c,SIZE(smat%data_c,1),info)
       IF (info.NE.0) THEN
          WRITE (6,*) 'Error in zhegst: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF

       ! --> solve a' * z' = eig * z' for eigenvalues eig between lb und ub

       nev = min(ne,hmat%matsize1)

       nex = min(max(nev/4, 45), hmat%matsize1-nev) !dimensioning for workspace

       ALLOCATE(eigenvectors_c(smat%matsize1,nev+nex))
       ALLOCATE(eigenvalues(nev+nex))
       eigenvectors_c = CMPLX(0.0,0.0)
       eigenvalues = 0.0

       do j = 1, hmat%matsize1
          do i = 1, j
             hmat%data_c(j,i) = conjg(hmat%data_c(i,j))
          end do
       end do

!       if(first_entry_franza) then
          call chase_c(hmat%data_c, hmat%matsize1, eigenvectors_c, eigenvalues, nev, nex, 25, 1e-10, 'R', 'S' )
!       else
!          call chase_c(hmat%data_c, hmat%matsize1, eigenvectors_c, eigenvalues, nev, nex, 25, 1e-10, 'A', 'S' )
!       end if

       ne = nev

       ! --> recover the generalized eigenvectors z by solving z' = l^t * z
       CALL ztrtrs('U','N','N',hmat%matsize1,nev,smat%data_c,eigenvectors_c,zmat%matsize1,info)
       IF (info.NE.0) THEN
          WRITE (6,*) 'Error in ztrtrs: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF

       DO i = 1, ne
          DO j = 1, hmat%matsize1
             zmat%data_c(j,i) = eigenvectors_c(j,i)
          END DO
          eig(i) = eigenvalues(i)
       END DO


       !TODO:  Store eigenvectors array to reuse it in next iteration



       DEALLOCATE(eigenvalues)
       DEALLOCATE(eigenvectors_c)
    ENDIF
    IF (info.NE.0) CALL judft_error("Diagonalization via ChASE failed", calledby = 'chase_diag')
  END SUBROUTINE chase_diag
#endif
END MODULE m_chase_diag
