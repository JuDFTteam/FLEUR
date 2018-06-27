!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!
! @authors: Miriam Hinzen, Gregor Michalicek
! Added MPI implementation, DW 2018
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

  !MPI
  INTERFACE 
     SUBROUTINE mpi_dchase_init( mpi_comm, n, nev, nex, xoff,yoff,xlen,ylen) BIND( c, name = 'dchase_init' )
       USE, INTRINSIC :: iso_c_binding
       INTEGER(c_int) :: mpi_comm, n, nev, nex, xoff,yoff,xlen,ylen
     END SUBROUTINE mpi_dchase_init
  END INTERFACE

  INTERFACE 
     SUBROUTINE mpi_zchase_init( mpi_comm, n, nev, nex, xoff,yoff,xlen,ylen) BIND( c, name = 'zchase_init' )
       USE, INTRINSIC :: iso_c_binding
       INTEGER(c_int) :: mpi_comm, n, nev, nex, xoff,yoff,xlen,ylen
     END SUBROUTINE mpi_zchase_init
  END INTERFACE
 
  INTERFACE 
     SUBROUTINE mpi_chase_r(h, v, ritzv, deg, tol, mode, opt ) BIND( c, name = 'dchase_solve_' )
       USE, INTRINSIC :: iso_c_binding
       REAL(c_double_complex)        :: h(*), v(*)
       INTEGER(c_int)                :: deg
       REAL(c_double)                :: ritzv(*), tol
       CHARACTER(len=1,kind=c_char)  :: mode, opt
     END SUBROUTINE mpi_chase_r
  END INTERFACE
 

  INTERFACE 
     SUBROUTINE mpi_chase_c(h, v, ritzv, deg, tol, mode, opt ) BIND( c, name = 'zchase_solve_' )
       USE, INTRINSIC :: iso_c_binding
       COMPLEX(c_double_complex)     :: h(*), v(*)
       INTEGER(c_int)                :: deg
       REAL(c_double)                :: ritzv(*), tol
       CHARACTER(len=1,kind=c_char)  :: mode, opt
     END SUBROUTINE mpi_chase_c
  END INTERFACE
  

  PRIVATE 

    INTEGER :: chase_eig_id

  PUBLIC init_chase, chase_diag

  CONTAINS

  SUBROUTINE init_chase(mpi,dimension,atoms,kpts,noco,l_real)
    USE m_types_mpimat
    USE m_types
    USE m_types_mpi
    USE m_judft
    USE m_eig66_io

    IMPLICIT NONE

    TYPE(t_mpi),               INTENT(IN)    :: mpi
    TYPE(t_dimension),         INTENT(IN)    :: dimension
    TYPE(t_atoms),             INTENT(IN)    :: atoms
    TYPE(t_kpts),              INTENT(IN)    :: kpts
    TYPE(t_noco),              INTENT(IN)    :: noco

    LOGICAL,                   INTENT(IN)    :: l_real

    INTEGER                                  :: nevd, nexd

    IF (juDFT_was_argument("-diag:chase")) THEN
       nevd = min(dimension%neigd,dimension%nvd+atoms%nlotot)
       nexd = min(max(nevd/4, 45),dimension%nvd+atoms%nlotot-nevd) !dimensioning for workspace
       chase_eig_id=open_eig(mpi%mpi_comm,DIMENSION%nbasfcn,nevd+nexd,kpts%nkpt,DIMENSION%jspd,&
                             noco%l_noco,.TRUE.,l_real,noco%l_soc,.FALSE.,mpi%n_size)
    END IF
  END SUBROUTINE init_chase

   SUBROUTINE chase_diag(hmat,smat,ikpt,jsp,iter,ne,eig,zmat)
     USE m_types_mpimat
    USE m_types
    USE m_judft
    USE iso_c_binding
    USE m_eig66_io

    !Simple driver to solve Generalized Eigenvalue Problem using the ChASE library
    IMPLICIT NONE

    CLASS(t_mat),              INTENT(INOUT) :: hmat,smat
    INTEGER,                   INTENT(IN)    :: ikpt
    INTEGER,                   INTENT(IN)    :: jsp
    INTEGER,                   INTENT(IN)    :: iter
    INTEGER,                   INTENT(INOUT) :: ne
    CLASS(t_mat), ALLOCATABLE, INTENT(OUT)   :: zmat
    REAL,                      INTENT(OUT)   :: eig(:)

    !Choose serial or parallel solver
    SELECT TYPE(hmat)
    CLASS is (t_mpimat)
       SELECT TYPE(smat)
       CLASS is (t_mpimat)
          CALL chase_diag_MPI(hmat,smat,ikpt,jsp,iter,ne,eig,zmat)
       CLASS default
          CALL judft_error("Inconsistent matrix setup")
       END SELECT
    CLASS is (t_mat)
       SELECT TYPE(smat)
       CLASS is (t_mat)
          CALL chase_diag_noMPI(hmat,smat,ikpt,jsp,iter,ne,eig,zmat)
       CLASS default
          CALL judft_error("Inconsistent matrix setup")
       END SELECT
    END SELECT
  END SUBROUTINE chase_diag
  
    SUBROUTINE chase_diag_noMPI(hmat,smat,ikpt,jsp,iter,ne,eig,zmat)

    USE m_types
    USE m_judft
    USE iso_c_binding
    USE m_eig66_io

    !Simple driver to solve Generalized Eigenvalue Problem using the ChASE library
    IMPLICIT NONE

    TYPE(t_mat),               INTENT(INOUT) :: hmat,smat
    INTEGER,                   INTENT(IN)    :: ikpt
    INTEGER,                   INTENT(IN)    :: jsp
    INTEGER,                   INTENT(IN)    :: iter
    INTEGER,                   INTENT(INOUT) :: ne
    CLASS(t_mat), ALLOCATABLE, INTENT(OUT)   :: zmat
    REAL,                      INTENT(OUT)   :: eig(:)

    INTEGER            :: i, j, nev, nex, nbands
    INTEGER            :: info

    CLASS(t_Mat),              ALLOCATABLE  :: zMatTemp
    REAL(c_double),            ALLOCATABLE  :: eigenvalues(:)

    ALLOCATE(t_mat::zmat)
    CALL zmat%alloc(hmat%l_real,hmat%matsize1,ne)

    nev = min(ne,hmat%matsize1)
    nex = min(max(nev/4, 45), hmat%matsize1-nev) !dimensioning for workspace

    ALLOCATE(eigenvalues(nev+nex))
    eigenvalues = 0.0

    ALLOCATE(t_mat::zmatTemp)
    CALL zMatTemp%alloc(hmat%l_real,hmat%matsize1,nev+nex)

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

       zMatTemp%data_r = 0.0

       do j = 1, hmat%matsize1
          do i = 1, j
             hmat%data_r(j,i) = hmat%data_r(i,j)
          end do
       end do
       if(iter.EQ.1) then
          call chase_r(hmat%data_r, hmat%matsize1, zMatTemp%data_r, eigenvalues, nev, nex, 25, 1e-6, 'R', 'S' )
       else
          CALL read_eig(chase_eig_id,ikpt,jsp,neig=nbands,eig=eigenvalues,zmat=zMatTemp)
          call chase_r(hmat%data_r, hmat%matsize1, zMatTemp%data_r, eigenvalues, nev, nex, 25, 1e-6, 'A', 'S' )
       end if

       ne = nev

       CALL write_eig(chase_eig_id,ikpt,jsp,nev+nex,nev+nex,&
                      eigenvalues(:(nev+nex)),zmat=zMatTemp)

       ! --> recover the generalized eigenvectors z by solving z' = l^t * z
       CALL dtrtrs('U','N','N',hmat%matsize1,nev,smat%data_r,smat%matsize1,zMatTemp%data_r,zmat%matsize1,info)
       IF (info.NE.0) THEN
          WRITE (6,*) 'Error in dtrtrs: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF

       DO i = 1, ne
          DO j = 1, hmat%matsize1
             zmat%data_r(j,i) = zMatTemp%data_r(j,i)
          END DO
          eig(i) = eigenvalues(i)
       END DO


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

       zMatTemp%data_c = CMPLX(0.0,0.0)

       do j = 1, hmat%matsize1
          do i = 1, j
             hmat%data_c(j,i) = conjg(hmat%data_c(i,j))
          end do
       end do

       if(iter.EQ.1) then
          call chase_c(hmat%data_c, hmat%matsize1, zMatTemp%data_c, eigenvalues, nev, nex, 25, 1e-6, 'R', 'S' )
       else
          CALL read_eig(chase_eig_id,ikpt,jsp,neig=nbands,eig=eigenvalues,zmat=zMatTemp)
          call chase_c(hmat%data_c, hmat%matsize1, zMatTemp%data_c, eigenvalues, nev, nex, 25, 1e-6, 'A', 'S' )
       end if

       ne = nev

       CALL write_eig(chase_eig_id,ikpt,jsp,nev+nex,nev+nex,&
                      eigenvalues(:(nev+nex)),zmat=zMatTemp)

       ! --> recover the generalized eigenvectors z by solving z' = l^t * z
       CALL ztrtrs('U','N','N',hmat%matsize1,nev,smat%data_c,smat%matsize1,zMatTemp%data_c,zmat%matsize1,info)
       IF (info.NE.0) THEN
          WRITE (6,*) 'Error in ztrtrs: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF

       DO i = 1, ne
          DO j = 1, hmat%matsize1
             zmat%data_c(j,i) = zMatTemp%data_c(j,i)
          END DO
          eig(i) = eigenvalues(i)
       END DO

    ENDIF
    IF (info.NE.0) CALL judft_error("Diagonalization via ChASE failed", calledby = 'chase_diag')
  END SUBROUTINE chase_diag_noMPI

  SUBROUTINE chase_diag_MPI(hmat,smat,ikpt,jsp,iter,ne,eig,zmat)
    use m_types_mpimat
    USE m_types
    USE m_judft
    USE iso_c_binding
    USE m_eig66_io
    
    !Simple driver to solve Generalized Eigenvalue Problem using the ChASE library
    IMPLICIT NONE

    TYPE(t_mpimat),               INTENT(INOUT) :: hmat,smat
    INTEGER,                   INTENT(IN)       :: ikpt
    INTEGER,                   INTENT(IN)       :: jsp
    INTEGER,                   INTENT(IN)       :: iter
    INTEGER,                   INTENT(INOUT)    :: ne
    CLASS(t_mat), ALLOCATABLE, INTENT(OUT)      :: zmat
    REAL,                      INTENT(OUT)      :: eig(:)

    INTEGER            :: i, j, nev, nex, nbands,xoff,yoff,xlen,ylen,ierr
    INTEGER            :: info

    CLASS(t_mat),                ALLOCATABLE  :: zMatTemp,chase_mat
    REAL,                        ALLOCATABLE  :: eigenvalues(:)
    include 'mpif.h'

    !Transform to standard problem using SCALAPACK
    IF (hmat%l_real) THEN
       CALL pdpotrf('U',smat%global_size1,smat%data_r,1,1,smat%blacs_desc,info)
       IF (info.NE.0) THEN
          WRITE (*,*) 'Error in pdpotrf: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF
        CALL pdsygst(1,'U',smat%global_size1,hmat%data_r,1,1,hmat%blacs_desc,smat%data_r,1,1,smat%blacs_desc,info)
        IF (info.NE.0) THEN
          WRITE (6,*) 'Error in pdsygst: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF
    ELSE
       CALL pzpotrf('U',smat%global_size1,smat%data_c,1,1,smat%blacs_desc,info)
       IF (info.NE.0) THEN
          WRITE (*,*) 'Error in pzpotrf: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF
        CALL pzhegst(1,'U',smat%global_size1,hmat%data_c,1,1,hmat%blacs_desc,smat%data_c,1,1,hmat%blacs_desc,info)
        IF (info.NE.0) THEN
          WRITE (6,*) 'Error in pzhegst: info =',info
          CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
       ENDIF
    END IF
    ! Now we are ready to set up chase
    
    nev = min(ne,hmat%global_size1)
    nex = min(max(nev/4, 45), hmat%global_size1-nev) !dimensioning for workspace
    IF (smat%l_real) THEN
       CALL mpi_dchase_init( hmat%mpi_com,hmat%global_size1, nev, nex, xoff,yoff,xlen,ylen)
    ELSE
       CALL mpi_zchase_init( hmat%mpi_com,hmat%global_size1, nev, nex, xoff,yoff,xlen,ylen)
    ENDIF
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,xlen,1,MPI_INTEGER,MPI_MAX,smat%mpi_com,ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,ylen,1,MPI_INTEGER,MPI_MAX,smat%mpi_com,ierr)
   
    ALLOCATE(t_mpimat::chase_mat)
    CALL chase_mat%init(smat%l_real,smat%global_size1,smat%global_size2,smat%mpi_com,.TRUE.,xlen,ylen)
    
    CALL chase_mat%copy(hmat,1,1) !redistribute matrix
    
    ALLOCATE(eigenvalues(nev+nex))
    eigenvalues = 0.0
    ALLOCATE(t_mpimat::zmatTemp)
    CALL zMatTemp%init(hmat%l_real,hmat%global_size1,nev+nex,MPI_COMM_SELF,.TRUE.,hmat%global_size1,nev+nex) !Generate a pseudo-distributed matrix


    IF (hmat%l_real) THEN
       IF(iter.EQ.1) THEN
          CALL chase_r(hmat%data_r, hmat%global_size1, zMatTemp%data_r, eigenvalues, nev, nex, 25, 1e-6, 'R', 'S' )
       ELSE
          CALL read_eig(chase_eig_id,ikpt,jsp,neig=nbands,eig=eigenvalues,zmat=zMatTemp)
          CALL chase_r(hmat%data_r, hmat%global_size1, zMatTemp%data_r, eigenvalues, nev, nex, 25, 1e-6, 'A', 'S' )
       END IF
    ELSE
       IF(iter.EQ.1) THEN
          CALL chase_c(hmat%data_c, hmat%global_size1, zMatTemp%data_c, eigenvalues, nev, nex, 25, 1e-6, 'R', 'S' )
       ELSE
          CALL read_eig(chase_eig_id,ikpt,jsp,neig=nbands,eig=eigenvalues,zmat=zMatTemp)
          CALL chase_c(hmat%data_c, hmat%global_size1, zMatTemp%data_c, eigenvalues, nev, nex, 25, 1e-6, 'A', 'S' )
       END IF
    ENDIF
    ne = nev
    CALL write_eig(chase_eig_id,ikpt,jsp,nev+nex,nev+nex,&
                      eigenvalues(:(nev+nex)),zmat=zMatTemp)

    CALL zmatTemp%copy(hmat,1,1) !Copy matrix into distributed form
    call zmatTemp%free()
    
    ! --> recover the generalized eigenvectors z by solving z' = l^t * z
    IF (smat%l_real) THEN
       CALL pdtrtrs('U','N','N',hmat%global_size1,nev,smat%data_r,1,1,smat%blacs_desc,&
            hmat%data_r,1,1,smat%blacs_desc,info)
    ELSE
       CALL pztrtrs('U','N','N',hmat%global_size1,nev,smat%data_c,1,1,smat%blacs_desc,&
            hmat%data_c,1,1,smat%blacs_desc,info)
    END IF
    IF (info.NE.0) THEN
       WRITE (6,*) 'Error in p?trtrs: info =',info
       CALL juDFT_error("Diagonalization failed",calledby="chase_diag")
    ENDIF

    !     Redistribute eigvec from ScaLAPACK distribution to each process
    !     having all eigenvectors corresponding to his eigenvalues as above
    !
    ALLOCATE(t_mpimat::zmat)
    CALL zmat%init(hmat%l_real,hmat%global_size1,hmat%global_size1,hmat%mpi_com,.FALSE.)
    CALL zmat%copy(hmat,1,1)

    
    DO i = 1, ne
       eig(i) = eigenvalues(i)
    END DO
  END SUBROUTINE chase_diag_MPI


#endif
END MODULE m_chase_diag
